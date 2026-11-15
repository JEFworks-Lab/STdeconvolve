// Variational EM for latent Dirichlet allocation.
//
// Fits the same model as the VEM method of the `topicmodels` package (LDA-C by
// David Blei): symmetric Dirichlet prior `alpha` on the per-document topic
// proportions, a free topic-word matrix `beta`, the same initialisation, the
// same M-step (including the Newton update of alpha) and the same EM stopping
// rule. The E-step differs in how the per-document optimum is computed:
//
//  * the word-level variational parameters phi are never stored. For a given
//    gamma the optimal phi is phi_nk = e_k * beta_kw / sum_j e_j * beta_jw with
//    e_k = exp(digamma(gamma_k) - digamma(sum gamma)), so one fixed-point
//    update of gamma only needs two passes over K per non-zero count and no
//    transcendental function calls,
//  * this fixed-point iteration is accelerated (SQUAREM), with safeguards
//    that keep it on the path of the plain iteration,
//  * documents are processed in parallel.
//
// All sums that end up in the model are accumulated in a fixed order, so the
// result does not depend on the number of threads.

#include <Rcpp.h>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <thread>
#include <vector>

#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

namespace {

// Smallest word probability; topicmodels stores log(beta) = -100 for words
// that received no mass in a topic.
const double BETA_FLOOR = 3.720075976020836e-44; // exp(-100)

// MT19937 seeded the way topicmodels does it (cokus.c), so that the random
// initialisation of beta for a given seed matches topicmodels::LDA.
class CokusMT {
public:
  explicit CokusMT(uint32_t seed) : pos_(N) {
    uint32_t x = seed | 1U;
    state_[0] = x;
    for (int j = 1; j < N; ++j) {
      x *= 69069U;
      state_[j] = x;
    }
  }
  double unif() { return static_cast<double>(next()) / 4294967296.0; }

private:
  static const int N = 624;
  static const int M = 397;
  uint32_t state_[N];
  int pos_;

  uint32_t next() {
    if (pos_ >= N) reload();
    uint32_t y = state_[pos_++];
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    return y ^ (y >> 18);
  }
  void reload() {
    for (int i = 0; i < N; ++i) {
      const uint32_t s1 = state_[(i + 1) % N];
      const uint32_t y = (state_[i] & 0x80000000U) | (s1 & 0x7FFFFFFFU);
      state_[i] = state_[(i + M) % N] ^ (y >> 1) ^ ((s1 & 1U) ? 0x9908B0DFU : 0U);
    }
    pos_ = 0;
  }
};

// digamma for x > 0: recurrence up to x >= 6, then the asymptotic series
// (absolute error below 1e-12). Safe to call from worker threads.
inline double digamma_pos(double x) {
  double r = 0.0;
  while (x < 6.0) {
    r -= 1.0 / x;
    x += 1.0;
  }
  const double f = 1.0 / (x * x);
  return r + std::log(x) - 0.5 / x -
    f * (1.0 / 12 - f * (1.0 / 120 - f * (1.0 / 252 - f * (1.0 / 240 -
    f * (1.0 / 132 - f * (691.0 / 32760))))));
}

// Dot product with a fixed association order (4 accumulators), so that the
// compiler can vectorise it without -ffast-math and results are reproducible.
inline double dot_k(const double* a, const double* b, int K) {
  double s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
  int k = 0;
  for (; k + 3 < K; k += 4) {
    s0 += a[k] * b[k];
    s1 += a[k + 1] * b[k + 1];
    s2 += a[k + 2] * b[k + 2];
    s3 += a[k + 3] * b[k + 3];
  }
  for (; k < K; ++k) s0 += a[k] * b[k];
  return (s0 + s1) + (s2 + s3);
}

// Run f(begin, end, thread_id) over [0, n) in chunks handed out dynamically.
template <class F>
void parallel_for(int n, int nthreads, int chunk, F f) {
  if (nthreads <= 1 || n <= chunk) {
    f(0, n, 0);
    return;
  }
  std::atomic<int> next(0);
  auto worker = [&](int tid) {
    for (;;) {
      const int b = next.fetch_add(chunk);
      if (b >= n) break;
      f(b, std::min(n, b + chunk), tid);
    }
  };
  std::vector<std::thread> threads;
  threads.reserve(nthreads - 1);
  for (int t = 1; t < nthreads; ++t) threads.emplace_back(worker, t);
  worker(0);
  for (auto& t : threads) t.join();
}

// Newton iterations for the symmetric Dirichlet parameter (LDA-C opt_alpha).
double opt_alpha(double ss, int D, int K) {
  double init_a = 100.0, log_a = std::log(init_a), df;
  int iter = 0;
  do {
    ++iter;
    double a = std::exp(log_a);
    if (std::isnan(a)) {
      init_a *= 10.0;
      a = init_a;
      log_a = std::log(a);
    }
    df = D * (K * R::digamma(K * a) - K * R::digamma(a)) + ss;
    const double d2f = D * (K * K * R::trigamma(K * a) - K * R::trigamma(a));
    log_a -= df / (d2f * a + df);
  } while (std::fabs(df) > 1e-5 && iter < 1000);
  return std::exp(log_a);
}

} // namespace

// [[Rcpp::export]]
Rcpp::List lda_vem_cpp(const Rcpp::IntegerVector& csc_p,
                       const Rcpp::IntegerVector& csc_i,
                       const Rcpp::NumericVector& csc_x,
                       int D, int V, int K,
                       double alpha, bool estimate_alpha,
                       Rcpp::Nullable<Rcpp::NumericMatrix> init_beta, int seed,
                       int em_max_iter, double em_tol,
                       int var_max_iter, double var_tol,
                       int nthreads, bool accelerate = true) {
  const int nnz = csc_x.size();
  const int* cp = csc_p.begin();
  const int* ci = csc_i.begin();
  const double* cx = csc_x.begin();
  if (nthreads < 1) nthreads = 1;
  if (var_max_iter < 0) var_max_iter = std::numeric_limits<int>::max();

  // ---- document-major (CSR) copy of the corpus; the input is word-major (CSC)
  std::vector<int> rp(D + 1, 0), rj(nnz);
  std::vector<double> rx(nnz), total(D, 0.0);
  for (int n = 0; n < nnz; ++n) ++rp[ci[n] + 1];
  for (int d = 0; d < D; ++d) rp[d + 1] += rp[d];
  {
    std::vector<int> fill(rp.begin(), rp.end() - 1);
    for (int w = 0; w < V; ++w) {
      for (int n = cp[w]; n < cp[w + 1]; ++n) {
        const int d = ci[n], at = fill[d]++;
        rj[at] = w;
        rx[at] = cx[n];
        total[d] += cx[n];
      }
    }
  }

  // ---- model. B[w * K + k] = p(word w | topic k)
  const size_t VK = static_cast<size_t>(V) * K, DK = static_cast<size_t>(D) * K;
  std::vector<double> B(VK), SS(VK);
  if (init_beta.isNotNull()) {
    Rcpp::NumericMatrix ib(init_beta);
    if (ib.nrow() != K || ib.ncol() != V) Rcpp::stop("`init_beta` must be a K x V matrix");
    for (size_t i = 0; i < VK; ++i) B[i] = std::max(ib[i], BETA_FLOOR);
  } else {
    CokusMT rng(static_cast<uint32_t>(seed));
    for (int k = 0; k < K; ++k) {
      double tot = 0.0;
      for (int w = 0; w < V; ++w) {
        const double v = 1.0 / V + rng.unif();
        B[static_cast<size_t>(w) * K + k] = v;
        tot += v;
      }
      for (int w = 0; w < V; ++w) B[static_cast<size_t>(w) * K + k] /= tot;
    }
  }

  // ---- variational parameters and per-document summaries
  std::vector<double> G(DK);        // gamma[d * K + k]
  std::vector<double> E(DK);        // exp(E_q[log theta_dk]) at the final gamma
  std::vector<double> part(D);      // part of the bound that needs no lgamma
  std::vector<double> elogsum(D);   // sum_k E_q[log theta_dk]
  std::vector<double> llh(D);
  int max_len = 0;
  for (int d = 0; d < D; ++d) max_len = std::max(max_len, rp[d + 1] - rp[d]);
  // per thread: two gammas, word weights and the beta columns of one document
  std::vector<std::vector<double> > work(
      nthreads, std::vector<double>(2 * static_cast<size_t>(K) + max_len +
                                    static_cast<size_t>(K) * max_len));
  std::vector<double> iters_by_thread(nthreads, 0.0);

  int vmax = var_max_iter;

  // E-step: optimise gamma for every document given beta and alpha
  auto e_step = [&]() {
    parallel_for(D, nthreads, 8, [&](int begin, int end, int tid) {
      double* g1 = work[tid].data();
      double* g2 = g1 + K;
      double* wt = g2 + K;
      double* Bd = wt + max_len;
      double niter = 0.0;
      for (int d = begin; d < end; ++d) {
        double* g = &G[static_cast<size_t>(d) * K];
        double* e = &E[static_cast<size_t>(d) * K];
        const int n0 = rp[d], len = rp[d + 1] - rp[d];
        const double* cnt = &rx[n0];
        // topic-word probabilities of the words in this document, one
        // contiguous column per topic, so that the loops below run over the
        // words of the document and vectorise
        for (int n = 0; n < len; ++n) {
          const double* b = &B[static_cast<size_t>(rj[n0 + n]) * K];
          for (int k = 0; k < K; ++k) Bd[static_cast<size_t>(k) * len + n] = b[k];
        }
        // wt[n] = sum_k e[k] * beta[k, word n]
        auto word_weights = [&](const double* ek) {
          for (int n = 0; n < len; ++n) wt[n] = 0.0;
          for (int k = 0; k < K; ++k) {
            const double v = ek[k];
            const double* RESTRICT col = Bd + static_cast<size_t>(k) * len;
            double* RESTRICT w = wt;
            for (int n = 0; n < len; ++n) w[n] += v * col[n];
          }
        };
        // Every document starts from the same uninformative gamma. Starting
        // from the gamma of the previous EM iteration would save updates, but
        // topics that were (nearly) removed from a document while beta was
        // still poor can then not come back, which leads to worse models.
        const double g0 = alpha + total[d] / K;
        for (int k = 0; k < K; ++k) g[k] = g0;
        const double gsum = K * alpha + total[d]; // invariant of the update
        const double dgsum = digamma_pos(gsum);
        // One fixed-point update of gamma: out = alpha + sum_n count_n * phi_n,
        // with phi_n the optimal word-topic responsibilities given `in`.
        // Returns the L1 distance between `in` and `out`.
        auto update = [&](const double* in, double* out) {
          for (int k = 0; k < K; ++k) e[k] = std::exp(digamma_pos(in[k]) - dgsum);
          word_weights(e);
          for (int n = 0; n < len; ++n) wt[n] = cnt[n] / std::max(wt[n], 1e-300);
          double change = 0.0;
          for (int k = 0; k < K; ++k) {
            out[k] = alpha + e[k] * dot_k(wt, Bd + static_cast<size_t>(k) * len, len);
            change += std::fabs(out[k] - in[k]);
          }
          niter += 1.0;
          return change;
        };
        // The plain iteration g <- update(g) converges linearly and often
        // slowly, so it is accelerated with SQUAREM (Varadhan & Roland 2008):
        // two updates give a step r and its change v, from which a longer step
        // is extrapolated, followed by one stabilising update. The objective
        // has several local optima when alpha < 1, and a long step can end up
        // in a different (in practice worse) one than the plain iteration, by
        // removing a topic from the document. Two safeguards keep the
        // accelerated iteration on the path of the plain one: the step length
        // is only allowed to grow gradually, and no step may take away more
        // than a fraction of the pseudo-counts of any topic. Convergence is
        // always judged by the size of a plain update.
        const double thresh = var_tol * gsum;
        double prev_res = std::numeric_limits<double>::infinity();
        double step_max = 1.0;
        const double step_growth = 4.0, shrink_limit = 0.5;
        for (int it = 0; it < vmax; it += 3) {
          const double res = update(g, g1);
          if (res < thresh) {
            for (int k = 0; k < K; ++k) g[k] = g1[k];
            break;
          }
          const double res2 = update(g1, g2);
          bool plain = res2 < thresh || !accelerate;
          if (res > prev_res) {
            // the last extrapolation did not help: be more careful
            step_max = std::max(1.0, step_max / step_growth);
            plain = true;
          }
          prev_res = res;
          double step = 1.0;
          if (!plain) {
            double rr = 0.0, vv = 0.0;
            for (int k = 0; k < K; ++k) {
              const double r = g1[k] - g[k], v = (g2[k] - g1[k]) - r;
              rr += r * r;
              vv += v * v;
            }
            step = (vv > 0.0) ? std::sqrt(rr / vv) : 1.0;
            if (step >= step_max) {
              step = step_max;
              step_max *= step_growth;
            }
            // shorten the step until no topic loses more than a fraction of
            // the pseudo-counts it has after the two plain updates
            while (step > 1.0) {
              bool ok = true;
              for (int k = 0; k < K; ++k) {
                const double r = g1[k] - g[k], v = (g2[k] - g1[k]) - r;
                const double x = g[k] + 2.0 * step * r + step * step * v;
                if (x < alpha + shrink_limit * (g2[k] - alpha)) { ok = false; break; }
              }
              if (ok) break;
              step = (step < 1.05) ? 1.0 : 0.5 * (step + 1.0);
            }
          }
          if (!(step > 1.0)) {
            // keep the result of the two plain updates
            for (int k = 0; k < K; ++k) g[k] = g2[k];
            if (res2 < thresh) break;
            prev_res = res2;
            continue;
          }
          for (int k = 0; k < K; ++k) {
            const double r = g1[k] - g[k], v = (g2[k] - g1[k]) - r;
            g1[k] = g[k] + 2.0 * step * r + step * step * v;
          }
          update(g1, g); // stabilising step
        }
        // bound at the final gamma (with the optimal phi for that gamma)
        double p = 0.0, es = 0.0;
        for (int k = 0; k < K; ++k) {
          const double el = digamma_pos(g[k]) - dgsum;
          e[k] = std::exp(el);
          es += el;
          p += (alpha - g[k]) * el;
        }
        word_weights(e);
        for (int n = 0; n < len; ++n) p += cnt[n] * std::log(std::max(wt[n], 1e-300));
        part[d] = p;
        elogsum[d] = es;
      }
      iters_by_thread[tid] += niter;
    });
    // lgamma terms on the main thread (lgamma is not guaranteed thread-safe)
    const double c0 = R::lgammafn(K * alpha) - K * R::lgammafn(alpha);
    double lik = 0.0;
    for (int d = 0; d < D; ++d) {
      const double* g = &G[static_cast<size_t>(d) * K];
      double v = part[d] + c0 - R::lgammafn(K * alpha + total[d]);
      for (int k = 0; k < K; ++k) v += R::lgammafn(g[k]);
      llh[d] = v;
      lik += v;
    }
    return lik;
  };

  // M-step for beta: expected topic-word counts, one word at a time
  auto m_step_beta = [&]() {
    parallel_for(V, nthreads, 4, [&](int begin, int end, int) {
      for (int w = begin; w < end; ++w) {
        const double* b = &B[static_cast<size_t>(w) * K];
        double* ss = &SS[static_cast<size_t>(w) * K];
        for (int k = 0; k < K; ++k) ss[k] = 0.0;
        for (int n = cp[w]; n < cp[w + 1]; ++n) {
          const double* e = &E[static_cast<size_t>(ci[n]) * K];
          const double s = cx[n] / std::max(dot_k(e, b, K), 1e-300);
          for (int k = 0; k < K; ++k) ss[k] += s * e[k];
        }
        for (int k = 0; k < K; ++k) ss[k] *= b[k];
      }
    });
    std::vector<double> tot(K, 0.0);
    for (int w = 0; w < V; ++w)
      for (int k = 0; k < K; ++k) tot[k] += SS[static_cast<size_t>(w) * K + k];
    for (int w = 0; w < V; ++w)
      for (int k = 0; k < K; ++k) {
        const size_t i = static_cast<size_t>(w) * K + k;
        const double v = SS[i] / tot[k];
        B[i] = (v > BETA_FLOOR) ? v : BETA_FLOOR; // also catches NaN
      }
  };

  // ---- EM
  std::vector<double> logLiks;
  int iter = 0;
  double converged = 1.0, lik_old = 0.0;
  while ((converged < 0 || converged > em_tol || iter <= 2) && iter < em_max_iter) {
    ++iter;
    const double lik = e_step();
    m_step_beta();
    if (estimate_alpha) {
      double ss = 0.0;
      for (int d = 0; d < D; ++d) ss += elogsum[d];
      alpha = opt_alpha(ss, D, K);
    }
    converged = (lik_old - lik) / lik_old;
    if (converged < 0 && vmax < std::numeric_limits<int>::max() / 2) vmax *= 2;
    lik_old = lik;
    logLiks.push_back(lik);
    Rcpp::checkUserInterrupt();
  }

  // ---- final inference with the fitted model
  e_step();

  Rcpp::NumericMatrix beta_out(K, V), gamma_out(D, K);
  for (size_t i = 0; i < VK; ++i) beta_out[i] = std::log(B[i]);
  for (int d = 0; d < D; ++d)
    for (int k = 0; k < K; ++k)
      gamma_out[d + static_cast<size_t>(D) * k] = G[static_cast<size_t>(d) * K + k];

  // most likely topic of every non-zero count (document-major order)
  Rcpp::IntegerVector wa_i(nnz), wa_j(nnz);
  Rcpp::NumericVector wa_v(nnz);
  for (int d = 0; d < D; ++d) {
    const double* e = &E[static_cast<size_t>(d) * K];
    for (int n = rp[d]; n < rp[d + 1]; ++n) {
      const double* b = &B[static_cast<size_t>(rj[n]) * K];
      int best = 0;
      double bv = e[0] * b[0];
      for (int k = 1; k < K; ++k) {
        const double v = e[k] * b[k];
        if (v > bv) { bv = v; best = k; }
      }
      wa_i[n] = d + 1;
      wa_j[n] = rj[n] + 1;
      wa_v[n] = best + 1;
    }
  }

  double var_iters = 0.0;
  for (int t = 0; t < nthreads; ++t) var_iters += iters_by_thread[t];

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta_out,
    Rcpp::Named("gamma") = gamma_out,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("loglikelihood") = Rcpp::wrap(llh),
    Rcpp::Named("iter") = iter,
    Rcpp::Named("logLiks") = Rcpp::wrap(logLiks),
    Rcpp::Named("wa_i") = wa_i,
    Rcpp::Named("wa_j") = wa_j,
    Rcpp::Named("wa_v") = wa_v,
    Rcpp::Named("var_iters") = var_iters);
}
