## how to find the local min quickly
x <- seq(from = -5, to = 10, by = 1)
f <- function(x) {x-5 + (x-5)^2}
plot(x, f(x), type='l')

## find the smallest y without seeing all ys
## can only ask for y for given x
optSearch <- function (f, interval,
                      lower = min(interval),
                      mid = round((min(interval) + max(interval))/2),
                      upper = max(interval),
                      minimum = TRUE,
                      tol = .Machine$double.eps^0.25)
{
  #if (minimum) {
  #

  if(lower == mid | mid == upper) {
    ## optimum found
    return (c(mid, f(mid)))
  }

  fl <- f(lower)
  fm <- f(mid)
  fu <- f(upper)

  print(lower)
  print(fl)
  print(mid)
  print(fm)
  print(upper)
  print(fu)

  if (fl < fu) {
    print('searching lower half')
    optSearch(f, interval = c(lower, mid))
  } else {
    print('searching upper half')
    optSearch(f, interval = c(mid, upper))
  }

  #
  #} else {
  #
  #}
}

results <- optSearch(f, interval = range(x))

plot(x, f(x), type='l')
points(results[1], results[2], col='red')
