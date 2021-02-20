times <- c(
  6.63,
  21.51,
  15.96,
  11.55,
  14.11,
  7.89,
  8.18,
  8.02,
  13.12,
  2.28,
  3.49,
  1.4,
  14.48,
  14.41,
  14.91,
  14.89,
  15.16,
  14.27
)

ndocs <- c(
  260,
  265,
  279,
  268,
  260,
  265,
  279,
  268,
  427,
  308,
  224,
  320,
  256,
  256,
  256,
  256,
  256,
  256
)

ngenes <- c(
  220,
  642,
  430,
  307,
  149,
  149,
  149,
  149,
  525,
  185,
  194,
  71,
  130,
  130,
  130,
  130,
  130,
  130
)

nterms <- c(
  185986,
  523998,
  220241,
  144161,
  242191,
  423128,
  173273,
  156392,
  190617,
  60858,
  50086,
  21052,
  1314351,
  1452816,
  1511123,
  1444423,
  1328152,
  1141268
)

largest_k <- c(
  35,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75
)

optimal_k <- c(
  12,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  75,
  20,
  75,
  10,
  75,
  75,
  75,
  75,
  75,
  75
)

samples <- c(
  "mob",
  "mob_rep1",
  "mob_rep2",
  "mob_rep3",
  "mob_common",
  "mob_rep1_common",
  "mob_rep2_common",
  "mob_rep3_common",
  "pdac_a1",
  "pdac_a2",
  "pdac_b1",
  "pdac_b2",
  "bregma04",
  "bregma09",
  "bregma14",
  "bregma19",
  "bregma24",
  "bregma29"
)

# number indicates max k that was fitted
# Common indicates common set of genes used in corpus
# the "75" was 10 to 75 by 5 for fitted K's
# the mob 35 was 2 to 35 by 1 for fitted K's
groups <- c(
  "mob35",
  "mob75",
  "mob75",
  "mob75",
  "mobCommon75",
  "mobCommon75",
  "mobCommon75",
  "mobCommon75",
  "pdac75",
  "pdac75",
  "pdac75",
  "pdac75",
  "merfish75",
  "merfish75",
  "merfish75",
  "merfish75",
  "merfish75",
  "merfish75"
)