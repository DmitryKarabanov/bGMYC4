test_that("adaptive defaults work on large and small trees", {
  data(east10)
  tr <- east10[[1]]
  n <- length(tr$tip.label)
  small <- ape::drop.tip(tr, tr$tip.label[6:n])
  
  # Большое дерево: t2 должен покрывать [2, Nnode]
  set.seed(1)
  r_big <- bgmyc.singlephy(tr, mcmc = 100, burnin = 10, thinning = 5)
  expect_true(all(r_big$par[, 3] >= 2 & r_big$par[, 3] <= 214))
  
  # Малое дерево: дефолты не падают
  set.seed(1)
  r_small <- bgmyc.singlephy(small, mcmc = 100, burnin = 10, thinning = 5)
  expect_true(all(r_small$par[, 3] >= 2 & r_small$par[, 3] <= 4))
})

test_that("explicit t2 and start are respected (back-compat)", {
  data(east10)
  tr <- east10[[1]]
  set.seed(1)
  r <- bgmyc.singlephy(tr, mcmc = 100, burnin = 10, thinning = 5,
                       t2 = 51, start = c(1, 0.5, 50))
  expect_true(all(r$par[, 3] <= 51))
})
