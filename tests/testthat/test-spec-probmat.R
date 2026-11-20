test_that("spec.probmat gives diagonal 1 and symmetric matrix", {
  data(east10)
  set.seed(1)
  r <- bgmyc.singlephy(east10[[1]], mcmc = 100, burnin = 10, thinning = 5)
  pm <- spec.probmat(r)
  expect_true(all(abs(diag(pm) - 1) < 1e-10))
  expect_true(all(abs(pm - t(pm)) < 1e-10))
  expect_true(all(pm >= -1e-12 & pm <= 1 + 1e-12))
})

test_that("spec.probmat on multibgmyc matches single-tree result for identical trees", {
  data(east10)
  tr <- east10[[1]]
  set.seed(1)
  r1 <- bgmyc.singlephy(tr, mcmc = 100, burnin = 10, thinning = 5)
  set.seed(1)
  r2 <- bgmyc.singlephy(tr, mcmc = 100, burnin = 10, thinning = 5)
  pm1 <- spec.probmat(r1)
  pm_multi <- spec.probmat(structure(list(r1, r2), class = "multibgmyc"))
  expect_true(max(abs(pm_multi - pm1)) < 1e-10)
})
