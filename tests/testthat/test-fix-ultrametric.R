test_that("fix_ultrametric repairs rounding noise and keeps topology", {
  data(east10)
  tr <- east10[[1]]
  fixed <- fix_ultrametric(tr)
  expect_true(ape::is.ultrametric(fixed))
  expect_true(ape::is.ultrametric(fixed, tol = 1e-6))
  expect_identical(fixed$tip.label, tr$tip.label)
  # dataprep принимает исправленное дерево
  expect_silent(dat <- bgmyc.dataprep(fixed))
})

test_that("fix_ultrametric rejects genuinely non-clock trees", {
  data(east10)
  bad <- east10[[1]]
  bad$edge.length[1] <- bad$edge.length[1] + 1
  expect_error(fix_ultrametric(bad), "non-clock")
})
