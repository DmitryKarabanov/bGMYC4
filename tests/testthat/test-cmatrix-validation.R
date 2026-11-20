test_that("bgmyc.spec validates cmatrix before use", {
  data(east10)
  set.seed(1)
  res <- bgmyc.singlephy(east10[[1]], mcmc = 100, burnin = 10, thinning = 5)
  tips <- east10[[1]]$tip.label
  
  # Без rownames — понятная ошибка
  bad <- matrix(1, nrow = length(tips), ncol = 2)
  expect_error(bgmyc.spec(res, cmatrix = bad), "rownames covering all tree tip labels")
  
  # С неполными rownames — понятная ошибка
  bad2 <- bad[1:5, , drop = FALSE]
  rownames(bad2) <- tips[1:5]
  expect_error(bgmyc.spec(res, cmatrix = bad2), "rownames covering all tree tip labels")
  
  # Хорошая cmatrix — работает и возвращает матрицы
  good <- matrix(1, nrow = length(tips), ncol = 2,
                 dimnames = list(tips, c("locA", "locB")))
  ok <- bgmyc.spec(res, cmatrix = good)
  expect_true(is.list(ok$matrices))
  expect_equal(length(ok$matrices), nrow(res$par))
})
