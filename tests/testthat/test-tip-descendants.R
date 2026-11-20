test_that(".tip_descendants is correct on trees with zero-length internal branches", {
  data(east10)
  tr <- east10[[1]]
  numtip <- length(tr$tip.label)
  
  # Создаём нулевую внутреннюю ветвь
  tr0 <- tr
  internal_edges <- which(tr0$edge[, 2] > numtip)
  set.seed(1)
  tr0$edge.length[sample(internal_edges, 1)] <- 0
  
  desc <- bGMYC4:::.tip_descendants(tr0)
  for (v in seq_len(tr0$Nnode)) {
    got <- sort(tr0$tip.label[desc[[v]]])
    want <- sort(ape::extract.clade(tr0, numtip + v)$tip.label)
    expect_identical(got, want, info = paste("node", v))
  }
})
