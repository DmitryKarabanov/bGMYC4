# bGMYC4 4.3.0

## Breaking changes

* `bgmyc.singlephy()`, `bgmyc.multiphylo()` and `bgmyc.chains()` now use
  adaptive defaults: `t2 = NULL` sets the threshold prior to the full range
  `[max(2, t1), Nnode]`, and `start = NULL` starts the chain at the midpoint
  of each prior interval. Explicit `t2` / `start` values behave as before.
  This replaces the old fixed defaults (`t2 = 51`, `start = c(1, 0.5, 50)`),
  which silently truncated the threshold prior on trees with more than ~52 nodes.

## Major improvements

* `fix_ultrametric()` now genuinely repairs rounding-level non-ultrametricity:
  branch lengths are rounded and all tips are then re-aligned to the deepest
  tip by extending only terminal branches, within a budget set by the rounding
  precision. Genuinely non-clock trees are still rejected with a clear message.

## Bug fixes

* `spec.probmat()` aligns tip order across trees before pooling, so pooled
  conspecificity probabilities are correct even when trees differ in tip order.
* `bgmyc.spec()` validates the `cmatrix` argument up front and returns a clear
  error instead of a cryptic `subscript out of bounds`.
* `.tip_descendants()` uses an order-robust traversal, so descendant sets are
  correct even when internal branches have zero length.
* `bgmyc.gibbs()`: removed duplicated threshold-prior warnings.
* `bgmyc.dataprep()`: trees with fewer than 3 tips are rejected with a clear error.
* `bgmyc.lik()`: returns `-Inf` instead of erroring when the threshold is out of range.
* `bgmyc.singlephy()`: the retained-samples message now reports an integer count.
* `bgmyc.convergence()`: the `flag` column is no longer `NA` when `mcmcse` is
  not installed.

## Documentation and tests

* Added `set.seed()` to all stochastic examples for reproducibility.
* Clarified the provenance of the bundled `east10` dataset.
* `CITATION` now cites the package (Karabanov) and the foundational method
  (Reid & Carstens 2012); Noah Reid is credited as contributor and copyright
  holder of the original bGMYC code.
* Removed unused `Suggests` (`future`, `future.apply`, `htmltools`).
* Added a `testthat` suite (230 checks) covering adaptive defaults,
  `fix_ultrametric()`, `spec.probmat()`, descendant traversal and `cmatrix` validation.

---

# bGMYC4 4.2.0 (2026-05-15)

*Initial release on CRAN (published as 4.1.6 due to version conflict).*

- Bayesian GMYC implementation for species delimitation.
- MCMC sampling of Yule/coalescent threshold with flexible priors.
- Multi-tree analysis with phylogenetic uncertainty pooling.
- Interactive diagnostics: trace plots, acceptance rates, ESS.
- Conspecificity probability matrices and point estimates.
- BEAST2-ready workflow with Nexus annotation handling.
