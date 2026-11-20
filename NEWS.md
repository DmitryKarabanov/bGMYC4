# bGMYC4 4.3.0 (2026-11-20)

## NEW FEATURES

- `fix_ultrametric()` is now an exported function. It repairs rounding-level
  noise from BEAST exports and refuses genuinely non-clock trees instead of
  silently stretching branches.
- New convergence diagnostics: `bgmyc.chains()` runs multiple independent
  chains on a single tree, `bgmyc.convergence()` computes split-R̂ and ESS,
  and `bgmyc.pool()` aggregates results across chains.
- `hpd_discrete()` computes highest posterior density intervals for discrete
  parameters (e.g., threshold counts).

## PERFORMANCE

- `spec.probmat()` is ~1700× faster (10 s → 0.006 s on 215 tips, 2000 samples).
  The species partition is a deterministic function of (tree, threshold), so
  tip descendants are computed once per tree and samples are tabulated by
  unique threshold instead of looping over every retained sample.
- `bgmyc.spec()` is ~250× faster (9.5 s → 0.04 s) via the same aggregation.
  Output is numerically equivalent to the previous implementation (max
  difference < 1e-14).

## IMPROVEMENTS

- `bgmyc.gibbs()` validates proposal widths and prior bounds, warns when
  `py2`/`pc2` exceed interpretable values, when the threshold prior is
  truncated by `Nnode`, when the chain is too short, or when burn-in eats
  almost the whole chain.
- `spec.probmat()` and `bgmyc.spec()` warn when trees contributed unequal
  numbers of samples (pooled probabilities then weight trees by sample count
  rather than equally).
- `bgmyc.spec()` validates `cmatrix`: row names must cover all tip labels,
  preventing silently wrong species-by-locality matrices.
- Numerical guards in `gmyc.5parameter.lik()`: stable handling of `0^0`,
  degenerate rates, and log-density underflow.
- `bgmyc.dataprep()` explicitly handles the root node.
- `plot.interactive.probmat()` protects graphics devices and removes unused
  arguments.

## BUG FIXES

- HPD type error in `summary.multibgmyc()`.
- Dead loop in `bgmyc.point()` when all clusters merge.
- Explicit `ape::` prefixes for functions called inside `bgmyc.dataprep()`,
  fixing "could not find function" errors in parallel workers.
- `bgmyc.multiphylo()` now forces its arguments before dispatching work to
  parallel workers. Previously, lazy argument promises (e.g. `mcmc = params$mcmc`
  passed from a calling script) were serialized to workers, where the caller's
  environment was unavailable, causing `object 'params' not found` errors when
  `cores > 1`. Sequential runs (`cores = 1`) were unaffected.
- `bgmyc.chains()` now forces its arguments before dispatching work to
  parallel workers. Previously, lazy argument promises (e.g. `mcmc = params$mcmc`
  passed from a calling script) were serialized to workers, where the caller's
  environment was unavailable, causing `object 'params' not found` errors when
  `cores > 1`. Sequential runs (`cores = 1`) were unaffected.

## DOCUMENTATION

- Full migration to roxygen2: all `.Rd` files are generated from `@export`
  tags, and `NAMESPACE` is rebuilt automatically.
- New package-level help page `?bGMYC4` with workflow overview.
- `inst/CITATION` now includes both the package and the foundational
  Reid & Carstens (2012) paper.
- README updated: correct bibliographic reference (BMC Evolutionary Biology
  12:196, not Mol. Ecol. Res.), honest description of parallelism (uses
  `parallel::makeCluster`, not `future`), and new section "Choosing ppcutoff"
  explaining single-linkage mechanics and stability plateaus.

## DEPENDENCIES

- `parallel` added to `Imports` (used for multi-tree parallelization).
- `future` and `future.apply` moved from `Imports` to `Suggests` (the package
  does not use them; they are only needed by the interactive script in
  `inst/scripts/`).

---

# bGMYC4 4.2.0 (2026-05-15)

*Initial release on CRAN (published as 4.1.6 due to version conflict).*

- Bayesian GMYC implementation for species delimitation.
- MCMC sampling of Yule/coalescent threshold with flexible priors.
- Multi-tree analysis with phylogenetic uncertainty pooling.
- Interactive diagnostics: trace plots, acceptance rates, ESS.
- Conspecificity probability matrices and point estimates.
- BEAST2-ready workflow with Nexus annotation handling.
