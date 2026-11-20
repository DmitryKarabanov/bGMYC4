## Submission 4.3.0

This is a major update to the existing CRAN package bGMYC4 (4.1.6).

### Changes since 4.1.6

- **Performance:** `spec.probmat()` is ~1700× faster and `bgmyc.spec()` ~250× faster
  via threshold-aggregated post-processing (output numerically equivalent, max diff < 1e-14).
- **New features:** `fix_ultrametric()` (now exported), `bgmyc.chains()`,
  `bgmyc.convergence()`, `bgmyc.pool()`, `hpd_discrete()`.
- **Robustness:** Input validation in `bgmyc.gibbs()`, numerical guards in
  `gmyc.5parameter.lik()`, explicit `ape::` prefixes for parallel workers,
  `cmatrix` validation in `bgmyc.spec()`.
- **Bug fixes:** HPD type error in `summary.multibgmyc()`; dead loop in `bgmyc.point()`;
  removed unused arguments from `plot.interactive.probmat()`; lazy argument promises
  in `bgmyc.multiphylo()` are now forced before dispatching work to parallel workers
  (fixes `object 'params' not found` errors with `cores > 1`).
- **Documentation:** Full migration to roxygen2; corrected bibliographic references
  (Reid & Carstens 2012 is BMC Evolutionary Biology 12:196, not Mol. Ecol. Res.);
  new section "Choosing ppcutoff" with empirical guidance.
- **Dependencies:** `parallel` added to Imports; `future`/`future.apply` moved to Suggests
  (the package itself does not use them; only the interactive script in `inst/scripts/` does).

### Test environments

- R 4.5.1, Windows 11 x64 (local)
- R CMD check --as-cran: 0 errors, 0 warnings, 0 notes