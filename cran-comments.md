## Test environments

* local: Windows, R 4.5.x
* win-builder: r-devel, r-release, r-oldrel

## R CMD check results

0 errors | 0 warnings | 1 note

    checking for future file timestamps ... NOTE
      unable to verify current time

This note is environmental (the check could not reach an external time
server) and does not reflect a problem with the package.

## Update notes

This is an update to the existing CRAN package bGMYC4 (currently 4.1.6).
Version 4.3.0 addresses reviewer comments: it introduces adaptive default
priors, repairs `fix_ultrametric()`, fixes tip-order handling in
`spec.probmat()`, strengthens input validation, and adds a `testthat` suite.
All changes are backward compatible except that the default values of `t2`
and `start` are now adaptive (documented under "Breaking changes" in NEWS.md).

## Downstream dependencies

There are no downstream dependencies.