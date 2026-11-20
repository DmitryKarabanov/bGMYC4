#' Fix rounding-level non-ultrametricity
#'
#' Rounds branch lengths to `digits` decimals; stops if the tree remains
#' non-ultrametric (re-export from BEAST or check the clock model instead).
#' Preprocessing only — the GMYC likelihood is unaffected.
#' @param tr A `"phylo"` object.
#' @param digits Number of decimals for rounding (default 8).
#' @return The (possibly rounded) tree.
#' @export
fix_ultrametric <- function(tr, digits = 8) {
  if (!ape::is.ultrametric(tr)) {
    tr$edge.length <- round(tr$edge.length, digits)
    if (!ape::is.ultrametric(tr))
      stop("Tree is still non-ultrametric after rounding; re-export from BEAST or check the clock model.")
  }
  tr
}
