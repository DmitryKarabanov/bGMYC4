\name{plot.interactive.probmat}
\alias{plot.interactive.probmat}
\title{Interactive Probability Matrix Plot for bGMYC Results}
\usage{
\method{plot}{interactive.probmat}(
  x,
  tree,
  palette = "green",
  show_tree = TRUE,
  tree_width = 0.3,
  width = 1400,
  height = 1000,
  save_html = NULL,
  custom_mode = TRUE,
  posterior_col = NULL,
  branch_gradient = NULL,
  tip_label_step = NULL,
  treeio_tree = NULL,
  ...
)
}
\arguments{
  \item{x}{A probability matrix (output from \code{spec.probmat}).}

  \item{tree}{A phylogenetic tree of class \code{phylo}.}

  \item{palette}{Color palette for heatmap: "green", "viridis", "RdYlBu", or "classic". Default: "green".}

  \item{show_tree}{Logical: whether to display phylogenetic tree alongside matrix. Default: TRUE.}

  \item{tree_width}{Relative width of tree panel (0.1–0.5) when \code{show_tree = TRUE}. Default: 0.3.}

  \item{width}{Plot width in pixels. Default: 1400.}

  \item{height}{Plot height in pixels. Default: 1000.}

  \item{save_html}{Optional path to save interactive plot as self-contained HTML file.}

  \item{custom_mode}{Logical: enable advanced vector-based visualization with posterior-colored branches? 
  Default: FALSE (backward compatible). TRUE for publication-quality output with interactive tree.}

  \item{posterior_col}{Optional name of column containing posterior probabilities in tree annotations.
  If NULL, auto-detects from c("posterior", "prob", "Posterior", "PROB").}

  \item{branch_gradient}{Optional list of color gradient control points for branch coloring.
  Format: list(c(position, R, G, B), ...). Default: red→orange→yellow→lightgreen→green ramp.}

  \item{tip_label_step}{Step for sparse tip labeling (e.g., 40 = label every 40th tip). 
  If NULL, auto-calculated as max(1, ceiling(n_tips / 40)).}

  \item{treeio_tree}{Optional pre-loaded \code{treedata} object (from \code{treeio::read.beast}).
  Avoids re-reading file when using custom_mode = TRUE.}

  \item{...}{Additional arguments passed to plotting functions.}
}
\value{
  An invisible \code{plotly} object. If \code{save_html} is specified, also saves to file.
}
\description{
  Creates an interactive HTML visualization of conspecificity probability matrix
  from bGMYC analysis, with optional posterior-colored phylogenetic tree.
  Supports two modes: standard (PNG-based tree overlay) and custom (vector-based subplot with interactive branches).
}
\examples{
\donttest{
data(east10)
set.seed(1)
result <- bgmyc.singlephy(east10[[1]], mcmc = 1000, burnin = 200, thinning = 10)
prob_mat <- spec.probmat(result)
p <- plot.interactive.probmat(prob_mat, east10[[1]])
}
}
\keyword{hplot}
\keyword{aplot}
\author{
  Qwen 3.6
  Dmitry Karabanov <dk@ibiw.ru>
}
