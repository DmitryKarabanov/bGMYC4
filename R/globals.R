# Suppress R CMD check NOTEs for variables used in NSE contexts
# (e.g., within dplyr pipelines in vignettes/examples).
# These are NOT written to .GlobalEnv.
utils::globalVariables(c(
  "bt", "sb", "numnod", "numtip", "numall", "nthresh", "internod",
  "nesting", "nested", "bt.ancs", "mrca.nodes", "nod.types", "n",
  "list.s.nod", "list.i.mat"
))
