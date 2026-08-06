#' @export
plot.bgmycprobmat <- function(x, tree,
                              palette = "classic",
                              save_pdf = NULL,
                              legend_cex = 0.5,
                              ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  # Data preparation
  probmat <- x
  ntip <- length(tree$tip.label)
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  edge <- tree$edge[which(tree$edge[, 2] < (ntip + 1)), ]
  tiplab <- tree$tip.label[edge[, 2]]
  orderedmat <- probmat[tiplab, ]
  orderedmat <- orderedmat[, tiplab]
  
  # PDF export setup
  if (!is.null(save_pdf)) {
    dir.create(dirname(save_pdf), showWarnings = FALSE, recursive = TRUE)
    grDevices::pdf(file = save_pdf, width = 10, height = 8, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  
  # Palette selection
  breaks <- c(0, 0.05, 0.49999999, 0.899999999, 0.9499999999, 1)
  if (palette == "viridis") {
    cols <- grDevices::colorRampPalette(c("#440154", "#31688E", "#35B779", "#FDE725"))(5)
  } else if (palette == "RdYlBu") {
    cols <- grDevices::colorRampPalette(c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B"))(5)
  } else {
    cols <- grDevices::heat.colors(5)
  }
  
  # Plot tree (left panel)
  par(fig = c(0, 0.5, 0, 1))
  ape::plot.phylo(tree, show.tip.label = FALSE, no.margin = TRUE, label.offset = 0.1)
  
  # Plot heatmap
  par(fig = c(0.47, 0.96, 0.035, 0.965), new = TRUE)
  image(x = c(1:ntip + 1), y = c(1:ntip), z = orderedmat,
        axes = FALSE, breaks = breaks, col = cols)
  
  # Plot legend (right strip)
  par(fig = c(0.968, 0.99, 0.035, 0.965), new = TRUE)
  leg <- matrix(nrow = 1, ncol = 20, data = seq(from = 0.025, to = 1, by = 0.05))
  image(y = seq(from = 0, to = 1, by = 0.05), x = 1, z = leg,
        col = cols, breaks = breaks, axes = FALSE)
  abline(h = c(0.95, 0.9, 0.5, 0.05), lwd = 0.5)
  text(x = 1, y = 0.975, labels = "p=0.95-1", cex = legend_cex, srt = 90)
  text(x = 1, y = 0.925, labels = "p=0.9-0.95", cex = legend_cex, srt = 90)
  text(x = 1, y = 0.7, labels = "p=0.5-0.9", cex = legend_cex, srt = 90)
  text(x = 1, y = 0.275, labels = "p=0.05-0.5", cex = legend_cex, srt = 90)
  text(x = 1, y = 0.025, labels = "p=0-0.05", cex = legend_cex, srt = 90)
}
