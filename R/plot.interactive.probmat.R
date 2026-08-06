#' @export
plot.interactive.probmat <- function(x, tree,
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
                                     ...) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' is required. Please install: install.packages('plotly')")
  }
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Please install: install.packages('base64enc')")
  }
  
  # === 1. Data preparation ===
  probmat <- x
  ntip <- length(tree$tip.label)
  tree <- ape::reorder.phylo(tree, order = "cladewise")
  edge <- tree$edge[which(tree$edge[, 2] < (ntip + 1)), ]
  tiplab <- tree$tip.label[edge[, 2]]
  orderedmat <- probmat[tiplab, ]
  orderedmat <- orderedmat[, tiplab]
  
  # === 2. Color palette ===
  if (palette == "green") {
    colors <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B",
                "#74C476", "#41AB5D", "#238B45", "#005A32")
  } else if (palette == "viridis") {
    colors <- c("#440154", "#31688E", "#35B779", "#FDE725")
  } else if (palette == "RdYlBu") {
    colors <- c("#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B")
  } else {
    colors <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A", "#E31A1C")
  }
  
  # === 3. Heatmap data ===
  df <- expand.grid(row = rownames(orderedmat), col = colnames(orderedmat))
  df$prob <- as.vector(orderedmat)
  df$hover_text <- paste0(
    "<b>", df$row, "</b> vs <b>", df$col, "</b><br>",
    "Probability: ", round(df$prob, 3),
    "<br>Interpretation: ",
    ifelse(df$prob >= 0.95, "Same species (p>=0.95)",
           ifelse(df$prob >= 0.90, "Likely same (p>=0.90)",
                  ifelse(df$prob >= 0.50, "Uncertain (0.5-0.9)",
                         ifelse(df$prob >= 0.05, "Likely different (p<0.5)",
                                "Different species (p<0.05)"))))
  )
  
  # === 4. Build heatmap ===
  p <- plotly::plot_ly(
    data = df, x = ~col, y = ~row, z = ~prob, type = "heatmap",
    colors = colors, text = ~hover_text, hoverinfo = "text", showscale = TRUE,
    width = width, height = height,
    hoverlabel = list(bgcolor = "white", font = list(color = "black"))
  )
  
  # === 5. Add tree on the left ===
  if (show_tree) {
    tmp_file <- tempfile(fileext = ".png")
    png(tmp_file, width = 400, height = height, bg = "transparent")
    par(mar = c(0, 0, 0, 0), mgp = c(0, 0, 0), oma = c(0, 0, 0, 0))
    ape::plot.phylo(tree, show.tip.label = FALSE, no.margin = TRUE,
                    label.offset = 0, direction = "rightwards", ...)
    dev.off()
    img_data <- base64enc::dataURI(file = tmp_file, mime = "image/png")
    unlink(tmp_file)
    p <- plotly::layout(p,
                        images = list(list(
                          source = img_data,
                          xref = "paper", yref = "paper",
                          x = 0, y = 1,
                          sizex = tree_width, sizey = 1,
                          xanchor = "left", yanchor = "top",
                          sizing = "stretch",
                          layer = "below"
                        )),
                        xaxis = list(
                          title = "", showticklabels = FALSE, zeroline = FALSE,
                          domain = c(tree_width, 1)
                        ),
                        yaxis = list(
                          title = "", showticklabels = FALSE, zeroline = FALSE
                        ),
                        margin = list(l = 0, r = 20, t = 40, b = 10),
                        plot_bgcolor = "rgba(0,0,0,0)",
                        paper_bgcolor = "rgba(0,0,0,0)",
                        title = list(text = "Conspecificity Probability Matrix", font = list(size = 16)),
                        showlegend = FALSE
    )
  } else {
    p <- plotly::layout(p,
                        xaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
                        yaxis = list(title = "", showticklabels = FALSE, zeroline = FALSE),
                        margin = list(l = 10, r = 20, t = 40, b = 10),
                        title = list(text = "Conspecificity Probability Matrix", font = list(size = 16)),
                        showlegend = FALSE
    )
  }
  
  # === 6. Save to HTML ===
  if (!is.null(save_html)) {
    dir.create(dirname(save_html), showWarnings = FALSE, recursive = TRUE)
    htmlwidgets::saveWidget(p, file = save_html, selfcontained = TRUE)
    message("Interactive plot saved to: ", normalizePath(save_html))
  }

  if (interactive() && is.null(save_html)) { 
    print(p)
  }
  
  invisible(p)
}
