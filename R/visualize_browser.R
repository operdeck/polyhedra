#!/usr/bin/env Rscript
# Browser-based rgl visualization without XQuartz.
# Usage: Rscript visualize_browser.R [polyname]

args <- commandArgs(trailingOnly = TRUE)
polyName <- if (length(args) >= 1) args[[1]] else "icosahedron"

options(rgl.useNULL = TRUE)
suppressPackageStartupMessages({
  library(rgl)
  library(htmlwidgets)
})

source("polyhedra.R")
source("draw.R")

if (!exists(polyName, mode = "any")) {
  stop(sprintf("Unknown polyhedron or object '%s'.\nAvailable names: %s",
               polyName, paste(sort(ls()), collapse = ", ")), call. = FALSE)
}

poly <- get(polyName)

if (!is.list(poly) || is.null(poly$coords) || is.null(poly$faces)) {
  cat("Note: the selected object is not a single polyhedron. Trying to draw it as a list of objects.\n")
}

open3d(useNULL = TRUE)
widget <- drawPoly(poly, label = polyName)

outfile <- file.path(getwd(), paste0("rgl_", gsub("[[:space:]]+", "_", polyName), ".html"))
saveWidget(widget, outfile, selfcontained = TRUE)
cat("Saved visualization to:", outfile, "\n")
browseURL(outfile)
