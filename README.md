# Polyhedra

Tools to create various polyhedra from ground principles — exploring construction, transformation, and visualization of regular and semi-regular solids.

## Repository structure

```
R/              ← Legacy R implementation (complete, functional)
wolfram/        ← New implementation: Wolfram Language + WLJS + Python (TODO)
```

## Wolfram / WLJS / Python (TODO)

The next-generation implementation will use:

- **Wolfram Language** — symbolic geometry, algebraic constructions, stellation
- **WLJS** — interactive browser-based 3D visualization
- **Python** — orchestration, data pipeline, notebook integration

> This is under development. See `wolfram/` for progress.

---

## R (legacy)

The original R implementation lives in `R/`. It constructs polyhedra from vertex coordinates and topology, supports transformations (dual, truncate, quasi, rhombic, stellation), discovery mode, and 2D layout generation.

### Quick start

```r
install.packages(c("data.table", "rgl", "testthat", "ggplot2",
                   "svglite", "colorspace", "htmlwidgets", "rmarkdown"))
```

```r
setwd("R")
source("polyhedra.R")
source("draw.R")

widget <- drawPolyWidget(list(tetrahedron, cube, icosahedron))
widget
```

### Render the notebook

```r
rmarkdown::render("R/polyhedra.Rmd")
```

### Browser visualization (no XQuartz)

```bash
Rscript R/visualize_browser.R icosahedron
```

### Key concepts

Only the coordinates of the basic polyhedra are given (tetrahedron, octahedron, icosahedron). Topology is discovered by the tool. Other polyhedra are derived via transformations:

```r
dodecahedron <- dual(icosahedron)
greatDodecahedron <- buildRegularPoly(coords = icosahedron$coords,
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6),
                                      name = "Great Dodecahedron")
```

Transformations available: `dual()`, `truncate()`, `quasi()`, `rhombic()`, `compose()`, and discovery mode via `discover()`.







