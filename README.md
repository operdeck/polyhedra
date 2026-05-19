# Polyhedra

R tools to create various polyhedra from ground principles.

![Compounds](snapshots/compound5octahedra.png)

## Getting started

This repository is script-based, so you can use it interactively from R or run the helper script to create browser-rendered polyhedron output.

### Install dependencies

In R:

```r
install.packages(c(
  "data.table",
  "rgl",
  "testthat",
  "ggplot2",
  "svglite",
  "colorspace",
  "htmlwidgets"
))
```

### Run browser visualization without XQuartz

From the repo root:

```bash
cd /Users/ottoperdeck/dev/polyhedra
Rscript visualize_browser.R icosahedron
```

This generates and opens an HTML file like `rgl_Icosahedron.html`.

### Use interactive R / VS Code

```r
setwd("/Users/ottoperdeck/dev/polyhedra")
source("polyhedra.R")
source("draw.R")

widget <- drawSinglePolyWidget(icosahedron)
widget
```

Or draw multiple objects:

```r
widget <- drawPolyWidget(list(tetrahedron, cube, icosahedron))
widget
```

Only the coordinates of the very basic polyhedra are given (tetrahedron, octahedron, icosahedron). The topology of those is not given but "discovered" by the tool. Other polyhedra are created as derived from these. E.g. the dodecahedron is created as

```r
dodecahedron <- dual(icosahedron)
```

## Browser-based visualization

This project can render `rgl` output as HTML widgets in a browser. That is the recommended mode on modern macOS.

Install the browser rendering dependency if needed:

```r
install.packages(c("data.table","rgl","testthat","ggplot2","svglite","colorspace","htmlwidgets"))
```

Then run:

```bash
cd /Users/ottoperdeck/dev/polyhedra
Rscript visualize_browser.R icosahedron
```

This creates and opens `rgl_Icosahedron.html` in your browser.

You can also use the new browser helpers from an R session:

```r
setwd("/Users/ottoperdeck/dev/polyhedra")
source("polyhedra.R")
source("draw.R")
widget <- drawSinglePolyWidget(icosahedron)
widget
```

Or for a list of shapes:

```r
widget <- drawPolyWidget(list(tetrahedron, cube, icosahedron))
widget
```

The basic stellated polyhedra ("Kepler-Poinsot polyhedra") have the same coordinates as the simpler ones and can be created from these

```r
greatDodecahedron <- buildRegularPoly(coords = icosahedron$coords, 
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6),
                                      name = "Great Dodecahedron")
smallStellatedDodecahedron <- buildRegularPoly(icosahedron$coords,
                                               polygonsize = 5,
                                               vertexsize = 5,
                                               exampleEdge = c(1,7),
                                               name = "Small Stellated Dodecahedron")
greatIcosahedron <- buildRegularPoly(icosahedron$coords,
                                     polygonsize = 3,
                                     vertexsize = 5,
                                     exampleEdge = c(2, 6),
                                     name = "Great Icosahedron")
greatStellatedDodecahedron <- dual(greatIcosahedron, name = "Great Stellated Dodecahedron", scaling = "vertex")
```

![Kepler Poinsots](snapshots/keplerpoinsots.png)

Basic descriptions of the polyhedra (according to Coxeter) are also generated along with the polyhedra.

Next to "dual", several other transformations like "truncate" and "rhombic" are available too, allowing to create a lot of different polyhedra

![Transformations](snapshots/rhombicregulars.png)

There also is a full discovery mode, where you only give a set of coordinates (from one of the other polyhedra) then it "discovers" a family of polyhedra sharing those vertex coordinates.

![Dodecahedron discoveries](snapshots/dodecahedron_discoveries.png)

The above is not a full list. Run the R notebook to create a more complete list.

Next up is a tool to automatically create a 2D layout of the polyhedra.







