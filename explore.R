source("polyhedra.R")
source("draw.R")

# Draw a single polygon. Offset is optional.
## Remaining Platonic solids
cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")
Platonics <- list(tetrahedron, octahedron, cube, icosahedron, dodecahedron)

## Regular star solids
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

KeplerPoinsots <- list(greatDodecahedron, smallStellatedDodecahedron, greatIcosahedron, greatStellatedDodecahedron)

Regulars <- c(Platonics, KeplerPoinsots)

Duals <- lapply(Regulars[seq(length(Regulars))%%2==1], dual)

# For printing
# clear3d()
# drawPoly(lapply(KeplerPoinsots, dual), start = c(0,2,0))
# drawPoly(lapply(Platonics, truncate), delta = c(2, 0, 0), label = "")
# snapshot3d("truncates.png", fmt="png", top=T)

clear3d()
drawAxes()

drawPoly(Regulars, delta = c(3, 0, 0), label = "Regular Polyhedra") # draw along the x-axis
drawPoly(Duals, start = c(0, 0, -3), delta = c(6, 0, 0), label = "Duals of Regulars")

combis <- lapply(Regulars[seq(length(Regulars))%%2==1], function(p) { return(compose(p, dual(p))) })
drawPoly(combis, start = c(0, 0, -6), delta = c(6, 0, 0), label = "Combined")

quasis <- lapply(Regulars, quasi)
drawPoly(quasis, start = c(0, 0, -9), delta = c(6, 0, 0), label = "Quasi")

# cubedirect <- buildRegularPoly(coords = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
#                                polygonsize = 4,
#                                vertexsize = 3,
#                                debug=T)

# this family (4) should be shown

compound5tetrahedra <- buildRegularPoly(dodecahedron$coords,
                                        polygonsize = 3,
                                        vertexsize = 3,
                                        exampleEdge = c(3, 8),
                                        name = "5 Tetrahedra")

# TODO!
# below 10 x {3,3} is still flickering but maybe this one has truly overlapping faces
# so two triangles really form a {6/2} together
dodecahedronStellations <- list(compound5tetrahedra,
                                dual(compound5tetrahedra, name = "Dual 5 Tetrahedra"),
                                quasi(compound5tetrahedra, name = "Quasi of 5 Tetrahedra"),
                                compose(compound5tetrahedra, dual(compound5tetrahedra), name = "10 Tetrahedra"))

drawPoly(dodecahedronStellations, start = c(3*4, 6, 0), delta = c(0, 0, -5), label = "Stellations Dodecahedron")

clear3d()
drawSinglePoly(dual(quasi(octahedron)))

# below does not look entirely OK, there are clashing faces but
# this could be due to the way {5/2} etc are trianglulized currently
# see https://en.wikipedia.org/wiki/Polytope_compound
# drawSinglePoly(compose(compound5tetrahedra, dual(compound5tetrahedra)), x = 0, y = -3, z = -6, 
#          label="Compound of 10 tetrahedra")

# When trying to construct directly, it tries to create two overlapping faces, which is
# not allowed. Lifting that restriction it creates exactly the same flickering polyhedron
# as when composing. Plus lifting that restriction used to fail a few other polyhedra 
# although that no longer seems to be the case. A possible way out could be to combine
# those faces
# compound10tetrahedra <- buildRegularPoly(dodecahedron$coords,
#                                         polygonsize = 3,
#                                         vertexsize = 6,
#                                         exampleEdge = c(3, 8), debug=T)

# Dodecahedron stellation family

compound5Cubes <- buildRegularPoly(dodecahedron$coords,
                                    polygonsize = 4,
                                    vertexsize = 6,
                                    exampleEdge = c(1, 8),
                                   name = "5 Cubes")
moreDodecahedronStellations <- list(compound5Cubes, dual(compound5Cubes), quasi(compound5Cubes))
drawPoly(moreDodecahedronStellations, start = c(3*5, 6, 0), delta = c(0, 0, -5))
# this one we saw before: dual(compound5Cubes)

clear3d()
# looks nice but faces overlap?
#drawPoly(quasi(compound10Cubes), x = 3, y = -3, label="many quasis") ### 5??

# saving to files:
# snapshot3d( "gallery.png", fmt = "png", top = TRUE )
# rgl.postscript("gallery.svg", fmt="svg")

# movie:
# movie3d(spin3d(axis = c(1, 1, 1)), duration = 3, dir = getwd())


rgl_init()
clear3d()
drawPoly(p, debug=T)
# orientation is clockwise when looking at face

# Dual is now very easy

# NB not sure the normalization should stay here - maybe we should apply a single factor to all


# Possibly the truncate also


# Try create a description of the polygon
# unique(sapply(vertexFigures, function(c) { return(paste(sapply(c$faces, function(f) {return(length(p$faces[[f]]))}), collapse = ",")) }))
# 
# Systematic list of all Archimedean solids and the operations to support:
#   
#   https://en.wikipedia.org/wiki/Archimedean_solid
# 
# Catalan solids, have a think
# 
# https://www.software3d.com/Archimedean.php
# 
# List of compounds I should be able to generate (currently not all e.g. 4 cubes):
#   
#   https://www.polyhedra.net/en/pictures.php?type=c
# 

# rgl.close()

