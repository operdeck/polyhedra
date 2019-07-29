library(rgl)

source("math.R")
source("polyhedra.R")

#open3d()
#shade3d( icosahedron3d() )
#writeOBJ(filename)

#rgl.spheres(vertices$x, vertices$y, vertices$z, r=0.2, color="yellow")
#rgl.texts(vertices$x, vertices$y, vertices$z, text = seq(nrow(vertices)), color="red")

# 3D tricks (movie, HTML interactive export, points highlighting)
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization#prepare-the-data

# Draw a single polygon. Offset is optional.
drawSinglePoly <- function(p, x=0, y=0, z=0, label="", debug=F)
{
  if (debug) {
    spacing <- 0.1
    alpha <- 0.6 # in debug make somewhat transparent
    spheres3d(x + p$vertices$x, y + p$vertices$y, z + p$vertices$z, color="grey", radius = 0.01)
    text3d(x + (1+spacing)*p$vertices$x, y + (1+spacing)*p$vertices$y, z + (1+spacing)*p$vertices$z, text = seq(nrow(p$vertices)), color="blue")
  } else {
    alpha <- 1
  }
  label <- paste(label, "(", description(p), ")")
  if (nchar(label) > 0) {
    text3d(x, y + min(p$vertices$y) - 1, z, text = label, color = "black", cex=0.7, pos = 1)
  }
  if (length(p$faces) > 0) { 
    if (!debug) {
      bodies <- findDistinctBodies(p) # this call is too heavy when in debug mode and may not work
      if (length(bodies) > 1) {
        bodyColors <- rainbow(length(bodies))
      }
    } else {
      bodies <- list()
    }
    faceType <- as.integer(factor(sapply(p$faces, length))) # faces considered same just by nr of edges
    if (max(faceType) > 1) {
      faceTypeColors <- rainbow(max(faceType))  
    }
    for (f in seq(length(p$faces))) {
      if (length(bodies) > 1) {
        faceColor <- bodyColors[which(sapply(bodies, function(b) { return(f %in% b)}))]
      } else {
        if (max(faceType) > 1) {
          faceColor <- faceTypeColors[faceType[f]]
        } else {
          faceColor <- rainbow(length(p$faces))[f]
        }
      }
      cx = mean(p$vertices$x[p$faces[[f]]])
      cy = mean(p$vertices$y[p$faces[[f]]])
      cz = mean(p$vertices$z[p$faces[[f]]])
      if (debug) {
        text3d(x + cx, y + cy, z + cz, text = paste0("F",f), color="black")
      }
      if (length(p$faces[[f]]) > 3) { 
        # for > 3 vertices triangulize: draw triangles between center of face and all edges
        
        # TODO this does not work well for {5/2} - just try plot in isolation
        
        rotatedFace <- shiftrotate(p$faces[[f]])
        for (t in seq(length(p$faces[[f]]))) {
          p1 <- p$faces[[f]][t]
          p2 <- rotatedFace[t]
          # NB not sure about the orientation of the triangle - may have to check on this
          triangles3d( x + c(p$vertices$x[c(p1,p2)],cx), 
                       y + c(p$vertices$y[c(p1,p2)],cy), 
                       z + c(p$vertices$z[c(p1,p2)],cz), 
                       col=faceColor, alpha=alpha) # "red"
        }
      } else {
        # NB not sure about the orientation of the triangle - may have to check on this
        triangles3d( x + p$vertices$x[p$faces[[f]]],
                     y + p$vertices$y[p$faces[[f]]], 
                     z + p$vertices$z[p$faces[[f]]], 
                     col=faceColor, alpha=alpha)
      }
    }
  }
}

drawPoly <- function(p, start = c(0, 0, 0), delta = c(1, 0, 0), label = "", debug=F)
{
  if (!is.null(names(p))) {
    drawSinglePoly(p, start[1], start[2], start[3], p$name, debug)
  } else {
    for (i in seq(length(p))) {
      drawSinglePoly(p[[i]], start[1] + (i-1)*delta[1], start[2] + (i-1)*delta[2], start[3] + (i-1)*delta[3], p[[i]]$name, debug)  
    }
    rgl.texts(start[1], start[2] + 2, start[3], text = label, color="blue", pos = 4, cex = 1)
  }
}
# rgl.close()

# The function rgl_init() will create a new RGL device if requested or if there is no opened device:
  
#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 640) {
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, width ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

## Gallery

# open3d()
# rgl.clear()
# par3d("cex" = 0.7)
rgl_init()

tetrahedron <- buildRegularPoly(vertices = rbind(data.frame(x=1, y=1, z=1), data.frame(x=1, y=-1, z=-1), data.frame(x=-1, y=1, z=-1), data.frame(x=-1, y=-1, z=1)),
                                polygonsize = 3,
                                vertexsize = 3,
                                name = "Tetrahedron")

octahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = c(-1,1), y = 0, z = 0), expand.grid(x = 0, y = c(-1,1), z = 0), expand.grid(x = 0, y = 0, z = c(-1,1))),
                               polygonsize = 3,
                               vertexsize = 4,
                               exampleEdge = c(1,3),
                               name = "Octahedron")

cube <- dual(octahedron, name = "Cube")

# all coords taken from https://en.wikipedia.org/wiki/Platonic_solid
phi <- (1+sqrt(5))/2
icosahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = 0, y = c(-1,1), z = c(-phi, phi)), 
                                                 expand.grid(x = c(-1,1), y = c(-phi, phi), z = 0), 
                                                 expand.grid(x = c(-phi, phi), y = 0, z = c(-1, 1))),
                                polygonsize = 3,
                                vertexsize = 5,
                                name = "Icosahedron",
                                debug = F)
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

Platonics <- list(tetrahedron, octahedron, cube, icosahedron, dodecahedron)

greatDodecahedron <- buildRegularPoly(vertices = icosahedron$vertices, 
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6),
                                      name = "Great Dodecahedron")
smallStellatedDodecahedron <- buildRegularPoly(icosahedron$vertices,
                                               polygonsize = 5,
                                               vertexsize = 5,
                                               exampleEdge = c(1,7),
                                               name = "Small Stellated Dodecahedron")
greatIcosahedron <- buildRegularPoly(icosahedron$vertices,
                                     polygonsize = 3,
                                     vertexsize = 5,
                                     exampleEdge = c(2, 6),
                                     name = "Great Icosahedron")
greatStellatedDodecahedron <- dual(greatIcosahedron, name = "Great Stellated Dodecahedron", scaling = "vertex")

KeplerPoinsots <- list(greatDodecahedron, smallStellatedDodecahedron, greatIcosahedron, greatStellatedDodecahedron)

Regulars <- c(Platonics, KeplerPoinsots)
stop()
drawPoly(Regulars, delta = c(3, 0, 0), label = "Regular Polyhedra") # draw along the x-axis

Duals <- lapply(Regulars[seq(length(Regulars))%%2==1], dual)
drawPoly(Duals, start = c(0, 0, -3), delta = c(6, 0, 0), label = "Duals of Regulars")

combis <- lapply(Regulars[seq(length(Regulars))%%2==1], function(p) { return(compose(p, dual(p))) })
drawPoly(combis, start = c(0, 0, -6), delta = c(6, 0, 0), label = "Combined")

archis <- lapply(Regulars[seq(length(Regulars))%%2==1], archi)
drawPoly(archis, start = c(0, 0, -9), delta = c(6, 0, 0), label = "Rhombic")

# cubedirect <- buildRegularPoly(vertices = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
#                                polygonsize = 4,
#                                vertexsize = 3,
#                                debug=T)

# this family (4) should be shown

compound5tetrahedra <- buildRegularPoly(dodecahedron$vertices,
                                        polygonsize = 3,
                                        vertexsize = 3,
                                        exampleEdge = c(3, 8),
                                        name = "5 Tetrahedra")

dodecahedronStellations <- list(compound5tetrahedra,
                                dual(compound5tetrahedra, name = "Dual 5 Tetrahedra"),
                                archi(compound5tetrahedra, name = "Rhombic of 5 Tetrahedra"),
                                compose(compound5tetrahedra, dual(compound5tetrahedra), name = "10 Tetrahedra"))

drawPoly(dodecahedronStellations, start = c(3*4, 6, 0), delta = c(0, 0, -5), label = "Stellations Dodecahedron")

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
# compound10tetrahedra <- buildRegularPoly(dodecahedron$vertices,
#                                         polygonsize = 3,
#                                         vertexsize = 6,
#                                         exampleEdge = c(3, 8), debug=T)

# Dodecahedron stellation family

compound5Cubes <- buildRegularPoly(dodecahedron$vertices,
                                    polygonsize = 4,
                                    vertexsize = 6,
                                    exampleEdge = c(1, 8),
                                   name = "5 Cubes")
moreDodecahedronStellations <- list(compound5Cubes, dual(compound5Cubes), archi(compound5Cubes))
drawPoly(moreDodecahedronStellations, start = c(3*5, 6, 0), delta = c(0, 0, -5))
# this one we saw before: dual(compound5Cubes)


# looks nice but faces overlap?
#drawPoly(archi(compound10Cubes), x = 3, y = -3, label="many archis") ### 5??

# saving to files:
# snapshot3d( "gallery.png", fmt = "png", top = TRUE )
# rgl.postscript("gallery.svg", fmt="svg")

# movie:
# movie3d(spin3d(axis = c(1, 1, 1)), duration = 3, dir = getwd())

stop("Brute force below")

# can we brute-force our way through the possibilities?
rgl_init()
p <- dodecahedron #icosahedron
vexPairs <- combn(x = seq(nrow(p$vertices)), m = 2)
vexPairs <- data.frame(t(vexPairs), dist=sapply(seq(ncol(vexPairs)), function(col) { return(distance(p$vertices[vexPairs[1,col],], p$vertices[vexPairs[2,col],]))}))
vexPairs <- vexPairs[order(vexPairs$dist),]
vexPairs$isUnique <- T
row.names(vexPairs) <- NULL
for (i in 2:nrow(vexPairs)) {
  if (deltaEquals(vexPairs$dist[i-1], vexPairs$dist[i])) { vexPairs$isUnique[i] <- F }
}
vexPairs <- vexPairs[(vexPairs$isUnique),] # this is now a short list of vertex index pairs

discovery <- expand.grid(polygonsize = 3:5, vertexsize = 3:20, vexPair = 1:nrow(vexPairs))
discovery$isPolygon <- F
for (i in seq(nrow(discovery))) {
  # NB it also creates polys with holes in it
  #some of the ones from dodecahedron are new!
  poly <- buildRegularPoly(p$vertices, discovery$polygonsize[i], discovery$vertexsize[i], c(vexPairs$X1[discovery$vexPair[i]], vexPairs$X2[discovery$vexPair[i]]))
  if (poly$n_faces > 0) {
    cat("Success!", discovery$polygonsize[i], discovery$vertexsize[i], vexPairs$X1[discovery$vexPair[i]], vexPairs$X2[discovery$vexPair[i]], fill=T)
    discovery$isPolygon[i] <- T
    drawPoly(poly, x=2*sum(discovery$isPolygon), debug=F)
  }
}

# Now that we have the orientation covered, try list the edges in a constructive way so we
# can build snub, chop, and maybe better duals as well
# We could test the misbehaving duals with even just a partial polyhedron perhaps

# figure out the topology

rgl_init
clear3d()
drawPoly(p, debug=T)
# orientation is clockwise when looking at face

# Dual is now very easy

# NB not sure the normalization should stay here - maybe we should apply a single factor to all


# Possibly the snub also


# Try create a description of the polygon
# unique(sapply(vexConnections, function(c) { return(paste(sapply(c$faces, function(f) {return(length(p$faces[[f]]))}), collapse = ",")) }))
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

