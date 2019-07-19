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
  if (nchar(label) > 0) {
    text3d(x, y + min(p$vertices$y) - 1, z, text = label, color = "black", cex=0.7, pos = 1)
  }
  if (length(p$faces) > 0 & !debug) { 
    bodies <- findDistinctBodies(p) # this call is too heavy when in debug mode and may not work
    if (length(bodies) > 1) {
      bodyColors <- rainbow(length(bodies))
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
        
        rotatedFace <- shift(p$faces[[f]])
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

drawPoly(Platonics, delta = c(3, 0, 0), label = "Platonic Polyhedra") # draw along the x-axis

Duals <- lapply(Platonics, dual)
drawPoly(Duals, start = c(0, 0, -3), delta = c(3, 0, 0), label = "Duals of Platonics")

KeplerPoinsots <- list() # etc


stop()


drawPoly(tetrahedron, label="Tetrahedron")

drawPoly(dual(tetrahedron), z = -3, label="Dual Tetrahedron")

drawPoly(compose(tetrahedron, dual(tetrahedron)), z = -6)

# cubedirect <- buildRegularPoly(vertices = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
#                                polygonsize = 4,
#                                vertexsize = 3,
#                                debug=T)

drawPoly(octahedron, x = 3, label="Octahedron")

cube <- dual(octahedron)
drawPoly(cube, x = 3, z = -3, label="Cube")

cubeWithDual <- compose(cube, octahedron)
drawPoly(cubeWithDual, x = 3, z = -6)

drawPoly(dodecahedron, x = 6, z = -3, label="Dodecahedron")

drawPoly(compose(icosahedron, dual(icosahedron)), x = 6, z = -6)

greatDodecahedron <- buildRegularPoly(vertices = icosahedron$vertices, 
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6))
drawPoly(greatDodecahedron, x = 9, label="Great Dodecahedron")

# the dual of this - doesnt work! strange - maybe because of {5/2} faces that occur, should be smallStellatedDodecahedron
# instead it creates a total mess
# drawPoly(dual(greatDodecahedron), debug=T) 

smallStellatedDodecahedron <- buildRegularPoly(icosahedron$vertices,
                                     polygonsize = 5,
                                     vertexsize = 5,
                                     exampleEdge = c(1,7))

drawPoly(smallStellatedDodecahedron, x = 9, z= -3, label="Small Stellated Dodecahedron") # its dual doesnt work either

# drawPoly(compose(smallStellatedDodecahedron, greatDodecahedron)) first is occluded completely by the other

greatIcosahedron <- buildRegularPoly(icosahedron$vertices,
                                     polygonsize = 3,
                                     vertexsize = 5,
                                     exampleEdge = c(2, 6))
drawPoly(greatIcosahedron, x = 12, label="Great Icosahedron")

greatStellatedDodecahedron <- dual(greatIcosahedron)
drawPoly(greatStellatedDodecahedron, x = 12, z = -3, label="Great Stellated Dodecahedron")

drawPoly(compose(greatStellatedDodecahedron, greatIcosahedron), x = 12, z = -6)

# Archi (are the same for the duals of regulars)
drawPoly(archi(cube), x = 3, y = 3, z = -1.5, label="{3,4,3,4}")
drawPoly(archi(icosahedron), x = 6, y = 3, z = -1.5, label="{3,5,3,5}")
drawPoly(archi(greatDodecahedron), x = 9, y = 3, z = -1.5, label="{5/2,5,5/2,5}")
drawPoly(archi(greatIcosahedron), x = 12, y = 3, z = -1.5, label="{3,5/2,3,5/2}")

compound5tetrahedra <- buildRegularPoly(dodecahedron$vertices,
                                        polygonsize = 3,
                                        vertexsize = 3,
                                        exampleEdge = c(3, 8))

drawPoly(compound5tetrahedra, x = 0, y = -3, label="Compound of 5 tetrahedra")
drawPoly(dual(compound5tetrahedra), x = 0, y = -3, z = -3, label="Dual of 5 tetrahedra")

drawPoly(archi(compound5tetrahedra), x = 0, y = -6, z = -1.5, label="5 octahedra")

# below does not look entirely OK, there are clashing faces
# see https://en.wikipedia.org/wiki/Polytope_compound
drawPoly(compose(compound5tetrahedra, dual(compound5tetrahedra)), x = 0, y = -3, z = -6, 
         label="Compound of 10 tetrahedra")

# When trying to construct directly, it tries to create two overlapping faces, which is
# not allowed. Lifting that restriction it creates exactly the same flickering polyhedron
# as when composing. Plus lifting that restriction used to fail a few other polyhedra 
# although that no longer seems to be the case. A possible way out could be to combine
# those faces
# compound10tetrahedra <- buildRegularPoly(dodecahedron$vertices,
#                                         polygonsize = 3,
#                                         vertexsize = 6,
#                                         exampleEdge = c(3, 8), debug=T)

# drawPoly(dodecahedron, debug=T)
compound10Cubes <- buildRegularPoly(dodecahedron$vertices,
                                    polygonsize = 4,
                                    vertexsize = 6,
                                    exampleEdge = c(1, 8))
drawPoly(compound10Cubes, x = 3, y = -3, label="Compound of 10 cubes") ### 5??

# looks nice but faces overlap?
#drawPoly(archi(compound10Cubes), x = 3, y = -3, label="many archis") ### 5??

# rgl.close()
# snapshot3d( "gallery.png", fmt = "png", top = TRUE )
# rgl.postscript("gallery.svg", fmt="svg")

# colouring
# if there multiple bodies -> each its own color
# else if there are multiple types of faces -> each its own color
# else rainbow


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
unique(sapply(vexConnections, function(c) { return(paste(sapply(c$faces, function(f) {return(length(p$faces[[f]]))}), collapse = ",")) }))

