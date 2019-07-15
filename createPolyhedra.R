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

# Draw a polygon. Offset is optional.
drawPoly <- function(p, x=0, y=0, z=0, label="", debug=F)
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
    rgl.texts(x, y, z + min(p$vertices$z) - 1, text = label, color = "black")
  }
  if (length(p$faces) > 0) {
    bodies <- findDistinctBodies(p)
    if (length(bodies) > 1) {
      bodyColors <- rainbow(length(bodies))
    }
    for (f in seq(length(p$faces))) {
      if (length(bodies) > 1) {
        faceColor <- bodyColors[which(sapply(bodies, function(b) { return(f %in% b)}))]
      } else {
        faceColor <- rainbow(length(p$faces))[f]
      }
      cx = mean(p$vertices$x[p$faces[[f]]])
      cy = mean(p$vertices$y[p$faces[[f]]])
      cz = mean(p$vertices$z[p$faces[[f]]])
      if (debug) {
        text3d(x + cx, y + cy, z + cz, text = paste0("F",f), color="black")
      }
      if (length(p$faces[[f]]) > 3) { 
        # for > 3 vertices triangulize: draw triangles between center of face and all edges
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

tetrahedron <- buildRegularPoly(vertices = rbind(data.frame(x=1, y=1, z=1), data.frame(x=1, y=-1, z=-1), data.frame(x=-1, y=1, z=-1), data.frame(x=-1, y=-1, z=1)),
                                polygonsize = 3,
                                vertexsize = 3,
                                debug = T)
drawPoly(tetrahedron, label="Tetrahedron")

drawPoly(dual(tetrahedron), z = -3, label="Dual Tetrahedron")

drawPoly(compose(tetrahedron, dual(tetrahedron)), z = -6)

# cubedirect <- buildRegularPoly(vertices = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
#                                polygonsize = 4,
#                                vertexsize = 3,
#                                debug=T)

octahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = c(-1,1), y = 0, z = 0), expand.grid(x = 0, y = c(-1,1), z = 0), expand.grid(x = 0, y = 0, z = c(-1,1))),
                               polygonsize = 3,
                               vertexsize = 4,
                               exampleEdge = c(1,3))
drawPoly(octahedron, x = 3, label="Octahedron")

cube <- dual(octahedron)
drawPoly(cube, x = 3, z = -3, label="Cube")

cubeWithDual <- compose(cube, octahedron)
drawPoly(cubeWithDual, x = 3, z = -6)

# all coords taken from https://en.wikipedia.org/wiki/Platonic_solid
phi <- (1+sqrt(5))/2

# dodecahedron seems to work even lacking same-plane condition
# dodecahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = c(-1,1), y = c(-1,1), z = c(-1,1)), 
#                                                   expand.grid(x = 0, y = c(-1/phi,1/phi), z = c(-phi, phi)), 
#                                                   expand.grid(x = c(-1/phi,1/phi), y = c(-phi, phi), z = 0),
#                                                   expand.grid(x = c(-phi, phi), y = 0, z = c(-1/phi,1/phi))),
#                                  polygonsize = 5,
#                                  vertexsize = 3,
#                                  exampleEdge = c(1,9),
#                                  debug = T)

icosahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = 0, y = c(-1,1), z = c(-phi, phi)), 
                                                 expand.grid(x = c(-1,1), y = c(-phi, phi), z = 0), 
                                                 expand.grid(x = c(-phi, phi), y = 0, z = c(-1, 1))),
                                polygonsize = 3,
                                vertexsize = 5,
                                debug = F)


drawPoly(icosahedron, x = 6, label="Icosahedron")

dodecahedron <- dual(icosahedron)
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


compound5tetrahedra <- buildRegularPoly(dodecahedron$vertices,
                                        polygonsize = 3,
                                        vertexsize = 3,
                                        exampleEdge = c(3, 8))

drawPoly(compound5tetrahedra, x = 0, y = -3, label="Compound of 5 tetrahedra")
drawPoly(dual(compound5tetrahedra), x = 0, y = -3, z = -3, label="Dual of 5 tetrahedra")

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
drawPoly(compound10Cubes, x = 3, y = -3, label="Compound of 10 cubes")

# rgl.close()


# colouring
# if there multiple bodies -> each its own color
# else if there are multiple types of faces -> each its own color
# else rainbow






