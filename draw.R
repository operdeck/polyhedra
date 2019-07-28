library(rgl)
library(data.table)
source("math.R")

#open3d()
#shade3d( icosahedron3d() )
#writeOBJ(filename)

#rgl.spheres(vertices$x, vertices$y, vertices$z, r=0.2, color="yellow")
#rgl.texts(vertices$x, vertices$y, vertices$z, text = seq(nrow(vertices)), color="red")

# 3D tricks (movie, HTML interactive export, points highlighting)
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization#prepare-the-data

drawPolygonTriangulate <- function(vertexIndices, vertexCoords, col="grey", alpha=1, offset=c(0,0,0))
{
  center <- apply(vertexCoords[vertexIndices], 2, mean)
  rotatedFace <- shift(vertexIndices)
  for (t in seq(length(vertexIndices))) {
    p1 <- vertexIndices[t]
    p2 <- rotatedFace[t]
    # NB not sure about the orientation of the triangle - may have to check on this
    triangles3d( offset[1] + c(vertexCoords$x[c(p1,p2)], center[1]), 
                 offset[2] + c(vertexCoords$y[c(p1,p2)], center[2]), 
                 offset[3] + c(vertexCoords$z[c(p1,p2)], center[3]), 
                 col=col, alpha=alpha) # "red"
  }
  
}

drawPolygon <- function(vertexIndices, vertexCoords, col="grey", alpha=1, offset=c(0,0,0))
{
  if (length(vertexIndices) > 3) { 
    drawPolygonTriangulate(vertexIndices, vertexCoords, col, alpha, offset)
  } else {
    # NB not sure about the orientation of the triangle - may have to check on this
    triangles3d( offset[1] + vertexCoords$x[vertexIndices],
                 offset[2] + vertexCoords$y[vertexIndices], 
                 offset[3] + vertexCoords$z[vertexIndices], 
                 col=col, alpha=alpha)
  }
  
}

# Draw a single polygon. Offset is optional.
drawSinglePoly <- function(p, x=0, y=0, z=0, label=ifelse(is.null(p$name),"",p$name), debug=F)
{
  if (debug) {
    spacing <- 0.1
    alpha <- 0.6 # in debug make somewhat transparent
    # origin
    spheres3d(x, y, z, color="red", radius = 0.02)
    # all vertices
    spheres3d(x + p$vertices$x, y + p$vertices$y, z + p$vertices$z, color="grey", radius = 0.01)
    text3d(x + (1+spacing)*p$vertices$x, y + (1+spacing)*p$vertices$y, z + (1+spacing)*p$vertices$z, text = seq(nrow(p$vertices)), color="blue")
  } else {
    alpha <- 1
  }
  if (!debug) {
    # avoid heavy description call in debug mode
    label <- paste(label, "(", description(p), ")")
  }
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
      if (debug) {
        cx = mean(p$vertices$x[p$faces[[f]]])
        cy = mean(p$vertices$y[p$faces[[f]]])
        cz = mean(p$vertices$z[p$faces[[f]]])
        text3d(x + cx, y + cy, z + cz, text = paste0("F",f), color="black")
      }
      drawPolygon(p$faces[[f]], p$vertices, faceColor, alpha, c(x, y, z)) # pass debug??
    }
  }
}

drawPoly <- function(p, start = c(0, 0, 0), delta = c(2, 0, 0), label = "", debug=F)
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

# draw a n/d polygon
testDrawPolygon <- function(n, d=1)
{
  clear3d()
  angles <- ((0:(n-1))/n)*2*pi
  coords <- data.table(x=sin(angles),y=cos(angles),z=0)
  face <- ((((1:n)-1)*d)%%n)+1
  drawPolygon(face, coords, col="red")
  # simplified <- list(vertices = coords, faces = list(face))
  # drawSinglePoly(simplified, debug=T)
}
