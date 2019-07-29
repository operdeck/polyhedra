library(rgl)
library(data.table)
source("math.R")

# The function rgl_init() will create a new RGL device if requested or if there is no opened device:

#' @param new.device a logical value. If TRUE, creates a new device
#' @param bg the background color of the device
#' @param width the width of the device
rgl_init <- function(new.device = FALSE, bg = "white", width = 640, height = width) {
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, height ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}


#open3d()
#shade3d( icosahedron3d() )
#writeOBJ(filename)

#rgl.spheres(vertices$x, vertices$y, vertices$z, r=0.2, color="yellow")
#rgl.texts(vertices$x, vertices$y, vertices$z, text = seq(nrow(vertices)), color="red")

# 3D tricks (movie, HTML interactive export, points highlighting)
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization#prepare-the-data

# TODO: obsolete
# drawPolygonTriangulate <- function(vertexIndices, vertexCoords, col="grey", alpha=1, offset=c(0,0,0))
# {
#   center <- colMeans(vertexCoords[vertexIndices])
#   for (t in seq(length(vertexIndices))) {
#     p1 <- vertexIndices[t]
#     p2 <- shiftrotate(vertexIndices)[t]
#     # NB not sure about the orientation of the triangle - may have to check on this
#     triangles3d( offset[1] + c(vertexCoords$x[c(p1,p2)], center[1]), 
#                  offset[2] + c(vertexCoords$y[c(p1,p2)], center[2]), 
#                  offset[3] + c(vertexCoords$z[c(p1,p2)], center[3]), 
#                  col=col, alpha=alpha)
#   }
#   
# }

drawStarPolygon <- function(vertexIndices, vertexCoords, col="grey", alpha=1, offset=c(0,0,0))
{
  # Find all line segment pairs of this face, excluding line segments that are already adjecent to eachother.
  segmentPairs <- as.data.table(t(combn(seq(length(vertexIndices)),2)))
  segmentPairs <- segmentPairs[abs(V1-V2)>1 & abs(V1-(V2%%length(vertexIndices)))>1]
  
  # A-B and C-D are the vertex indices of the segment pairs that will be checked for intersection
  A <- vertexIndices[unlist(segmentPairs[,1])]
  B <- shiftrotate(vertexIndices)[unlist(segmentPairs[,1])]
  C <- vertexIndices[unlist(segmentPairs[,2])]
  D <- shiftrotate(vertexIndices)[unlist(segmentPairs[,2])]
  
  # The segment intersections will be calculated from the coordinates of the 4 points
  intersections <- data.table( 
    seg1From = A, seg1To = B, seg2From = C, seg2To = D, # keep the vertex indices
    Ax = vertexCoords$x[A], Ay = vertexCoords$y[A], Az = vertexCoords$z[A],
    Bx = vertexCoords$x[B], By = vertexCoords$y[B], Bz = vertexCoords$z[B],
    Cx = vertexCoords$x[C], Cy = vertexCoords$y[C], Cz = vertexCoords$z[C],
    Dx = vertexCoords$x[D], Dy = vertexCoords$y[D], Dz = vertexCoords$z[D])
  
  # The math to calculate the intersection is easily derived. The intersections are at A + alpha BA = C + beta DC. We want
  # both alpha and beta and both in a numerical stable way because the very compact list of segment pairs we started with,
  # so we may need to swap AB and CD and thus alpha and beta etc.
  intersections[, beta  := ((Cx-Ax)*(By-Ay) - (Cy-Ay)*(Bx-Ax)) / ((Dy-Cy)*(Bx-Ax) - (Dx-Cx)*(By-Ay))]
  intersections[, alpha := ((Ax-Cx)*(Dy-Cy) - (Ay-Cy)*(Dx-Cx)) / ((By-Ay)*(Dx-Cx) - (Bx-Ax)*(Dy-Cy))]
  
  # Intersections are only the ones for which alpha and beta are in 0-1 range. The intersections are not in a particular
  # order although we would want that. Ordering them by inner angle in 3D turns out to be mathematically impossible although
  # potentially we could use the direction of the inner vectors wrt the direction of the inner vector of the whole face.
  intersections <- intersections[alpha<1 & alpha>0 & beta<1 & beta>0]
  intersections[, c("Ix","Iy","Iz") := list(Ax + alpha*(Bx-Ax), Ay + alpha*(By-Ay), Az + alpha*(Bz-Az))]
  intersections[, intersectionIdx := seq(.N)]
  
  # for debugging only:
  spheres3d(intersections$Ix, intersections$Iy, intersections$Iz, color="green", radius = 0.02)
  text3d(intersections$Ix+0.05, intersections$Iy, intersections$Iz, color="black", text=intersections$intersectionIdx)
  
  # Intersections is a very compact representation, expand it to get a list from all face vertices to all intersection points.
  allIntersections <- rbind(data.table(from = intersections$seg1From,
                                       to = intersections$seg1To,
                                       alpha = intersections$alpha,
                                       Ix = intersections$Ix, Iy = intersections$Iy, Iz = intersections$Iz, intersectionIdx = intersections$intersectionIdx),
                            data.table(from = intersections$seg1To,
                                       to = intersections$seg1From,
                                       alpha = (1-intersections$alpha),
                                       Ix = intersections$Ix, Iy = intersections$Iy, Iz = intersections$Iz, intersectionIdx = intersections$intersectionIdx),
                            data.table(from = intersections$seg2From,
                                       to = intersections$seg2To,
                                       alpha = intersections$beta,
                                       Ix = intersections$Ix, Iy = intersections$Iy, Iz = intersections$Iz, intersectionIdx = intersections$intersectionIdx),
                            data.table(from = intersections$seg2To,
                                       to = intersections$seg2From,
                                       alpha = (1-intersections$beta),
                                       Ix = intersections$Ix, Iy = intersections$Iy, Iz = intersections$Iz, intersectionIdx = intersections$intersectionIdx))
  
  # Then, for each of the line segments of the face, cut it at the nearest intersection point. 
  allIntersections <- allIntersections[, isClosestIntersection := (seq(.N) == which.min(alpha)), by=c("from", "to")]
  allIntersections <- allIntersections[(isClosestIntersection)][order(from, to)]
  
  # Now, we can draw the stars by drawing triangles from each of the face vertices to the intersection points connected to it. 
  # Then draw by providing triplets with the vertices of the triangles.
  triangles <- dcast(allIntersections[, triangleside:=seq(.N), by=from], from~triangleside, value.var=c("Ix","Iy","Iz"))
  triangles3d(offset[1] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Ix_1[t], triangles$Ix_2[t], vertexCoords$x[triangles$from[t]])})),
              offset[2] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Iy_1[t], triangles$Iy_2[t], vertexCoords$y[triangles$from[t]])})),
              offset[3] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Iz_1[t], triangles$Iz_2[t], vertexCoords$z[triangles$from[t]])})),
              col=col, alpha=alpha)

  # To draw the middle piece we need the intersections in order. Sorting by inner angle doesnt work in 3D so we do it
  # simply by going round and each time finding the closest intersection that is not used yet. Might be somewhat slow.
  middlepolygon <- unique(allIntersections[, c("intersectionIdx", "Ix", "Iy", "Iz")])
  middlepolygon[, rank := 0]
  currentPt <- 3 # start somewhere randomly
  nAssignedOrder <- 0
  repeat {
    middlepolygon[(currentPt), rank := 1+nAssignedOrder]
    nAssignedOrder <- nAssignedOrder+1
    if (nAssignedOrder == nrow(middlepolygon)) {
      break
    } else {
      distances <- sqrt(rowSums((t(t(middlepolygon[(rank==0),2:4]) - as.numeric(middlepolygon[currentPt,2:4])))^2))  
      nextPt <- which(middlepolygon$rank==0)[which.min(distances)[1]]
    }
    currentPt <- nextPt
  }
  polygon3d( offset[1] + middlepolygon$Ix[order(middlepolygon$rank)],
             offset[2] + middlepolygon$Iy[order(middlepolygon$rank)], 
             offset[3] + middlepolygon$Iz[order(middlepolygon$rank)], 
             col="yellow", alpha=alpha)
}

drawAxes <- function()
{
  arrow3d(p0 = c(0,0,0), p1=c(1.8,0,0), width=0.2, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,1.8,0), width=0.2, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,0,1.8), width=0.2, type="rotation")
  texts3d(c(2,0,0), c(0,2,0), c(0,0,2), text = c("x","y","z"), color="black")
}

drawPolygon <- function(vertexIndices, vertexCoords, col="grey", alpha=1, offset=c(0,0,0), label=NULL, drawlines=F)
{
  if (drawlines) {
    lines3d(vertexCoords[c(vertexIndices, vertexIndices[1])]$x, vertexCoords[c(vertexIndices, vertexIndices[1])]$y, vertexCoords[c(vertexIndices, vertexIndices[1])]$z, color="blue")
  }
  if (!is.null(label)) {
    center <- colMeans(vertexCoords[vertexIndices])
    text3d(offset[1] + center[1], offset[2] + center[2], offset[3] + center[3], text = label, color="black")
  }
  
  # NB face orientation is not normalized
  if (length(vertexIndices) == 3) {
    triangles3d( offset[1] + vertexCoords$x[vertexIndices],
                 offset[2] + vertexCoords$y[vertexIndices],
                 offset[3] + vertexCoords$z[vertexIndices],
                 col=col, alpha=alpha)
  } else if (length(vertexIndices) > 3) {
    ang <- innerAngles(vertexCoords[vertexIndices])
    if(sum(ang) > 2*pi) {
      drawStarPolygon(vertexIndices, vertexCoords, col, alpha, offset)
      # drawPolygonTriangulate(vertexIndices, vertexCoords, col, alpha, offset)
    } else {
      polygon3d( offset[1] + vertexCoords$x[vertexIndices],
                 offset[2] + vertexCoords$y[vertexIndices], 
                 offset[3] + vertexCoords$z[vertexIndices], 
                 col=col, alpha=alpha)
    }
  }
  
}

# Draw a single polygon. Offset is optional.
drawSinglePoly <- function(p, x=0, y=0, z=0, label=ifelse(is.null(p$name),"",p$name), debug=F)
{
  if (debug) {
    drawAxes()
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
      drawPolygon(p$faces[[f]], p$vertices, faceColor, alpha, c(x, y, z), debug, label=paste0("F",f)) 
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



# draw a n/d polygon
testDrawPolygon <- function(n, d=1, label="")
{
  rgl_init()
  clear3d()
  drawAxes()
  angles <- ((0:(n-1))/n)*2*pi
  coords <- data.table(x=sin(angles),y=cos(angles),z=0)
  face <- ((((1:n)-1)*d)%%n)+1 #+ as.numeric(sapply(0:(d-1), function(x){return(rep(x,n%/%d))}))
  
  #innerAngles(coords[face])
  drawPolygon(face, coords, label=paste(n,d,sep="/"), drawlines = T)
  # simplified <- list(vertices = coords, faces = list(face))
  # drawSinglePoly(simplified, debug=T)
}

testDrawPolygons <- function()
{
  testDrawPolygon(3)
  testDrawPolygon(4) # triangulization
  testDrawPolygon(5) 
  #testDrawPolygon(6,2) # interesting case but would create multiple faces really
  testDrawPolygon(5,2) # simple star
  testDrawPolygon(8,3) # more complex star
  #testDrawPolygon(8,2) # also a star, generation of this is slightly more complex
}
