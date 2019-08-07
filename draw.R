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

#rgl.spheres(coords$x, coords$y, coords$z, r=0.2, color="yellow")
#rgl.texts(coords$x, coords$y, coords$z, text = seq(nrow(coords)), color="red")

# 3D tricks (movie, HTML interactive export, points highlighting)
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization#prepare-the-data


# the rgl method is not reliable so still doing it myself
drawPolygonTriangulate <- function(face, coords, col="grey", alpha=1, offset=c(0,0,0))
{
  center <- apply(coords[face,], 2, mean)
  for (t in seq(length(face))) {
    p1 <- face[t]
    p2 <- shiftrotate(face)[t]
    # NB not sure about the orientation of the triangle - may have to check on this
    triangles3d( offset[1] + c(coords[c(p1,p2),1], center[1]),
                 offset[2] + c(coords[c(p1,p2),2], center[2]),
                 offset[3] + c(coords[c(p1,p2),3], center[3]),
                 col=col, alpha=alpha)
  }

}


# Draws a star polygon with intersecting edges

# NB perhaps could be a lot simpler by getting the vertices in actual order, finding
# the intersections of the edges of neighbouring vertices, then defining a new face 
# with 2*N vertices. Perhaps polygon3d can then even draw that directly. Also perhaps
# this function should just be a "decompose", returning a list of simple faces that can
# be passed on to actual drawing functions.

drawStarPolygon <- function(face, coords, col="grey", alpha=1, offset=c(0,0,0), debug=F)
{
  # Find all line segment pairs of this face, excluding line segments that are already adjecent to eachother.
  segmentPairs <- as.data.table(t(combn(seq(length(face)),2)))
  segmentPairs <- segmentPairs[abs(V1-V2)>1 & abs(V1-(V2%%length(face)))>1]
  
  # A-B and C-D are the vertex indices of the segment pairs that will be checked for intersection
  A <- face[unlist(segmentPairs[,1])]
  B <- shiftrotate(face)[unlist(segmentPairs[,1])]
  C <- face[unlist(segmentPairs[,2])]
  D <- shiftrotate(face)[unlist(segmentPairs[,2])]
  
  # The segment intersections will be calculated from the coordinates of the 4 points
  intersections <- data.table( 
    seg1From = A, seg1To = B, seg2From = C, seg2To = D, # keep the vertex indices
    Ax = coords[A,1], Ay = coords[A,2], Az = coords[A,3],
    Bx = coords[B,1], By = coords[B,2], Bz = coords[B,3],
    Cx = coords[C,1], Cy = coords[C,2], Cz = coords[C,3],
    Dx = coords[D,1], Dy = coords[D,2], Dz = coords[D,3])
  
  # The math to calculate the intersection is easily derived. The intersections are at A + alpha BA = C + beta DC. We want
  # both alpha and beta and both in a numerical stable way because the very compact list of segment pairs we started with,
  # so we may need to swap AB and CD and thus alpha and beta etc. But can still be singular!
  
  intersections[, beta  := ((Cx-Ax)*(By-Ay) - (Cy-Ay)*(Bx-Ax)) / ((Dy-Cy)*(Bx-Ax) - (Dx-Cx)*(By-Ay))]
  intersections[, alpha := ((Ax-Cx)*(Dy-Cy) - (Ay-Cy)*(Dx-Cx)) / ((By-Ay)*(Dx-Cx) - (Bx-Ax)*(Dy-Cy))]
  if (any(is.na(intersections$beta) | is.na(intersections$alpha))) {
    intersections[, beta  := ((Cz-Az)*(By-Ay) - (Cy-Ay)*(Bz-Az)) / ((Dy-Cy)*(Bz-Az) - (Dz-Cz)*(By-Ay))]
    intersections[, alpha := ((Az-Cz)*(Dy-Cy) - (Ay-Cy)*(Dz-Cz)) / ((By-Ay)*(Dz-Cz) - (Bz-Az)*(Dy-Cy))]
    if (any(is.na(intersections$beta) | is.na(intersections$alpha))) {
      intersections[, beta  := ((Cz-Az)*(Bx-Ax) - (Cx-Ax)*(Bz-Az)) / ((Dx-Cx)*(Bz-Az) - (Dz-Cz)*(Bx-Ax))]
      intersections[, alpha := ((Az-Cz)*(Dx-Cx) - (Ax-Cx)*(Dz-Cz)) / ((Bx-Ax)*(Dz-Cz) - (Bz-Az)*(Dx-Cx))]
      if (any(is.na(intersections$beta) | is.na(intersections$alpha))) {
        stop("No intersections - still singular")
      }
    }
  }
  
  # Intersections are only the ones for which alpha and beta are in 0-1 range. The intersections are not in a particular
  # order although we would want that. Ordering them by inner angle in 3D turns out to be mathematically impossible although
  # potentially we could use the direction of the inner vectors wrt the direction of the inner vector of the whole face.
  intersections <- intersections[alpha<1 & alpha>0 & beta<1 & beta>0]
  intersections[, c("Ix","Iy","Iz") := list(Ax + alpha*(Bx-Ax), Ay + alpha*(By-Ay), Az + alpha*(Bz-Az))]
  intersections[, intersectionIdx := seq(.N)]
  
  # for debugging only:
  if (debug) {
    spheres3d(intersections$Ix, intersections$Iy, intersections$Iz, color="green", radius = 0.02)
    text3d(intersections$Ix+0.05, intersections$Iy, intersections$Iz, color="black", text=intersections$intersectionIdx)
  }
  
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
  triangles3d(offset[1] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Ix_1[t], triangles$Ix_2[t], coords[triangles$from[t],1])})),
              offset[2] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Iy_1[t], triangles$Iy_2[t], coords[triangles$from[t],2])})),
              offset[3] + as.numeric(sapply(seq(nrow(triangles)), function(t) {c(triangles$Iz_1[t], triangles$Iz_2[t], coords[triangles$from[t],3])})),
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
  drawPolygonTriangulate(order(middlepolygon$rank), as.matrix(middlepolygon[,c("Ix","Iy","Iz"),with=F]), col=col, offset=offset)
}

drawAxes <- function()
{
  arrow3d(p0 = c(0,0,0), p1=c(1.8,0,0), width=0.2, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,1.8,0), width=0.2, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,0,1.8), width=0.2, type="rotation")
  texts3d(c(2,0,0), c(0,2,0), c(0,0,2), text = c("x","y","z"), color="black")
}

drawPolygon <- function(face, coords, col="grey", alpha=1, offset=c(0,0,0), label=NULL, drawlines=F)
{
  if (drawlines) {
    lines3d(coords[c(face, face[1]),1] + offset[1],
            coords[c(face, face[1]),2] + offset[2],
            coords[c(face, face[1]),3] + offset[3], color="blue")
  }
  if (!is.null(label)) {
    center <- apply(coords[face,], 2, mean)
    text3d(offset[1] + center[1], offset[2] + center[2], offset[3] + center[3], text = label, color="black")
  }
  
  # NB face orientation is not normalized
  if (length(face) == 3) {
    triangles3d( offset[1] + coords[face,1],
                 offset[2] + coords[face,2],
                 offset[3] + coords[face,3],
                 col=col, alpha=alpha)
  } else if (length(face) > 3) {
    ang <- innerAngles(coords[face,])
    
    if((sum(ang) > 2*pi) & !deltaEquals(sum(ang), 2*pi)) {
      drawStarPolygon(face, coords, col, alpha, offset)
      # drawPolygonTriangulate(face, coords, col, alpha, offset)
    } else {
      drawPolygonTriangulate(face, coords, col, alpha, offset)
      # polygon3d( offset[1] + coords$x[face],
      #            offset[2] + coords$y[face], 
      #            offset[3] + coords$z[face], 
      #            col=col, alpha=alpha)
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
    
    # coords
    spheres3d(x + p$coords[,1], y + p$coords[,2], z + p$coords[,3], color="green", radius = 0.02)
    text3d(x + (1+spacing)*p$coords[,1], y + (1+spacing)*p$coords[,2], z + (1+spacing)*p$coords[,3], text = seq(nrow(p$coords)), color="blue")
    
    # edges
    if("coordPairToEdge" %in% names(p)) {
      for (i in 1:(nrow(p$coordPairToEdge)-1)) {
        for (j in (i+1):nrow(p$coordPairToEdge)) {
          if (p$coordPairToEdge[i,j] != 0) {
            mid <- apply(p$coords[c(i,j),],2,mean)
            text3d(mid[1], mid[2], mid[3], text = paste0("e",p$coordPairToEdge[i,j]), color="black")
          }
        }
      }
    }    
  } else {
    alpha <- 1
  }
  if (!debug) {
    # avoid heavy description call in debug mode
    label <- paste(label, "(", description(p), ")")
  }
  if (nchar(label) > 0) {
    text3d(x, y + min(p$coords[,2]) - 1, z, text = label, color = "black", cex=0.7, pos = 1)
  }
  if (length(p$faces) > 0) { 
    if (length(p$bodies) > 1) {
      bodyColors <- rainbow(length(p$bodies))
    }
    faceType <- as.integer(factor(sapply(p$faces, length))) # faces considered same just by nr of edges
    if (max(faceType) > 1) {
      faceTypeColors <- rainbow(max(faceType))  
    }
    for (f in seq(length(p$faces))) {
      if (length(p$bodies) > 1) {
        faceColor <- bodyColors[which(sapply(p$bodies, function(b) { return(f %in% b)}))]
      } else {
        if (max(faceType) > 1) {
          faceColor <- faceTypeColors[faceType[f]]
        } else {
          faceColor <- rainbow(length(p$faces))[f]
        }
      }
      
      drawPolygon(p$faces[[f]], p$coords, faceColor, alpha, c(x, y, z), label=ifelse(debug,paste0("F",f),""), drawlines=debug) 
    }
  }
}

drawPoly <- function(p, start = c(0, 0, 0), delta = c(2, 0, 0), label = "", debug=F)
{
  if (!is.null(names(p))) { # not testing whether there is a name, testing whether this is a list with poly's or not
    drawSinglePoly(p, start[1], start[2], start[3], ifelse(is.null(p$name), "", p$name), debug)
  } else {
    for (i in seq(length(p))) {
      drawSinglePoly(p[[i]], start[1] + (i-1)*delta[1], start[2] + (i-1)*delta[2], start[3] + (i-1)*delta[3], p[[i]]$name, debug)  
    }
    rgl.texts(start[1], start[2] + 2, start[3], text = label, color="blue", pos = 4, cex = 1)
  }
}
# rgl.close()


# draw a n/d polygon
testDrawPolygon <- function(n, d=1, label="")
{
  rgl_init()
  clear3d()
  drawAxes()
  angles <- ((0:(n-1))/n)*2*pi
  coords <- as.matrix(data.table(x=sin(angles),y=cos(angles),z=0))
  face <- ((((1:n)-1)*d)%%n)+1 #+ as.numeric(sapply(0:(d-1), function(x){return(rep(x,n%/%d))}))
  
  #innerAngles(coords[face])
  drawPolygon(face, coords, label=paste(n,d,sep="/"), drawlines = T)
  # simplified <- list(coords = coords, faces = list(face))
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

## Debugging / Testing

testDrawPoly <- function()
{
  clear3d()
  drawAxes()
  drawSinglePoly(tetrahedron)  
  drawSinglePoly(octahedron)
  
  drawSinglePoly(tetrahedron, debug=T)  
  
  drawSinglePoly(cube)
  
  drawPoly(smallStellatedDodecahedron)
}

testDrawGallery <- function()
{
  clear3d()
  drawAxes()
  drawPoly(Platonics, start = c(1, 1, 1), delta = c(1, 1, 1))
  
  drawPoly(KeplerPoinsots, start = c(1, 1, 3), delta = 1.5*c(1, 1, 1))
}


