library(rgl)
library(data.table)
source("geometry.R")

# NB gradient colors possible e.g.
# triangles3d(x=c(1,3,2),y=c(2,3,1),z=c(0,0,0),color=c("blue","yellow"))
# args to 3d can be matrix e.g. points3d(cube$coords)
# lines3d(matrix(c(1,2,3,1,2,6),ncol=3,byrow=T),color="grey")

eToStr <- function(poly, e)
{
  return (paste0(e, " (", paste(sort(poly$edgeToVertices[e,]),collapse="-"), ")"))  
}

fToStr <- function(f)
{
  return(paste0("F",f))
}

drawInit <- function(new.device = FALSE, bg = "white", width = 640, height = width) {
  if( new.device | rgl.cur() == 0 ) {
    rgl.open()
    par3d(windowRect = 50 + c( 0, 0, width, height ) )
    rgl.bg(color = bg )
  }
  rgl.clear(type = c("shapes", "bboxdeco"))
  rgl.viewpoint(theta = 15, phi = 20, zoom = 0.7)
}

drawSegments <- function(coordsFrom, coordsTo, ...)
{
  if (!is.matrix(coordsFrom) & !is.matrix(coordsTo)) {
    segments3d(c(coordsFrom[1], coordsTo[1]), 
               c(coordsFrom[2], coordsTo[2]), 
               c(coordsFrom[3], coordsTo[3]), ...)
  } else {
    if (is.matrix(coordsFrom)) {
      if (dim(coordsFrom)[2]==2) coordsFrom <- as.matrix(data.table(coordsFrom, z=0)) # provisional 2D support
    }
    if (is.matrix(coordsTo)) {
      if (dim(coordsTo)[2]==2) coordsTo <- as.matrix(data.table(coordsTo, z=0)) # provisional 2D support
    }
    if (!is.matrix(coordsFrom)) coordsFrom <- matrix(coordsFrom, ncol=3, byrow=T)
    if (!is.matrix(coordsTo)) coordsTo <- matrix(coordsTo, ncol=3, byrow=T)
    coordsFrom <- data.table(coordsFrom)[,row:=seq(.N)]
    coordsTo <- data.table(coordsTo)[,row:=seq(.N)]
    segments3d(rbind(coordsFrom, coordsTo)[order(row)], ...)
  }
}

drawLines <- function(coords, closed=F, ...)
{
  if (is.matrix(coords)) {
    if (dim(coords)[2]==2) coords <- as.matrix(data.table(coords, z=0)) # provisional 2D support  
  } 
  if (!is.matrix(coords)) coords <- matrix(coords, ncol=3, byrow=T)
  lines3d(coords, ...)
  if (closed) lines3d(coords[c(nrow(coords),1),], ...)
}

drawTexts <- function(coords, text, ...)
{
  if (is.matrix(coords)) {
    if (dim(coords)[2]==2) coords <- as.matrix(data.table(coords, z=0)) # provisional 2D support  
  } 
  if (!is.matrix(coords)) coords <- matrix(coords, ncol=3, byrow=T)
  texts3d(coords, text=text, ...)
}

drawDots <- function(coords, label=NULL, radius=0.01, ...)
{
  if (is.matrix(coords)) {
    if (dim(coords)[2]==2) coords <- as.matrix(data.table(coords, z=0)) # provisional 2D support  
  } 
  if (!is.matrix(coords)) coords <- matrix(coords, ncol=3, byrow=T)
  spheres3d(coords, radius=radius, ...)
  if (!is.null(label)) texts3d(1.1*coords, text=label, color="black")
}

#open3d()
#shade3d( icosahedron3d() )
#writeOBJ(filename)
# todo figure out how to update that obj file with a mtl
# so to set the face colrs. which app shows these colors
# properly anyway?

#rgl.spheres(coords$x, coords$y, coords$z, r=0.2, color="yellow")
#rgl.texts(coords$x, coords$y, coords$z, text = seq(nrow(coords)), color="red")

# 3D tricks (movie, HTML interactive export, points highlighting)
# http://www.sthda.com/english/wiki/a-complete-guide-to-3d-visualization-device-system-in-r-r-software-and-data-visualization#prepare-the-data


# the rgl method is not reliable so still doing it myself
drawConvexPolygon <- function(coords, ...)
{
  if (nrow(coords) < 3) stop("Face to triangulate should have >= 3 points.")

  for (t in (seq(nrow(coords)-2)+1)) {
    # NB not sure about the orientation of the triangle - may have to check on this
    triangles3d(coords[c(1, t, t+1),], ...)
  }
}

# Draws a star polygon with intersecting edges
drawStarPolygon <- function(coords, ...)
{
  #drawTexts(coords,text=seq(nrow(coords)))
  #drawLines(coords, closed=T)
  intersections <- as.data.table(t(combn(seq(nrow(coords)),2)))[ V2 != (V1+1) & (V1-1) != (V2 %% nrow(coords))]
  n<-nrow(coords)
  intersectionPts <- rbindlist(lapply(seq(nrow(intersections)), function(i) {
    intersection <- intersect_2Segments(coords[intersections[i]$V1,], coords[(intersections[i]$V1 %% n)+1,], 
                                        coords[intersections[i]$V2,], coords[(intersections[i]$V2 %% n)+1,])
    if (intersection$status != "intersect" || is.na(intersection$alpha) || is.na(intersection$beta)) return()
    return(data.table(seg1a = intersections[i]$V1, seg1b = (intersections[i]$V1 %% n)+1,
                      seg2a = intersections[i]$V2, seg2b = (intersections[i]$V2 %% n)+1,
                      f1 = intersection$alpha, f2 = intersection$beta,
                      x = intersection$I0[1], y = intersection$I0[2], z = intersection$I0[3]))
  }))
  if (nrow(intersectionPts) == 0) {
    stop("No intersections") # delegate to convex?
  }
  
  currentSegStart <- 1
  center <- apply(coords, 2, mean)
  cnt <- 0
  repeat {
    currentSegEnd <- (currentSegStart %% n)+1
    alphas <- sapply(seq(nrow(intersectionPts)), function(i) {
      if (intersectionPts[i]$seg1a==currentSegStart & intersectionPts[i]$seg1b==currentSegEnd) return(intersectionPts[i]$f1)
      if (intersectionPts[i]$seg1b==currentSegStart & intersectionPts[i]$seg1a==currentSegEnd) return(1-intersectionPts[i]$f1)
      if (intersectionPts[i]$seg2a==currentSegStart & intersectionPts[i]$seg2b==currentSegEnd) return(intersectionPts[i]$f2)
      if (intersectionPts[i]$seg2b==currentSegStart & intersectionPts[i]$seg2a==currentSegEnd) return(1-intersectionPts[i]$f2)
      return(Inf)
    })
    l <- which.min(alphas)
    #drawDots(intersectionPts[l, 7:9], radius=0.03, color="purple")
    triangle1 <- matrix(c(coords[currentSegStart, ], center, intersectionPts[l, 7:9]), ncol = 3, byrow=T)
    # print(triangle1)
    # print(class(triangle1))
    # if (!isNormalOutwardFacing(triangle1)) {
    #   triangle1 <- triangle1[c(3,2,1),]  
    # }
    triangles3d(triangle1, ...)
    if (intersectionPts[l]$seg1a==currentSegStart | intersectionPts[l]$seg1b==currentSegStart) {
      pNext <- ifelse(intersectionPts[l]$f2 < 0.5, intersectionPts[l]$seg2a, intersectionPts[l]$seg2b)   
    } else {
      pNext <- ifelse(intersectionPts[l]$f1 < 0.5, intersectionPts[l]$seg1a, intersectionPts[l]$seg1b)   
    }
    triangle2 <- matrix(c(center, intersectionPts[l, 7:9], coords[pNext, ]), ncol = 3, byrow=T)
    # if (!isNormalOutwardFacing(triangle2)) {
    #   triangle2 <- triangle1[c(3,2,1),]  
    # }
    triangles3d(triangle2, ...)
    
    cnt <- cnt+1
    if (cnt > nrow(coords)) break;
    if (pNext == 1) break
    currentSegStart <- pNext
  }
  if (cnt != nrow(coords)) warning(paste("Different # of intersection points than expected. Found", 
                                      cnt, "intersections from", nrow(intersectionPts), "candidates, expected", nrow(coords)))
}

drawAxes <- function()
{
  arrow3d(p0 = c(0,0,0), p1=c(1.8,0,0), width=0.1, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,1.8,0), width=0.1, type="rotation")
  arrow3d(p0 = c(0,0,0), p1=c(0,0,1.8), width=0.1, type="rotation")
  texts3d(c(2,0,0), c(0,2,0), c(0,0,2), text = c("x","y","z"), color="black")
}

drawPolygon <- function(face, coords, col="grey", alpha=1, offset=c(0,0,0), label=NULL, drawlines=F, drawvertices=F)
{
  spacing <- 0.1
  center <- apply(coords[face,], 2, mean)
  if (drawlines) {
    drawLines(t(t(coords[face,]) + offset), closed=T, color="orange")
  }
  if (drawvertices) {
    drawDots(t(offset + t(coords[face,])), color="green", radius = 0.02)
    drawTexts(t(offset + (1+spacing)*t(coords[face,])), text = face, color="blue")
  }
  if (!is.null(label)) {
    drawTexts(offset + (1+spacing)*center, text = label, color="black")
  }
  
  # NB face orientation is not normalized
  if (length(face) == 3) {
    drawConvexPolygon(t(t(coords[face,])+offset), col=col, alpha=alpha)
  } else if (length(face) > 3) {
    ang <- innerAngles(coords[face,])
    
    if((sum(ang) > 2*pi) & !deltaEquals(sum(ang), 2*pi)) {
      if (!isFlatFace(coords[face,])) {
        drawConvexPolygon(t(t(coords[face,])+offset), col=col, alpha=alpha)
        drawTexts(offset + center, text = "!", color="white")
      } else {
        drawStarPolygon(t(t(coords[face,])+offset), col=col, alpha=alpha)
      }
    } else {
      drawConvexPolygon(t(t(coords[face,])+offset), col=col, alpha=alpha)
    }
  }
}

# Assigns colors to faces. 
# If there is one body and all faces are the same, applies the (first) color provider to the faces
# If there is one body but different types of faces, then applys the (first) color provider to the face types
# If there are multiple bodies, but just a single color provider, that provider is applied per body
# If there are multiple bodies and multiple providers, the first/second logic is applied per body

assignColors <- function(p, colorProvider = rainbow)
{
  isSingleBody <- (length(p$bodies) == 1)
  isSingleProvider <- (!is.list(colorProvider))

  colors <- sapply(seq(length(p$faces)), function(f) {
    fbody <- which(sapply(p$bodies, function(b) { return(f %in% b)}))
    faceTypeIndexInCurrentBody <- as.integer(factor(sapply(p$faces [p$bodies[[fbody]]], length)))

    if (!isSingleBody & !isSingleProvider) {
      # multiple bodies and multiple providers - one per body, recycling
      bodyColorProvider <- colorProvider[[ ((fbody-1) %% (length(colorProvider))) + 1 ]]
    } else {
      bodyColorProvider <- ifelse(isSingleProvider, colorProvider, colorProvider[[1]])
    }
    
    if (!isSingleBody & isSingleProvider) {
      # color per body
      faceColor <- bodyColorProvider(length(p$bodies))[fbody]    
    } else {
      if (max(faceTypeIndexInCurrentBody) == 1) {
        # color per face
        faceColor <- bodyColorProvider(length(p$bodies[[fbody]])) [ which(p$bodies[[fbody]] == f) ]  
      } else {
        # color per face type
        faceColor <- bodyColorProvider(max(faceTypeIndexInCurrentBody)) [faceTypeIndexInCurrentBody[ which(p$bodies[[fbody]] == f) ]]  
      } 
    }    
    return(faceColor)
  })    
  
  return(colors)
}

# Draw a single polygon. Offset is optional.
drawSinglePoly <- function(p, offset=c(0,0,0), label=ifelse(is.null(p$name),"",p$name), debug=F, colorProvider = rainbow)
{
  if (debug) {
    drawAxes()
    spacing <- 0.1
    alpha <- 0.6 # in debug make somewhat transparent
    
    # coords
    drawDots(t(offset + t(p$coords)), color="green", radius = 0.02)
    drawTexts(t(offset + (1+spacing)*t(p$coords)), text = seq(nrow(p$coords)), color="blue")
    
    # edges
    if("coordPairToEdge" %in% names(p)) {
      for (i in 1:(nrow(p$coordPairToEdge)-1)) {
        for (j in (i+1):nrow(p$coordPairToEdge)) {
          if (p$coordPairToEdge[i,j] != 0) {
            mid <- apply(p$coords[c(i,j),],2,mean)
            drawTexts(mid, text = paste0("e",p$coordPairToEdge[i,j]), color="black")
          }
        }
      }
    }    
  } else {
    alpha <- 1
  }
  if (!debug) {
    # avoid heavy description call in debug mode
    if (label != "") {
      label <- paste(label, "(", description(p), ")")
    }
  }
  if (nchar(label) > 0) {
    drawTexts( c(offset[1], offset[2] + min(p$coords[,2]) - 1, offset[3]), text = label, color = "black", cex=0.7, pos = 1)
  }
  if (length(p$faces) > 0) { 
    colors <- assignColors(p, colorProvider)
    for (f in seq(length(p$faces))) {
      drawPolygon(p$faces[[f]], p$coords, colors[f], alpha, offset, label=ifelse(debug,fToStr(f),""), drawlines=debug) 
    }
  }
}

drawPoly <- function(p, start = c(0, 0, 0), delta = c(2, 0, 0), label = "", debug=F, colorProvider = rainbow)
{
  if (!is.null(names(p))) { # not testing whether there is a name, testing whether this is a list with poly's or not
    drawSinglePoly(p, offset=start, ifelse(is.null(p$name), "", p$name), debug, colorProvider)
  } else {
    for (i in seq(length(p))) {
      drawSinglePoly(p[[i]], offset=start+(i-1)*delta, p[[i]]$name, debug, colorProvider)  
    }
    # overall label
    title3d(main=label, color="blue")
    #rgl.texts(start[1], start[2] + 2, start[3], text = label, color="blue", pos = 4, cex = 1)
  }
}
# rgl.close()

drawSkeleton <- function(p)
{
  clear3d()
  p<-cube
  p <- buildRegularPoly(dodecahedron$coords,
                                          polygonsize = 3,
                                          vertexsize = 3,
                                          exampleEdge = c(3, 8),
                                          name = "5 Tetrahedra")
  cols <- assignColors(p)
  edges <- rbindlist(lapply(which(upper.tri(p$coordPairToEdge) & (p$coordPairToEdge != 0)), function(edge) {
    dim <- nrow(p$coordPairToEdge)
    col <- (edge-1)%%dim+1
    row <- (edge-1)%/%dim+1
    return (data.table(vex1=row, vex2=col, F1=p$coordPairToFaces[row,col], F2=p$coordPairToFaces[col,row]))
  }))
  
  #drawSegments(p$coords[edges$vex1,], p$coords[edges$vex2,], color="green")
  
  d <- 0.1
  center <- apply(p$coords, 2, mean)
  for (i in seq_len(nrow(edges))) {
    p1 <- p$coords[edges$vex1[i],]
    p2 <- p$coords[edges$vex2[i],]
    
    d1 <- crossproduct(center-p1, p2-p1)
    d2 <- crossproduct(center-p2, p1-p2)
    
    x1a <- p1 + d*d1/vectorlength(d1)
    x2a <- p2 + d*d2/vectorlength(d2)
    x1b <- p1 - d*d1/vectorlength(d1)
    x2b <- p2 - d*d2/vectorlength(d2)
    
    drawConvexPolygon(matrix(c(x1a, x2a, x2b, x1b), ncol = 3, byrow = T), 
                      color=c(cols[edges$F1[i]], cols[edges$F2[i]]))
  }
}

# draw a n/d polygon
testDrawPolygon <- function(n, d=1)
{
  drawInit()
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
  testDrawPolygon(7,2) 
  testDrawPolygon(7,3) 
  testDrawPolygon(8,3) # more complex star
  #testDrawPolygon(8,2) # also a star, generation of this is slightly more complex
  
  
  p <- dodecahedron
  f1 <- p$coords[ p$faces[[1]], ]
  clear3d()
  polygon3d( f1[,1], f1[,2], f1[,3] ) # fine for simple planes
  
  p <- smallStellatedDodecahedron
  f1 <- p$coords[ p$faces[[1]], ]
  clear3d()
  polygon3d( f1[,1], f1[,2], f1[,3] ) # triangulation fails for star polyhedra
  drawStarPolygon(p$faces[[1]], p$coords) # works (TODO: why not interface with just coords?)
  drawPolygon(p$faces[[1]], p$coords, label="Hello", drawlines=T, drawvertices=T)
  
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
  drawSinglePoly(quasi(octahedron))
  
  drawPoly(smallStellatedDodecahedron)
  drawPoly(greatIcosahedron, colorProvider = function(n) { return(sample(c("green", "red"), n, replace = T))})
  
  drawSinglePoly( octahedron, colorProvider = heat.colors )
  
  drawSinglePoly(compose(icosahedron, dual(icosahedron)),
                 colorProvider = list( heat.colors, rainbow ))
  
  drawSinglePoly(compose(icosahedron, dual(icosahedron)),
                 colorProvider = list( function(n) { return( tail(rainbow(32),n)) }, function(n) { return( head(rainbow(32), n)) } ))
  
  library(colorspace)
  drawSinglePoly(icosahedron, colorProvider=qualitative_hcl)
}

testDrawGallery <- function()
{
  clear3d()
  drawAxes()
  drawPoly(Platonics, start = c(1, 1, 1), delta = c(1, 1, 1))
  
  drawPoly(KeplerPoinsots, start = c(1, 1, 3), delta = 1.5*c(1, 1, 1), colorProvider = heat.colors)
}

# # play3d(spin3d())

# drawConvexPolygon( cube$coords[cube$faces[[1]], ], alpha=0.5, texture="arabic.png", textype="rgb", shininess=60, color="red")
# bg3d(sphere = TRUE, texture = system.file("textures/sunsleep.png", package = "rgl"), back = "filled" )

# Replace DrawPoly with mfrow/title
#
# > clear3d()
# > mfrow3d(ceiling(sqrt(length(KeplerPoinsots))), ceiling(sqrt(length(KeplerPoinsots))), sharedMouse = T)
# > for (i in seq(length(KeplerPoinsots))) {
#   +   if (i > 1) next3d()
#   +   drawSinglePoly(KeplerPoinsots[[i]])
#   + }
# > rglwidget(elementId = "KeplerPoinsots")
#
# title3d
# 
