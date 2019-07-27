# Discover new polyhedra

source("polyhedra.R")

discover <- function(p, debug=F)
{
  clear3d()
  
  # find unique distance pairs
  vexPairs <- combn(x = seq(nrow(p$vertices)), m = 2)
  vexPairs <- data.table(t(vexPairs), 
                         dist=sapply(seq(ncol(vexPairs)), function(col) { return(distance(p$vertices[vexPairs[1,col],], p$vertices[vexPairs[2,col],]))}))
  vexPairs <- vexPairs[order(vexPairs$dist),]
  vexPairs$isUnique <- T
  row.names(vexPairs) <- NULL
  for (i in 2:nrow(vexPairs)) {
    if (deltaEquals(vexPairs$dist[i-1], vexPairs$dist[i])) { vexPairs$isUnique[i] <- F }
  }
  vexPairs <- vexPairs[(vexPairs$isUnique),] # this is now a short list of vertex index pairs
  
  searchgrid <- expand.grid(facedim = 3:5, vertexdim = 3:20, vexPair = 1:nrow(vexPairs))
  searchgrid$isPolygon <- F
  searchgrid$vertex1 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$X1[searchgrid$vexPair[i]]})
  searchgrid$vertex2 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$X2[searchgrid$vexPair[i]]})
  for (i in seq(nrow(searchgrid))) {
    # NB it also creates polys with holes in it
    # some of the ones from dodecahedron are new!
    poly <- buildRegularPoly(p$vertices, searchgrid$facedim[i], searchgrid$vertexdim[i], c(searchgrid$vertex1[i], searchgrid$vertex2[i]), debug=debug)
    if (length(poly$faces) > 0) {
      print("Success!")
      searchgrid$isPolygon[i] <- T
      print(searchgrid[i, ])
      
      #print(poly)
      drawSinglePoly(poly, debug=T)
      
      #descr <- description(poly)
      descr<-"issues w discover(dodecahedron) 3rd"
      print(descr)
      searchgrid$isPolygon[i] <- T
      
      drawSinglePoly(poly, x=2*sum(searchgrid$isPolygon), label=paste("Discovery",sum(searchgrid$isPolygon),"from",p$name), debug=F)
    }
  }
}

# if (debug) print("All done!")
#discover(archi(dodecahedron), debug=F)

stop("issue:")
#discover(dodecahedron) 
poly <- buildRegularPoly(dodecahedron$vertices, 5, 6, c(1, 4))
clear3d()
drawSinglePoly(poly, debug = T) # this shows the error drawing {5/2}
topo <- getTopology(poly$faces, debug=T) # this errors out on edge 4-1 / into a loop

# face is "below" origin
faces <- poly$faces
debug <- T
simplified <- list(vertices=poly$vertices, faces=list(poly$faces[[3]], poly$faces[[4]]), name="debug")
print(simplified$faces[[1]])
print(simplified$faces[[2]])
drawSinglePoly(simplified, debug = T) # shows drawing error and orientation issue when using [[1]],[[2]]
getTopology(simplified$faces) # loops

# with a fix we get into a loop in 
getTopology(simplified$faces, debug=T)


# intersections of two linesegments A-B and C-D
face <- simplified$faces[[2]]
# A <- face[1]
# B <- face[2]
# C <- face[4] # skip 2-3 and 1-5 as these intersections are covered by definition
# D <- face[5]

# drawSinglePoly(simplified, debug=T)

# find all segment pairs of this face, excluding segments that are already adjecent to eachother
segmentpairstartindices<-as.data.table(t(combn(seq(length(face)),2)))
segmentpairstartindices <- segmentpairstartindices[abs(V1-V2)>1 & abs(V1-(V2%%length(face)))>1]

# vertex indices of segments AB and CD
A <- face[unlist(segmentpairstartindices[,1])]
B <- shift(face)[unlist(segmentpairstartindices[,1])]
C <- face[unlist(segmentpairstartindices[,2])]
D <- shift(face)[unlist(segmentpairstartindices[,2])]

# coordinates of the segments
segments <- data.table( 
  seg1From = A, seg1To = B, seg2From = C, seg2To = D,
  Ax = simplified$vertices$x[A],
  Ay = simplified$vertices$y[A],
  Az = simplified$vertices$z[A],
  Bx = simplified$vertices$x[B],
  By = simplified$vertices$y[B],
  Bz = simplified$vertices$z[B],
  Cx = simplified$vertices$x[C],
  Cy = simplified$vertices$y[C],
  Cz = simplified$vertices$z[C],
  Dx = simplified$vertices$x[D],
  Dy = simplified$vertices$y[D],
  Dz = simplified$vertices$z[D])
# intersect = A + alpha BA = C + beta DC
segments[, beta  := ((Cx-Ax)*(By-Ay) - (Cy-Ay)*(Bx-Ax)) / ((Dy-Cy)*(Bx-Ax) - (Dx-Cx)*(By-Ay))]
segments[, alpha := ((Ax-Cx)*(Dy-Cy) - (Ay-Cy)*(Dx-Cx)) / ((By-Ay)*(Dx-Cx) - (Bx-Ax)*(Dy-Cy))]

segments <- segments[alpha<1 & alpha>0 & beta<1 & beta>0]

segments[, c("Ix","Iy","Iz") := list(Ax + alpha*(Bx-Ax), Ay + alpha*(By-Ay), Az + alpha*(Bz-Az))]
center <- apply(segments[, c("Ix","Iy","Iz")], 2, mean)
# sapply(seq(nrow(segments)), function(i) {
#   return(vectorAngle(t(t(segments[1, c("Ix","Iy","Iz")]) - center),
#                      t(t(segments[i, c("Ix","Iy","Iz")]) - center)))
# })*360/(2*pi)
segments[, intersectionIdx := seq(.N)]
# to draw, sort by angle is not possible in 3D
# start somewhere, find the closest that is not used etc 

spheres3d(segments$Ix, segments$Iy, segments$Iz, color="green", radius = 0.02)
text3d(segments$Ix+0.05, segments$Iy, segments$Iz, color="black", text=segments$intersectionIdx)

# intersect if beta > 0 and beta < 1 (delta)
# gives a minimal list of new intersections

# segments is a very compact representation, expand to easier to use listing all combinations
allsegs <- rbind(data.table(from = segments$seg1From,
                            to = segments$seg1To,
                            alpha = segments$alpha,
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz, intersectionIdx = segments$intersectionIdx),
                 data.table(from = segments$seg1To,
                            to = segments$seg1From,
                            alpha = (1-segments$alpha),
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz, intersectionIdx = segments$intersectionIdx),
                 data.table(from = segments$seg2From,
                            to = segments$seg2To,
                            alpha = segments$beta,
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz, intersectionIdx = segments$intersectionIdx),
                 data.table(from = segments$seg2To,
                            to = segments$seg2From,
                            alpha = (1-segments$beta),
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz, intersectionIdx = segments$intersectionIdx))
allsegs <- allsegs[, isClosestIntersection := (seq(.N) == which.min(alpha)), by=c("from", "to")][(isClosestIntersection)][order(from, to)]

# stars
triangles <- dcast(allsegs[, triangleside:=seq(.N), by=from], from~triangleside, value.var=c("Ix","Iy","Iz"))
for (t in seq(nrow(triangles))) {
  triangles3d(c(triangles$Ix_1[t], triangles$Ix_2[t], simplified$vertices$x[triangles$from[t]]),
              c(triangles$Iy_1[t], triangles$Iy_2[t], simplified$vertices$y[triangles$from[t]]),
              c(triangles$Iz_1[t], triangles$Iz_2[t], simplified$vertices$z[triangles$from[t]]),
              col="yellow")
}

# middle
middlepolygon <- unique(allsegs[, c("intersectionIdx", "Ix", "Iy", "Iz")])
middlepolygon[, isUsed := F]
currentPt <- 1
repeat {
  middlepolygon[(currentPt), isUsed := T]
  if (sum(!middlepolygon$isUsed) == 0) {
    nextPt <- 1
  } else {
    distances <- sqrt(rowSums((t(t(middlepolygon[!(isUsed),2:4]) - as.numeric(middlepolygon[currentPt,2:4])))^2))  
    nextPt <- which(!middlepolygon$isUsed)[which.min(distances)[1]]
  }
  triangles3d(c(middlepolygon$Ix[currentPt], middlepolygon$Ix[nextPt], center[1]),
              c(middlepolygon$Iy[currentPt], middlepolygon$Iy[nextPt], center[2]),
              c(middlepolygon$Iz[currentPt], middlepolygon$Iz[nextPt], center[3]),
              col="green")
  if (sum(!middlepolygon$isUsed) == 0) break
  currentPt <- nextPt
}

# todo check w {8/3} also
clear3d()
n <- 7
d <- 3
angles <- ((0:(n-1))/n)*2*pi
coords <- data.table(x=sin(angles),y=cos(angles),z=0)
face <- ((((1:n)-1)*d)%%n)+1
simplified <- list(vertices = coords, faces = list(face))
drawSinglePoly(simplified, debug=T)
