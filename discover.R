# Discover new polyhedra

source("polyhedra.R")

discover <- function(p, debug=F)
{
  clear3d()
  
  # find unique distance pairs
  vexPairs <- combn(x = seq(nrow(p$vertices)), m = 2)
  vexPairs <- data.frame(t(vexPairs), 
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
simplified <- list(vertices=poly$vertices, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
print(simplified$faces[[1]])
print(simplified$faces[[2]])
drawSinglePoly(simplified, debug = T) # shows drawing error and orientation issue 
getTopology(simplified$faces) # loops

# with a fix we get into a loop in 
getTopology(simplified$faces, debug=T)


# intersections of two linesegments A-B and C-D
face <- simplified$faces[[1]]
# A <- face[1]
# B <- face[2]
# C <- face[4] # skip 2-3 and 1-5 as these intersections are covered by definition
# D <- face[5]

# find all segment pairs of this face, excluding segments that are already adjecent to eachother
segmentpairstartindices<-as.data.frame(combn(seq(length(face)),2))
#segmentpairstartindices<- expand.grid(seq(length(face)), seq(length(face)))
segmentpairstartindices <- segmentpairstartindices[,
                                                   apply(segmentpairstartindices,2,function(pair) {
                                                     abs(pair[1]-pair[2])>1 &
                                                       abs(pair[1]-(pair[2]%%length(face)))>1})]

# vertex indices of segments AB and CD
A <- face[as.numeric(segmentpairstartindices[1,])]
B <- shift(face)[as.numeric(segmentpairstartindices[1,])]
C <- face[as.numeric(segmentpairstartindices[2,])]
D <- shift(face)[as.numeric(segmentpairstartindices[2,])]

# coordinates of the segments
segments <- data.frame( seg1From = A, seg1To = B, seg2From = C, seg2To = D,
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
segments$beta <- 
  ((segments$Cx-segments$Ax)*(segments$By-segments$Ay) - (segments$Cy-segments$Ay)*(segments$Bx-segments$Ax)) /
  ((segments$Dy-segments$Cy)*(segments$Bx-segments$Ax) - (segments$Dx-segments$Cx)*(segments$By-segments$Ay))
segments$alpha <-
  ((segments$Ax-segments$Cx)*(segments$Dy-segments$Cy) - (segments$Ay-segments$Cy)*(segments$Dx-segments$Cx)) /
  ((segments$By-segments$Ay)*(segments$Dx-segments$Cx) - (segments$Bx-segments$Ax)*(segments$Dy-segments$Cy))

segments$Ix <- segments$Ax + segments$alpha*(segments$Bx-segments$Ax)
segments$Iy <- segments$Ay + segments$alpha*(segments$By-segments$Ay)
segments$Iz <- segments$Az + segments$alpha*(segments$Bz-segments$Az)

# # just for checking
# segments$Jx <- segments$Cx + segments$beta*(segments$Dx-segments$Cx)
# segments$Jy <- segments$Cy + segments$beta*(segments$Dy-segments$Cy)
# segments$Jz <- segments$Cz + segments$beta*(segments$Dz-segments$Cz)

spheres3d(segments$Ix, segments$Iy, segments$Iz, color="green", radius = 0.02)

# intersect if beta > 0 and beta < 1 (delta)
# gives a minimal list of new intersections

# segments is a very compact representation, expand to easier to use listing all combinations
allsegs <- rbind(segments[, c(1:2, 18:21)], 
                 data.frame(seg1From = segments$seg1To,
                            seg1To = segments$seg1From,
                            alpha = (1-segments$alpha),
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz),
                 data.frame(seg1From = segments$seg2From,
                            seg1To = segments$seg2To,
                            alpha = segments$beta,
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz),
                 data.frame(seg1From = segments$seg2To,
                            seg1To = segments$seg2From,
                            alpha = (1-segments$beta),
                            Ix = segments$Ix, Iy = segments$Iy, Iz = segments$Iz))
