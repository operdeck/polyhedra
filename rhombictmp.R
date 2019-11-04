# rhombic

p <- dodecahedron
p <- icosahedron
p<-greatDodecahedron
clear3d()
drawPoly(p, debug = T)

rgl_init(new.device = T)
clear3d()
drawAxes()

# construct new coordinates from the edges of the current solid
newFace <- list()
for (e in seq(nrow(p$edgeToFaces))) {
  print(e)
  F1 <- p$edgeToFaces[e,1]
  F2 <- p$edgeToFaces[e,2]
  pts <- ((which(e == p$coordPairToEdge) - 1) %% ncol(p$coordPairToEdge))+1
  P1 <- pts[1]
  P2 <- pts[2]
  
  drawPolygon(p$faces[[F1]], p$coords, drawlines = T, drawvertices = T, label=paste0("F",F1))
  drawPolygon(p$faces[[F2]], p$coords, drawlines = T, drawvertices = T, label=paste0("F",F2))
  
  n1 <- normal(p$coords[p$faces[[F1]][1],],
               p$coords[p$faces[[F1]][2],],
               p$coords[p$faces[[F1]][3],])
  
  n2 <- normal(p$coords[p$faces[[F2]][1],],
               p$coords[p$faces[[F2]][2],],
               p$coords[p$faces[[F2]][3],])
  
  center1 <- apply(p$coords[p$faces[[F1]],], 2, mean)
  center2 <- apply(p$coords[p$faces[[F2]],], 2, mean)
  
  polyLines(c(center1, center1+n1), color = "red")
  polyLines(c(center2, center2+n2), color = "red")

  dp1p2 <- distance(p$coords[P1,], p$coords[P2,])
  dn1n2 <- distance(n1, n2)
  alpha <- dp1p2/dn1n2
  
  p1_a <- p$coords[P1,] + alpha*n1
  p1_b <- p$coords[P1,] + alpha*n2
  p2_a <- p$coords[P2,] + alpha*n2
  p2_b <- p$coords[P2,] + alpha*n1
  
  polyLines(c(p1_a, p1_b, p2_a, p2_b, p1_a), color="green")
  
  newFace[[1+length(newFace)]] <- 
    data.table(x=c(p1_a[1], p1_b[1], p2_a[1], p2_b[1]),
               y=c(p1_a[2], p1_b[2], p2_a[2], p2_b[2]),
               z=c(p1_a[3], p1_b[3], p2_a[3], p2_b[3]),
               edge=e,
               frompoint=c(P1, P1, P2, P2), 
               fromface=c(F1, F2, F2, F1))
}

nc <- rbindlist(newFace) # but: contains duplicate coordinates
nc[, idx := seq(.N)]
nc[, idx := idx[1], by=c("fromface","frompoint")]
nc[, idx := as.integer(factor(idx))] # same coordinates now have same index

# faces from the edges
squaresFromEdges <- nc[,.(f = list(idx)),by=edge]

# de-duplicated coordinates
newCoords <- as.matrix( nc[, .(x = x[1], y = y[1], z = z[1]), by=idx][order(idx)][, 2:4] )

# find new coords for existing faces
newTopo <- unique(nc[, c("frompoint","fromface","idx")])

# construct new faces mostly by lookup
facesFromOldFaces <-
  lapply(seq(length(p$faces)), 
         function(oldFace) {
           return(sapply(p$faces[[oldFace]], 
                         function(i) { return(newTopo[fromface==oldFace & frompoint==i]$idx)}))})

facesFromOldVertices <-
  lapply(seq(length(p$vertexFigures)),
         function(oldVertex) { 
           v <- p$vertexFigures[[oldVertex]]
           return(sapply(v$faces, function(i) {return(newTopo[fromface==i & frompoint==v$center]$idx)}))})

newRhombic <-
  setPoly(coords = newCoords/vectorlength(newCoords), # unit length
          faces = c(squaresFromEdges$f, facesFromOldFaces, facesFromOldVertices),
          name = "rhombic")

clear3d()
drawPoly(newRhombic)

