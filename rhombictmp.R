# rhombic

p <- dodecahedron

clear3d()
drawPoly(p, debug = T)

rgl_init(new.device = T)
clear3d()
drawAxes()

e <- 1 # 27 gives error

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
  lines3d( c(center1[1], center1[1]+n1[1]), c(center1[2], center1[2]+n1[2]), c(center1[3], center1[3]+n1[3]), color="red")
  lines3d( c(center2[1], center2[1]+n2[1]), c(center2[2], center2[2]+n2[2]), c(center2[3], center2[3]+n2[3]), color="red")
  
  dp1p2 <- distance(p$coords[P1,], p$coords[P2,])
  dn1n2 <- distance(n1, n2)
  alpha <- dp1p2/dn1n2
  
  p1_a <- p$coords[P1,] + alpha*n1
  p1_b <- p$coords[P1,] + alpha*n2
  p2_a <- p$coords[P2,] + alpha*n2
  p2_b <- p$coords[P2,] + alpha*n1
  
  lines3d( c(p1_a[1], p1_b[1], p2_a[1], p2_b[1], p1_a[1]),
           c(p1_a[2], p1_b[2], p2_a[2], p2_b[2], p1_a[2]),
           c(p1_a[3], p1_b[3], p2_a[3], p2_b[3], p1_a[3]), color="green")
}

