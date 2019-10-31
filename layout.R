# automatic layout gen

source("polyhedra.R")
source("draw.R")

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- quasi(cube)
drawPoly(xx, debug=T)

# for layout use a new device
rgl_init(new.device = T)
clear3d()
drawAxes()

# draw 2D face
drawFace2D <- function(face2d)
{
  points3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, color="black" )
  text3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, 
         text=face2d$points3D, color="blue" )
  lines3d(x = face2d$coords2D[,1][c(seq(nrow(face2d$coords2D)),1)], 
          y = face2d$coords2D[,2][c(seq(nrow(face2d$coords2D)),1)], z = 0, color="red" )
  center <- apply(face2d$coords2D, 2, mean)
  text3d(x = center[1], y = center[2], z = 0, text=paste0("F",face2d$face), color="blue" )
}

# Start with a layout with just the first face
layout <- list()

# place it on 2D plane
currentFace <- 1
coords3D <- xx$coords[xx$faces[[currentFace]],]
angles <- innerAngles(coords3D)
radii <- distance(coords3D)
offset <- c(1,1)
layout[[currentFace]] <- list( face = currentFace, 
                               points3D = xx$faces[[currentFace]],
                               coords2D = matrix(c(offset[1] + cos(cumsum(angles)) * radii, offset[2] + sin(cumsum(angles)) * radii),nrow = nrow(coords3D)))

drawFace2D(layout[[currentFace]])

# repeat from here
repeat {
  # find a free edge
  # one of xx$edgeToFaces that has [,1] in layout and [,2] not (or vice versa)
  placedFaces <- sapply(layout, function(f) {return(f$face)})
  candidateEdges <- which((xx$edgeToFaces[,1] %in% placedFaces & !(xx$edgeToFaces[,2] %in% placedFaces)) |
                            (xx$edgeToFaces[,2] %in% placedFaces & !(xx$edgeToFaces[,1] %in% placedFaces)))
  if (length(candidateEdges) == 0) stop("Done! (or perhaps need to place unconnected face)")
  candidateEdge <- candidateEdges[1] # TODO iterate/recurse here and fail if there are none
  if (xx$edgeToFaces[candidateEdge,1] %in% placedFaces) {
    currentFace <- xx$edgeToFaces[candidateEdge,1]
    otherFace <- xx$edgeToFaces[candidateEdge,2]
  } else {
    otherFace <- xx$edgeToFaces[candidateEdge,1]
    currentFace <- xx$edgeToFaces[candidateEdge,2]
  }
  vertices <- ((which(xx$coordPairToEdge==candidateEdge)-1) %/% nrow(xx$coordPairToEdge))+1
  cat("Placed:", currentFace, ", adding:", otherFace, fill = T)
  cat("vertices:", paste(vertices, collapse=","), fill=T)
  
  # place it on 2D plane
  coords3D <- xx$coords[xx$faces[[otherFace]],]
  angles <- innerAngles(coords3D)
  radii <- distance(coords3D)
  offset <- c(5,2)
  layout[[otherFace]] <- list( face = otherFace, 
                               points3D = xx$faces[[otherFace]],
                               coords2D = matrix(c(offset[1] + cos(cumsum(angles)) * radii, y = offset[2] + sin(cumsum(angles)) * radii),nrow = nrow(coords3D)))
  
  drawFace2D(layout[[otherFace]])
  
  # identify the vertices of the edge in both faces
  if (length(vertices) != 2) stop("An edge must have 2 vertices")
  p1current <- which(layout[[currentFace]]$points3D == vertices[1])
  p1other <- which(layout[[otherFace]]$points3D == vertices[1])
  p2current <- which(layout[[currentFace]]$points3D == vertices[2])
  p2other <- which(layout[[otherFace]]$points3D == vertices[2])
  if (length(c(p1current,p1other,p2current,p2other)) != 4) stop("All vertices must occur in both faces")
  
  # rotate and translate the connected face
  vecCurrent <- as.numeric(layout[[currentFace]]$coords2D[p2current,]) - as.numeric(layout[[currentFace]]$coords2D[p1current,])
  vecOther <- as.numeric(layout[[otherFace]]$coords2D[p2other,]) - as.numeric(layout[[otherFace]]$coords2D[p1other,])
  theta <-  vectorAngle2D(vecCurrent, vecOther)
  rotation <- rotationMatrix2D(theta)
  
  # transform points P of other face
  # as: p1current + rotation (P - p1other)
  moved <- sweep(layout[[otherFace]]$coords2D, 2, as.numeric(layout[[otherFace]]$coords2D[p1other,]))
  rotated <- t(rotation %*% t(moved))
  translated <- sweep(rotated, 2, as.numeric(layout[[currentFace]]$coords2D[p1current,]), "+")
  layout[[otherFace]]$coords2D <- translated
  drawFace2D(layout[[otherFace]])
  
  cat("Placed:", otherFace, fill = T)
}

