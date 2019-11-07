# automatic layout gen

source("polyhedra.R")
source("draw.R")

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- quasi(cube)
xx <- rhombic(dodecahedron)
xx <- tetrahedron
clear3d()
drawPoly(xx, debug=T)

# draw 2D face
drawFace2D <- function(face2d, ...)
{
  #points3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, color="black" )
  text3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, 
         text=face2d$points3D, color="black", cex=0.5 )
  lines3d(x = face2d$coords2D[,1][c(seq(nrow(face2d$coords2D)),1)], 
          y = face2d$coords2D[,2][c(seq(nrow(face2d$coords2D)),1)], z = 0, ...)
  center <- apply(face2d$coords2D, 2, mean)
  text3d(x = center[1], y = center[2], z = 0, text=paste0("F",face2d$face), color="blue", cex=0.5 )
}

# evaluate size of layout
evalLayout  <- function(l)
{
  maxs <- sapply(l, function(f) { return(apply(f$coords2D,2,max))})
  mins <- sapply(l, function(f) { return(apply(f$coords2D,2,min))})
  max_x <- max(maxs[1,])
  max_y <- max(maxs[2,])
  min_x <- min(mins[1,])
  min_y <- min(mins[2,])
  layoutsize <- (max_x-min_x)*(max_y-min_y)
  return(layoutsize)
}

# project 3D face onto 2D
projectFace <- function(poly, faceno)
{
  coords3D <- poly$coords[poly$faces[[faceno]],]
  angles <- innerAngles(coords3D)
  radii <- distance(coords3D)
  return(list( face = faceno, 
               points3D = xx$faces[[faceno]],
               coords2D = matrix(c(cos(cumsum(angles)) * radii, 
                                   sin(cumsum(angles)) * radii),
                                 nrow = nrow(coords3D))))
}

# put a new face into position, given a global 2D layout
positionNextFace <- function(polyhedron3D, edge, placedFaces) 
{
  if (polyhedron3D$edgeToFaces[edge,1] %in% placedFaces) {
    placedFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
    newFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
  } else {
    newFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
    placedFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
  }
  vertices <- ((which(polyhedron3D$coordPairToEdge==edge)-1) %/% nrow(polyhedron3D$coordPairToEdge))+1
  cat("Placed:", placedFaceIdx, ", adding:", newFaceIdx, fill = T)
  cat("vertices:", paste(vertices, collapse=","), fill=T)
  
  # place the to be added face on 2D plane
  face2D <- projectFace(polyhedron3D, newFaceIdx)
  drawFace2D(face2D, color="green")
  
  # lookup the one that is placed already
  currentFace2D <- layout[[which(sapply(layout, function(x){return(x$face)})==placedFaceIdx)]]
  
  # identify the vertices of the edge in both faces
  if (length(vertices) != 2) stop("An edge must have 2 vertices")
  p1current <- which(currentFace2D$points3D == vertices[1])
  p1other <- which(face2D$points3D == vertices[1])
  p2current <- which(currentFace2D$points3D == vertices[2])
  p2other <- which(face2D$points3D == vertices[2])
  if (length(c(p1current,p1other,p2current,p2other)) != 4) stop("All vertices must occur in both faces")
  
  # rotate and translate the connected face
  vecCurrent <- as.numeric(currentFace2D$coords2D[p2current,]) - as.numeric(currentFace2D$coords2D[p1current,])
  vecOther <- as.numeric(face2D$coords2D[p2other,]) - as.numeric(face2D$coords2D[p1other,])
  theta <-  vectorAngle2D(vecCurrent, vecOther)
  rotation <- rotationMatrix2D(theta)
  
  # transform points P of other face
  # as: p1current + rotation (P - p1other)
  moved <- sweep(face2D$coords2D, 2, as.numeric(face2D$coords2D[p1other,]))
  rotated <- t(rotation %*% t(moved))
  translated <- sweep(rotated, 2, as.numeric(currentFace2D$coords2D[p1current,]), "+")
  face2D$coords2D <- translated
  
  return(face2D)
}

bestResult <- -Inf
stop()

##### todo refactor until clean then make recursive

# for layout use a new device
rgl_init(new.device = T)
clear3d()
drawAxes()

# Start with a layout with just the first face
layout <- list()

# Place first face on 2D plane - does not matter where we start
face2D <- projectFace(xx, 1) # just to put off the algo - set to 1 later again


# repeat from here
repeat {
  
  # place face2D at this level which will become another function argument
  layout[[1+length(layout)]] <- face2D
  drawFace2D(face2D, color="red")
  level <- length(layout) # prep for search version, this will become a function argument (?)
  
  cat("Level", level, "placed face:", face2D$face, fill = T)
  
  # update the list of available edges, and store it with the layout itself
  shiftFacePoints <- shiftrotate(face2D$points3D)
  newEdges <- sapply(seq(length(face2D$points3D)), function(i) {return(xx$coordPairToEdge[face2D$points3D[i],shiftFacePoints[i]])})
  if (level > 1) {
    layout[[level]]$candidateEdges <-
      setdiff(c(newEdges, layout[[level-1]]$candidateEdges), intersect(layout[[level-1]]$candidateEdges, newEdges))
  } else {
    layout[[level]]$candidateEdges <- newEdges
  }
  
  # keep track of the layout min/max
  newMinCoords <- apply(face2D$coords2D,2,min)
  newMaxCoords <- apply(face2D$coords2D,2,max)
  if (level > 1) {
    layout[[level]]$minCoords <- apply(matrix(c(newMinCoords, layout[[level-1]]$minCoords), ncol=2, byrow = T), 2, min)
    layout[[level]]$maxCoords <- apply(matrix(c(newMaxCoords, layout[[level-1]]$maxCoords), ncol=2, byrow = T), 2, max)
  } else {
    layout[[level]]$minCoords <- newMinCoords
    layout[[level]]$maxCoords <- newMaxCoords
  }
  lines3d( x = c(layout[[level]]$minCoords[1], layout[[level]]$maxCoords[1], layout[[level]]$maxCoords[1], layout[[level]]$minCoords[1], layout[[level]]$minCoords[1]),
           y = c(layout[[level]]$minCoords[2], layout[[level]]$minCoords[2], layout[[level]]$maxCoords[2], layout[[level]]$maxCoords[2], layout[[level]]$minCoords[2]),
           z = 0,
           color="yellow")
  
  # get list of all faces placed
  placedFaces <- sapply(layout[1:level], function(f) {return(f$face)})
  
  # check if we're done
  if (length(layout[[level]]$candidateEdges) == 0) {
    # Done! (or perhaps need to place unconnected face)
    # Evaluate it here
    break
  }
  
  # this will become a loop over all candidates
  edge <- rev(layout[[level]]$candidateEdges)[1] # TODO iterate/recurse here and fail if there are none
  
  # place it into position, given the current layout (this involves a lot of linear algebra)
  face2D <- positionNextFace(xx, edge, placedFaces) 
  
  #
  # now recurse, calling yourself with (newFace2D)
  # 
}

evalLayout(layout)

