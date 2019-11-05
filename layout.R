# automatic layout gen

source("polyhedra.R")
source("draw.R")

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- quasi(cube)
xx <- rhombic(dodecahedron)
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


# for layout use a new device
rgl_init(new.device = T)
clear3d()
drawAxes()

#drawFace2D(layout[[currentFace]])

bestResult <- -Inf

findBestLayout <- function(p, j)
{
  # Faces in the layout 1..j
  placedFaces <- sapply(layout[1:j], function(f) {return(f$face)})
  
  # Candidate edges have one face in the placed faces and one that is not
  candidateEdges <- which((p$edgeToFaces[,1] %in% placedFaces & !(p$edgeToFaces[,2] %in% placedFaces)) |
                            (p$edgeToFaces[,2] %in% placedFaces & !(p$edgeToFaces[,1] %in% placedFaces)))
  
  # No candidates? We may be done!
  if (length(candidateEdges) == 0) {
    cat("No candidate edges at level", j, fill=T)
    return
  }
  
  cat("Level", j, "candidate edges:", candidateEdges, fill=T)
  
  # Try all of them
  for (e in candidateEdges) {
    cat("Level", j, "trying edge", e, fill=T)
    
    # Determine faces and vertices around this edge
    if (p$edgeToFaces[e,1] %in% placedFaces) {
      currentFace <- p$edgeToFaces[e,1]
      otherFace <- p$edgeToFaces[e,2]
    } else {
      otherFace <- p$edgeToFaces[e,1]
      currentFace <- p$edgeToFaces[e,2]
    }
    vertices <- ((which(p$coordPairToEdge==e)-1) %/% nrow(p$coordPairToEdge))+1
    if (length(vertices) != 2) stop("An edge must have 2 vertices")
    cat("Level", j, "placed:", currentFace, ", new face:", otherFace, fill = T)
    cat("Level", j, "vertices:", vertices, fill=T)
    
    # Project the new face to 2D
    newFace2D <- projectFace(p, otherFace) 
    drawFace2D(newFace2D, color="blue")
    
    ## TODO below
    ## layout[[j]] not good
    ## instead find which in layout[[1..j]] matches "currentFace" and use that index
    ## however with just one j doesnt matter
  
    # identify where the edge is in both faces
    p1current <- which(layout[[j]]$points3D == vertices[1])
    p1other <- which(layout[[j]]$points3D == vertices[1])
    p2current <- which(newFace2D$points3D == vertices[2])
    p2other <- which(newFace2D$points3D == vertices[2])
    if (length(c(p1current,p1other,p2current,p2other)) != 4) stop("All vertices must occur in both faces")
    
    # determine the rotation of the new face to put it next to the edge
    vecCurrent <- as.numeric(layout[[j]]$coords2D[p2current,]) - as.numeric(layout[[j]]$coords2D[p1current,])
    vecOther <- as.numeric(newFace2D$coords2D[p2other,]) - as.numeric(newFace2D$coords2D[p1other,])
    theta <-  vectorAngle2D(vecCurrent, vecOther)
    rotation <- rotationMatrix2D(theta)
    
    # do the actual rotation and translation
    # as: p1current + rotation (P - p1other)
    moved <- sweep(newFace2D$coords2D, 2, as.numeric(newFace2D$coords2D[p1other,]))
    rotated <- t(rotation %*% t(moved))
    translated <- sweep(rotated, 2, as.numeric(layout[[j]]$coords2D[p1current,]), "+")
    newFace2D$coords2D <- translated
    drawFace2D(newFace2D, color="green") # NOT GOOD, places inside by accident...
    
    # Now, evaluate the new layout with newFace2D added
    # if already worse then go to next
    
    # If we have now placed all faces then evaluate
    # and store the metric as best, keeping this layout somewhere
    # and go to next
    
    # extend layout with newFace2D
    # and call findBestLayout with j+1
  }
}

# Start with a layout with just the first face
layout <- list()

# Place first face on 2D plane - does not matter where we start
face2D <- projectFace(xx, 4) # just to put off the algo - set to 1 later again
layout[[1+length(layout)]] <- face2D

findBestLayout(xx, 1)


# repeat from here
repeat {
  # find a free edge
  # one of xx$edgeToFaces that has [,1] in layout and [,2] not (or vice versa)
  placedFaces <- sapply(layout, function(f) {return(f$face)})
  candidateEdges <- which((xx$edgeToFaces[,1] %in% placedFaces & !(xx$edgeToFaces[,2] %in% placedFaces)) |
                            (xx$edgeToFaces[,2] %in% placedFaces & !(xx$edgeToFaces[,1] %in% placedFaces)))
  if (length(candidateEdges) == 0) {
    # Done! (or perhaps need to place unconnected face)
    break
  }
  
  
  
  candidateEdge <- rev(candidateEdges)[1] # TODO iterate/recurse here and fail if there are none
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

evalLayout(layout)

