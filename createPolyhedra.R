library(rgl)

#open3d()
#shade3d( icosahedron3d() )
#writeOBJ(filename)

distance <- function(a, b)
{
  return (sqrt((a$x-b$x)^2 + (a$y-b$y)^2 + (a$z-b$z)^2 ))
}

deltaEquals <- function(x, y, delta = 1e-6)
{
  return(abs(x-y)<delta)
}

#rgl.spheres(vertices$x, vertices$y, vertices$z, r=0.2, color="yellow")
#rgl.texts(vertices$x, vertices$y, vertices$z, text = seq(nrow(vertices)), color="red")

# Identifies connected faces in given polyhedron (sharing an edge) and puts the face numbers
# of distinct bodies into separate lists. Useful for color assigment.
findDistinctBodies <- function(p)
{
  bodies <- list()
  bodyIndex <- 1
  bodies[[bodyIndex]] <- c(1) # start with face 1
  repeat 
  {
    i <- 1 # index in list of faces in current body, keeps incrementing while new faces may be added at end of this list
    repeat
    {
      unassigned <- setdiff(1:length(p$faces), unlist(bodies))
      if (length(unassigned) == 0) break # done, all faces assigned to a body
      # find faces connect to the face at position i
      connectedfaces <- unassigned[sapply(unassigned, function(f) {return(length(intersect(p$faces[[f]], 
                                                                                           p$faces[[bodies[[bodyIndex]][i]]])) > 1)})]
      if (length(connectedfaces) > 0) {
        bodies[[bodyIndex]] <- c(bodies[[bodyIndex]], connectedfaces)
      }
      i <- i+1
      if (i > length(bodies[[bodyIndex]])) break # done with all faces in the current body
    }
    if (length(unassigned) == 0) break # done, all faces assigned to a body
    
    bodyIndex <- bodyIndex + 1
    bodies[[bodyIndex]] <- unassigned[1]
  }
  return(bodies)
}

# Draw a polygon. Offset is optional.
drawPoly <- function(p, x=0, y=0, z=0, label="")
{
  rgl.texts(x + p$vertices$x, y + p$vertices$y, z + p$vertices$z, text = seq(nrow(p$vertices)), color="blue")
  if (nchar(label) > 0) {
    rgl.texts(x, y, z + min(p$vertices$z) - 1, text = label, color = "black")
  }
  if (length(p$faces) > 0) {
    bodies <- findDistinctBodies(p)
    if (length(bodies) > 1) {
      bodyColors <- rainbow(length(bodies))
    }
    for (f in seq(length(p$faces))) {
      if (length(bodies) > 1) {
        faceColor <- bodyColors[which(sapply(bodies, function(b) { return(f %in% b)}))]
      } else {
        faceColor <- rainbow(length(p$faces))[f]
      }
      if (length(p$faces[[f]]) > 3) { 
        # for > 3 vertices triangulize: draw triangles between center of face and all edges
        cx = mean(p$vertices$x[p$faces[[f]]])
        cy = mean(p$vertices$y[p$faces[[f]]])
        cz = mean(p$vertices$z[p$faces[[f]]])
        for (t in seq(length(p$faces[[f]]))) {
          p1 <- p$faces[[f]][t]
          p2 <- p$faces[[f]][(t %% length(p$faces[[f]])) + 1]
          # NB not sure about the orientation of the triangle - may have to check on this
          triangles3d( x + c(p$vertices$x[c(p1,p2)],cx), 
                       y + c(p$vertices$y[c(p1,p2)],cy), 
                       z + c(p$vertices$z[c(p1,p2)],cz), 
                       col=faceColor) # "red"
        }
      } else {
        # NB not sure about the orientation of the triangle - may have to check on this
        triangles3d( x + p$vertices$x[p$faces[[f]]],
                     y + p$vertices$y[p$faces[[f]]], 
                     z + p$vertices$z[p$faces[[f]]], 
                     col=faceColor)
      }
    }
  }
}
# rgl.close()

# given a distance d and a schaffly (?) description eg {4, 3}





# Build up a face from existing vertices in given polygon, within constraints passed in as
# the number of edges/vertices of this face, number of faces per vertex and the length of
# an edge. The face does not need to be regular but will have all edges of the same length.
buildFace <- function(p, polygonsize, vertexsize, edgelength, aFace = c(), debug = F)
{
  if (debug) print(aFace)
  
  # If we're done, verify that the face is not already in the polygon
  if (length(aFace) == polygonsize) {
    if (length(p$faces) > 0) {
      c1 <- sort(aFace)
      for (f in seq(length(p$faces))) {
        if (identical(c1, sort(p$faces[[f]]))) {
          if (debug) cat("Duplicate", fill=T)
          return(NULL)
        }
      }
    }
    if (debug) cat("Face complete", fill=T)
    # consider turning it upside down here to normalize the normal
    return(aFace)
  }
  
  # Check vertex not fully occupied
  candidates <- setdiff(which( sapply(seq(nrow(p$vertices)), function(i) { sum(unlist(p$faces) == i)}) < vertexsize ), aFace)
  if (debug) cat("Free vertices:", candidates, fill = T)
  
  if (length(candidates) == 0) { return(FALSE) }
  
  # Check vertex right distance to previous point
  if (length(aFace) > 0) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[length(aFace)],], p$vertices[candidates,]), edgelength))]  
    if (debug) cat("Right distance to prev:", candidates, fill = T)
  }
  
  # Check last vertex right distance to first point
  if (length(aFace) == (polygonsize-1)) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[1],], p$vertices[candidates,]), edgelength))]  
    if (debug) cat("Right distance to first:", candidates, fill = T)
  }
  
  # TODO:  
  # for polysize > 3 check they're all in the same plane
  # any further points need to be in the same plane
  # plane can be formed after 3 pts https://mathinsight.org/forming_planes
  # dist to plane: https://mathinsight.org/distance_point_plane
  
  # Loop over the remaining candidates and try them out
  for (c in candidates) {
    f <- buildFace(p, polygonsize, vertexsize, edgelength, c(aFace, c), debug)
    if (!is.null(f)) { return(f) }
  }
  
  return(NULL)
}

# Build a polygon given a set of vertices (full x, y, z coordinates), the two vertices
# of an example edge (used to determine the global edge size of this polyhedron)
buildRegularPoly <- function(vertices, polygonsize, vertexsize, exampleEdge = c(1,2), debug=F)
{
  # Scale the vertices to unit length
  vertices <- vertices / apply(vertices, 1, function(r) { sqrt(sum(r^2)) })
  
  # Determine global edge size from given example edge
  edgelength <- distance(vertices[exampleEdge[1],], vertices[exampleEdge[2],])
  
  # Add new faces one by one until all vertices have the specified number of faces
  poly <- list(vertices = vertices, faces = list())
  repeat {
    # Check if there's any not fully occupied vertices
    freeVertices <- which( sapply(seq(nrow(poly$vertices)), function(i) { sum(unlist(poly$faces) == i)}) < vertexsize )
    if (length(freeVertices) == 0) {
      if (debug) print("All done!")
      break
    }
    
    if (debug) cat("Building face:", length(poly$faces) + 1, fill=T)
    
    f <- buildFace(poly, polygonsize, vertexsize, edgelength, debug = debug)
    if (debug) print(f)
    
    if (!is.null(f)) {
      poly$faces[[length(poly$faces) + 1]] <- f
      if (debug) drawPoly(poly)
    } else {
      if (debug) print("Can't construct polygon!")
      poly$faces <- list()
      break
    }
  } 
  
  poly[["n_faces"]] <- length(poly$faces)
  poly[["n_edges"]] <- length(poly$faces) + nrow(poly$vertices) - 2 # Euler's formula
  poly[["n_vertices"]] <- nrow(poly$vertices)
  
  return(poly)
}


# Create dual polygon from given. Vertices become faces and vice versa.
dual <- function(p)
{
  # for every face, calculate the midpoint
  newVertices <- as.data.frame(t(sapply(seq(length(p$faces)), function(f) { apply(p$vertices[p$faces[[f]],], 2, mean) })))
  # rescale to unit lenth (NB may not apply to all inverses - TODO reconsider)
  newVertices <- newVertices / apply(newVertices, 1, function(r) { sqrt(sum(r^2)) })
  # vertices become the new faces
  newFaces <- list()
  for (v in seq(nrow(p$vertices))) {
    # which faces contain vertex v  
    connectedFaces <- which(sapply(p$faces, function(f) { return (v %in% f) }))
    # but they can be in an arbitrary order
    f1 <- connectedFaces[1]
    
    aNewFace <- c()
    neighbourVertices <- c()
    for (i in seq(length(connectedFaces))) {
      # an edge from v along f1 to next point w (so the edge v-w is an edge of f1) - there should be two depending on direction chosen,
      # then choose the one not already used
      w_up <- p$faces[[f1]] [(which(p$faces[[f1]] == v) %% 3) + 1]
      w_down <- p$faces[[f1]] [((which(p$faces[[f1]] == v) + 1 + length(p$faces[[f1]])) %% 3) + 1]
      w <- ifelse(w_up %in% neighbourVertices, w_down, w_up)
      # next face is the other face (there can be only 1) that contains v-w or w-v
      f2 <- setdiff(which(sapply(p$faces, function(f) { return ((w %in% f) & (v %in% f)) })), f1)
      aNewFace <- c(aNewFace, f1)
      neighbourVertices <- c(neighbourVertices, w)
      f1 <- f2
    }
    newFaces[[v]] <- aNewFace
  }
  poly <- list(vertices = newVertices, faces = newFaces)
  
  poly[["n_faces"]] <- length(poly$faces)
  poly[["n_edges"]] <- length(poly$faces) + nrow(poly$vertices) - 2 # Euler's formula
  poly[["n_vertices"]] <- nrow(poly$vertices)
  
  return(poly)
}

# Combine two polyhedra into one
compose <- function(p1, p2)
{
  return(list(vertices = rbind(p1$vertices, p2$vertices),
              faces = c(p1$faces, lapply(p2$faces, function(f) { return(f+nrow(p1$vertices)) })),
              n_faces = p1$n_faces + p2$n_faces,
              n_edges = p1$n_edges + p2$n_edges,
              n_vertices = p1$n_vertices + p2$n_vertices))
}

## Gallery

open3d()
rgl.clear()
par3d("cex" = 0.7)

tetrahedron <- buildRegularPoly(vertices = rbind(data.frame(x=1, y=1, z=1), data.frame(x=1, y=-1, z=-1), data.frame(x=-1, y=1, z=-1), data.frame(x=-1, y=-1, z=1)),
                                polygonsize = 3,
                                vertexsize = 3,
                                debug = T)
drawPoly(tetrahedron, label="Tetrahedron")

tetrahedron2 <- dual(tetrahedron)
drawPoly(tetrahedron2, z = -3, label="Dual Tetrahedron")

drawPoly(compose(tetrahedron, dual(tetrahedron)), z = -6)

# seems to work despite not checking faces
cubedirect <- buildRegularPoly(vertices = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
                               polygonsize = 4,
                               vertexsize = 3)

octahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = c(-1,1), y = 0, z = 0), expand.grid(x = 0, y = c(-1,1), z = 0), expand.grid(x = 0, y = 0, z = c(-1,1))),
                               polygonsize = 3,
                               vertexsize = 4,
                               exampleEdge = c(1,3))
drawPoly(octahedron, x = 3, label="Octahedron")

cube <- dual(octahedron)
drawPoly(cube, x = 3, z = -3, label="Cube")

cubeWithDual <- compose(cube, octahedron)
drawPoly(cubeWithDual, x = 3, z = -6)

# all coords taken from https://en.wikipedia.org/wiki/Platonic_solid
phi <- (1+sqrt(5))/2

# dodecahedron seems to work even lacking same-plane condition
# dodecahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = c(-1,1), y = c(-1,1), z = c(-1,1)), 
#                                                   expand.grid(x = 0, y = c(-1/phi,1/phi), z = c(-phi, phi)), 
#                                                   expand.grid(x = c(-1/phi,1/phi), y = c(-phi, phi), z = 0),
#                                                   expand.grid(x = c(-phi, phi), y = 0, z = c(-1/phi,1/phi))),
#                                  polygonsize = 5,
#                                  vertexsize = 3,
#                                  exampleEdge = c(1,9),
#                                  debug = T)

icosahedron <- buildRegularPoly(vertices = rbind(expand.grid(x = 0, y = c(-1,1), z = c(-phi, phi)), 
                                                 expand.grid(x = c(-1,1), y = c(-phi, phi), z = 0), 
                                                 expand.grid(x = c(-phi, phi), y = 0, z = c(-1, 1))),
                                polygonsize = 3,
                                vertexsize = 5,
                                debug = F)


drawPoly(icosahedron, x = 6, label="Icosahedron")

dodecahedron <- dual(icosahedron)
drawPoly(dodecahedron, x = 6, z = -3, label="Dodecahedron")

drawPoly(compose(icosahedron, dual(icosahedron)), x = 6, z = -6)


greatIcosahedron <- buildRegularPoly(icosahedron$vertices,
                                     polygonsize = 3,
                                     vertexsize = 5,
                                     exampleEdge = c(2, 6))
drawPoly(greatIcosahedron, x = 9, label="Great Icosahedron")

greatStellatedDodecahedron <- dual(greatIcosahedron)
drawPoly(greatStellatedDodecahedron, x = 9, z = -3, label="Great Stellated Dodecahedron")

drawPoly(compose(greatStellatedDodecahedron, greatIcosahedron), x = 9, z = -6)

# NB below one not working yet
stop()

# this one doesnt work because of violation of face constraint
greatDodecahedron <- buildRegularPoly(vertices = icosahedron$vertices, 
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6), debug=T)

compound5tetrahedra <- buildRegularPoly(dodecahedron$vertices,
                                        polygonsize = 3,
                                        vertexsize = 3,
                                        exampleEdge = c(8, 15), debug = T)
# rgl.close()

# colouring
# if there multiple bodies -> each its own color
# else if there are multiple types of faces -> each its own color
# else rainbow


crossProduct <- function(ab, ac){
  abci = ab$y * ac$z - ac$y * ab$z
  abcj = ac$x * ab$z - ab$x * ac$z
  abck = ab$x * ac$y - ac$x * ab$y
  return (c(x=abci, y=abcj, z=abck))
}

# create cross product from the two vectors from three of the points to get the normal to the plane described by these 3 points
n <- crossProduct(cubedirect$vertices[3,]-cubedirect$vertices[1,], cubedirect$vertices[5,]-cubedirect$vertices[1,])

# then the inner product of n with the vector of any of those to a given point is 0 when in the same plane
n %*% unlist(cubedirect$vertices[7,]-cubedirect$vertices[1,]) # 7 is in the same plane
n %*% unlist(cubedirect$vertices[6,]-cubedirect$vertices[1,]) # but this one is not
n %*% unlist(cubedirect$vertices[8,]-cubedirect$vertices[1,]) # but this one is not

# or apply directly to a set of vertices that should be tested
n %*% apply(cubedirect$vertices[1:8,], 1, function(x) {return(x-unlist(cubedirect$vertices[1,]))})




