# Functions specifically to constructing polyhedra

# shift elements from f one up or down
# NB could be generalized to shift more than 1 but not needed here
shift <- function(f, direction = 1)
{
  if (length(f) < 2) return(f)
  if (direction == 1) return (f[c(2:length(f),1)])
  return(f[c(length(f),1:(length(f)-1))])
}


# Identifies connected faces in given polyhedron (i.e. polygons sharing an edge) and puts 
# the face numbers of distinct bodies into separate lists. Useful for color assigment.
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

# Build up a face from existing vertices in given polygon, within constraints passed in as
# the number of edges/vertices of this face, number of faces per vertex and the length of
# an edge. The face does not need to be regular but will have all edges of the same length.
buildFace <- function(p, polygonsize, vertexsize, edgelength, aFace = c(), debug = F)
{
  if (debug) cat("Face", aFace, fill=T)
  
  # Bail out if face has 3 vertices and those already exist in one of the polygon's faces
  if (length(aFace) == 3 & length(p$faces) > 0) {
    if (any( sapply( seq(length(p$faces)), function(f) { return (length(intersect(aFace, p$faces[[f]])))}) == 3 )) {
      if (debug) cat("Same plane as an existing face", fill=T)
      return(NULL)
    }
  }
  
  # Otherwise, we may be done
  if (length(aFace) == polygonsize) {
    
    # Final check: if (first) 3 points of a face are in aFace then not good
    if (length(p$faces) > 0) {
      # calculate normal from the first three of the face points (the rest is guaranteed to be in the same plane)
      normal <- normal(p$vertices[aFace[1],], p$vertices[aFace[2],], p$vertices[aFace[3],])
      
      for (f in p$faces) {
        verticesInNewPlane <- deltaEquals(normal %*% apply(p$vertices[f,], 1, function(x) {return(x-unlist(p$vertices[aFace[1],]))}), 0)
        if (sum(verticesInNewPlane) >= 3) {
          if (debug) {
            cat("Existing face", f, "has 3 or more points in same plane as new face", aFace, fill=T)
          }
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
  
  # Check vertex right distance to previous point
  if (length(candidates) > 0 & length(aFace) > 0) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[length(aFace)],], p$vertices[candidates,]), edgelength))]  
    if (debug) cat("Right distance to prev:", candidates, fill = T)
  }
  
  # Check last vertex right distance to first point
  if (length(candidates) > 0 & length(aFace) == (polygonsize-1)) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[1],], p$vertices[candidates,]), edgelength))]  
    if (debug) cat("Right distance to first:", candidates, fill = T)
  }
  
  # Check if the 4th and further points are in the same plane as the first 3
  if (length(candidates) > 0 & length(aFace) >= 3) {
    # calculate normal from three of the face points
    normal <- normal(p$vertices[aFace[1],], p$vertices[aFace[2],], p$vertices[aFace[3],])
    
    # then the inner product of n with the vector of any of those to a given point is 0 when in the same plane
    candidates <- candidates[which(deltaEquals(normal %*% apply(p$vertices[candidates,], 1, function(x) {return(x-unlist(p$vertices[aFace[1],]))}), 0))]
    if (debug) cat("In the same plane as the first three:", candidates, fill = T)
  }
  
  # Loop over the remaining candidates and try them out
  for (c in candidates) {
    f <- buildFace(p, polygonsize, vertexsize, edgelength, c(aFace, c), debug)
    if (!is.null(f)) { return(f) }
  }
  
  return(NULL)
}

# TODO consider building polygon from just symbol like {3,5}

# Build a polygon given a set of vertices (full x, y, z coordinates), the two vertices
# of an example edge (used to determine the global edge size of this polyhedron)
buildRegularPoly <- function(vertices, polygonsize, vertexsize, exampleEdge = c(1,2), debug=F)
{
  # Scale the vertices to unit length
  vertices <- normalizedistances(vertices)
  
  # Determine global edge size from given example edge
  edgelength <- distance(vertices[exampleEdge[1],], vertices[exampleEdge[2],])
  
  # Add new faces one by one until all vertices have the specified number of faces
  poly <- list(vertices = vertices, faces = list())
  edges <- matrix(data = 0, nrow = nrow(vertices), ncol = nrow(vertices))
  repeat {
    # Check if there's any not fully occupied vertices
    freeVertices <- which( sapply(seq(nrow(poly$vertices)), function(i) { sum(unlist(poly$faces) == i)}) < vertexsize )
    if (length(freeVertices) == 0) {
      if (debug) print("All done!")
      break
    }
    
    if (debug) cat("Building face:", length(poly$faces) + 1, fill=T)
    
    # Find existing edges in the poly that do not have two faces attached to it yet
    # pass in the first of such as the seed for the new face. This ensures that we first
    # finish bodies that are partially constructed.
    
    startEdgesX <- which(apply(edges, 1, function(nConnectedFaces) {return(any(nConnectedFaces==1))} ))
    if (length(startEdgesX) > 0) {
      e1 <- startEdgesX[1]
      e2 <- which(edges[e1,] == 1)[1] # must be one
      if (is.na(e2)) {
        stop("Not good")
      }
      if (debug) cat("Starting with edge", e1, "-", e2, fill=T)
      f <- buildFace(poly, polygonsize, vertexsize, edgelength, aFace = c(e1, e2), debug = debug)
    } else {
      f <- buildFace(poly, polygonsize, vertexsize, edgelength, debug = debug)
    }
    
    if (debug) print(f)
    
    if (!is.null(f)) {
      poly$faces[[length(poly$faces) + 1]] <- f
      for (i in seq(length(f))) {
        edges[f[i], shift(f)[i]] <- 1 + edges[f[i], shift(f)[i]]
        edges[shift(f)[i], f[i]] <- 1 + edges[shift(f)[i], f[i]]
      }
      if (debug) drawPoly(poly, debug=T)
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


# TODO: dual of some doesnt work, maybe the ones with pentagons
# the dual of this - doesnt work! strange - maybe because of {5/2} faces that occur, should be smallStellatedDodecahedron
#drawPoly(dual(greatDodecahedron)) 

# Create dual polygon from given. Vertices become faces and vice versa.
dual <- function(p)
{
  # for every face, calculate the midpoint
  newVertices <- as.data.frame(t(sapply(seq(length(p$faces)), function(f) { apply(p$vertices[p$faces[[f]],], 2, mean) })))
  # rescale to unit lenth (NB may not apply to all inverses - TODO reconsider, maybe they should
  # be scaled while building the new faces, to have the same length as the length of the current points around it.
  newVertices <- normalizedistances(newVertices)
  
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

