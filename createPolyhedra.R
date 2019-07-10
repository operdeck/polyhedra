library(rgl)

open3d()
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

cube.vertices <- expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1))
cube.d <- distance(cube.vertices[1,], cube.vertices[2, ])
cube.facesperedge <- 3
cube.polygonsize <- 4

octahedron.vertices <- rbind(expand.grid(x = c(-1,1), y = 0, z = 0), expand.grid(x = 0, y = c(-1,1), z = 0), expand.grid(x = 0, y = 0, z = c(-1,1)))
octahedron.d <- distance(octahedron.vertices[1,], octahedron.vertices[3, ])
octahedron.facesperedge <- 4
octahedron.polygonsize <- 3

#rgl.spheres(vertices$x, vertices$y, vertices$z, r=0.2, color="yellow")
#rgl.texts(vertices$x, vertices$y, vertices$z, text = seq(nrow(vertices)), color="red")

drawPoly <- function(p)
{
  rgl.clear()
  rgl.texts(p$vertices$x, p$vertices$y, p$vertices$z, text = seq(nrow(p$vertices)), color="blue")
  if (length(p$faces) > 0) {
    for (f in seq(length(p$faces))) {
      # triangulize, draw triangles between center of face and all edges
      cx = mean(p$vertices$x[p$faces[[f]]])
      cy = mean(p$vertices$y[p$faces[[f]]])
      cz = mean(p$vertices$z[p$faces[[f]]])
      for (t in seq(length(p$faces[[f]]))) {
        p1 <- p$faces[[f]][t]
        p2 <- p$faces[[f]][(t %% length(p$faces[[f]])) + 1]
        # NB not sure about the orientation of the triangle - may have to check on this
        triangles3d( c(p$vertices$x[c(p1,p2)],cx), c(p$vertices$y[c(p1,p2)],cy), c(p$vertices$z[c(p1,p2)],cz))
      }
    }
  }
}
# rgl.close()

# given a distance d and a schaffly (?) description eg {4, 3}

# to build a new regular face with E edges
# 1. choose a first vertex that is not fully taken yet (less than F faces containing this)
# 2. choose a second vertex that has distance d to the first one (and has not been taken)
# 3. choose a third vertex that has distance d, has not been taken, and builds a face {1,2,3} with a normal pointing outwards
# (if we need more)
# 4. for any next vertex, choose so that distance fullfilled, hasnt been taken, and is in the same plane as 1-2-3

# start with initial just for testing
# poly$faces[[1]] <- c(1, 2, 4, 3)
# poly$faces[[2]] <- c(1, 3, 7, 5)
# poly$faces[[3]] <- c(3, 4, 8, 7)

# the last one also needs to have desired distance to the first vertex so at this point
# v1 should be in the 2nd half of above expression, if not --> fail

# any further points need to be in the same plane
# plane can be formed after 3 pts https://mathinsight.org/forming_planes
# dist to plane: https://mathinsight.org/distance_point_plane

buildFace <- function(p, aFace = c())
{
  print(aFace)

  # If we're done, verify that the face is not already in the polygon
  if (length(aFace) == p$polygonsize) {
    if (length(p$faces) > 0) {
      c1 <- sort(aFace)
      for (f in seq(length(p$faces))) {
        if (identical(c1, sort(p$faces[[f]]))) {
          cat("Duplicate", fill=T)
          return(NULL)
        }
      }
    }
    cat("We're OK", fill=T)
    # consider turning it upside down here to normalize the normal
    return(aFace)
  }
  
  # Check vertex not fully occupied
  candidates <- setdiff(which( sapply(seq(nrow(p$vertices)), function(i) { sum(unlist(p$faces) == i)}) < p$facesperedge ), aFace)
  cat("Free vertices:", candidates, fill = T)
  
  if (length(candidates) == 0) { return(FALSE) }
  
  # Check vertex right distance to previous point
  if (length(aFace) > 0) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[length(aFace)],], p$vertices[candidates,]), p$edgesize))]  
    cat("Right distance to prev:", candidates, fill = T)
  }
  
  # Check last vertex right distance to first point
  if (length(aFace) == (p$polygonsize-1)) {
    candidates <- candidates[which(deltaEquals(distance(p$vertices[aFace[1],], p$vertices[candidates,]), p$edgesize))]  
    cat("Right distance to first:", candidates, fill = T)
  }
    
  # for polysize > 3 check they're all in the same plane
  
  for (c in candidates) {
    f <- buildFace(p, c(aFace, c))
    if (!is.null(f)) { return(f) }
  }
  
  return(NULL)
}

poly <- list(vertices = octahedron.vertices,
             faces = list(),
             edgesize = octahedron.d,
             facesperedge = octahedron.facesperedge,
             polygonsize = octahedron.polygonsize)

repeat {
  # Check if there's any not fully occupied vertices
  freeVertices <- which( sapply(seq(nrow(poly$vertices)), function(i) { sum(unlist(poly$faces) == i)}) < poly$facesperedge )
  if (length(freeVertices) == 0) {
    print("All done!")
    cat("The polyhedron has", length(poly$faces), "faces,", length(poly$faces)*poly$polygonsize/2, "edges,", nrow(poly$vertices), "vertices", fill = T)
    break
  }

  cat("Building face:", length(poly$faces) + 1, fill=T)
  
  f <- buildFace(poly)
  print(f)
  
  if (!is.null(f)) {
    poly$faces[[length(poly$faces) + 1]] <- f
    drawPoly(poly)
  } else {
    print("Can't construct polygon!")
    break
  }
} 


