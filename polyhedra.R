# Functions specifically to constructing polyhedra

source("math.R")

# shift elements from f one up or down
# NB could be generalized to shift more than 1 but not needed here
shift <- function(f, direction = 1)
{
  if (length(f) < 2) return(f)
  if (direction == 1) return (f[c(2:length(f),1)])
  return(f[c(length(f),1:(length(f)-1))])
}

# Returns T if outward facing (= rotating anti-clockwise when looking down towards face in direction of origin)
isNormalOutwardFacing <- function(p, f)
{
  n <- normal( p$vertices[f[1],], p$vertices[f[2],], p$vertices[f[3],] )
  mid <- apply(p$vertices[f,],2,mean)
  return (vectorlength(mid+n) > vectorlength(mid-n))
}

# # Identifies connected faces in given polyhedron (i.e. polygons sharing an edge) and puts 
# # the face numbers of distinct bodies into separate lists. Useful for color assigment.
# findDistinctBodies <- function(p)
# {
#   bodies <- list()
#   bodyIndex <- 1
#   bodies[[bodyIndex]] <- c(1) # start with face 1
#   repeat 
#   {
#     i <- 1 # index in list of faces in current body, keeps incrementing while new faces may be added at end of this list
#     repeat
#     {
#       unassigned <- setdiff(1:length(p$faces), unlist(bodies))
#       if (length(unassigned) == 0) break # done, all faces assigned to a body
#       # find faces connect to the face at position i
#       connectedfaces <- unassigned[sapply(unassigned, function(f) {return(length(intersect(p$faces[[f]], 
#                                                                                            p$faces[[bodies[[bodyIndex]][i]]])) > 1)})]
#       if (length(connectedfaces) > 0) {
#         bodies[[bodyIndex]] <- c(bodies[[bodyIndex]], connectedfaces)
#       }
#       i <- i+1
#       if (i > length(bodies[[bodyIndex]])) break # done with all faces in the current body
#     }
#     if (length(unassigned) == 0) break # done, all faces assigned to a body
#     
#     bodyIndex <- bodyIndex + 1
#     bodies[[bodyIndex]] <- unassigned[1]
#   }
#   return(bodies)
# }

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
    
    # We're done, but just make sure to have the normal face outward
    if (!isNormalOutwardFacing(p, aFace)) {
      if (debug) cat("Flip face to make normal outward facing", fill=T)
      aFace <- rev(aFace)
    }
    
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
buildRegularPoly <- function(vertices, polygonsize, vertexsize, exampleEdge = c(1,2), name = "", debug=F)
{
  # Scale the vertices to unit length
  vertices <- normalizedistances(vertices)
  
  # Determine global edge size from given example edge
  edgelength <- distance(vertices[exampleEdge[1],], vertices[exampleEdge[2],])
  
  # Add new faces one by one until all vertices have the specified number of faces
  poly <- list(vertices = vertices, faces = list(), name = name)
  edges <- matrix(data = 0, nrow = nrow(vertices), ncol = nrow(vertices))
  repeat {
    # Check if there's any not fully occupied vertices
    freeVertices <- which( sapply(seq(nrow(poly$vertices)), function(i) { sum(unlist(poly$faces) == i)}) < vertexsize )
    if (length(freeVertices) == 0) {
      if (any(edges == 1)) {
        if (debug) print("Done but with gaps")  
        poly$faces <- list()
      } else {
        if (debug) print("All done!")
      }
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
  
  return(poly)
}

# Derives the topology of a polyhedron based on just a list of faces. Each face consists of a list
# of references to vertices as everywhere else. The returned structure contains
#
# - a square matrix vexToRFace mapping vertex pairs to the index of the face on the right-hand side of it
# - a square and symetrical matrix vexToEdge mapping vertex pairs to the number (index) of the edge connecting both
# - an unordered list vexConnections with for every vertex of the polygon (NB in compounds there could be multiple of these for a vertex):
#       - the index of the center vertex
#       - the indices of the faces surrounding the vertex (in order)
#       - the indices of the vertices connected to this vertex (in order)
getTopology <- function(faces, debug=F)
{
  n_vertices <- max(sapply(faces, max))
  foundAnyGaps <- F
  
  getFaceForDebug <- function(faceidx)
  {
    if (faceidx < 1 | faceidx > length(faces)) {
      return (paste0("F",faceidx,"**invalid**"))
    }
    return (paste(paste0("F",faceidx,":"),paste(faces[[faceidx]],collapse=";")))
  }
  vexToEdge <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  vexToRFace <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  for (currentFaceNr in seq(length(faces))) {
    currentFace <- faces[[currentFaceNr]]
    currentFaceS <- shift(currentFace)
    for (i in seq(length(currentFace))) {
      # cat("Edge", currentFace[i], "-", currentFaceS[i], "has face", currentFaceNr, "on the right", fill = T)
      
      # vexToRFace for a pair [a,b] gives the face number on the right-hand side of edge connection a-b (when looked from above)
      
      if (0 == vexToRFace[currentFace[i], currentFaceS[i]]) {
        # normal situation
        vexToRFace[currentFace[i], currentFaceS[i]] <- currentFaceNr
      } else {
        if (0 == vexToRFace[currentFaceS[i], currentFace[i]]) {
          # face is turned the other way around - can happen when below origin
          vexToRFace[currentFaceS[i], currentFace[i]] <- currentFaceNr
          # stop(paste("Right-hand face for edge", currentFace[i], "-", currentFaceS[i], 
          #            "already present:", vexToRFace[currentFace[i], currentFaceS[i]], 
          #            "while wanting to set to", currentFaceNr,
          #            "is the orientation consistent?"))
        } else {
          # edge is fully occupied, this is an error
          stop(paste("Both edges", currentFace[i], "-", currentFaceS[i], 
                     "already present:", vexToRFace[currentFace[i], currentFaceS[i]], 
                     "and:", vexToRFace[currentFaceS[i], currentFace[i]], 
                     "while wanting to set to", currentFaceNr))
        }
      }
      
      # vexToEdge for a pair [a,b] gives the index number of the edge a-b. vexToEdge[a,b] = vexToEdge[b,a]
      if (0 == vexToEdge[currentFace[i], currentFaceS[i]]) {
        vexToEdge[currentFaceS[i], currentFace[i]] <- 1+max(vexToEdge)
        vexToEdge[currentFace[i], currentFaceS[i]] <- vexToEdge[currentFaceS[i], currentFace[i]]
      }
    }
  }
  
  # Knowing which edges connect to which faces we can now build up oriented lists of 
  # faces and eges around any vertex
  
  isConnectedToPreviousStartFace <- rep(F, length(faces))
  isUsedAsStartFace <- rep(F, length(faces))
  vexConnections <- list()
  startingFace <- NULL
  repeat {
    if (is.null(startingFace)) {
      # start just somewhere
      startingFace <- 1
    } else {
      if (any(isConnectedToPreviousStartFace & !isUsedAsStartFace)) {
        # prefer to continue on a face connected to any of the faces previously used
        startingFace <- which(isConnectedToPreviousStartFace & !isUsedAsStartFace)[1]
      } else {
        # otherwise there may be a disjunct set of faces (multiple bodies) so start there
        if (any(!isUsedAsStartFace)) {
          startingFace <- which(!isUsedAsStartFace)[1]
        } else {
          # all faces used, we're done
          break;
        }
      }
    }
    isUsedAsStartFace[startingFace] <- T

  # for (startingFace in seq(length(faces))) 
    if (debug) cat("start face:", getFaceForDebug(startingFace), fill=T)
    
    # check which vex we already have with this face in it
    coveredPoints <- c()
    if (length(vexConnections) > 0) {
      vexCoveringThisStartingFace <- which(sapply(vexConnections, function(con) { return(startingFace %in% con$faces)}))
      if (length(vexCoveringThisStartingFace) > 0) {
        coveredPoints <- sapply(vexCoveringThisStartingFace, function(i) {return(vexConnections[[i]]$center)})
        if (debug) cat("skipping already covered points for this face:", paste0(coveredPoints, collapse=";"), fill=T)
      }
    }
    
    for (centerPoint in setdiff(faces[[startingFace]], coveredPoints))
    {
      if (debug) cat("start building vertex around", centerPoint, fill=T)
      vexToFaces <- c()
      vexToVertex <- c()
      f <- startingFace
      
      repeat
      {
        vexToFaces <- c(vexToFaces, f) # keep track of an oriented list of connected faces around centerPoint
        isConnectedToPreviousStartFace[f] <- T
        if (debug) cat("center",centerPoint, "face",getFaceForDebug(f), fill=T)
        i <- which(faces[[f]] == centerPoint) # index in face f of centerPoint
        prevPoint <- shift(faces[[f]],-1)[i] # point on an edge of face f that connects to centerPoint
        vexToVertex <- c(vexToVertex, prevPoint) # keep track of an oriented list of connected points around centerPoint
        
        nextPoint <- shift(faces[[f]])[i] # find next face f that connects at edge nextPoint - centerPoint
        f <- vexToRFace[ nextPoint, centerPoint ]
        if (f == 0) {
          foundAnyGaps <- T
          cat("WARNING: there's a gap, no face to the right of the edge from", centerPoint, "to", nextPoint, fill=T)
          # don't add this vex although in theory we could try adding vertices with gaps by first trying going
          # the other way from the start - becomes messy and what if there are multiple gaps
          break
        }

        if (f == startingFace | f == 0) {
          # check if we don't have this already - not just checking the same center point but also
          # the set of faces as we could have multiple bodies just sharing a vertex but no faces
          if (!any(sapply(vexConnections, function(con) { return( centerPoint == con$center &
                                                                  setequal(vexToFaces, con$faces) &
                                                                  setequal(vexToVertex, con$vex))})))
          {
            if (debug) cat("adding vex center",centerPoint,"faces:",paste(paste0("F",vexToFaces), collapse = ";"),"vex:",paste(vexToVertex, collapse=";"), fill=T)
            vexConnections[[1+length(vexConnections)]] <- list(center = centerPoint,
                                                               faces = vexToFaces,
                                                               vex = vexToVertex)
          } else {
            # not sure this can even happen with the more careful choice of starting vertices in place now
            if (debug) cat("skipping adding vex center",centerPoint,"faces:",paste(paste0("F",vexToFaces), collapse = ";"),"vex:",paste(vexToVertex, collapse=";"), fill=T)
          }
          break # next vertex
        }
      }    
    }
  }
  
  return (list(vexToEdge=vexToEdge, vexToRFace=vexToRFace, vexConnections=vexConnections, hasGaps=foundAnyGaps))
}

# Return list of list of distinct bodies in p. List contains the face indices.
findDistinctBodies <- function(p, topo = getTopology(p$faces), debug=F)
{
  getConnectedFaces <- function(conn, body, debug=F)
  {
    bodycount <- 1
    repeat {
      f <- body[bodycount]
      # vertices around this face
      vex <- which(sapply(conn, function(c) { return(f %in% c$faces)}))
      # faces connected to any of these vertices
      connectedFaces <- unique(as.vector(sapply(vex, function(x) {return(conn[[x]]$faces)})))
      if (debug) cat("connected faces to", f, ":", paste(connectedFaces, collapse = ", "), fill=T)
      newFaces <- setdiff(connectedFaces, body)
      if (debug) cat("new faces:", paste(newFaces, collapse = ", "), fill=T)
      #if (length(newFaces) == 0) break
      body <- c(body, newFaces)
      bodycount <- bodycount + 1
      if (bodycount > length(body)) break
    }
    return(body)
  }
  
  bodies <- list()
  for (faceIdx in seq(length(p$faces))) {
    b <- getConnectedFaces(topo$vexConnections, faceIdx)
    if (!any(sapply(bodies, function(body) {return(setequal(body, b))}))) {
      bodies[[1+length(bodies)]] <- b
    }
  }
  return(bodies)
}

dual <- function(p, name=paste("dual", p$name), scaling = "edge", debug=F)
{
  topo <- getTopology(p$faces)
  
  newVertexCoords <- as.data.frame(t(sapply(p$faces, function(f) { return(apply(p$vertices[f,],2,mean))})))
  newFaces <- lapply(topo$vexConnections, function(c) {return(c$faces)})
  
  # scale so that mid of a new vertex is at same distance from origin as mid
  # of an old vertex - works at least for regulars, otherwise somewhat arbitrary
  if (scaling == "edge") {
    scale <- vectorlength(apply(p$vertices[p$faces[[1]][1:2],],2,mean)) / 
      vectorlength(apply(newVertexCoords[newFaces[[1]][1:2],],2,mean))
  } else if (scaling == "vertex") {
    scale <- vectorlength(p$vertices[p$faces[[1]][1],])/
      vectorlength(newVertexCoords[newFaces[[1]][1],])
  } else if (scaling == "none") {
    scale <- 1
  } else {
    stop(paste("Wrong scaling argument:", scaling))
  }
  
  pDual <- list( vertices = newVertexCoords*scale, faces = newFaces, name = name)
  
  # TODO it is possible we end up with uneven faces e.g.
  # dual(archi(dodecahedron)) --> break the face up into triangles?
  
  # make sure all faces are oriented consistently
  for (i in seq(length(pDual$faces))) {
    if (!isNormalOutwardFacing(pDual, pDual$faces[[i]])) {
      if (debug) cat("Flip face to make normal outward facing", fill=T)
      pDual$faces[[i]] <- rev(pDual$faces[[i]])
    }
  }
  
  return(pDual)
}

# Create an derived polyhedron by truncating all vertices to the mid of the faces
archi <- function(p, name=paste("archi", p$name), debug=F)
{
  topo <- getTopology(p$faces)
  
  # new points are midpoints of all edges ; the index of each edge is obtained using lookup in the vexToEdge matrix
  archiPoints <- as.data.frame(normalizedistances(t(sapply(seq(max(topo$vexToEdge)), function(e) { 
    w <- which(topo$vexToEdge==e & upper.tri(topo$vexToEdge)) # shall be one and only one
    return(apply(p$vertices[c((w - 1) %/% nrow(topo$vexToEdge) + 1, 
                              (w - 1) %% nrow(topo$vexToEdge) + 1),],2,mean))
  }))))
  # one set of faces comes from the vertices neighbouring each point
  archiFaces1 <- lapply(topo$vexConnections, function(vex) { return(topo$vexToEdge[vex$vex, vex$center])})
  # the other set comes from the faces, using midpoints of their edges
  archiFaces2 <- lapply(p$faces, function(f) { sapply(seq(length(f)), function(j) {return(topo$vexToEdge[f[j], shift(f)[j]])})})
  
  pArchi <- list(vertices = archiPoints, faces = c(archiFaces1, archiFaces2), name=name)
  
  # make sure all faces are oriented consistently
  for (i in seq(length(pArchi$faces))) {
    if (!isNormalOutwardFacing(pArchi, pArchi$faces[[i]])) {
      if (debug) cat("Flip face to make normal outward facing", fill=T)
      pArchi$faces[[i]] <- rev(pArchi$faces[[i]])
    }
  }
  
  return(pArchi)
}


# Combine two polyhedra into one
compose <- function(p1, p2, name=paste("compose", paste(p1$name, p2$name, sep=",")), debug=F)
{
  # TODO unify vertices!!
  # TODO finish
  
  # start with giving the p2 vertices the identity reference
  p2NewReference <- nrow(p1$vertices) + seq(nrow(p2$vertices)) 
  
  # then see which ones are identical to p1 vertices and track that index
  for (v1 in seq(nrow(p1$vertices))) { 
    samePoints <- which(deltaEquals(distance(p1$vertices[v1,], p2$vertices[seq(nrow(p2$vertices)),]),0))
    if (length(samePoints) > 0) {
      p2NewReference[ samePoints ] <- v1  
    }
  }
  
  # TODO now relabel indices the faces of p2 using the mapping p2NewReference
  
  return(list(vertices = rbind(p1$vertices, p2$vertices),
              faces = c(p1$faces, lapply(p2$faces, function(f) { return(f+nrow(p1$vertices)) })),
              name = name))
}

# Create new polyhedron by chopping off the vertices replacing each by a new face
# TODO lin alg to find new points is not ok yet
snub <- function(p, name = paste("snub", p$name), debug=F)
{
  topo <- getTopology(p$faces)
  
  # every vertex becomes a new face with all new points
  allPoints <- NULL
  allFaces <- list()
  for (v in topo$vexConnections) {
    # create new points close to vertex center C in direction of the connected points P
    # new point = C + alpha*(PC)
    angles <- innerAngles(p$vertices[v$vex,], center=p$vertices[v$center,])
    # find alpha such that the sides of the newly created faces are equal
    alpha <- 0.5 - 0.5*(sin(angles/2)/(1+sin(angles/2))) # not trivial but easily derived
    
    newPoints <- p$vertices[v$vex,]*alpha + (1-alpha)*data.frame(x=rep(p$vertices$x[v$center], length(v$vex)), 
                                                                 y=rep(p$vertices$y[v$center], length(v$vex)),
                                                                 z=rep(p$vertices$z[v$center], length(v$vex)))
    newPoints$from <- v$center
    newPoints$to <- v$vex
    if (is.null(allPoints)) {
      allPoints <- newPoints
    } else {
      allPoints <- rbind(allPoints, newPoints)
    }
    allFaces[[v$center]] <- (nrow(allPoints)-nrow(newPoints)+1):nrow(allPoints)
  }
  # now transform the old faces
  newPointsLookup <- matrix(data = NA, nrow = nrow(allPoints), ncol = nrow(allPoints)) # sparse?
  for (i in seq(nrow(allPoints))) {
    # TODO maybe check for inconsistency if there's a value already
    newPointsLookup[allPoints$from[i], allPoints$to[i]] <- i
  }
  for (f in p$faces)
  {
    snubbedFace <- as.vector(sapply(seq(length(f)), function(idx) {
      return(c(newPointsLookup[f[idx], shift(f)[idx]], 
               newPointsLookup[shift(f)[idx], f[idx]]))}))
    allFaces[[length(allFaces)+1]] <- snubbedFace
  }
  
  pSnub <- list(vertices = allPoints[,1:3], faces = allFaces, name=name)
  
  # make sure all faces are oriented consistently
  for (i in seq(length(pSnub$faces))) {
    if (!isNormalOutwardFacing(pSnub, pSnub$faces[[i]])) {
      if (debug) cat("Flip face to make normal outward facing", fill=T)
      pSnub$faces[[i]] <- rev(pSnub$faces[[i]])
    }
  }
  
  return(pSnub)
}

# clear3d()
# drawSinglePoly(snub(cube), debug=T)

description <- function(p, debug=F)
{
  combineDescriptions <- function(descrs, suffixIdentical = "")
  {
    descrs <- data.frame(table(unlist(descrs)), stringsAsFactors = F)
    if (nrow(descrs) == 1) {
      # all are identical
      if (descrs$Freq[1] == 1) {
        # there is just one
        return(as.character(descrs$Var1[1]))
      } else {
        # there are multiple but all are the same
        return(paste0(descrs$Var1[1], suffixIdentical))
      }
    } else {
      # returning them as a frequency table (12 x a + 5 x b)
      return(paste(sapply(seq(nrow(descrs)), function(i) {return(paste0(descrs$Freq[i], "x", descrs$Var1[i]))}), collapse=" + "))
    }
  }
  
  getFaceDescription <- function(f, addlengths=debug, addangles=debug)
  {
    sidelengths <- round(distance(p$vertices[f,], p$vertices[shift(f),]),6)
    baseDescr <- as.character(length(f))
    if (addlengths) {
      baseDescr <- paste0(baseDescr, " [length: ", combineDescriptions(sidelengths), "]")
    }
    
    # TODO check for regularity of angles and face lengths 
    angles <- innerAngles(p$vertices[f,])
    if (deltaEquals(sum(angles), 2*pi)) {
      # simple polygon
      fDescr <- baseDescr
    } else {
      # polygon with multiple rounds
      fDescr <- paste(baseDescr, sum(angles)/(2*pi), sep="/") # TODO maybe round
    }
    #cat("Face:",f,fDescr,fill=T)
    return(fDescr)
  }
  
  getVertexDescription <- function(aVertex)
  {
    faceDescriptions <- sapply(aVertex$faces, function(fi) { return(getFaceDescription(p$faces[[fi]])) })
    vexDescriptions <- getFaceDescription(aVertex$vex, addlengths=F)
    
    if (length(unique(faceDescriptions)) == 1 & length(unique(vexDescriptions)) == 1) {
      vexDescription <- paste0("{",faceDescriptions[1],",",vexDescriptions[1],"}")
    } else {
      vexDescription <- paste0("{",paste(faceDescriptions, collapse=","),"}") # TODO put in best order
    }
    
    return(vexDescription)
  }
  
  getBodyDescription <- function(bodyVertices)
  {
    vertexDescriptions <- sapply(bodyVertices, getVertexDescription)
    # TODO if just a single one return that
    # otherwise switch over to counting the faces
    bodyDescription <- combineDescriptions(vertexDescriptions)
    return(bodyDescription)
  }
  
  topo <- getTopology(p$faces)  
  bodies <- findDistinctBodies(p, topo)
  
  bodyDescriptions <- lapply(bodies, function(b) {
    # subselect vertices that have faces in current body
    getBodyDescription(topo$vexConnections[sapply(topo$vexConnections, function(t) {return(length(intersect(t$faces, b)) > 0 )})])
  })
  polyDescription <- combineDescriptions(bodyDescriptions, paste0("X", length(bodyDescriptions)))
  
  return(polyDescription)
}

# description(octahedron)
# description(cube)
# description(smallStellatedDodecahedron)
# description(archi(cube))
# description(compose(cube, octahedron))
# 
# sapply(lapply(Regulars, dual), description)

# p <- octahedron #smallStellatedDodecahedron
# f <- p$faces[[1]]
# innerAngles(p$vertices[f,])*360/2/pi

# to test the above
# p1<-smallStellatedDodecahedron
# p2<-dual(p1, scaling="vertex")

# # checking the normals of the face
# p <- octahedron
# drawPoly(p, debug=T)

