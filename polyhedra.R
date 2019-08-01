# Functions specifically to constructing polyhedra

library(data.table)

source("math.R")
source("draw.R")

# Returns T if outward facing (= rotating anti-clockwise when looking down towards face in direction of origin)
isNormalOutwardFacing <- function(p, f)
{
  n <- normal( p$coords[f[1],], p$coords[f[2],], p$coords[f[3],] )
  mid <- apply(p$coords[f,],2,mean)
  return (vectorlength(mid+n) > vectorlength(mid-n))
}

# Build up a face from existing coords in given polygon, within constraints passed in as
# the number of edges/vertices of this face, number of faces per vertex and the length of
# an edge. The face does not need to be regular but will have all edges of the same length.
buildFace <- function(p, polygonsize, vertexsize, edgelength, aFace = c(), debug = F)
{
  if (debug) cat("Face candidate:", paste0("{",paste0(aFace, collapse = "-"),"}"), fill=T)
  
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
      normal <- normal(p$coords[aFace[1],], p$coords[aFace[2],], p$coords[aFace[3],])
      
      for (f in p$faces) {
        verticesInNewPlane <- deltaEquals(normal %*% apply(p$coords[f,], 1, function(x) {return(x-unlist(p$coords[aFace[1],]))}), 0)
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
  candidates <- setdiff(which( sapply(seq(nrow(p$coords)), function(i) { sum(unlist(p$faces) == i)}) < vertexsize ), aFace)
  if (debug) cat("Free vertices:", candidates, fill = T)
  
  # Check vertex right distance to previous point
  if (length(candidates) > 0 & length(aFace) > 0) {
    candidates <- candidates[which(deltaEquals(distance(p$coords[aFace[length(aFace)],], p$coords[candidates,]), edgelength))]  
    if (debug) cat("Right distance to prev:", candidates, fill = T)
  }
  
  # Check last vertex right distance to first point
  if (length(candidates) > 0 & length(aFace) == (polygonsize-1)) {
    candidates <- candidates[which(deltaEquals(distance(p$coords[aFace[1],], p$coords[candidates,]), edgelength))]  
    if (debug) cat("Right distance to first:", candidates, fill = T)
  }
  
  # Check if the 4th and further points are in the same plane as the first 3
  if (length(candidates) > 0 & length(aFace) >= 3) {
    # calculate normal from three of the face points
    normal <- normal(p$coords[aFace[1],], p$coords[aFace[2],], p$coords[aFace[3],])
    
    # then the inner product of the normal with the vector of face point 1 to the candidates is 0 when in the same plane
    candidatevectors <- t(apply(matrix(p$coords[candidates,], ncol=3), 1, function(x) {return(x-p$coords[aFace[1],])}))
    candidates <- candidates[which(deltaEquals(apply(candidatevectors, 1, function(x) {return(normal %*% x)}), 0))]
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

# Build a polyhedron given a set of coords (full x, y, z coordinates), the two vertices
# of an example edge (used to determine the global edge size of this polyhedron)
buildRegularPoly <- function(coords, polygonsize, vertexsize, exampleEdge = c(1,2), name = "", debug=F)
{
  if (!is.matrix(coords)) coords <- as.matrix(coords)
  
  # Scale the coords to unit length
  coords <- normalizedistances(coords)
  
  if (debug) {
    spacing <- 0.1
    spheres3d(coords[,1], coords[,2], coords[,3], color="green", radius = 0.02)
    text3d((1+spacing)*coords[,1], (1+spacing)*coords[,2], (1+spacing)*coords[,3], text = seq(nrow(coords)), color="blue")
  }
  
  # Determine global edge size from given example edge
  edgelength <- distance(coords[exampleEdge[1],], coords[exampleEdge[2],])
  
  # Add new faces one by one until all vertices have the specified number of faces
  poly <- list(coords = coords, faces = list(), name = name)
  edges <- matrix(data = 0, nrow = nrow(coords), ncol = nrow(coords))
  repeat {
    # Check if there's any coordinates not having the desired number of faces connected
    freeVertices <- which( sapply(seq(nrow(poly$coords)), function(i) { sum(unlist(poly$faces) == i)}) < vertexsize )
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
        edges[f[i], shiftrotate(f)[i]] <- 1 + edges[f[i], shiftrotate(f)[i]]
        edges[shiftrotate(f)[i], f[i]] <- 1 + edges[shiftrotate(f)[i], f[i]]
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

# Figures out the topology of given polyhedron. Returns a list with multiple
# elements, depending on flags.
#
# $edges: list of list of
#         faces = c(F1,F2) : the two neighbouring faces (index into p$faces), not in any order and one of them can be empty if there are gaps
#         vertices = c(V1, V2) : the two vertices (index into $vertices), also not in order but neither can be empty
# $vertices: list of list of
#         faces = c(...) : ordered list of faces (index into p$faces) around this vertex
#         edges = c(...) : ordered list of edges (index into $edges) around this vertex
#         points = c(...) : ordered list of points (index into p$coords) around this vertex
#         center = P : center of this vertex (index into p$coords) (vertices != coordinates, e.g. with mulitple bodies)
# $bodies: list of
#         c(...) : faces in this body (index into p$faces)

topology <- function(p, edges=T, vertexFigures=T, bodies=T)
{
  
}

# Derives the topology of a polyhedron based on just a list of faces. Each face consists of a list
# of references to coords as everywhere else. The returned structure contains
#
# - a square matrix vexToRFace mapping vertex pairs to the index of the face on the right-hand side of it
# - a square and symetrical matrix vexToEdge mapping vertex pairs to the number (index) of the edge connecting both
# - an unordered list vexConnections with for every vertex of the polygon (NB in compounds there could be multiple of these for a vertex):
#       - the index of the center vertex
#       - the indices of the faces surrounding the vertex (in order)
#       - the indices of the coords connected to this vertex (in order)

# TODO make more robust
# edge list vexToEdge[i,j] gives (newly created) edge number for vertex i to vertex j by following along all faces
#           vexToEdge[j,i] = vexToEdge[i,j] 
# face list vexToFace[i,j] gives the RHS face if face is oriented normal etc
#           in case of collission a face can be turned upside down - previous face assignments to be corrected!
#           but a face can only be turned once 
# face orientation = array for all faces starts with NA then T if follows face orientation, F otherwise
#
# with that we know which faces are neighbouring - the symmetry axis of the VexToFace 
# so this should solve the bodies problem immediately
# gaps exist if there is no symmetry
# there can also be orphaned vertices
#
# we also know which faces occur at any given vertex -->
# build vexFigure by orienting the faces by ordering the row elements of vexToFace using the symmetry
# vexFigure indices correspond to coordinate indices

getTopology <- function(faces, debug=F)
{
  if (!is.null(names(faces))) {
    if (is.list(faces) & "faces" %in% names(faces)) {
      faces <- faces$faces # this method is commonly called with a polygon instead of just the faces
    }
  }
  n_vertices <- max(sapply(faces, max))
  foundAnyGaps <- F
  
  getFaceForDebug <- function(faceidx)
  {
    if (faceidx < 1 | faceidx > length(faces)) {
      return (paste0("F",faceidx,"**invalid**"))
    }
    return (paste(paste0("F",faceidx,":"),paste(faces[[faceidx]],collapse=";")))
  }
  
  # Build a list of all edges by going round all faces, storing
  # the edge index in a 2D matrix indexed by vertex indices. Not
  # currently keeping a list of edges->faces directly but could be
  # done trivially if needed.
  
  vexToEdge <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  vexToRFace <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  for (currentFaceNr in seq(length(faces))) {
    currentFace <- faces[[currentFaceNr]]
    currentFaceS <- shiftrotate(currentFace)
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
          # TODO perhaps rotate ? or keep the fact that this face is rotated in another
          # matrix of i/j to face, with rotated faces as -1, then for those faces use
          # vexToRFace the other way around?
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
        prevPoint <- shiftrotate(faces[[f]],-1)[i] # point on an edge of face f that connects to centerPoint
        vexToVertex <- c(vexToVertex, prevPoint) # keep track of an oriented list of connected points around centerPoint
        
        nextPoint <- shiftrotate(faces[[f]])[i] # find next face f that connects at edge nextPoint - centerPoint
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
      isConnectedToFace <- sapply(conn, function(c) { return(f %in% c$faces)})
      if (length(isConnectedToFace) == 0) return(body) # f not connected, treat as separate body
      vex <- which(isConnectedToFace)
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
  
  newVertexCoords <- t(sapply(p$faces, function(f) { return(apply(p$coords[f,],2,mean))}))
  newFaces <- lapply(topo$vexConnections, function(c) {return(c$faces)})
  
  # scale so that mid of a new vertex is at same distance from origin as mid
  # of an old vertex - works at least for regulars, otherwise somewhat arbitrary
  if (scaling == "edge") {
    scale <- vectorlength(apply(p$coords[p$faces[[1]][1:2],],2,mean)) / 
      vectorlength(apply(newVertexCoords[newFaces[[1]][1:2],],2,mean))
  } else if (scaling == "vertex") {
    scale <- vectorlength(p$coords[p$faces[[1]][1],])/
      vectorlength(newVertexCoords[newFaces[[1]][1],])
  } else if (scaling == "none") {
    scale <- 1
  } else {
    stop(paste("Wrong scaling argument:", scaling))
  }
  
  pDual <- list( coords = newVertexCoords*scale, faces = newFaces, name = name)
  
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

# Create an derived polyhedron by truncating all coords to the mid of the faces
archi <- function(p, name=paste("archi", p$name), debug=F)
{
  topo <- getTopology(p$faces)
  
  # new points are midpoints of all edges ; the index of each edge is obtained using lookup in the vexToEdge matrix
  archiPoints <- normalizedistances(t(sapply(seq(max(topo$vexToEdge)), function(e) { 
    w <- which(topo$vexToEdge==e & upper.tri(topo$vexToEdge)) # shall be one and only one
    return(apply(p$coords[c((w - 1) %/% nrow(topo$vexToEdge) + 1, 
                            (w - 1) %% nrow(topo$vexToEdge) + 1),],2,mean))
  })))
  # one set of faces comes from the vertices neighbouring each point
  archiFaces1 <- lapply(topo$vexConnections, function(vex) { return(topo$vexToEdge[vex$vex, vex$center])})
  # the other set comes from the faces, using midpoints of their edges
  archiFaces2 <- lapply(p$faces, function(f) { sapply(seq(length(f)), function(j) {return(topo$vexToEdge[f[j], shiftrotate(f)[j]])})})
  
  pArchi <- list(coords = archiPoints, faces = c(archiFaces1, archiFaces2), name=name)
  
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
  # TODO unify coords!!
  # TODO finish
  
  # start with giving the p2 coords the identity reference
  p2NewReference <- nrow(p1$coords) + seq(nrow(p2$coords)) 
  
  # then see which ones are identical to p1 coords and track that index
  for (v1 in seq(nrow(p1$coords))) { 
    samePoints <- which(deltaEquals(distance(p1$coords[v1,], p2$coords[seq(nrow(p2$coords)),]),0))
    if (length(samePoints) > 0) {
      p2NewReference[ samePoints ] <- v1  
    }
  }
  
  # TODO now relabel indices the faces of p2 using the mapping p2NewReference
  
  return(list(coords = rbind(p1$coords, p2$coords),
              faces = c(p1$faces, lapply(p2$faces, function(f) { return(f+nrow(p1$coords)) })),
              name = name))
}

# Create new polyhedron by chopping off the coords replacing each by a new face
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
    angles <- innerAngles(p$coords[v$vex,], center=p$coords[v$center,])
    # find alpha such that the sides of the newly created faces are equal
    alpha <- 0.5 - 0.5*(sin(angles/2)/(1+sin(angles/2))) # not trivial but easily derived
    
    newPoints <- p$coords[v$vex,]*alpha + (1-alpha)*data.table(x=rep(p$coords[v$center,1], length(v$vex)), 
                                                               y=rep(p$coords[v$center,2], length(v$vex)),
                                                               z=rep(p$coords[v$center,3], length(v$vex)))
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
      return(c(newPointsLookup[f[idx], shiftrotate(f)[idx]], 
               newPointsLookup[shiftrotate(f)[idx], f[idx]]))}))
    allFaces[[length(allFaces)+1]] <- snubbedFace
  }
  
  pSnub <- list(coords = allPoints[,1:3], faces = allFaces, name=name)
  
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
    descrFreqs <- data.table(table(unlist(descrs)), stringsAsFactors = F)
    if (nrow(descrFreqs) == 1) {
      # all are identical
      if (descrFreqs$N[1] == 1) {
        # there is just one
        return(as.character(descrFreqs$V1[1]))
      } else {
        # there are multiple but all are the same
        return(paste0(descrFreqs$V1[1], suffixIdentical))
      }
    } else {
      # returning them as a frequency table (12 x a + 5 x b)
      return(paste(sapply(seq(nrow(descrFreqs)), function(i) {return(paste0(descrFreqs$N[i], "x", descrFreqs$V1[i]))}), collapse=" + "))
    }
  }
  
  getFaceDescription <- function(f, addlengths=debug, addangles=debug)
  {
    sidelengths <- round(distance(p$coords[f,], p$coords[shiftrotate(f),]),6)
    baseDescr <- as.character(length(f))
    if (addlengths) {
      baseDescr <- paste0(baseDescr, " [length: ", combineDescriptions(sidelengths), "]")
    }
    
    # TODO check for regularity of angles and face lengths 
    angles <- innerAngles(p$coords[f,])
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

# Platonic solids
# all coords taken from https://en.wikipedia.org/wiki/Platonic_solid

tetrahedron <- buildRegularPoly(coords = rbind(data.table(x=1, y=1, z=1), data.table(x=1, y=-1, z=-1), data.table(x=-1, y=1, z=-1), data.table(x=-1, y=-1, z=1)),
                                polygonsize = 3,
                                vertexsize = 3,
                                name = "Tetrahedron")

# constructing the cube directly will help as this is the first one with faces with > 3 sides
cube <- buildRegularPoly(coords = expand.grid(x = c(-1, 1), y = c(-1, 1), z = c(-1, 1)),
                         polygonsize = 4,
                         vertexsize = 3,
                         name = "Cube")

octahedron <- buildRegularPoly(coords = rbind(expand.grid(x = c(-1,1), y = 0, z = 0), expand.grid(x = 0, y = c(-1,1), z = 0), expand.grid(x = 0, y = 0, z = c(-1,1))),
                               polygonsize = 3,
                               vertexsize = 4,
                               exampleEdge = c(1,3),
                               name = "Octahedron")

icosahedron <- buildRegularPoly(coords = rbind(expand.grid(x = 0, y = c(-1,1), z = c(-phi, phi)), 
                                               expand.grid(x = c(-1,1), y = c(-phi, phi), z = 0), 
                                               expand.grid(x = c(-phi, phi), y = 0, z = c(-1, 1))),
                                polygonsize = 3,
                                vertexsize = 5,
                                name = "Icosahedron",
                                debug = F)


# description(octahedron)
# description(cube)
# description(smallStellatedDodecahedron)
# description(archi(cube))
# description(compose(cube, octahedron))
# 
# sapply(lapply(Regulars, dual), description)

# p <- octahedron #smallStellatedDodecahedron
# f <- p$faces[[1]]
# innerAngles(p$coords[f,])*360/2/pi

# to test the above
# p1<-smallStellatedDodecahedron
# p2<-dual(p1, scaling="vertex")

# # checking the normals of the face
# p <- octahedron
# drawPoly(p, debug=T)

