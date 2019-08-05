library(data.table)

# Things that do not depend on coordinates


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
# - a square matrix coordPairToFaces mapping vertex pairs to the index of the face on the right-hand side of it
# - a square and symetrical matrix coordPairToEdge mapping vertex pairs to the number (index) of the edge connecting both
# - an unordered list vexConnections with for every vertex of the polygon (NB in compounds there could be multiple of these for a vertex):
#       - the index of the center vertex
#       - the indices of the faces surrounding the vertex (in order)
#       - the indices of the coords connected to this vertex (in order)

# TODO make more robust
# edge list coordPairToEdge[i,j] gives (newly created) edge number for vertex i to vertex j by following along all faces
#           coordPairToEdge[j,i] = coordPairToEdge[i,j] 
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

# TODO make structures part of p, not a seperate structure

getTopology <- function(p, debug=F)
{
  faces <- p$faces # possibility to cache the topology or parts of it
  
  n_vertices <- max(sapply(faces, max))
  foundAnyGaps <- F # TODO derive from face matrix below (a-symmetries)
  
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
  
  coordPairToEdge <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  coordPairToFaces <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  for (f in seq(length(faces))) {
    currentFace <- faces[[f]]
    currentFaceS <- shiftrotate(currentFace)
    for (i in seq(length(currentFace))) {
      if (0 == coordPairToFaces[currentFace[i], currentFaceS[i]]) {
        coordPairToFaces[currentFace[i], currentFaceS[i]] <- f
      } else {
        if (0 == coordPairToFaces[currentFaceS[i], currentFace[i]]) {
          coordPairToFaces[currentFaceS[i], currentFace[i]] <- f
        } else {
          stop(paste("The faces next to coords", currentFace[i], "-", currentFaceS[i], 
                     "already set:", coordPairToFaces[currentFace[i], currentFaceS[i]], 
                     "and:", coordPairToFaces[currentFaceS[i], currentFace[i]], 
                     "while wanting to set to", f))
        }
      }
      if (0 == coordPairToEdge[currentFace[i], currentFaceS[i]]) {
        coordPairToEdge[currentFaceS[i], currentFace[i]] <- 1+max(coordPairToEdge)
        coordPairToEdge[currentFace[i], currentFaceS[i]] <- coordPairToEdge[currentFaceS[i], currentFace[i]]
      }
    }
  }
  
  # Build edgeToFaces from those two matrices. coordPairToEdge[n, 1:2] gives the two faces next to n.
  edgeToFaces <- matrix(nrow=max(coordPairToEdge), ncol=2, data=0)
  for (i in seq(nrow(coordPairToEdge))) {
    for (j in seq(ncol(coordPairToEdge))) {
      edge <- coordPairToEdge[i,j]
      if (edge != 0) {
        edgeToFaces[edge, 1] <- coordPairToFaces[i,j]  
        edgeToFaces[edge, 2] <- coordPairToFaces[j,i]  
      }
    }
  }
  
  # Finding bodies using above structure. Two faces are in the same body if they share an edge.
  addConnectedFaces <- function (face, edges, currbody = c())
  {
    if (face==0) return( c() ) # not a valid face
    if (face %in% currbody) return( c() ) # already in the current body
    b <- c(face)
    nextEdges <- which(edgeToFaces[edges,1] == face | edgeToFaces[edges,2] == face)
    if (length(nextEdges) == 0) return(b)
    for (nextEdge in edges[nextEdges]) {
      b <- c(b, addConnectedFaces(edgeToFaces[nextEdge,1], setdiff(edges, edges[nextEdges]), b))
      b <- c(b, addConnectedFaces(edgeToFaces[nextEdge,2], setdiff(edges, edges[nextEdges]), b))
    }
    return(b)
  }
  
  bodies <- list()
  repeat {
    remainingEdges <- which(apply(edgeToFaces, 1, function(e2f) {(e2f[1] != 0 & !(e2f[1] %in% unlist(bodies))) & (e2f[2] != 0 & !(e2f[2] %in% unlist(bodies)))}))
    if (length(remainingEdges) == 0) break
    startEdge <- remainingEdges[1]
    currentBody <- unique(c(addConnectedFaces(edgeToFaces[startEdge,1], setdiff(remainingEdges, startEdge)), 
                            addConnectedFaces(edgeToFaces[startEdge,2], setdiff(remainingEdges, startEdge))))
    bodies[[1+length(bodies)]] <- currentBody
  }
  
  # Vertex figures: start with a raw list of all faces connected to a coord, then process each of these
  # to set order (via edges) and split (if there are multiple solids sharing a coord)
  vexFigures <- list()
  rawVertexFigures <- lapply(seq(nrow(coordPairToFaces)), function(i) {return(setdiff(unique(c(coordPairToFaces[i,],
                                                                                               coordPairToFaces[,i])),0))})
  for (vi in seq(length(rawVertexFigures))) {
    # all the edges that connect pairs of the faces in this raw vertex
    v <- rawVertexFigures[[vi]]
    if (length(v) == 0) break # in case of gaps
    # TODO more restrictive, use coords to vex
    remainingRawEdges <- which(sapply(seq(nrow(edgeToFaces)), function(i) {return((edgeToFaces[i,1] %in% c(0,v)) & (edgeToFaces[i,2] %in% c(0, v)))}))
    remainingFaces <- v
    repeat {
      currentFace <- remainingFaces[1] # start
      vertexFaceFig <- currentFace
      vertexVexFig <- c() # TODO = vertexCoordFig
      repeat {
        # the edges neighbouring current face
        nbredges <- remainingRawEdges[edgeToFaces[remainingRawEdges,1]==currentFace | edgeToFaces[remainingRawEdges,2]==currentFace]
        if (length(nbredges) == 0) break # no neighbour edges remaining, must be done with this vertex
        nextFace <- ifelse(edgeToFaces[nbredges[1],1] == currentFace, edgeToFaces[nbredges[1],2], edgeToFaces[nbredges[1],1])
        vertexVexFig <- c(vertexVexFig, which(coordPairToEdge[vi,]==nbredges[1])) # other coord connected by edge
        remainingRawEdges <- setdiff(remainingRawEdges, nbredges[1])
        if (length(remainingRawEdges) == 0) break # really done, breaks second repeat also
        if (nextFace == remainingFaces[1]) break # back to start face
        vertexFaceFig <- c(vertexFaceFig, nextFace)
        currentFace <- nextFace
      }
      #print(vertexFaceFig)
      #print(vertexVexFig)
      vexFigures[[1+length(vexFigures)]] <- list(center = vi,
                                                         faces = vertexFaceFig,
                                                         vex = vertexVexFig)
      
      remainingFaces <- setdiff(remainingFaces, vertexFaceFig)
      if (length(remainingFaces) == 0) break
    }
  }
  


  # # TODO find vertex figures by going rowwise through coordPairToFaces. Every row is one or
  # # more vertex figures. Order by connecting the faces found in one row that share 2 coords.
  # 
  # # Knowing which edges connect to which faces we can now build up oriented lists of
  # # faces and eges around any vertex
  # 
  # isConnectedToPreviousStartFace <- rep(F, length(faces))
  # isUsedAsStartFace <- rep(F, length(faces))
  # vexConnections <- list()
  # startingFace <- NULL
  # repeat {
  #   if (is.null(startingFace)) {
  #     # start just somewhere
  #     startingFace <- 1
  #   } else {
  #     if (any(isConnectedToPreviousStartFace & !isUsedAsStartFace)) {
  #       # prefer to continue on a face connected to any of the faces previously used
  #       startingFace <- which(isConnectedToPreviousStartFace & !isUsedAsStartFace)[1]
  #     } else {
  #       # otherwise there may be a disjunct set of faces (multiple bodies) so start there
  #       if (any(!isUsedAsStartFace)) {
  #         startingFace <- which(!isUsedAsStartFace)[1]
  #       } else {
  #         # all faces used, we're done
  #         break;
  #       }
  #     }
  #   }
  #   isUsedAsStartFace[startingFace] <- T
  # 
  #   # for (startingFace in seq(length(faces)))
  #   if (debug) cat("start face:", getFaceForDebug(startingFace), fill=T)
  # 
  #   # check which vex we already have with this face in it
  #   coveredPoints <- c()
  #   if (length(vexConnections) > 0) {
  #     vexCoveringThisStartingFace <- which(sapply(vexConnections, function(con) { return(startingFace %in% con$faces)}))
  #     if (length(vexCoveringThisStartingFace) > 0) {
  #       coveredPoints <- sapply(vexCoveringThisStartingFace, function(i) {return(vexConnections[[i]]$center)})
  #       if (debug) cat("skipping already covered points for this face:", paste0(coveredPoints, collapse=";"), fill=T)
  #     }
  #   }
  # 
  #   for (centerPoint in setdiff(faces[[startingFace]], coveredPoints))
  #   {
  #     if (debug) cat("start building vertex around", centerPoint, fill=T)
  #     vexToFaces <- c()
  #     vexToVertex <- c()
  #     f <- startingFace
  # 
  #     repeat
  #     {
  #       vexToFaces <- c(vexToFaces, f) # keep track of an oriented list of connected faces around centerPoint
  #       isConnectedToPreviousStartFace[f] <- T
  #       if (debug) cat("center",centerPoint, "face",getFaceForDebug(f), fill=T)
  #       i <- which(faces[[f]] == centerPoint) # index in face f of centerPoint
  #       prevPoint <- shiftrotate(faces[[f]],-1)[i] # point on an edge of face f that connects to centerPoint
  #       vexToVertex <- c(vexToVertex, prevPoint) # keep track of an oriented list of connected points around centerPoint
  # 
  #       nextPoint <- shiftrotate(faces[[f]])[i] # find next face f that connects at edge nextPoint - centerPoint
  #       f <- coordPairToFaces[ nextPoint, centerPoint ]
  #       if (f == 0) {
  #         foundAnyGaps <- T
  #         cat("WARNING: there's a gap, no face to the right of the edge from", centerPoint, "to", nextPoint, fill=T)
  #         # don't add this vex although in theory we could try adding vertices with gaps by first trying going
  #         # the other way from the start - becomes messy and what if there are multiple gaps
  #         break
  #       }
  # 
  #       if (f == startingFace | f == 0) {
  #         # check if we don't have this already - not just checking the same center point but also
  #         # the set of faces as we could have multiple bodies just sharing a vertex but no faces
  #         if (!any(sapply(vexConnections, function(con) { return( centerPoint == con$center &
  #                                                                 setequal(vexToFaces, con$faces) &
  #                                                                 setequal(vexToVertex, con$vex))})))
  #         {
  #           if (debug) cat("adding vex center",centerPoint,"faces:",paste(paste0("F",vexToFaces), collapse = ";"),"vex:",paste(vexToVertex, collapse=";"), fill=T)
  #           vexConnections[[1+length(vexConnections)]] <- list(center = centerPoint,
  #                                                              faces = vexToFaces,
  #                                                              vex = vexToVertex)
  #         } else {
  #           # not sure this can even happen with the more careful choice of starting vertices in place now
  #           if (debug) cat("skipping adding vex center",centerPoint,"faces:",paste(paste0("F",vexToFaces), collapse = ";"),"vex:",paste(vexToVertex, collapse=";"), fill=T)
  #         }
  #         break # next vertex
  #       }
  #     }
  #   }
  # }
  
  return (list(coordPairToEdge=coordPairToEdge, 
               coordPairToFaces=coordPairToFaces, 
               #vexConnections=vexConnections, 
               vexConnections=vexFigures,
               bodies=bodies,
               hasGaps=foundAnyGaps))
}

# Return list of list of distinct bodies in p. List contains the face indices.
findDistinctBodies <- function(p, topo = getTopology(p), debug=F)
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

testTopology <- function()
{
  t <- getTopology(cube)
  t <- getTopology(compose(cube, dual(cube))) # two bodies
  
  # two strange pyramids sharing a coord
  coords <- as.matrix(rbind(data.table(x=cos(2*pi*(0:3)/4), y=sin(2*pi*(0:3)/4), z=0),
                            data.table(x=cos(2*pi*(1/8+(0:3)/4)), y=sin(2*pi*(1/8+(0:3)/4)), z=0),
                            data.table(x=0,y=0,z=1)))
  strangedualpyramid <- list(coords=coords, faces=list(c(1:4), c(5:8), 
                                                       c(1,2,9), c(2,3,9), c(3,4,9), c(4,1,9),
                                                       c(5,6,9), c(6,7,9), c(7,8,9), c(8,5,9)))
  t <- getTopology(strangedualpyramid)
  drawPoly(strangedualpyramid)
  
  # todo similar with gaps
  
  poly <- buildRegularPoly(dodecahedron$coords, 5, 6, c(1, 4))
  drawSinglePoly(poly)
  
  simplified <- list(coords=poly$coords, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
  t <- getTopology(simplified) # errors
  drawSinglePoly(simplified, debug=T)
}
# poly <- buildRegularPoly(dodecahedron$coords, 5, 6, c(1, 4))
# clear3d()
# drawSinglePoly(poly, debug = T) # this shows the error drawing {5/2} # but now goes into a loop?
# # huh, now it works but not in non-debug
# topo <- getTopology(poly, debug=T) # this errors out on edge 4-1 / into a loop
# 
# # face is "below" origin
# faces <- poly$faces
# debug <- T
# simplified <- list(coords=poly$coords, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
# print(simplified$faces[[1]])
# print(simplified$faces[[2]])
# drawSinglePoly(simplified, debug = T) # shows drawing error and orientation issue when using [[1]],[[2]]
# getTopology(simplified) # loops
# 
# # with a fix we get into a loop in 
# getTopology(simplified, debug=T)


