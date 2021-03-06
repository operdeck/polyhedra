library(data.table)

# Things that do not depend on coordinates

dropZeros <- function(a) { return(a[a!=0]) }

# Derives the topology of a polyhedron based on just a list of faces. Each face consists of a list
# of references to coords as everywhere else. The returned structure contains
#
# - a square matrix coordPairToFaces mapping point (coordinate) pairs to the index of the faces next to it - no order implied
# - a square and symetrical matrix coordPairToEdge mapping point (coordinate) pairs to the number (index) of the edge connecting both
# - a matrix edgeToFaces with the two (or less in case of gaps) faces connected to an edge
# - a matrix edgeToVertices with the two vertices that form an edge, reverse of coordPairToEdge
# - a list bodies with each element a distinct set of faces belonging to a seperate body
# - an unordered list vertexFigures with for every vertex of the polygon (NB in compounds there could be multiple of these per point):
#       - the coordinate index of the center point
#       - the indices of the faces surrounding the vertex (in order)
#       - the indices of the coords connected to this vertex (in order)

topology <- function(p, debug=F)
{
  faces <- p$faces # possibility to cache the topology or parts of it
  
  n_vertices <- max(c(0,unlist(sapply(faces,max))))

  getFaceForDebug <- function(faceidx)
  {
    if (faceidx < 1 | faceidx > length(faces)) {
      return (paste0("F",faceidx,"**invalid**"))
    }
    return (paste(paste0("F",faceidx,":"),paste(faces[[faceidx]],collapse=";")))
  }
  
  # Map pairs of coordinates to face and edge indices by going round all faces

  coordPairToEdge <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  coordPairToFaces <- matrix(nrow=n_vertices, ncol=n_vertices, data=0)
  for (f in safeseq(length(faces))) {
    currentFace <- faces[[f]]
    currentFaceS <- shiftrotate(currentFace)
    for (i in safeseq(length(currentFace))) {
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
  
  # Build edgeToFaces and edgeToVertices from those two matrices.
  
  edgeToFaces <- matrix(nrow=max(c(0,coordPairToEdge)), ncol=2, data=0)
  edgeToVertices <- matrix(nrow=max(c(0,coordPairToEdge)), ncol=2, data=0)
  for (i in safeseq(nrow(coordPairToEdge))) {
    for (j in safeseq(ncol(coordPairToEdge))) {
      edge <- coordPairToEdge[i,j]
      if (edge != 0) {
        edgeToFaces[edge, 1] <- coordPairToFaces[i,j]  
        edgeToFaces[edge, 2] <- coordPairToFaces[j,i]  
        
        edgeToVertices[edge, 1] <- i
        edgeToVertices[edge, 2] <- j
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
    # remainingEdges <- which(apply(edgeToFaces, 1, function(e2f) {(e2f[1] != 0 & !(e2f[1] %in% unlist(bodies))) & (e2f[2] != 0 & !(e2f[2] %in% unlist(bodies)))}))
    remainingEdges <- which(apply(edgeToFaces, 1, function(e2f) {!(e2f[1] %in% unlist(bodies)) & !(e2f[2] %in% unlist(bodies))}))
    if (length(remainingEdges) == 0) break
    startEdge <- remainingEdges[1]
    currentBody <- unique(c(addConnectedFaces(edgeToFaces[startEdge,1], setdiff(remainingEdges, startEdge)), 
                            addConnectedFaces(edgeToFaces[startEdge,2], setdiff(remainingEdges, startEdge))))
    bodies[[1+length(bodies)]] <- currentBody
  }
  
  # Vertex figures: start with a raw list of all faces connected to a coord, then process each of these
  # to set order (via edges) and split (if there are multiple solids sharing a coord)
  vexFigures <- list()
  rawVertexFigures <- lapply(safeseq(nrow(coordPairToFaces)), function(i) {return(setdiff(unique(c(coordPairToFaces[i,],
                                                                                               coordPairToFaces[,i])),0))})
  for (vertexCenterPoint in safeseq(length(rawVertexFigures))) {
    # all the edges that connect pairs of the faces in this raw vertex
    remainingFaces <- rawVertexFigures[[vertexCenterPoint]]
    if (length(remainingFaces) == 0) break # in case of gaps
    remainingRawEdges <- dropZeros(coordPairToEdge[vertexCenterPoint,])
    repeat {
      currentFace <- remainingFaces[1] # start
      vertexFaceFig <- currentFace
      vertexCoordFig <- c() 
      repeat {
        # the edges neighbouring current face
        nbredges <- remainingRawEdges[edgeToFaces[remainingRawEdges,1]==currentFace | edgeToFaces[remainingRawEdges,2]==currentFace]
        if (length(nbredges) == 0) break # no neighbour edges remaining, must be done with this vertex
        nextFace <- ifelse(edgeToFaces[nbredges[1],1] == currentFace, edgeToFaces[nbredges[1],2], edgeToFaces[nbredges[1],1])
        vertexCoordFig <- c(vertexCoordFig, which(coordPairToEdge[vertexCenterPoint,]==nbredges[1])) # other coord connected by edge
        remainingRawEdges <- setdiff(remainingRawEdges, nbredges[1])
        if (length(remainingRawEdges) == 0) break # really done, breaks second repeat also
        if (nextFace == remainingFaces[1]) break # back to start face
        vertexFaceFig <- c(vertexFaceFig, nextFace)
        currentFace <- nextFace
      }
      #print(vertexFaceFig)
      #print(vertexCoordFig)
      vexFigures[[1+length(vexFigures)]] <- list(center = vertexCenterPoint,
                                                         faces = vertexFaceFig,
                                                         vex = vertexCoordFig)
      
      remainingFaces <- setdiff(remainingFaces, vertexFaceFig)
      if (length(remainingFaces) == 0) break
    }
  }

  return (list(coordPairToEdge=coordPairToEdge, 
               coordPairToFaces=coordPairToFaces,
               edgeToFaces=edgeToFaces,
               edgeToVertices=edgeToVertices,
               vertexFigures=vexFigures,
               bodies=bodies))
}

testTopology <- function()
{
  t <- topology(cube)
  t <- topology(compose(cube, dual(cube))) # two bodies
  
  # two strange pyramids sharing a coord
  coords <- as.matrix(rbind(data.table(x=cos(2*pi*(0:3)/4), y=sin(2*pi*(0:3)/4), z=0),
                            data.table(x=cos(2*pi*(1/8+(0:3)/4)), y=sin(2*pi*(1/8+(0:3)/4)), z=0),
                            data.table(x=0,y=0,z=1)))
  strangedualpyramid <- list(coords=coords, faces=list(c(1:4), c(5:8), 
                                                       c(1,2,9), c(2,3,9), c(3,4,9), c(4,1,9),
                                                       c(5,6,9), c(6,7,9), c(7,8,9), c(8,5,9)))
  t <- topology(strangedualpyramid)
  drawPoly(strangedualpyramid)
  
  cubewfewgaps <- list(coords=cube$coords, faces=cube$faces[c(T,F,T,T,F,T)])
  drawSinglePoly(cubewfewgaps, debug=F)
  t <- topology(cubewfewgaps)
  
  cubewithmanygaps <- list(coords=cube$coords, faces=cube$faces[c(F,T,F,F,T,F)])
  drawSinglePoly(cubewithmanygaps, debug=F)
  
  poly <- buildRegularPoly(dodecahedron$coords, 5, 6, c(1, 4))
  drawSinglePoly(poly)
  
  simplified <- list(coords=poly$coords, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
  t <- topology(simplified) # errors
  drawSinglePoly(simplified, debug=F)
}
# poly <- buildRegularPoly(dodecahedron$coords, 5, 6, c(1, 4))
# clear3d()
# drawSinglePoly(poly, debug = T) # this shows the error drawing {5/2} # but now goes into a loop?
# # huh, now it works but not in non-debug
# topo <- topology(poly, debug=T) # this errors out on edge 4-1 / into a loop
# 
# # face is "below" origin
# faces <- poly$faces
# debug <- T
# simplified <- list(coords=poly$coords, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
# print(simplified$faces[[1]])
# print(simplified$faces[[2]])
# drawSinglePoly(simplified, debug = T) # shows drawing error and orientation issue when using [[1]],[[2]]
# topology(simplified) # loops
# 
# # with a fix we get into a loop in 
# topology(simplified, debug=T)


