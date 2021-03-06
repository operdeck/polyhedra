# Functions specifically to constructing polyhedra

library(data.table)

source("geometry.R")
source("draw.R")
source("topology.R")

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
    if (!isNormalOutwardFacing(p$coords[aFace[1:3], ])) {
      if (debug) cat("Flip new face to make normal outward facing", fill=T)
      aFace <- rev(aFace)
    }
    
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

# Internal call to set everything after coords and faces are defined
setPoly <- function(coords, faces, name, original=NULL, originalFaceMap=NULL, debug=F)
{
  # TODO consider making sure the list of coordinates is not superfluous
  # ie only contains references from faces, nothing more, nothing less
  
  if (ncol(coords) != 3) stop("Coords must be 3D")
  coords <- as.matrix(coords)
  
  if (length(faces) < 1) return (NULL)
  
  # make sure all faces are oriented consistently
  for (i in seq(length(faces))) {
    if (!is.null(faces[[i]])) {
      if (!isNormalOutwardFacing(coords[faces[[i]][1:3], ])) {
        if (debug) cat("Flip face", i, "to make normal outward facing", fill=T)
        faces[[i]] <- rev(faces[[i]])
      }
    }
  }
  
  newPoly <- list(coords=coords, faces=faces, name=name, original=original, originalFaceMap=originalFaceMap)
  t <- topology(newPoly)
  for (x in names(t)) {
    newPoly[[x]] <- t[[x]] # copy over the topology attributes
  }
  
  return (newPoly)
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
    drawDots(coords, color="green", radius = 0.02)
    drawTexts((1+spacing)*coords, text = seq(nrow(coords)), color="blue")
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
  
  newPoly <- setPoly(poly$coords, poly$faces, poly$name, debug=debug)
  
  return(newPoly)
}

dual <- function(p, name=paste("dual", p$name), scaling = "edge", debug=F)
{
  newVertexCoords <- t(sapply(p$faces, function(f) { return(apply(p$coords[f,],2,mean))}))
  newFaces <- lapply(p$vertexFigures, function(c) {return(c$faces)})
  
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
  
  pDual <- setPoly(newVertexCoords*scale, newFaces, name, debug=debug)
  
  return(pDual)
}

# Create new faces at the vertices of the old solid effectively
# keeping existing faces but rotating them
quasi <- function(p, name=paste("quasi", p$name), debug=F)
{
  # new points are midpoints of all edges ; the index of each edge is obtained using lookup in the coordPairToEdge matrix
  archiPoints <- normalizedistances(t(sapply(seq(max(p$coordPairToEdge)), function(e) { 
    w <- which(p$coordPairToEdge==e & upper.tri(p$coordPairToEdge)) # shall be one and only one
    return(apply(p$coords[c((w - 1) %/% nrow(p$coordPairToEdge) + 1, 
                            (w - 1) %% nrow(p$coordPairToEdge) + 1),],2,mean))
  })))
  # one set of faces comes from the vertices neighbouring each point
  archiFaces1 <- lapply(p$vertexFigures, function(vex) { return(p$coordPairToEdge[vex$vex, vex$center])})
  # the other set comes from the faces, using midpoints of their edges
  archiFaces2 <- lapply(p$faces, function(f) { sapply(seq(length(f)), function(j) {return(p$coordPairToEdge[f[j], shiftrotate(f)[j]])})})
  
  pArchi <- setPoly(archiPoints, c(archiFaces1, archiFaces2), name, debug=debug)
  
  return(pArchi)
}

# creates new faces at the vertices of the old solid, changing 
# existing faces to have twice as many points
truncate <- function(p, name = paste("truncated", p$name), debug=F)
{
  newpts <- list()
  for(v in p$vertexFigures) 
  {
    if(debug){
      drawTexts(p$coords[v$center,], text = v$center, color="red")
      for (pt in v$vex) 
      {
        drawLines(p$coords[c(v$center, pt),], color="blue")
        drawTexts(p$coords[pt,], text = pt, color="red")
      }
    }
    # new point from center of vertex figure is on line to connected point at relative distance alpha
    alpha <- sapply(seq(length(v$vex)), 
                    function(i) {return(1/(2+vectorlength(p$coords[v$vex[i],]-p$coords[shiftrotate(v$vex)[i],])/
                                             vectorlength(p$coords[v$vex[i],]-p$coords[v$center,])))})
    newpts[[v$center]] <- as.data.table(t(sapply(seq(length(v$vex)),
                                                 function(i) {return((1-alpha[i])*p$coords[v$center,] + alpha[i]*p$coords[v$vex[i], ])})))
    newpts[[v$center]]$center <- v$center
    newpts[[v$center]]$connected <- v$vex
    if (debug) {
      drawDots(as.matrix(newpts[[v$center]][,c("x","y","z")]), color="darkgreen", radius=0.01)
      drawTexts(as.matrix(newpts[[v$center]][,c("x","y","z")]),
                text = newpts[[v$center]]$connected, color="darkgreen")
    }
  }
  
  newcoords <- rbindlist(newpts)
  newcoords[,newcoordidx := seq(.N)]
  
  # the new faces at the old vertices
  newfaces <- newcoords[, .( connected = list(connected), coords = list(newcoordidx)) , by=center]
  
  # the old faces require some careful indexing
  oldfaces <- list()
  for (f in p$faces) {
    prv <- shiftrotate(f,-1)
    nxt <- shiftrotate(f)
    oldfaces[[1+length(oldfaces)]] <- 
      as.numeric(sapply(seq(length(f)), function(i) {return(c(newcoords[center==f[i]&connected==prv[i]]$newcoordidx,
                                                              newcoords[center==f[i]&connected==nxt[i]]$newcoordidx))}))
  }
  
  # done!
  trunc2 <- setPoly(coords = as.matrix(newcoords[,1:3]),
                    faces = c(newfaces$coords, oldfaces), name=name, debug=debug)
  return (trunc2)
}

# Edges become squares, existing faces stay but are elevated, and
# vertices become new faces
# TODO: rhombic(smallStellatedDodecahedron) fails!
rhombic <- function(p, name=paste("rhombic", p$name), debug=F)
{
  # construct new coordinates from the edges of the current solid
  newFace <- list()
  for (e in seq(nrow(p$edgeToFaces))) {
    if (debug) print(e)
    F1 <- p$edgeToFaces[e,1]
    F2 <- p$edgeToFaces[e,2]
    pts <- ((which(e == p$coordPairToEdge) - 1) %% ncol(p$coordPairToEdge))+1
    P1 <- pts[1]
    P2 <- pts[2]
    
    if (debug) {
      drawPolygon(p$faces[[F1]], p$coords, drawlines = T, drawvertices = T, label=paste0("F",F1))
      drawPolygon(p$faces[[F2]], p$coords, drawlines = T, drawvertices = T, label=paste0("F",F2))
    }
    
    n1 <- normal(p$coords[p$faces[[F1]][1],],
                 p$coords[p$faces[[F1]][2],],
                 p$coords[p$faces[[F1]][3],])
    
    n2 <- normal(p$coords[p$faces[[F2]][1],],
                 p$coords[p$faces[[F2]][2],],
                 p$coords[p$faces[[F2]][3],])
    
    if (debug) {
      center1 <- apply(p$coords[p$faces[[F1]],], 2, mean)
      center2 <- apply(p$coords[p$faces[[F2]],], 2, mean)
      drawLines(c(center1, center1+n1), color = "red")
      drawLines(c(center2, center2+n2), color = "red")
    }
    
    dp1p2 <- distance(p$coords[P1,], p$coords[P2,])
    dn1n2 <- distance(n1, n2)
    alpha <- dp1p2/dn1n2
    
    p1_a <- p$coords[P1,] + alpha*n1
    p1_b <- p$coords[P1,] + alpha*n2
    p2_a <- p$coords[P2,] + alpha*n2
    p2_b <- p$coords[P2,] + alpha*n1
    
    if (debug) {
      drawLines(c(p1_a, p1_b, p2_a, p2_b, p1_a), color="green")
    }
    
    newFace[[1+length(newFace)]] <- 
      data.table(x=c(p1_a[1], p1_b[1], p2_a[1], p2_b[1]),
                 y=c(p1_a[2], p1_b[2], p2_a[2], p2_b[2]),
                 z=c(p1_a[3], p1_b[3], p2_a[3], p2_b[3]),
                 edge=e,
                 frompoint=c(P1, P1, P2, P2), 
                 fromface=c(F1, F2, F2, F1))
  }
  
  nc <- rbindlist(newFace) # but: contains duplicate coordinates
  nc[, idx := seq(.N)]
  nc[, idx := idx[1], by=c("fromface","frompoint")]
  nc[, idx := as.integer(factor(idx))] # same coordinates now have same index
  
  # faces from the edges
  squaresFromEdges <- nc[,.(f = list(idx)),by=edge]
  
  # de-duplicated coordinates
  newCoords <- as.matrix( nc[, .(x = x[1], y = y[1], z = z[1]), by=idx][order(idx)][, 2:4] )
  
  # find new coords for existing faces
  newTopo <- unique(nc[, c("frompoint","fromface","idx")])
  
  # construct new faces mostly by lookup
  facesFromOldFaces <-
    lapply(seq(length(p$faces)), 
           function(oldFace) {
             return(sapply(p$faces[[oldFace]], 
                           function(i) { return(newTopo[fromface==oldFace & frompoint==i]$idx)}))})
  
  facesFromOldVertices <-
    lapply(seq(length(p$vertexFigures)),
           function(oldVertex) { 
             v <- p$vertexFigures[[oldVertex]]
             return(sapply(v$faces, function(i) {return(newTopo[fromface==i & frompoint==v$center]$idx)}))})
  
  newRhombic <-
    setPoly(coords = newCoords/max(vectorlength(newCoords)), # unit length
            faces = c(squaresFromEdges$f, facesFromOldFaces, facesFromOldVertices),
            name = name, debug=debug)
  
  return(newRhombic)
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
  
  return(setPoly(coords = rbind(p1$coords, p2$coords),
                 faces = c(p1$faces, lapply(p2$faces, function(f) { return(f+nrow(p1$coords)) })),
                 name,
                 debug=debug))
}

# Create new polyhedron by chopping off the vertices replacing each by a new face
# but keeping existing faces intact, in effect doubling their number of vertices
# TODO lin alg to find new points is not ok yet
# TODO make use of topology
# truncate(quasi(cube)) errors out but should give something even if not regular
# truncate <- function(p, name = paste("truncate", p$name), debug=F)
# {
#   # every vertex becomes a new face with all new points
#   allPoints <- NULL
#   allFaces <- list()
#   for (v in p$vertexFigures) {
#     # create new points close to vertex center C in direction of the connected points P
#     # new point = C + alpha*(PC)
#     angles <- innerAngles(p$coords[v$vex,], center=p$coords[v$center,])
#     # find alpha such that the sides of the newly created faces are equal
#     alpha <- 0.5 - 0.5*(sin(angles/2)/(1+sin(angles/2))) # not trivial but easily derived
#     
#     newPoints <- p$coords[v$vex,]*alpha + (1-alpha)*data.table(x=rep(p$coords[v$center,1], length(v$vex)), 
#                                                                y=rep(p$coords[v$center,2], length(v$vex)),
#                                                                z=rep(p$coords[v$center,3], length(v$vex)))
#     newPoints$from <- v$center
#     newPoints$to <- v$vex
#     if (is.null(allPoints)) {
#       allPoints <- newPoints
#     } else {
#       allPoints <- rbind(allPoints, newPoints)
#     }
#     allFaces[[v$center]] <- (nrow(allPoints)-nrow(newPoints)+1):nrow(allPoints)
#   }
#   # now transform the old faces
#   newPointsLookup <- matrix(data = NA, nrow = nrow(allPoints), ncol = nrow(allPoints)) # sparse?
#   for (i in seq(nrow(allPoints))) {
#     # TODO maybe check for inconsistency if there's a value already
#     newPointsLookup[allPoints$from[i], allPoints$to[i]] <- i
#   }
#   for (f in p$faces)
#   {
#     truncatebedFace <- as.vector(sapply(seq(length(f)), function(idx) {
#       return(c(newPointsLookup[f[idx], shiftrotate(f)[idx]], 
#                newPointsLookup[shiftrotate(f)[idx], f[idx]]))}))
#     allFaces[[length(allFaces)+1]] <- truncatebedFace
#   }
#   
#   pTruncate <- setPoly(coords = allPoints[,1:3], faces = allFaces, name=name, debug=debug)
#   
#   return(pTruncate)
# }

# clear3d()
# drawSinglePoly(truncate(cube), debug=T)

name <- function(p)
{
  return(p$name)
}

description <- function(p, debug=F)
{
  if ("original" %in% names(p)) p <- p$original
  
  combineDescriptions <- function(descrs, suffixIdentical = "")
  {
    descrFreqs <- data.table(table(unlist(descrs)), stringsAsFactors = F)[order(-N,V1)]
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
      #print(descrFreqs)
      return(paste(sapply(seq(nrow(descrFreqs)), 
                          function(i) {return(paste0(descrFreqs$N[i], "x", descrFreqs$V1[i]))}), collapse=" + "))
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
      if (!isFlatFace(p$coords[f,])) {
        # not a real polygon, not flat!
        fDescr <- paste0("!",baseDescr)
      } else {
        # star polygon
        fDescr <- paste(baseDescr, sum(angles)/(2*pi), sep="/") # TODO maybe round
      }
    }
    #cat("Face:",f,fDescr,fill=T)
    return(fDescr)
  }
  
  getVertexDescription <- function(aVertex)
  {
    faceDescriptions <- sapply(aVertex$faces, function(fi) { 
      if (fi == 0) return("GAP")
      else return(getFaceDescription(p$faces[[fi]])) })
    vexDescriptions <- getFaceDescription(aVertex$vex, addlengths=F)
    
    if (length(unique(faceDescriptions)) == 1 & length(unique(vexDescriptions)) == 1) {
      vexDescription <- paste0("{",faceDescriptions[1],",",vexDescriptions[1],"}")
    } else {
      # create rotate variations then find the one with lowest sort order
      faceDescriptionRotations <- sapply(1:length(faceDescriptions), function(i) {
        return(paste(shiftrotate(faceDescriptions,i-1),collapse=","))
        })
      vexDescription <- paste0("{",(sort(faceDescriptionRotations))[1],"}") 
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
  
  bodyDescriptions <- lapply(p$bodies, function(b) {
    # subselect vertices that have faces in current body
    getBodyDescription(p$vertexFigures[sapply(p$vertexFigures, function(t) {return(length(intersect(t$faces, b)) > 0 )})])
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


dodecahedron <- dual(icosahedron, name = "Dodecahedron")

greatDodecahedron <- buildRegularPoly(coords = icosahedron$coords, 
                                      polygonsize = 5, vertexsize = 5, exampleEdge = c(1,6),
                                      name = "Great Dodecahedron")
smallStellatedDodecahedron <- buildRegularPoly(icosahedron$coords,
                                               polygonsize = 5,
                                               vertexsize = 5,
                                               exampleEdge = c(1,7),
                                               name = "Small Stellated Dodecahedron")
greatIcosahedron <- buildRegularPoly(icosahedron$coords,
                                     polygonsize = 3,
                                     vertexsize = 5,
                                     exampleEdge = c(2, 6),
                                     name = "Great Icosahedron")
greatStellatedDodecahedron <- dual(greatIcosahedron, name = "Great Stellated Dodecahedron", 
                                   scaling = "vertex")


Platonics <- list(tetrahedron, octahedron, cube, icosahedron, dodecahedron)
KeplerPoinsots <- list(greatDodecahedron, smallStellatedDodecahedron, greatIcosahedron, greatStellatedDodecahedron)
Regulars <- c(Platonics, KeplerPoinsots)


testDescription <- function()
{
  description(octahedron)
  description(cube)
  description(smallStellatedDodecahedron)
  description(quasi(cube))
  description(compose(cube, octahedron))
  
  sapply(lapply(Regulars, dual), description)
}

testRhombic <- function()
{
  p <- dodecahedron
  p <-greatDodecahedron
  p <- icosahedron
  clear3d()
  drawPoly(p, debug = T)
  
  drawInit(new.device = T)
  clear3d()
  drawAxes()
  r <- rhombic(p, debug=F)
  drawPoly(r)
}

kerst2019 <- function()
{
  drawInit(new.device = T, width = 1536, height = 1536/sqrt(2))
  rgl.viewpoint(zoom=0.2)
  compound5tetrahedra <- buildRegularPoly(dodecahedron$coords,
                                          polygonsize = 3,
                                          vertexsize = 3,
                                          exampleEdge = c(3, 8),
                                          name = "5 Tetrahedra")
  stars <- list(greatStellatedDodecahedron, greatIcosahedron,
                 smallStellatedDodecahedron,
                 quasi(greatStellatedDodecahedron),
                 compose(greatStellatedDodecahedron, dual(greatStellatedDodecahedron)),
                 compound5tetrahedra, 
                 quasi(compound5tetrahedra))
  
  xmasCols <- function(n, alpha=1)
  {
    return(sample(c("darkgreen","red","gold"), n, replace=T))
  }
  rgl.bg(back = "filled", color="black")
  clear3d()
  set.seed(1)
  sizeOfUniverse <- 40
  nStarsInUniverse <- 5
  drawTexts(c(0,sizeOfUniverse,0), "Een goed 2020!", color="white")
  for (i in seq(nStarsInUniverse)) {
    print(i)
    drawSinglePoly(stars[[round(runif(1,1,length(stars)))]],
                   offset = runif(3,-sizeOfUniverse,sizeOfUniverse), 
                   colorProvider = xmasCols, label = "")
  }
  snapshot3d("kerst2019.png", fmt = "png", top = TRUE )
  rgl.postscript("kerst2019.svg", "svg", drawText = T) 
  rgl.postscript("kerst2019.pdf", "pdf", drawText = T) 
}
