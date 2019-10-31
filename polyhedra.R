# Functions specifically to constructing polyhedra

library(data.table)

source("math.R")
source("draw.R")
source("topology.R")

# Returns T if outward facing (= rotating anti-clockwise when looking down towards face in direction of origin)
isNormalOutwardFacing <- function(coords, f)
{
  n <- normal( coords[f[1],], coords[f[2],], coords[f[3],] )
  mid <- apply(coords[f,],2,mean)
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
    if (!isNormalOutwardFacing(p$coords, aFace)) {
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
setPoly <- function(coords, faces, name, debug=F)
{
  # make sure all faces are oriented consistently
  for (i in seq(length(faces))) {
    if (!isNormalOutwardFacing(coords, faces[[i]])) {
      if (debug) cat("Flip face", i, "to make normal outward facing", fill=T)
      faces[[i]] <- rev(faces[[i]])
    }
  }

  newPoly <- list(coords=coords, faces=faces, name=name)
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
  
  newPoly <- setPoly(poly$coords, poly$faces, poly$name, debug)
  
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
  
  pDual <- setPoly(newVertexCoords*scale, newFaces, name, debug)

  return(pDual)
}

# Create an derived polyhedron by truncating all coords to the mid of the faces
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
  
  pArchi <- setPoly(archiPoints, c(archiFaces1, archiFaces2), name, debug)

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
  
  return(setPoly(coords = rbind(p1$coords, p2$coords),
              faces = c(p1$faces, lapply(p2$faces, function(f) { return(f+nrow(p1$coords)) })),
              name,
              debug))
}

# Create new polyhedron by chopping off the vertices replacing each by a new face
# but keeping existing faces intact, in effect doubling their number of vertices
# TODO lin alg to find new points is not ok yet
# TODO make use of topology
# truncate(quasi(cube)) errors out but should give something even if not regular
truncate <- function(p, name = paste("truncate", p$name), debug=F)
{
  # every vertex becomes a new face with all new points
  allPoints <- NULL
  allFaces <- list()
  for (v in p$vertexFigures) {
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
    truncatebedFace <- as.vector(sapply(seq(length(f)), function(idx) {
      return(c(newPointsLookup[f[idx], shiftrotate(f)[idx]], 
               newPointsLookup[shiftrotate(f)[idx], f[idx]]))}))
    allFaces[[length(allFaces)+1]] <- truncatebedFace
  }
  
  pTruncate <- setPoly(coords = allPoints[,1:3], faces = allFaces, name=name, debug)

  return(pTruncate)
}

# clear3d()
# drawSinglePoly(truncate(cube), debug=T)

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
    faceDescriptions <- sapply(aVertex$faces, function(fi) { 
      if (fi == 0) return("GAP")
      else return(getFaceDescription(p$faces[[fi]])) })
    vexDescriptions <- getFaceDescription(aVertex$vex, addlengths=F)
    
    if (length(unique(faceDescriptions)) == 1 & length(unique(vexDescriptions)) == 1) {
      vexDescription <- paste0("{",faceDescriptions[1],",",vexDescriptions[1],"}")
    } else {
      # shiftrotate them so we consistently take the same one
      descrs <- sapply(1:length(faceDescriptions), function(i) {return(paste(shiftrotate(faceDescriptions,i+1),collapse=","))})
      vexDescription <- paste0("{",(sort(descrs))[1],"}") 
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

# description(octahedron)
# description(cube)
# description(smallStellatedDodecahedron)
# description(quasi(cube))
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


stop("new trunc")

rgl_init()
clear3d()
drawAxes()
xx <- quasi(cube) # pooff
xx <- greatDodecahedron
drawPoly(xx, debug=T)

rgl_init(new.device = T)
clear3d()
drawAxes()

newpts <- list()
for(v in xx$vertexFigures) 
{
  text3d( x=xx$coords[v$center,1], 
          y=xx$coords[v$center,2], 
          z=xx$coords[v$center,3], color="red", texts = v$center)
  for (p in v$vex) 
  {
    lines3d( x=xx$coords[c(v$center, p),1], 
             y=xx$coords[c(v$center, p),2], 
             z=xx$coords[c(v$center, p),3], color="blue")
    text3d( x=xx$coords[p,1], 
            y=xx$coords[p,2], 
            z=xx$coords[p,3], color="red", texts = p)
  }
  # new point from center of vertex figure is on line to connected point at relative distance alpha
  alpha <- sapply(seq(length(v$vex)), 
                  function(i) {return(1/(2+vectorlength(xx$coords[v$vex[i],]-xx$coords[shiftrotate(v$vex)[i],])/
                                           vectorlength(xx$coords[v$vex[i],]-xx$coords[v$center,])))})
  newpts[[v$center]] <- as.data.table(t(sapply(seq(length(v$vex)),
                   function(i) {return((1-alpha[i])*xx$coords[v$center,] + alpha[i]*xx$coords[v$vex[i], ])})))
  newpts[[v$center]]$center <- v$center
  newpts[[v$center]]$connected <- v$vex
  text3d( x=newpts[[v$center]]$x, 
          y=newpts[[v$center]]$y, 
          z=newpts[[v$center]]$z, color="darkgreen", texts = newpts[[v$center]]$connected)
}

newcoords <- rbindlist(newpts)
newcoords[,newcoordidx := seq(.N)]

# the new faces at the old vertices
newfaces <- newcoords[, .( connected = list(connected), coords = list(newcoordidx)) , by=center]

# the old faces require some careful indexing
oldfaces <- list()
for (f in xx$faces) {
  prv <- shiftrotate(f,-1)
  nxt <- shiftrotate(f)
  oldfaces[[1+length(oldfaces)]] <- 
    as.numeric(sapply(seq(length(f)), function(i) {return(c(newcoords[center==f[i]&connected==prv[i]]$newcoordidx,
                                                            newcoords[center==f[i]&connected==nxt[i]]$newcoordidx))}))
}

# done!
trunc_xx <- setPoly(coords = as.matrix(newcoords[,1:3]),
        faces = c(newfaces$coords, oldfaces), name="truncate2")
drawPoly(trunc_xx)

