# automatic layout gen

# TODO
# build up layout in batches? eg optimize first 6 then next etc hopefully gives quicker convergence
# store 2d projection with the 3D solid (??)
# but do add edgeToVertex Nx2 matrix - there already is ege to faces

library(ggplot2)
library(svglite)
library(data.table)

source("polyhedra.R")
source("draw.R")

# draw 2D face
## TODO reasonably obsolete now that the layout is plotted w ggplot
drawFace2D <- function(face2d, outline = TRUE, color = "green")
{
  if (outline) {
    #points3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, color="black" )
    text3d(x = face2d$coords2D[,1], y = face2d$coords2D[,2], z = 0, 
           text=face2d$points3D, color="black", cex=0.5 )
    # if (!is.null(edges)) {
    #   text3d(x = (face2d$coords2D[,1] + face2d$coords2D[shiftrotate(seq(nrow(face2d$coords2D))),1])/2, 
    #          y = (face2d$coords2D[,2] + face2d$coords2D[shiftrotate(seq(nrow(face2d$coords2D))),2])/2, z = 0, 
    #          text=edges, color="black", cex=0.7 )
    #   
    # }
    lines3d(x = face2d$coords2D[,1][c(seq(nrow(face2d$coords2D)),1)], 
            y = face2d$coords2D[,2][c(seq(nrow(face2d$coords2D)),1)], z = 0, color = color)
    center <- apply(face2d$coords2D, 2, mean)
    text3d(x = center[1], y = center[2], z = 0, text=paste0("F",face2d$face), color="blue", cex=0.5 )
  } else {
    if (nrow(face2d$coords2D) < 3) stop("Face to triangulate should have >= 3 points.")
    
    for (t in (seq(nrow(face2d$coords2D)-2)+1)) {
      triangle <- face2d$coords2D[c(1, t, t+1),]
      triangles3d( triangle[,1], triangle[,2], 0, col=color)
    }
    lines3d(x = face2d$coords2D[,1][c(seq(nrow(face2d$coords2D)),1)], 
            y = face2d$coords2D[,2][c(seq(nrow(face2d$coords2D)),1)], z = 0, color = "black")
  }
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

#' Put a new face into position, given a global 2D layout and and edge to connect to
#'
#' @param polyhedron3D The 3D solid with full topology
#' @param edge Index of the edge to connect to
#' @param placedFaces Array of all faces from the 3D solid that are currently in the layout already
#' @param level Current level in the layout (to index the global layout2D structure)
#' @param trace Debug flag
#'
#' @return Fully defined 2D face
positionNextFace <- function(polyhedron3D, edge, placedFaces, level, trace=T) 
{
  if (polyhedron3D$edgeToFaces[edge,1] %in% placedFaces) {
    placedFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
    newFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
  } else {
    newFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
    placedFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
  }
  vertices <- ((which(polyhedron3D$coordPairToEdge==edge)-1) %/% nrow(polyhedron3D$coordPairToEdge))+1
  if (trace) {cat(paste0(rep(" ", level),collapse=""), "finding next face for edge", edge, "and existing face", placedFaceIdx, ": adding face", newFaceIdx, "with shared vertices:", paste(vertices, collapse=","), fill = T)}
  
  # place the to be added face on 2D plane
  face2D <- projectFace(polyhedron3D, newFaceIdx)
  #drawFace2D(face2D, color="green")
  
  # lookup the one that is placed already
  currentFace2D <- layout2D[[which(sapply(layout2D[1:level], function(x){return(x$face)})==placedFaceIdx)]]
  
  # identify the vertices of the edge in both faces
  if (length(vertices) != 2) stop("An edge must have 2 vertices")
  p1current <- which(currentFace2D$points3D == vertices[1])
  p1other <- which(face2D$points3D == vertices[1])
  p2current <- which(currentFace2D$points3D == vertices[2])
  p2other <- which(face2D$points3D == vertices[2])
  if (length(c(p1current,p1other,p2current,p2other)) != 4) stop("All vertices must occur in both faces")
  
  # rotate and translate the connected face
  vecCurrent <- as.numeric(currentFace2D$coords2D[p2current,]) - as.numeric(currentFace2D$coords2D[p1current,])
  vecOther <- as.numeric(face2D$coords2D[p2other,]) - as.numeric(face2D$coords2D[p1other,])
  theta <-  vectorAngle2D(vecCurrent, vecOther)
  rotation <- rotationMatrix2D(theta)
  
  # transform points P of other face
  # as: p1current + rotation (P - p1other)
  moved <- sweep(face2D$coords2D, 2, as.numeric(face2D$coords2D[p1other,]))
  rotated <- t(rotation %*% t(moved))
  translated <- sweep(rotated, 2, as.numeric(currentFace2D$coords2D[p1current,]), "+")
  face2D$coords2D <- translated
  
  # TODO: enrich the face with more attributes
  # - coords2D = matrix
  # - min/max cooordinates
  # - center
  # - connectionEdge (was sourceEdge)
  # - connectionPoints
  # - faceReference (was face)
  # - vertexReference (was "points3D")
  # so distinguish from layout attributes like running min/max and candidateEdges
  return(face2D)
}

drawLayout <- function(level, debug=T, original)
{
  colors <- assignColors(original)
  
  plotdata <- rbindlist(lapply(layout2D[1:level], function(l) {
    data.table(x=l$coords2D[,1], y=l$coords2D[,2], face=l$face, color=colors[l$face], vex=xx$faces[[l$face]])
  }))
  colorMapIdentity <- unique(colors)
  names(colorMapIdentity) <- colorMapIdentity
  
  # should be part of the data already
  centers <- rbindlist(lapply(seq(level), function(l) {
    center <- apply(layout2D[[l]]$coords2D,2,mean)
    return(data.table(x=center[1],y=center[2],level=l,face=layout2D[[l]]$face))
  }))
  
  # TODO consider adding min/max lines linetype dotted
  p <-ggplot(plotdata, aes(x,y)) + 
    geom_polygon(aes(fill = color, group = face), color="black", size=0.1, show.legend = FALSE) +
    scale_fill_manual(values = colorMapIdentity)+
    scale_x_continuous(limits = c(min(layout2D[[level]]$minCoords), max(layout2D[[level]]$maxCoords)))+
    scale_y_continuous(limits = c(min(layout2D[[level]]$minCoords), max(layout2D[[level]]$maxCoords)))+
    #theme_void()+
    ggtitle(paste("Layout of", original$name), 
            subtitle = paste("round",ncalls,"eval=",round(bestEval,5),round(difftime(Sys.time(), startTime, units = "mins"),2),"mins"))
  
  if (debug) {
    p <- p + geom_label(mapping=aes(label=vex), size=3, color="black", alpha=0.6)+
      geom_text(data=centers, mapping=aes(x,y,label=paste0(level," (F",face,")")), inherit.aes = F, size=3)
  } else {
    p <- p + theme_void()
  }
  
  print(p)
  
  ggsave(file=file.path("layouts",paste0("layout ", original$name, ".svg")), plot=p, width=10, height=10)
  # 
  # rgl_init(new.device = T)
  # rgl.viewpoint(theta = 0, phi = 0, zoom = 1)
  # for (l in layout2D[1:level]) { 
  #   if (!is.null(colors)) {
  #     drawFace2D(l, outline, colors[l$face])    
  #   } else {
  #     drawFace2D(l, outline) 
  #   }
  # }
  # if (outline) {
  #   lines3d( x = c(layout2D[[level]]$minCoords[1], layout2D[[level]]$maxCoords[1], layout2D[[level]]$maxCoords[1], layout2D[[level]]$minCoords[1], layout2D[[level]]$minCoords[1]),
  #            y = c(layout2D[[level]]$minCoords[2], layout2D[[level]]$minCoords[2], layout2D[[level]]$maxCoords[2], layout2D[[level]]$maxCoords[2], layout2D[[level]]$minCoords[2]),
  #            z = 0,
  #            color="yellow")
  # }
}

getLayoutDigest <- function(level, candidate=NULL)
{
  if (is.null(candidate)) {
    digest <- paste(sort(sapply(layout2D[1:level], function(f) {return(f$sourceEdge)})),collapse="-")
  } else {
    digest <- paste(sort(c(sapply(layout2D[1:level], function(f) {return(f$sourceEdge)}),candidate)),collapse="-")
  }
  return(digest)
}

edgeDescr <- function(polyhedron3D, e)
{
  return (paste0("(", paste0(sort(((which(polyhedron3D$coordPairToEdge==e)-1) %% nrow(polyhedron3D$coordPairToEdge))+1),collapse="-"), ")"))  
}

# Find the best 2D layout given global layout
best2DLayout <- function(polyhedron3D, level = 1, evaluator, trace=T, face2D = projectFace(polyhedron3D, 1), sourceEdge = NA, maxrounds = 0)
{
  if (level == 1) 
  {
    layout2D <<- list()
    allLayoutDigests <<- list()
    bestEval <<- Inf
    bestDigest <<- ""
    ncalls <<- 0
    nevalskips <<- 0
    ndigestskips <<- 0
  }
  
  if (maxrounds < 0 & bestDigest != "") return() # stop when there is 1 solution
  if (((ncalls+1) > maxrounds) & (maxrounds > 0)) {
    return()
  }
  ncalls <<- ncalls+1
  
  #
  # Update layout with new face2D
  #
  
  face2D$sourceEdge <- sourceEdge # strange to do here - should be done when creating face2d
  layout2D[[level]] <<- face2D
  allLayoutDigests[[1+length(allLayoutDigests)]] <<- getLayoutDigest(level)
  
  # update the list of available edges, and store it with the layout itself
  shiftFacePoints <- shiftrotate(face2D$points3D)
  newEdges <- sapply(seq(length(face2D$points3D)), function(i) {return(polyhedron3D$coordPairToEdge[face2D$points3D[i],shiftFacePoints[i]])})
  if (level > 1) {
    layout2D[[level]]$candidateEdges <<-
      setdiff(c(newEdges, layout2D[[level-1]]$candidateEdges), intersect(layout2D[[level-1]]$candidateEdges, newEdges))
  } else {
    layout2D[[level]]$candidateEdges <<- newEdges
  }
  
  # keep track of the layout min/max
  newMinCoords <- apply(face2D$coords2D,2,min)
  newMaxCoords <- apply(face2D$coords2D,2,max)
  if (level > 1) {
    layout2D[[level]]$minCoords <<- apply(matrix(c(newMinCoords, layout2D[[level-1]]$minCoords), ncol=2, byrow = T), 2, min)
    layout2D[[level]]$maxCoords <<- apply(matrix(c(newMaxCoords, layout2D[[level-1]]$maxCoords), ncol=2, byrow = T), 2, max)
  } else {
    layout2D[[level]]$minCoords <<- newMinCoords
    layout2D[[level]]$maxCoords <<- newMaxCoords
  }
  
  if (trace) {
    cat(paste0(rep(" ", level),collapse=""),
        "Entered level", level, "placed face:", face2D$face, 
        "along edge", sourceEdge, edgeDescr(polyhedron3D, sourceEdge), 
        "eval =", evaluator(level), fill = T)
    drawLayout(level=level, debug = T, polyhedron3D)
    invisible(readline(prompt="Press [enter] to continue"))
  }
  
  #
  # Are we done?
  #
  
  if (evaluator(level) > bestEval) {
    if (trace) {cat(paste0(rep(" ", level),collapse=""),"*** early exit **** at level", level,"eval=", evaluator(level), "but best is", bestEval, fill = T)}
    nevalskips <<- nevalskips+1
    return()
  }
  
  #
  # Get candidates for next "move"
  #
  
  # get list of all faces placed
  placedFaces <- sapply(layout2D[1:level], function(f) {return(f$face)})
  
  if (level == 1) {
    # for first level don't go through all the edges
    candidates <- layout2D[[level]]$candidateEdges[1]
  } else {
    candidates <- layout2D[[level]]$candidateEdges
  }
  if (length(candidates) == 0) {
    if (evaluator(level) < bestEval & !deltaEquals(evaluator(level), bestEval)) {
      bestEval <<- evaluator(level)
      bestDigest <<- getLayoutDigest(level)
      if (trace) {cat(paste0(rep(" ", level),collapse=""),"Done with better evaluation! At level", level, "eval=", evaluator(level), fill = T)  }
      if (!trace) {cat("Found layout, round",ncalls,"eval=",round(bestEval,5),round(difftime(Sys.time(), startTime, units = "mins"),5),"mins", fill = T)  }
      drawLayout(level, debug=F, original = polyhedron3D)
    } else {
      if (trace) {cat(paste0(rep(" ", level),collapse=""),"Done but no improvement. At level", level, "eval=", evaluator(level), fill = T)  }
    }
  } else {
    candidates <- candidates[sample.int(length(candidates))] # shuffle to speed up search
    
    # first find the 2D faces
    candidateFaces <- list()
    for (edge in candidates) {
      layoutDigest <- getLayoutDigest(level, edge)
      if (layoutDigest %in% allLayoutDigests) {
        if (trace) {cat(paste0(rep(" ", level),collapse=""),"*** skip processing edge", edge, " **** at level", level,"layout done already:", layoutDigest, fill = T)}
        ndigestskips <<- ndigestskips+1
        next
      }
      
      if (trace) {cat(paste0(rep(" ", level),collapse=""),"select edge:", edge, edgeDescr(polyhedron3D, edge), "from", candidates, fill = T)}
      
      # place it into position
      newFace2D <- positionNextFace(polyhedron3D, edge, placedFaces, level, trace) 
      candidateFaces[[1+length(candidateFaces)]] <- list(edge=edge, face=newFace2D) # should be no need for sep struct
    }
    
    # then order them by (early) evaluation
    if (length(candidateFaces) > 0) {
      for (c in seq(length(candidateFaces))) {
        candidateMinCoords <- apply(candidateFaces[[c]]$face$coords2D,2,min) # should be part of the face already
        candidateMaxCoords <- apply(candidateFaces[[c]]$face$coords2D,2,max)
        
        currentLayoutWithCandidateMin <- apply(matrix(c(candidateMinCoords, layout2D[[level]]$minCoords), ncol=2, byrow = T), 2, min)
        currentLayoutWithCandidateMax <- apply(matrix(c(candidateMaxCoords, layout2D[[level]]$maxCoords), ncol=2, byrow = T), 2, max)
        
        candidateFaces[[c]]$eval <- evaluator(level, currentLayoutWithCandidateMin, currentLayoutWithCandidateMax)
      }
      candidateOrder <- order(sapply(candidateFaces,function(x){return(x$eval)}))
      
      # then recurse, starting with the one with smallest local eval result
      for (c in candidateOrder) {
        best2DLayout(polyhedron3D, level+1, evaluator, trace=trace, candidateFaces[[c]]$face, candidateFaces[[c]]$edge, maxrounds = maxrounds)
      }
    }
  }
  
  if (trace) {cat(paste0(rep(" ", level),collapse=""),"Returning from level", level, fill = T)}
}

layoutAreaEvaluator <- function(level, min=layout2D[[level]]$minCoords, max=layout2D[[level]]$maxCoords)
{
  # area - the smaller the better
  rslt <- ((max[1]-min[1])*(max[2]-min[2]))
  # # prefer horizontal by penalizing vertical orientation 
  # if ((max[2]-min[2]) > (max[1]-min[1])) rslt <- rslt*(1+1/(4*level))
  # this really slows down convergence - better swap at the end
  return(rslt)
}

layoutSquarenessEvaluator <- function(level, min=layout2D[[level]]$minCoords, max=layout2D[[level]]$maxCoords)
{
  # abs diff x y - the closer the better
  return (abs((max[1]-min[1]) - (max[2]-min[2])))
}

########

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- icosahedron
xx <- tetrahedron
xx <- rhombic(cube)
xx <- truncate(icosahedron)
xx <- quasi(cube)
xx <- rhombic(dodecahedron)
xx <- dodecahedron
#xx <- cube
#clear3d()
#drawPoly(xx, debug=T)

# for layout use a new device
# rgl_init(new.device = T)
# clear3d()
# drawAxes()

# Start with a layout with just the first face
startTime <- Sys.time()

best2DLayout(xx, evaluator = layoutAreaEvaluator, trace=F, maxrounds=15000)

endTime <- Sys.time()
cat("#calls:", ncalls, fill=T)
cat("#early skips because of repeated layout", ndigestskips, fill=T)
cat("#early skips because of early evaluation", nevalskips, fill=T)
print(bestEval)
print(bestDigest) # layout can be reconstructed from this
cat("Elapsed:", as.double(difftime(endTime, startTime, units = "mins")), "mins", fill=T)

stop()

# test overlap
set.seed(125)
xx <- rhombic(dodecahedron)
xx<-quasi(octahedron)
xx<- truncate(dodecahedron)
best2DLayout(xx, evaluator = layoutSquarenessEvaluator, trace=T, maxrounds=-1)
# overlap of #45 and #54, and #15 and #38

isBetween <- function(m, x1, x2)
{
  return ((m > x1) & (m < x2) & !deltaEquals(m ,x1) & !deltaEquals(m, x2))
}
isBoundingBoxOverlap <- function(fig1, fig2)
{
  # TODO keep these min/max and also centers with the placed faces
  fig1Min <- apply(fig1$coords2D,2,min)
  fig1Max <- apply(fig1$coords2D,2,max)
  fig2Min <- apply(fig2$coords2D,2,min)
  fig2Max <- apply(fig2$coords2D,2,max)
  
  if (fig1Max[1] < fig2Min[1] | deltaEquals(fig1Max[1], fig2Min[1])) return(F)
  if (fig1Max[2] < fig2Min[2] | deltaEquals(fig1Max[2], fig2Min[2])) return(F)
  if (fig1Min[1] > fig2Max[1] | deltaEquals(fig1Min[1], fig2Max[1])) return(F)
  if (fig1Min[2] > fig2Max[2] | deltaEquals(fig1Min[2], fig2Max[2])) return(F)
  
  # if (isBetween(fig2Min[1], fig1Min[1], fig1Max[1]) & isBetween(fig2Min[2], fig1Min[2], fig1Max[2])) return(T)
  # if (isBetween(fig2Min[1], fig1Min[1], fig1Max[1]) & isBetween(fig2Max[2], fig1Min[2], fig1Max[2])) return(T)
  # if (isBetween(fig2Max[1], fig1Min[1], fig1Max[1]) & isBetween(fig2Min[2], fig1Min[2], fig1Max[2])) return(T)
  # if (isBetween(fig2Max[1], fig1Min[1], fig1Max[1]) & isBetween(fig2Max[2], fig1Min[2], fig1Max[2])) return(T)
  
  return(T)
}
isEdgeConnectedInLayout <- function(fig1, fig2)
{
  edges1 <- which(xx$edgeToFaces[,2]==fig1$face | xx$edgeToFaces[,1]==fig1$face)
  edges2 <- which(xx$edgeToFaces[,2]==fig2$face | xx$edgeToFaces[,1]==fig2$face)
  if ((fig1$sourceEdge %in% edges2) | (fig2$sourceEdge %in% edges1)) return(T)
  
  return(F)
}
isBoundingCircleOverlap <- function(fig1, fig2)
{
  center1 <- apply(fig1$coords2D,2,mean)
  center2 <- apply(fig2$coords2D,2,mean)
  maxRadius1 <- max(distance(fig1$coords2D,center1))
  maxRadius2 <- max(distance(fig2$coords2D,center2))
  d <- distance(center1, center2)
  if ((maxRadius1+maxRadius2 > d) & !deltaEquals(maxRadius1+maxRadius2, 2)) return(T)
  return(F)
}
# TODO this may not belong here
getSingleSharedPoint <- function(fig1, fig2)
{
  # is here a single shared point P between them?
  p1 <- -1
  p2 <- -1
  for (i in seq(nrow(fig1$coords2D))) {
    for (j in seq(nrow(fig2$coords2D))) {
      if (deltaEquals(0, distance(fig1$coords2D[i,], fig2$coords2D[j,]))) {
        if (p1 != -1) return(NULL) # more than one shared point
        p1 <- i
        p2 <- j
      }
    }
  }
  if (p1 == -1) return(NULL) # no shared point
  
  originalP <- xx$faces[[fig1$face]][p1]
  if (originalP != xx$faces[[fig2$face]][p2]) stop("Topology error")
  return(originalP)
}
# TODO this does not belong here. When connecting a face, a check should be done
# that the angles of all placed connected faces of the vertex figure are < 2pi.
isSingleSharedPointOK <- function(fig1, fig2, originalP)
{
  # Figure out which faces share that point in 2D
  p1 <- which(xx$faces[[fig1$face]] == originalP)
  vertex <- which( sapply(xx$vertexFigures, function(v) {return(v$center == originalP)}))
  if (length(vertex) != 1) stop("Topology error")
  connectedFacesIn2D <- which(sapply(layout2D, function(l) {
    if (!(l$face %in% xx$vertexFigures[[vertex]]$faces)) return(F)
    c <- which(xx$faces[[l$face]] == originalP)
    if (length(c) != 1) return(F)
    return(deltaEquals(0, distance(fig1$coords2D[p1,], l$coords2D[c,])))
  }))
  
  # is their total angle around P > 2pi?
  # sum angles connectedFacesIn2D around originalP
  angles <- sapply(layout2D[connectedFacesIn2D], function(l) {
    prevF <- shiftrotate(xx$faces[[l$face]], -1)
    nextF <- shiftrotate(xx$faces[[l$face]])
    c <- which(xx$faces[[l$face]] == originalP)
    return(vectorAngle(xx$coords[xx$faces[[l$face]][c],] - xx$coords[prevF[c],], 
                       xx$coords[xx$faces[[l$face]][c],] - xx$coords[nextF[c],]))
  })
  if (sum(angles) < 2*pi) return(T)
  
  return(F) # no overlap
}

layout2D <<- list()
allLayoutDigests <<- list()

combis <- combn(length(layout2D),2,simplify=F)
n<-0
for (facePair in combis) {
  if (isBoundingBoxOverlap(layout2D[[facePair[1]]], layout2D[[facePair[2]]])) {
    if (!isEdgeConnectedInLayout(layout2D[[facePair[1]]], layout2D[[facePair[2]]])) {
      if (isBoundingCircleOverlap(layout2D[[facePair[1]]], layout2D[[facePair[2]]])) {
        p <- getSingleSharedPoint(layout2D[[facePair[1]]], layout2D[[facePair[2]]])
        if (is.null(p)) {
          cat("Possible overlap:", facePair[1], "and", facePair[2], fill=T)
          n<-n+1
        } else {
          # TODO this does not belong here
          if (!isSingleSharedPointOK(layout2D[[facePair[1]]], layout2D[[facePair[2]]], p)) {
            cat("Overlap:", facePair[1], "and", facePair[2], fill=T)
            n<-n+1
          } else {
            # no overlap
            cat("NO Overlap:", facePair[1], "and", facePair[2], fill=T)
          }
        }
      }
    }
  }
}
print(n)

# nog niet goed - truncate octahedron mist overlap van bv 22/29
# toch snijlijnen check ? single shared p is erg ingewikkeld wel nodig voor deuken

