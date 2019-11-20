# automatic layout gen

# TODO
# build up layout in batches? eg optimize first 6 then next etc hopefully gives quicker convergence
# store 2d projection with the 3D solid (??)
# but do add edgeToVertex Nx2 matrix - there already is ege to faces

library(ggplot2)
library(svglite)
library(data.table)

source("polyhedra.R")
source("geometry2D.R")

#' Put a new face into position, given a global 2D layout and and edge to connect to
#'
#' @param polyhedron3D The 3D solid with full topology
#' @param edge Index of the edge to connect to
#' @param placedFaces Array of all faces from the 3D solid that are currently in the layout already
#' @param level Current level in the layout (to index the global layout2D structure)
#' @param trace Debug flag
#'
#' @return Fully defined 2D face
positionNextFace <- function(polyhedron3D, edge=NA, placedFaces=c(), level=1, debug=T) 
{
  if (length(placedFaces) > 0) {
    if (polyhedron3D$edgeToFaces[edge,1] %in% placedFaces) {
      placedFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
      newFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
    } else {
      newFaceIdx <- polyhedron3D$edgeToFaces[edge,1]
      placedFaceIdx <- polyhedron3D$edgeToFaces[edge,2]
    }
    vertices <- ((which(polyhedron3D$coordPairToEdge==edge)-1) %/% nrow(polyhedron3D$coordPairToEdge))+1
    if (debug) {cat(paste0(rep(" ", level),collapse=""), 
                    "finding next face for edge", edge, 
                    "and existing face", paste0("F",placedFaceIdx), 
                    ":", paste0("F",newFaceIdx), 
                    "with shared vertices:", paste(vertices, collapse=","), fill = T)}
  } else {
    newFaceIdx <- 1 # start at face 1
    placedFaceIdx <- NA
  }  
  
  # project 3D face on 2D using a cache so not need to do this over and over again
  face2D <- NULL
  if (newFaceIdx <= length(face3Dto2DprojectionCache)) {
    face2D <- face3Dto2DprojectionCache[[newFaceIdx]]
  }
  if (is.null(face2D)) {
    face2D <- list( faceReference = newFaceIdx, 
                    vexReferences = polyhedron3D$faces[[newFaceIdx]],
                    coords2D = projectFace(polyhedron3D$coords[polyhedron3D$faces[[newFaceIdx]],]))
    face3Dto2DprojectionCache[[newFaceIdx]] <<- face2D
  }
  #drawFace2D(face2D, color="green")
  
  if (length(placedFaces) > 0) {
    # lookup the one that is placed already
    currentFace2D <- layout2D[[which(sapply(layout2D[1:level], function(x){return(x$faceReference)})==placedFaceIdx)]]
    
    # identify the vertices of the edge in both faces
    if (length(vertices) != 2) stop("An edge must have 2 vertices")
    p1current <- which(currentFace2D$vexReferences == vertices[1])
    p1other <- which(face2D$vexReferences == vertices[1])
    p2current <- which(currentFace2D$vexReferences == vertices[2])
    p2other <- which(face2D$vexReferences == vertices[2])
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
  }
  
  # enrich the face with meta attributes to speed up things later
  face2D$minCoords <- apply(face2D$coords2D,2,min)
  face2D$maxCoords <- apply(face2D$coords2D,2,max)
  face2D$center2D <- apply(face2D$coords2D,2,mean)
  face2D$connectionEdge <- edge
  face2D$connectedToFaceReference <- placedFaceIdx
  
  # TODO
  # - connectionPoints
  # - faceReference (was face) (done)
  return(face2D)
}

constructFlapje <- function(level, original)
{
  allConnectionEdges <- sapply(layout2D[1:level], function(l) {return(l$connectionEdge)})
  
  flapForOneFace <- function(face, original)
  {
    vex1 <- face$vexReferences
    vex2 <- shiftrotate(vex1)
    
    # Index of P1 and P2 for edges that are not used to connect 2D projections
    edge <- sapply(seq(length(vex1)), function(i){return(original$coordPairToEdge[vex1[i],vex2[i]])})
    idx1 <- which(!(edge %in% allConnectionEdges))
    idx2 <- (idx1%%length(vex1))+1
    
    wideFmt <- data.table(P1x = face$coords2D[idx1,1], 
                          P1y = face$coords2D[idx1,2], 
                          P2x = face$coords2D[idx2, 1],
                          P2y = face$coords2D[idx2, 2],
                          edge = edge[idx1])
    wideFmt[, d := distance(c(P1x, P1y), c(P2x, P2y))]
    wideFmt[, h := d/3] # height
    wideFmt[, alpha := h/tan(60*(2*pi/360))] # angle fixed to 60 degr
    wideFmt[, c("P3x", "P3y") := list(P1x + (alpha/d)*(P2x-P1x) - (h/d)*(P1y-P2y), P1y + (alpha/d)*(P2y-P1y) - (h/d)*(P2x-P1x))]
    wideFmt[, c("P4x", "P4y") := list(P2x - (alpha/d)*(P2x-P1x) - (h/d)*(P1y-P2y), P2y - (alpha/d)*(P2y-P1y) - (h/d)*(P2x-P1x))]
    
    return(wideFmt)
  }
  
  allWide <- rbindlist(lapply(layout2D[1:level],flapForOneFace,original))
  allWide[, n:=seq(.N), by=edge]
  allWide <- allWide[n==1]
  allWide[, eseq := seq(.N)]
  
  return(allWide)
}

drawLayout <- function(level, debug=T, original)
{
  colors <- assignColors(original)
  
  plotdata <- rbindlist(lapply(layout2D[1:level], function(l) {
    data.table(x=l$coords2D[,1], 
               y=l$coords2D[,2], 
               faceReference=l$faceReference, 
               color=colors[l$faceReference], 
               vex=original$faces[[l$faceReference]])
  }))
  colorMapIdentity <- unique(colors)
  names(colorMapIdentity) <- colorMapIdentity
  
  centers <- rbindlist(lapply(seq(level), function(l) {
    return(data.table(x=layout2D[[l]]$center2D[1],
                      y=layout2D[[l]]$center2D[2],
                      level=l,
                      faceReference=layout2D[[l]]$faceReference))
  }))
  
  fWide <- constructFlapje(level, original)
  fLong <- data.table( x = c(fWide$P1x, fWide$P3x, fWide$P4x, fWide$P2x),
                       y = c(fWide$P1y, fWide$P3y, fWide$P4y, fWide$P2y),
                       edge = rep(fWide$eseq, 4))
  fCenters <- fLong[, .(x=mean(x), y=mean(y)), by=edge]

  size <- 1.1*max(layout2D[[level]]$layoutMaxCoords[2]-layout2D[[level]]$layoutMinCoords[2],
                  layout2D[[level]]$layoutMaxCoords[1]-layout2D[[level]]$layoutMinCoords[1])
  xcenter <- mean(c(layout2D[[level]]$layoutMinCoords[1], layout2D[[level]]$layoutMaxCoords[1]))
  ycenter <- mean(c(layout2D[[level]]$layoutMinCoords[2], layout2D[[level]]$layoutMaxCoords[2]))
  
  p <-ggplot(plotdata, aes(x,y)) + 
    geom_polygon(data=fLong, mapping=aes(x,y,group=edge), fill="lightgrey", color="black", linetype="dotted", size=0.1, inherit.aes = F)+
    geom_text(data=fCenters, mapping=aes(x,y,label=edge), inherit.aes = F, size=2) +
    geom_polygon(aes(fill = color, group = faceReference), color="black", size=0.1, show.legend = FALSE) +
    scale_fill_manual(values = colorMapIdentity)+
    scale_x_continuous(limits = c(xcenter-size/2, xcenter+size/2))+
    scale_y_continuous(limits = c(ycenter-size/2, ycenter+size/2))+
    # scale_x_continuous(limits = c(min(layout2D[[level]]$layoutMinCoords), max(layout2D[[level]]$layoutMaxCoords)))+
    # scale_y_continuous(limits = c(min(layout2D[[level]]$layoutMinCoords), max(layout2D[[level]]$layoutMaxCoords)))+
    ggtitle(paste("Layout of", original$name), 
            subtitle = paste("round",ncalls,"eval=",round(bestEval,5),round(difftime(Sys.time(), startTime, units = "mins"),2),"mins"))
  
  if (debug) {
    p <- p + geom_label(mapping=aes(label=vex), size=3, color="black", alpha=0.6)+
      geom_text(data=centers, mapping=aes(x,y,label=paste0(level," (F",faceReference,")")), inherit.aes = F, size=3) +
      geom_hline(yintercept = layout2D[[level]]$layoutMinCoords[2], colour="yellow", linetype="dashed") + 
      geom_hline(yintercept = layout2D[[level]]$layoutMaxCoords[2], colour="yellow", linetype="dashed") + 
      geom_vline(xintercept = layout2D[[level]]$layoutMinCoords[1], colour="yellow", linetype="dashed") + 
      geom_vline(xintercept = layout2D[[level]]$layoutMaxCoords[1], colour="yellow", linetype="dashed") +
      theme_minimal()
  } else {
    p <- p + theme_void()
  }
  p <- p+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  print(p)
  
  if(stopAtFirstLayout) stop("Stopping after first layout, for debugging")
  
  ggsave(file=file.path("layouts",paste0("layout ", original$name, ".svg")), plot=p, width=10, height=10)
}

getLayoutDigest <- function(level, candidate=NULL)
{
  if (is.null(candidate)) {
    digest <- paste(sort(sapply(layout2D[1:level], function(f) {return(f$connectionEdge)})),collapse="-")
  } else {
    digest <- paste(sort(c(sapply(layout2D[1:level], function(f) {return(f$connectionEdge)}),candidate)),collapse="-")
  }
  return(digest)
}

edgeDescr <- function(polyhedron3D, e)
{
  return (paste0("(", paste(sort(polyhedron3D$edgeToVertices[e,]),collapse="-"), ")"))  
}

isBoundingBoxOverlap <- function(fig1, fig2)
{
  if (fig1$maxCoords[1] < fig2$minCoords[1] | deltaEquals(fig1$maxCoords[1], fig2$minCoords[1])) return(F)
  if (fig1$maxCoords[2] < fig2$minCoords[2] | deltaEquals(fig1$maxCoords[2], fig2$minCoords[2])) return(F)
  if (fig1$minCoords[1] > fig2$maxCoords[1] | deltaEquals(fig1$minCoords[1], fig2$maxCoords[1])) return(F)
  if (fig1$minCoords[2] > fig2$maxCoords[2] | deltaEquals(fig1$minCoords[2], fig2$maxCoords[2])) return(F)
  
  return(T)
}

isBoundingCircleOverlap <- function(fig1, fig2)
{
  maxRadius1 <- max(distance(fig1$coords2D, fig1$center2D))
  maxRadius2 <- max(distance(fig2$coords2D, fig2$center2D))
  d <- distance(fig1$center2D, fig2$center2D)
  if ((maxRadius1+maxRadius2 < d) | deltaEquals(maxRadius1+maxRadius2, 2)) return(F)
  
  return(T)
}

# Find the best 2D layout given global layout
best2DLayout <- function(polyhedron3D, level = 1, evaluator, trace=T, debug=T, face2D = NULL, maxrounds = 0)
{
  spacing <- paste0(rep(" ", level),collapse = "")
  if (level == 1) 
  {
    layout2D <<- list()
    allLayoutDigests <<- list()
    face3Dto2DprojectionCache <<- list()
    bestEval <<- Inf
    bestDigest <<- ""
    ncalls <<- 0
    nevalskips <<- 0
    ndigestskips <<- 0
    noverlapdetects <<- 0
    stopAtFirstLayout <<- F
  }
  
  if (maxrounds < 0 & bestDigest != "") return() # stop when there is 1 solution
  if (((ncalls+1) > maxrounds) & (maxrounds > 0)) {
    return()
  }
  ncalls <<- ncalls+1
  
  #
  # Update layout with new face2D
  #
  
  if (is.null(face2D)) {
    face2D <- positionNextFace(polyhedron3D, debug=debug)
  }
  layout2D[[level]] <<- face2D
  allLayoutDigests[[1+length(allLayoutDigests)]] <<- getLayoutDigest(level)
  
  # update the list of available edges, and store it with the layout itself
  shiftFacePoints <- shiftrotate(face2D$vexReferences)
  newEdges <- sapply(seq(length(face2D$vexReferences)), function(i) {return(polyhedron3D$coordPairToEdge[face2D$vexReferences[i],shiftFacePoints[i]])})
  if (level > 1) {
    layout2D[[level]]$layoutCandidateEdges <<-
      setdiff(c(newEdges, layout2D[[level-1]]$layoutCandidateEdges), intersect(layout2D[[level-1]]$layoutCandidateEdges, newEdges))
  } else {
    layout2D[[level]]$layoutCandidateEdges <<- newEdges
  }
  
  # keep track of the layout min/max
  if (level > 1) {
    layout2D[[level]]$layoutMinCoords <<- apply(matrix(c(face2D$minCoords, layout2D[[level-1]]$layoutMinCoords), ncol=2, byrow = T), 2, min)
    layout2D[[level]]$layoutMaxCoords <<- apply(matrix(c(face2D$maxCoords, layout2D[[level-1]]$layoutMaxCoords), ncol=2, byrow = T), 2, max)
  } else {
    layout2D[[level]]$layoutMinCoords <<- face2D$minCoords
    layout2D[[level]]$layoutMaxCoords <<- face2D$maxCoords
  }
  
  if (debug) {
    cat(spacing,
        "Entered level", level, "placed face:", face2D$faceReference, 
        "along edge", face2D$connectionEdge, edgeDescr(polyhedron3D, face2D$connectionEdge), 
        "eval =", evaluator(level), fill = T)
    if (trace) {
      drawLayout(level=level, debug = T, polyhedron3D)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
  
  #
  # Are we done?
  #
  if (evaluator(level) > bestEval) {
    if (debug) {cat(spacing,"*** early exit **** at level", level,"eval=", evaluator(level), "but best is", bestEval, fill = T)}
    nevalskips <<- nevalskips+1
    return()
  }
  
  #
  # Get candidates for next "move"
  #
  
  # get list of all faces placed
  placedFaces <- sapply(layout2D[1:level], function(f) {return(f$faceReference)})
  
  if (level == 1) {
    # for first level don't go through all the edges
    candidates <- layout2D[[level]]$layoutCandidateEdges[1]
  } else {
    candidates <- layout2D[[level]]$layoutCandidateEdges
  }
  if (length(candidates) == 0) {
    if (evaluator(level) < bestEval & !deltaEquals(evaluator(level), bestEval)) {
      bestEval <<- evaluator(level)
      bestDigest <<- getLayoutDigest(level)
      if (debug) {cat(spacing,"Done with better evaluation! At level", level, "eval=", evaluator(level), fill = T)  }
      if (!debug) {cat("Found layout, round",ncalls,"eval=",round(bestEval,5),round(difftime(Sys.time(), startTime, units = "mins"),5),"mins", fill = T)  }
      drawLayout(level, debug=debug, original = polyhedron3D)
    } else {
      if (debug) {cat(spacing,"Done but no improvement. At level", level, "eval=", evaluator(level), fill = T)  }
    }
  } else {
    candidates <- candidates[sample.int(length(candidates))] # shuffle to speed up search
    
    # first find the 2D faces
    candidateFaces <- list()
    for (edge in candidates) {
      layoutDigest <- getLayoutDigest(level, edge)
      if (layoutDigest %in% allLayoutDigests) {
        if (debug) {cat(spacing,"*** skip processing edge", edge, " **** at level", level,"layout done already:", layoutDigest, fill = T)}
        ndigestskips <<- ndigestskips+1
        next
      }
      
      if (debug) {cat(spacing,"select edge:", edge, edgeDescr(polyhedron3D, edge), "from", candidates, fill = T)}
      
      # place it into position
      newFace2D <- positionNextFace(polyhedron3D, edge, placedFaces, level, debug) 
      
      candidateFaces[[1+length(candidateFaces)]] <- newFace2D
    }
    
    # then order them by (early) evaluation
    if (length(candidateFaces) > 0) {
      
      candidateEvaluation <- sapply(candidateFaces, function(c) {
        currentLayoutWithCandidateMin <- apply(matrix(c(c$minCoords, layout2D[[level]]$layoutMinCoords), ncol=2, byrow = T), 2, min)
        currentLayoutWithCandidateMax <- apply(matrix(c(c$maxCoords, layout2D[[level]]$layoutMaxCoords), ncol=2, byrow = T), 2, max)
        
        return (evaluator(level, currentLayoutWithCandidateMin, currentLayoutWithCandidateMax))
      })
      
      # then recurse, starting with the one with smallest local eval result
      for (c in order(candidateEvaluation)) {
        nextFace <- candidateFaces[[c]]
        
        ## TODO call
        # checkOverlap(nextFace, level)
        
        # last check if all is well - check if there is no overlap with existing faces
        hasOverlap <- NULL
        for (l in safeseq(level)) {
          if (layout2D[[l]]$faceReference != nextFace$connectedToFaceReference) {
            if (isBoundingBoxOverlap(layout2D[[l]], nextFace)) {
              
              # optionally bypass the connection pts
              # possibleOverlappingPoints <-
              #   which(!(nextFace$vexReferences %in% polyhedron3D$edgeToVertices[nextFace$connectionEdge,]))
              
              ## Any of l's points in nextFace or vice versa...?
              if (nextFace$connectedToFaceReference == layout2D[[l]]$faceReference) {
                # exclude points along the edge that directly connects these two faces
                connectedIdx <- layout2D[[l]]$vexReferences %in% polyhedron3D$edgeToVertices[nextFace$connectionEdge,]
                pip1 <- apply(layout2D[[l]]$coords2D[!connectedIdx], 1, isPointInFace, nextFace$coords2D)
              } else {
                pip1 <- apply(layout2D[[l]]$coords2D, 1, isPointInFace, nextFace$coords2D)
              }
              if (any(pip1)) {
                hasOverlap <- list(f=layout2D[[l]]$faceReference, p=layout2D[[l]]$vexReferences[pip1])
                break
              }
              pip2 <- apply(nextFace$coords2D, 1, isPointInFace, layout2D[[l]]$coords2D)
              if (any(pip2)) {
                hasOverlap <- list(f=layout2D[[l]]$faceReference, p=nextFace$vexReferences[pip2])
                break
              }
            }
          }
        }
        if (!is.null(hasOverlap)) {
          noverlapdetects <<- noverlapdetects+1
          if (debug) {
            cat(spacing,"overlap", 
                paste0("F",nextFace$faceReference), 
                "with", paste0("F",hasOverlap$f), 
                "at", paste(paste0("P", hasOverlap$p),collapse=","), fill=T)
          }
          #stopAtFirstLayout<<-T
          next
        }
        
        best2DLayout(polyhedron3D, level+1, evaluator, trace=trace, debug=debug, candidateFaces[[c]], maxrounds = maxrounds)
      }
    }
  }
  
  if (debug) {cat(spacing,"Returning from level", level, fill = T)}
}

layoutAreaEvaluator <- function(level, min=layout2D[[level]]$layoutMinCoords, max=layout2D[[level]]$layoutMaxCoords)
{
  # area - the smaller the better
  rslt <- ((max[1]-min[1])*(max[2]-min[2]))
  # # prefer horizontal by penalizing vertical orientation 
  # if ((max[2]-min[2]) > (max[1]-min[1])) rslt <- rslt*(1+1/(4*level))
  # this really slows down convergence - better swap at the end
  return(rslt)
}

layoutSquarenessEvaluator <- function(level, min=layout2D[[level]]$layoutMinCoords, max=layout2D[[level]]$layoutMaxCoords)
{
  # abs diff x y - the closer the better
  return (abs((max[1]-min[1]) - (max[2]-min[2])))
}

########

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- icosahedron
xx <- tetrahedron
xx <- truncate(cube)
xx <- dodecahedron
xx <- quasi(cube)
xx <- rhombic(cube)
xx <- truncate(cube)
xx <- truncate(octahedron)
xx <- cube
xx <- rhombic(dodecahedron)
#clear3d()
#drawPoly(xx, debug=T)

# for layout use a new device
# rgl_init(new.device = T)
# clear3d()
# drawAxes()

# Start with a layout with just the first face
startTime <- Sys.time()

# best2DLayout(xx, evaluator = layoutAreaEvaluator, trace=F, debug=F, maxrounds=50000)
best2DLayout(xx, evaluator = layoutAreaEvaluator, trace=F, debug=F)

here <- function() {}

endTime <- Sys.time()
cat("#calls:", ncalls, fill=T)
cat("#early skips because of repeated layout", ndigestskips, fill=T)
cat("#early skips because of eval results", nevalskips, fill=T)
cat("#overlaps detected", noverlapdetects, fill=T)
print(bestEval)
print(bestDigest) # layout can be reconstructed from this
cat("Elapsed:", as.double(difftime(endTime, startTime, units = "mins")), "mins", fill=T)

