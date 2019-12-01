# automatic layout gen

# TODO
# build up layout in batches? eg optimize first 6 then next etc hopefully gives quicker convergence
# store 2d projection with the 3D solid (??)

library(ggplot2)
library(svglite)
library(data.table)

source("geometry.R")
source("polyhedra.R") # only for testing

set.seed(1234)

eToStr <- function(polyhedron3D, e)
{
  return (paste0(e, " (", paste(sort(polyhedron3D$edgeToVertices[e,]),collapse="-"), ")"))  
}

fToStr <- function(f)
{
  return(paste0("F",f))
}

elapsedTimeToStr <- function()
{
  secs <- as.double(difftime(Sys.time(), startTime, units = "secs"))
  hrs <- secs %/% 3600
  mins <- (secs - hrs*3600) %/% 60
  return(paste0(hrs,":",mins,":",round(secs-hrs*3600-mins*60,2)))
}

#' Put a new face into position, given a global 2D layout and and edge to connect to
#'
#' @param polyhedron3D The 3D solid with full topology
#' @param edge Index of the edge to connect to
#' @param placedFaces Array of all faces from the 3D solid that are currently in the layout already
#' @param level Current level in the layout (to index the global layout2D structure)
#' @param debug Debug flag
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
    #vertices <- ((which(polyhedron3D$coordPairToEdge==edge)-1) %/% nrow(polyhedron3D$coordPairToEdge))+1
    vertices <- polyhedron3D$edgeToVertices[edge,]
    if (debug) {cat(paste0(rep(" ", level),collapse=""), 
                    "positioning next face at edge", eToStr(polyhedron3D, edge), 
                    "existing", fToStr(placedFaceIdx), 
                    "new", fToStr(newFaceIdx), fill = T)}
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
    currentFace2D <- layout2D[[which(placedFaces==placedFaceIdx)]]
    
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
  
  return(face2D)
}

drawLayout <- function(original, level=length(layout2D), debug=T)
{
  constructAllFlapjes <- function(level, original)
  {
    allConnectionEdges <- sapply(layout2D[1:level], function(l) {return(l$connectionEdge)})
    
    flapForOneFace <- function(face, original)
    {
      relativeSize <- 1/4 # size relative to edge
      angleDegrees <- 50 # angle of flapje in degrees
      
      vex1 <- face$vexReferences
      vex2 <- shiftrotate(vex1)
      
      # Index of P1 and P2 for edges that are not used to connect 2D projections
      edge <- sapply(seq(length(vex1)), function(i){return(original$coordPairToEdge[vex1[i],vex2[i]])})
      idx1 <- which(!(edge %in% allConnectionEdges))
      idx2 <- (idx1%%length(vex1))+1
      
      flap <- data.table(P1x = face$coords2D[idx1,1], 
                         P1y = face$coords2D[idx1,2], 
                         P2x = face$coords2D[idx2, 1],
                         P2y = face$coords2D[idx2, 2],
                         edge = edge[idx1])
      flap[, d := distance(c(P1x, P1y), c(P2x, P2y))]
      flap[, h := d*relativeSize] # height
      flap[, alpha := h/tan(angleDegrees*(2*pi/360))] # angle fixed to 60 degr
      flap[, c("P3x", "P3y") := list(P1x + (alpha/d)*(P2x-P1x) - (h/d)*(P1y-P2y), P1y + (alpha/d)*(P2y-P1y) - (h/d)*(P2x-P1x))]
      flap[, c("P4x", "P4y") := list(P2x - (alpha/d)*(P2x-P1x) - (h/d)*(P1y-P2y), P2y - (alpha/d)*(P2y-P1y) - (h/d)*(P2x-P1x))]
      
      return(flap)
    }
    
    all <- rbindlist(lapply(layout2D[1:level],flapForOneFace,original))
    all[, n:=seq(.N), by=edge] # only keep the first for all pairs of flapjes for the same edge
    all <- all[n==1]
    all[, eseq := seq(.N)] # sequential numbering of the edges (not their actual indices)
    
    return(all)
  }
  
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
  
  size <- 1.1*max(layout2D[[level]]$layoutMaxCoords[2]-layout2D[[level]]$layoutMinCoords[2],
                  layout2D[[level]]$layoutMaxCoords[1]-layout2D[[level]]$layoutMinCoords[1])
  xcenter <- mean(c(layout2D[[level]]$layoutMinCoords[1], layout2D[[level]]$layoutMaxCoords[1]))
  ycenter <- mean(c(layout2D[[level]]$layoutMinCoords[2], layout2D[[level]]$layoutMaxCoords[2]))
  
  p <-ggplot(plotdata, aes(x,y)) 
  
  if (!debug) {
    flapjes <- constructAllFlapjes(level, original)
    flapjesLongFmt <- data.table( x = c(flapjes$P1x, flapjes$P3x, flapjes$P4x, flapjes$P2x),
                                  y = c(flapjes$P1y, flapjes$P3y, flapjes$P4y, flapjes$P2y),
                                  edge = rep(flapjes$eseq, 4))
    fCenters <- flapjesLongFmt[, .(x=mean(x), y=mean(y)), by=edge]
    
    p <- p +
      geom_polygon(data=flapjesLongFmt, mapping=aes(group=edge), fill="lightgrey", color="black", linetype="dotted", size=0.1)+
      geom_text(data=fCenters, mapping=aes(label=edge), size=2)
  }
  
  p <- p + 
    geom_polygon(aes(fill = color, group = faceReference), color="black", size=0.1, show.legend = FALSE) +
    scale_fill_manual(values = colorMapIdentity)+
    scale_x_continuous(limits = c(xcenter-size/2, xcenter+size/2))+
    scale_y_continuous(limits = c(ycenter-size/2, ycenter+size/2))
  
  if (debug) {
    p <- p + geom_label(mapping=aes(label=vex), size=3, color="black", alpha=0.6)+
      geom_text(data=centers, mapping=aes(x,y,label=paste0(level," (",fToStr(faceReference),")")), inherit.aes = F, size=3) +
      geom_hline(yintercept = layout2D[[level]]$layoutMinCoords[2], colour="yellow", linetype="dashed") + 
      geom_hline(yintercept = layout2D[[level]]$layoutMaxCoords[2], colour="yellow", linetype="dashed") + 
      geom_vline(xintercept = layout2D[[level]]$layoutMinCoords[1], colour="yellow", linetype="dashed") + 
      geom_vline(xintercept = layout2D[[level]]$layoutMaxCoords[1], colour="yellow", linetype="dashed") +
      theme_minimal() +
      ggtitle(paste("Layout of", original$name), 
              subtitle = paste("round",polyStatus[["ncalls"]],"eval=",round(bestEval,5),elapsedTimeToStr()))
  } else {
    p <- p + theme_void() +
      ggtitle(paste("Layout of", original$name), subtitle = description(original))
  }
  p <- p+theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
  #ggsave(file=file.path("layouts",paste0("layout ", original$name, ".svg")), plot=p, width=10, height=10)
}


addFaceToLayout <- function(face2D, polyhedron3D, level) 
{
  layout2D[[level]] <<- face2D
  
  # update the list of available edges, and store it with the layout itself
  shiftFacePoints <- shiftrotate(face2D$vexReferences)
  newEdges <- sapply(seq(length(face2D$vexReferences)), function(i) {
    return(polyhedron3D$coordPairToEdge[face2D$vexReferences[i],shiftFacePoints[i]])
  })
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
}

#' Turns a layout into a digest string so that geometrically identical but different
#' layouts have the same digest string. This is to avoid doing duplicate work when just
#' the order of processing certain edges is different.
#'
#' @param level Current depth of the layout
#' @param candidate Candidate edge to include, or not - used for early exit
layoutToDigest <- function(level, candidate=NULL)
{
  connections <- sapply(layout2D[1:level], function(f) {return(f$connectionEdge)})
  digest <- paste(sort(c(connections,candidate)),collapse="-")
  return(digest)
}

#' Reconstruct global layout2D object from just a digest string. 
#'
#' @param digest Digest string, defaults to the best digest
#' @param poly Original 3D solid
digestToLayout <- function(poly, digest = bestDigest)
{
  layout2D <<- list()
  faces <- list()
  
  # first one is implicit
  level <- 1
  face <- positionNextFace(poly, debug=F)
  addFaceToLayout(face, poly, level) 
  faces[[level]] <- face$faceReference
  
  edges <- as.numeric(unlist(strsplit(digest,"-",fixed = T)))
  repeat {
    #print(unlist(faces))
    e <- edges[which(apply(poly$edgeToFaces[edges,],1,function(x) {return(1==length(intersect(x,unlist(faces))))}))[1]]
    if (is.na(e)) break
    #print(e)
    level <- level + 1
    face <- positionNextFace(poly, edge = e, placedFaces = unlist(faces), level = level, debug=F)
    faces[[level]] <- face$faceReference
    addFaceToLayout(face, poly, level) 
  }
}

logPoly <- function(elem = NULL)
{
  if (!("polyStatus" %in% ls(.GlobalEnv))) {
    polyStatus <<- list()  
  }
  if (is.null(elem)) {
    print(data.table(counter=names(unlist(polyStatus)), value=unlist(polyStatus))[order(-value)])
  } else {
    if (!(elem %in% names(polyStatus))) {
      polyStatus[[elem]] <<- 1L  
    } else {
      polyStatus[[elem]] <<- polyStatus[[elem]]+1L
    }
  }
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


# TODO skip checking points that are in both if faces are connected
# TODO bounding box check for points before going down PIP routine
# TODO consider going reverse through levels and bail out if whole layout fails bb test
# TODO PIP check can be avoided for points on the connection edges
checkOverlap <- function(nextFace, level, debug=T)
{
  spacing <- paste0(rep(" ", level),collapse = "")
  # last check if all is well - check if there is no overlap with existing faces
  hasOverlap <- NULL
  for (l in safeseq(level)) {
    currentFace <- layout2D[[l]]
    if (currentFace$faceReference != nextFace$connectedToFaceReference) {
      if (isBoundingBoxOverlap(currentFace, nextFace)) {
        logPoly("face bounding box overlap")
        # for all P in F1:
        #   if P outside of bounding box F2 then skip rest of loop
        #   if P is equal to any point of F2 (distance = 0) then skip rest of loop
        #   if isPointInFace(P, F2) then return TRUE (is overlap) 
        
        # and vice versa for F1, F1
        for (i in seq(nrow(currentFace$coords2D))) {
          if (currentFace$coords2D[i,1] < nextFace$minCoords[1]) next
          if (currentFace$coords2D[i,1] > nextFace$maxCoords[1]) next
          if (currentFace$coords2D[i,2] < nextFace$minCoords[2]) next
          if (currentFace$coords2D[i,2] > nextFace$maxCoords[2]) next
          logPoly("pip check")
          if (isPointInFace(currentFace$coords2D[i,], nextFace$coords2D)) {
            hasOverlap <- list(f=currentFace$faceReference, p=currentFace$vexReferences[i])  
            break;
          }
        }
        if (is.null(hasOverlap)) {
          for (i in seq(nrow(nextFace$coords2D))) {
            if (nextFace$coords2D[i,1] < currentFace$minCoords[1]) next
            if (nextFace$coords2D[i,1] > currentFace$maxCoords[1]) next
            if (nextFace$coords2D[i,2] < currentFace$minCoords[2]) next
            if (nextFace$coords2D[i,2] > currentFace$maxCoords[2]) next
            logPoly("pip check")
            if (isPointInFace(nextFace$coords2D[i,], currentFace$coords2D)) {
              hasOverlap <- list(f=currentFace$faceReference, p=nextFace$vexReferences[i])  
              break;
            }
          }
        }
      }
    }
  }
  if (!is.null(hasOverlap)) {
    logPoly("overlap detected")
    if (debug) {
      cat(spacing, "overlap", fToStr(nextFace$faceReference), 
          "with", fToStr(hasOverlap$f), 
          "at", paste(paste0("P", hasOverlap$p),collapse=","), fill=T)
    }
    # stopAtFirstLayout<<-T
  }
  
  return(!is.null(hasOverlap))
}

# area - the smaller the better
layoutAreaEvaluator <- function(level, min=NULL, max=NULL)
{
  if (is.null(min)) min <- layout2D[[level]]$layoutMinCoords
  if (is.null(max)) max <- layout2D[[level]]$layoutMaxCoords
  width <- max[1]-min[1]
  height <- max[2]-min[2]
  return (width*height)
}

# abs diff x y - the closer the better
layoutSquarenessEvaluator <- function(level, min=NULL, max=NULL)
{
  if (is.null(min)) min <- layout2D[[level]]$layoutMinCoords
  if (is.null(max)) max <- layout2D[[level]]$layoutMaxCoords
  width <- max[1]-min[1]
  height <- max[2]-min[2]
  return (abs(width - height))
}

# least amount of remaining space on portrait A4
layoutA4Evaluator <- function(level, min=NULL, max=NULL)
{
  if (!is.null(min) | !is.null(max)) return(1) # does not support early exit because it is relative
  if (is.null(min)) min <- layout2D[[level]]$layoutMinCoords
  if (is.null(max)) max <- layout2D[[level]]$layoutMaxCoords
  width <- max[1]-min[1]
  height <- max[2]-min[2]
  ratio <- sqrt(2)
  if (height > ratio*width) {
    return (height*abs(width - height/ratio))  
  } else {
    return (width*abs(height - ratio*width))
  }
}

#' Resursive function to find the best 2D layout for a solid. Uses a global layout list
#' and other globals to reduce stack use.
#'
#' @param polyhedron3D Solid to map to 2D 
#' @param face2D Next 2D face to add. Default to NULL in first call.
#' @param level Resursive level, defaults to 1.
#' @param evaluator How to evaluate the layout. Maximizing to smallest value.
#' @param debugLevel 0=quiet, 1=print when new layout found, 2=print all steps, 3=stop after each iteration
#' @param maxrounds -1=stop at first layout, 0=keep searching, n = stop after n iterations
best2DLayout <- function(polyhedron3D, face2D = NULL, level = 1, evaluator = layoutAreaEvaluator, debugLevel=1, maxrounds = 0)
{
  spacing <- paste0(rep(" ", level),collapse = "")
  if (level == 1) 
  {
    layout2D <<- list()
    allLayoutDigests <<- list()
    face3Dto2DprojectionCache <<- list()
    bestEval <<- Inf
    bestDigest <<- ""
    stopAtFirstLayout <<- F
    polyStatus <<- list()
    polyStatus[["ncalls"]] <- 0L
    startTime <<- Sys.time()
  }
  
  if (maxrounds < 0 & bestDigest != "") return() # stop when there is 1 solution
  if (((polyStatus[["ncalls"]]+1) > maxrounds) & (maxrounds > 0)) {
    return()
  }
  logPoly("ncalls")
  
  #
  # Update layout with new face2D
  #
  
  if (is.null(face2D)) {
    face2D <- positionNextFace(polyhedron3D, debug=(debugLevel>=2))
  }
  
  addFaceToLayout(face2D, polyhedron3D, level) 
  allLayoutDigests[[1+length(allLayoutDigests)]] <<- layoutToDigest(level)
  
  if (debugLevel >= 2) {
    cat(spacing,
        "Entered level", level, "placed face:", face2D$faceReference, 
        "along edge", eToStr(polyhedron3D, face2D$connectionEdge), 
        "eval =", evaluator(level), fill = T)
    if (debugLevel == 3) {
      plot <- drawLayout(polyhedron3D, level=level, debug = T)
      print(plot)
      invisible(readline(prompt="Press [enter] to continue"))
    }
  }
  
  #
  # Are we done?
  #
  if (evaluator(level) > bestEval) {
    if (debugLevel >= 2) {
      cat(spacing,"*** early exit **** at level", level,"eval=", evaluator(level), "but best is", bestEval, fill = T)
    }
    logPoly("early skips because of eval results")
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
      bestDigest <<- layoutToDigest(level)
      if (debugLevel >= 2) {
        cat(spacing,"Done with better evaluation! At level", level, "eval=", evaluator(level), fill = T)  
      }
      if (debugLevel >= 1) {
        cat("Found layout, round",polyStatus[["ncalls"]],"eval=",round(bestEval,5),elapsedTimeToStr(), fill = T)  
        plot <- drawLayout(original = polyhedron3D, level, debug=T)
        print(plot)
      }
      if (debugLevel >= 3) {
        invisible(readline(prompt="Press [enter] to continue"))
      }
    } else {
      if (debugLevel >= 2) {
        cat(spacing,"Done but no improvement. At level", level, "eval=", evaluator(level), fill = T)  
      }
    }
  } else {
    candidates <- candidates[sample.int(length(candidates))] # shuffle to speed up search
    
    # first find the 2D faces
    candidateFaces <- list()
    for (edge in candidates) {
      layoutDigest <- layoutToDigest(level, edge)
      if (layoutDigest %in% allLayoutDigests) {
        if (debugLevel >= 2) {
          cat(spacing,"*** skip processing edge", eToStr(polyhedron3D, edge), "at level", level,"layout done already:", layoutDigest, fill = T)
        }
        logPoly("early skips because of repeated layout")
        next
      }
      
      if (debugLevel >= 2) {
        cat(spacing,"select edge:", eToStr(polyhedron3D, edge), "from", candidates, fill = T)
      }
      
      # place it into position
      newFace2D <- positionNextFace(polyhedron3D, edge, placedFaces, level, debug=(debugLevel >= 2)) 
      if (checkOverlap(newFace2D, level, debug=(debugLevel >= 2))) { 
        next
      }
      
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
        
        best2DLayout(polyhedron3D, candidateFaces[[c]], level+1, evaluator, debugLevel, maxrounds)
      }
    }
  }
  
  if (debugLevel >= 2) {
    cat(spacing,"Returning from level", level, fill = T)
  }
}

#' Gets 2D layout of a solid, returning as a ggplot object that can be printed or saved to e.g. SVG.
#'
#' @param poly The 3D solid
#'
#' @return A ggplot object
get2DLayout <- function(poly)
{
  best2DLayout(poly, evaluator = layoutA4Evaluator, debugLevel=0, maxrounds = 5000)
  digestToLayout(poly, digest = bestDigest)
  return(drawLayout(original=poly, debug=F))
}

testLayout <- function()
{
  cube <- dual(octahedron, name = "Cube")
  dodecahedron <- dual(icosahedron, name = "Dodecahedron")
  
  xx <- tetrahedron
  xx <- truncate(cube)
  xx <- quasi(cube)
  xx <- truncate(cube)
  xx <- truncate(octahedron)
  xx <- rhombic(dodecahedron)
  xx <- dodecahedron
  xx <- icosahedron
  xx <- rhombic(dodecahedron)
  #clear3d()
  #drawPoly(xx, debug=T)
  
  # for layout use a new device
  # rgl_init(new.device = T)
  # clear3d()
  # drawAxes()
  
  # Start with a layout with just the first face
  
  #best2DLayout(xx, evaluator = layoutAreaEvaluator, debugLevel=1, maxrounds=5000)
  best2DLayout(xx, evaluator = layoutAreaEvaluator, debugLevel=1, maxrounds=0)
  
  
  #best2DLayout(xx, evaluator = layoutA4Evaluator, debug=T, maxrounds = -1)
  
  logPoly()
  
  print(bestEval)
  print(bestDigest) # layout can be reconstructed from this
  cat("Elapsed:", elapsedTimeToStr(), fill=T)
  
  digestToLayout(xx)
  drawLayout(xx, debug=F)
}

#testLayout()

testComplexLayout <- function()
{
  # here segments intersection w faces would work
  xx <- compose(tetrahedron, dual(tetrahedron))
  F1 <- xx$faces[[1]]
  F2 <- xx$faces[[5]]
  F3 <- xx$faces[[6]]
  
  # for this one, segments interect faces not going to work,
  # should get lines from faces through faces then check if pts
  # if really are in the polygons and check whether they
  # split existing segments
  xx <- greatDodecahedron
  F1 <- xx$faces[[4]]
  F2 <- xx$faces[[1]]
  F3 <- xx$faces[[9]]
  
  clear3d()
  drawPoly(xx, debug = T)
  
  rgl_init(new.device = T)
  #drawAxes()
  drawPolygon(F1, xx$coords, alpha=0.8, label=paste0("F",which(sapply(xx$faces, function(x){return(identical(x,F1))}))), drawlines=T, drawvertices=T)
  drawPolygon(F2, xx$coords, alpha=0.8, label=paste0("F",which(sapply(xx$faces, function(x){return(identical(x,F2))}))), drawlines=T, drawvertices=T)
  drawPolygon(F3, xx$coords, alpha=0.8, label=paste0("F",which(sapply(xx$faces, function(x){return(identical(x,F3))}))), drawlines=T, drawvertices=T)
  
  # get intersection line
  # http://geomalgorithms.com/a05-_intersect-1.html#intersect3D_2Planes()
  #
  # then intersect this with all segments
  
  allIntersectionSegments <- list()
  for (F2 in xx$faces) {
    # find intersection lines of plane of F1 with all other planes
    i <- intersect3D_2Planes(planeToNormalForm(F1, xx$coords), planeToNormalForm(F2, xx$coords))
    if (i$status == "intersect") {
      spheres3d(i$P0[1], i$P0[2], i$P0[3], color="red", radius = 0.01)
      spheres3d(i$P1[1], i$P1[2], i$P1[3], color="red", radius = 0.01)
      lines3d(c(i$P0[1],i$P1[1]), c(i$P0[2],i$P1[2]), c(i$P0[3], i$P1[3]), color="red")
      
      # Find all intersections of this line with the primary face
      intersectionPoints <- list()
      for (n in seq(length(F1))) {
        S0 <- xx$coords[F1[n],]
        S1 <- xx$coords[shiftrotate(F1)[n],]
        
        i_seg <- intersect_2Segments(i$P0, i$P1, S0, S1, firstIsLine = T)
        if (i_seg$status == "intersect") {
          # not list coords twice
          isNewPoint <- T
          if (length(intersectionPoints) >= 1) {
            for (p in seq(length(intersectionPoints))) {
              if (deltaEquals(i_seg$I0[1], intersectionPoints[[p]]$x)) {
                if (deltaEquals(i_seg$I0[2], intersectionPoints[[p]]$y)) {
                  if (deltaEquals(i_seg$I0[3], intersectionPoints[[p]]$z)) {
                    isNewPoint <- F
                    break
                  }  
                }  
              }  
            }
          }
          if (isNewPoint) {
            intersectionPoints[[1+length(intersectionPoints)]] <-
              data.table(x = i_seg$I0[1],
                         y = i_seg$I0[2],
                         z = i_seg$I0[3])
          }
          spheres3d(i_seg$I0[1], i_seg$I0[2], i_seg$I0[3], color="blue", radius = 0.07)
        }
      }
      # For now assume there are just 2 intersection points of the line with the face
      # TODO this is not true in general. We should check the other face(s) as well and 
      # define the intersection segments as those that are inside both polygons. There may
      # even be multiple segments from one line (at least in theory).
      
      # do the same for the other face
      # if there are only 2 intersections assume a line through them
      # otherwise counting magic
      
      if (length(intersectionPoints) != 2) stop("Assuming exactly 2 intersections with primary face")
      intersectionPointsTable <- rbindlist(intersectionPoints)
      # TODO avoid duplicate segments
      allIntersectionSegments[[1+length(allIntersectionSegments)]] <- 
        data.table( primaryFace = which(sapply(xx$faces, identical, F1)),
                    intersectFace = which(sapply(xx$faces, identical, F2)),
                    P0_x = intersectionPointsTable$x[1],
                    P0_y = intersectionPointsTable$y[1],
                    P0_z = intersectionPointsTable$z[1],
                    P1_x = intersectionPointsTable$x[2],
                    P1_y = intersectionPointsTable$y[2],
                    P1_z = intersectionPointsTable$z[2])
    }
  }
  
  apply(rbindlist(allIntersectionSegments), 1, function(seg) {
    lines3d(c(seg[3],seg[6]), c(seg[4],seg[7]), c(seg[5], seg[8]), color="blue")  
  })
  
  # then, intersect all those new face segments with eachother
  # and test if these intersections are inside the polygons
  
}
