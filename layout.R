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
  # drawInit(new.device = T)
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

# TODO change to segmentation returning list of coords and segments
# global hullCoords/Edges are returns of this segmentation function

#' Get intersection lines of all faces of a polygon.
#'
#' @param poly 
segmentation <- function(poly)
{
  # Start with current polygon 
  # hullCoords is a list of unique coordinates of the hull
  # TODO not global but return from this func
  
  poly <- compose(tetrahedron, dual(tetrahedron))
  poly <- greatDodecahedron
  
  spacing <- 1.1
  primaryFace<-1
  
  hullCoords <<- as.list(as.data.table(t(poly$coords))) 
  vexCells <- which(upper.tri(poly$coordPairToEdge) & (poly$coordPairToEdge != 0))
  hullEdges <<- lapply(vexCells, function(edge) {
    dim <- nrow(poly$coordPairToEdge)
    col <- (edge-1)%%dim+1
    row <- (edge-1)%/%dim+1
    #
    # for debugging only face 1
    #
    if (poly$coordPairToFaces[row,col] != primaryFace & poly$coordPairToFaces[col,row] != primaryFace) return()
    return (data.table(vex1=row, vex2=col, srcVex1=row, srcVex2=col, 
                       F1=poly$coordPairToFaces[row,col], F2=poly$coordPairToFaces[col,row]))
  })
  
  clear3d()
  drawPolygon(poly$faces[[primaryFace]], poly$coords, alpha=0.5, col = "yellow", label=fToStr(primaryFace), drawlines=T, drawvertices=T)

  getHullCoordIdx <- function(p)
  {
    for (i in safeseq(length(hullCoords))) {
      if (deltaEquals(0, distance(p, hullCoords[[i]]))) {
        return(i)
      }
    } 
    hullCoords[[1+length(hullCoords)]] <<- p
    return(length(hullCoords))
  }

  findHullEdge <- function(a, b)
  {
    idx <- which(sapply(hullEdges, function(edge) {
      return( (edge$vex1==a && edge$vex2==b) || (edge$vex1==b && edge$vex2==a))
    }))   
    if (length(idx) == 0) return(NULL)
    return(idx)
  }
  
  addHullEdge <- function(a, b, F1, F2=F1, srca=NA, srcb=NA)
  {
    # should we check first if it already exists?
    hullEdges[[1+length(hullEdges)]] <<- 
      data.table(vex1=a, vex2=b, srcVex1=srca, srcVex2=srcb, F1=F1, F2=F2)
  }
  
  splitHullEdge <- function(a, b, mid)
  {
    if (b != mid & a != mid) {
      old <- findHullEdge(a,b)
      if (is.null(findHullEdge(a,mid))) {
        if (is.null(old)) stop("Splitting old edge that cannot be found")
        addHullEdge(a,mid,hullEdges[[old]]$F1,hullEdges[[old]]$F2,hullEdges[[old]]$srcVex1,hullEdges[[old]]$srcVex2)
      }
      if (is.null(findHullEdge(mid,b))) {
        if (is.null(old)) stop("Splitting old edge that cannot be found")
        addHullEdge(mid,b,hullEdges[[old]]$F1,hullEdges[[old]]$F2,hullEdges[[old]]$srcVex1,hullEdges[[old]]$srcVex2)
      }
      if (!is.null(old)) hullEdges[[old]] <<- NULL
    }
  }
  
  
  # All face intersections
  facePairs <- combn(length(poly$faces),2)
  
  # Exclude face pairs that are connected anyway (poly$edgeToFaces)
  isEdgeConnected <- apply(facePairs,2,function(fp) { 
    return( any(apply(poly$edgeToFaces, 1, function(edge) { 
      return (2 == length(intersect(fp,edge)))
    })))
  })
  
  if (all(isEdgeConnected)) return(data.table())
  
  # TODO we could (perhaps) further optimize to exclude face pairs with a too large distance to eachother
  
  # Find intersection lines of the face pairs - these are not segments yet but unbounded lines
  intersectionLines <- rbindlist(apply(facePairs[,!isEdgeConnected], 2, function(facepair) {
    i <- intersect3D_2Planes(planeToNormalForm(poly$faces[[facepair[1]]], poly$coords), 
                             planeToNormalForm(poly$faces[[facepair[2]]], poly$coords))
    if (i$status == "intersect") {
        return(data.table(F1=facepair[1], F2=facepair[2], 
                        P0_x=i$P0[1], P0_y=i$P0[2], P0_z=i$P0[3], 
                        P1_x=i$P1[1], P1_y=i$P1[2], P1_z=i$P1[3]))
    }
  }))

  # for debugging restrict to one face
  intersectionLines <- intersectionLines[F1==primaryFace | F2==primaryFace]
    
  for (i in unique(unlist(intersectionLines[,c("F1","F2")]))) {
    if (i != primaryFace) {
      drawPolygon(poly$faces[[i]], poly$coords, alpha=0.3, col = "grey", label=fToStr(i), drawlines=T, drawvertices=T)
    }
  }
  
  # Turn lines into segments by truncating to the edges of the face
  # TODO: for now assuming either face returns the same result. This is not true for complex cases
  # so we should do both then intersect (overlap) the resulting segments.
  dummy <- lapply(safeseq(nrow(intersectionLines)), function(i) {
    faceVertices <- poly$faces[[intersectionLines[i]$F1]]
    faceCoords <- poly$coords[faceVertices, ]
    #next_idx <- shiftrotate(seq(nrow(faceCoords))) # maybe just replace by doing j%% + 1
    #prev_idx <- shiftrotate(seq(nrow(faceCoords)), -1)
    skipNext <- F
    intersectionSegmentP0 <- NULL # really .. intersection1/2
    intersectionSegmentP1 <- NULL
    # would probably become a 4x3 matrix in the general case
    for (j in seq(nrow(faceCoords))) {
      if (skipNext) { 
        skipNext <- F
        next
      }
      i_seg <- intersect_2Segments(c(intersectionLines[i]$P0_x, intersectionLines[i]$P0_y, intersectionLines[i]$P0_z), 
                                   c(intersectionLines[i]$P1_x, intersectionLines[i]$P1_y, intersectionLines[i]$P1_z), 
                                   faceCoords[j,], 
                                   faceCoords[(j%%nrow(faceCoords))+1,], 
                                   firstIsLine = T)
      # if (intersectionLines[i]$F2==4) {
      #   print(i_seg)
      # }
      if (i_seg$status == "intersect") {
        intersectionVertexIndex <- getHullCoordIdx(i_seg$I0)
        if (is.null(intersectionSegmentP0)) {
          intersectionSegmentP0 <- intersectionVertexIndex
        } else if (is.null(intersectionSegmentP1)) {
          intersectionSegmentP1 <- intersectionVertexIndex
          
          cat("Add edge:", intersectionSegmentP0, "-", intersectionSegmentP1, fill=T)
          addHullEdge(intersectionSegmentP0, intersectionSegmentP1,
                      intersectionLines[i]$F1, intersectionLines[i]$F2)
        } else {
          stop("Expected no more than 2 intersections of two faces")
        }
        
        # if the intersection is the segment end, skip next round as this would start with the same point
        if (deltaEquals(0, distance(i_seg$I0, faceCoords[(j%%nrow(faceCoords))+1,]))) {
          skipNext <- T  
        } else {
          if (!deltaEquals(0, distance(i_seg$I0, faceCoords[j,]))) {
            # somewhere in the segment
            cat("Split edge:", 
                faceVertices[j], 
                intersectionVertexIndex, 
                faceVertices[(j%%nrow(faceCoords))+1], fill=T)
            splitHullEdge(faceVertices[j], faceVertices[(j%%nrow(faceCoords))+1], 
                          intersectionVertexIndex)
          }
        }
        
        if (!is.null(intersectionSegmentP0) & !is.null(intersectionSegmentP1)) {
          break # both intersections found - assuming there are only 2! maybe we should not and assert?!
        }
      }
    }
    if (is.null(intersectionSegmentP0) | is.null(intersectionSegmentP1)) stop("Expected two points")
  })
  
  # show new points
  if (length(hullCoords) > nrow(poly$coords)) {
    for (i in (nrow(poly$coords)+1):length(hullCoords)) {
      drawDots(hullCoords[[i]], color="red", radius = 0.03)
      drawTexts(spacing*hullCoords[[i]], text=i, color="red")
    }
  }
  # show face-face intersection segments
  print(rbindlist(hullEdges))
  faceSegments <- rbindlist(hullEdges)[is.na(srcVex1) & is.na(srcVex2)]
  for (i in safeseq(nrow(faceSegments)))  {
    #drawDots(c(hullCoords[[segments[i]$S0]], hullCoords[[segments[i]$S1]]), color="red", radius = 0.03)
    #drawTexts(spacing*hullCoords[[h[i]$S0]], text=h[i]$S0, color="red")
    #drawTexts(spacing*hullCoords[[h[i]$S1]], text=h[i]$S1, color="red")
    drawSegment(spacing*hullCoords[[faceSegments[i]$vex1]], hullCoords[[faceSegments[i]$vex2]], color="red")
  }
  
  # Intersect all those new segments amongst eachother
  segmentPairs <- combn(nrow(faceSegments),2)
  for (i in safeseq(ncol(segmentPairs))) {
    seg1_start <- hullCoords[[faceSegments[segmentPairs[1,i]]$vex1]]
    seg1_end <- hullCoords[[faceSegments[segmentPairs[1,i]]$vex2]]
    seg2_start <- hullCoords[[faceSegments[segmentPairs[2,i]]$vex1]]
    seg2_end <- hullCoords[[faceSegments[segmentPairs[2,i]]$vex2]]
    
    i_seg <- intersect_2Segments(seg1_start, seg1_end, seg2_start, seg2_end)
    
    #print(segmentPairs[,i])
    # intersect and one of the two flags t then skip? no new pt
    if (i_seg$status == "intersect") {
      #print(i_seg)
      # only if genuinly new point
      if (!deltaEquals(0, distance(i_seg$I0, seg1_start)) & !deltaEquals(0, distance(i_seg$I0, seg1_end)) &
          !deltaEquals(0, distance(i_seg$I0, seg2_start)) & !deltaEquals(0, distance(i_seg$I0, seg2_end))) 
      {
        print("New pt") # !! TODO: register but also SPLIT the segments and repeat??
        idx <- getHullCoordIdx(i_seg$I0)
        drawDots(i_seg$I0, color="purple", radius = 0.02)
        drawTexts(spacing*i_seg$I0, text=idx, color="purple")
      }
    }
  }
  
  #return(merge(intersectionSegments, intersectionLines, by=c("F1","F2"))) # merge only for debug
  #return(intersectionSegments)
}

# Find segment end points for the intersection of a line with one polygon
hull <- function(poly)
{
  spacing <- 1.1
  
  poly <- greatDodecahedron
  poly <- compose(tetrahedron, dual(tetrahedron))
  #poly <- greatIcosahedron
  
  segments <- segmentation(poly)
  
  clear3d()
  primaryFace <- 1
  for (i in unique(unlist(segments[F1 == primaryFace | F2 == primaryFace,1:2]))) {
    if (i == primaryFace) {
      drawPolygon(poly$faces[[i]], poly$coords, alpha=0.5, col = "yellow", label=fToStr(i), drawlines=T, drawvertices=T)
    } else {
      drawPolygon(poly$faces[[i]], poly$coords, alpha=0.3, col = "grey", label=fToStr(i), drawlines=T, drawvertices=T)
    }
  }
  for (i in which(segments$F1 == primaryFace | segments$F2 == primaryFace))  {
    # for debugging only:  
    drawDots(c(hullCoords[[segments[i]$S0]], hullCoords[[segments[i]$S1]]), color="red", radius = 0.03)
    drawTexts(spacing*hullCoords[[segments[i]$S0]], text=segments[i]$S0, color="red")
    drawTexts(spacing*hullCoords[[segments[i]$S1]], text=segments[i]$S1, color="red")
  }
  
  # below belongs in segmentation routine too
  

  
  # Points:
  # We have old face vertices, new segment end points (possibly 0) and segment intersections (possibly 0)
  # Label them and give them a unique index (eventually: globally)
  
  # poly$faces[[primaryFace]]
  # start at 1
  # traverse first edge (1-2) until you encounter an intersect
  
  # Split segments:
  # We have old face edges and new segments (possibly 0)
  # Now see if any of the points are in (not end points) and if so split into multiple segments (could be > 2)
  
  # Traverse:
  # Start at old face vertex V1 and old face edge E
  #   go from V1 to E until you encounter a (new) vertex V2
  #   at that vertex take the sharpest (left) possible 
}

