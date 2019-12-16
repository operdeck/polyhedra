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
segmentation <- function(poly, debug=F)
{
  # Start with current polygon 
  # hullCoords is a list of unique coordinates of the hull
  # TODO not global but return from this func
  
  #poly <- compose(tetrahedron, dual(tetrahedron))
  #poly <- greatDodecahedron
  #poly <- greatIcosahedron
  
  spacing <- 1.1
  #debug<-T
  debugPrimaryFace<-NA

  hullCoords <<- as.list(as.data.table(t(poly$coords))) 
  vexCells <- which(upper.tri(poly$coordPairToEdge) & (poly$coordPairToEdge != 0))
  hullEdges <<- lapply(vexCells, function(edge) {
    dim <- nrow(poly$coordPairToEdge)
    col <- (edge-1)%%dim+1
    row <- (edge-1)%/%dim+1
    return (data.table(vex1=row, vex2=col, srcVex1=row, srcVex2=col, 
                       F1=poly$coordPairToFaces[row,col], F2=poly$coordPairToFaces[col,row],
                       onOriginal=T))
  })
  if (debug) {
    clear3d()
    if (!is.na(debugPrimaryFace)) {
      hullEdges <- hullEdges[sapply(hullEdges, function(e) { return(e$F1 == debugPrimaryFace || e$F2 == debugPrimaryFace)})]
      drawPolygon(poly$faces[[debugPrimaryFace]], poly$coords, alpha=0.5, col = "yellow", label=fToStr(debugPrimaryFace), drawlines=T, drawvertices=T)
    }
  }

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

  findHullSourceEdge <- function(a, b)
  {
    idx <- which(sapply(hullEdges, function(edge) {
      return( (edge$srcVex1==a && edge$srcVex2==b) || (edge$srcVex1==b && edge$srcVex2==a))
    }))   
    if (length(idx) == 0) return(NULL)
    return(idx)
  }
  
  addHullEdge <- function(a, b, F1, F2=F1, srca=a, srcb=b, onOriginalPoly=FALSE, debug=F)
  {
    if (is.null(findHullEdge(a,b))) {
      hullEdges[[1+length(hullEdges)]] <<- 
        data.table(vex1=a, vex2=b, srcVex1=srca, srcVex2=srcb, F1=F1, F2=F2, onOriginal=onOriginalPoly)
    } else {
      if (debug) cat("Segment",a,"-",b,"already exists so skipping adding",fill=T)
    }
  }
  
  splitHullEdge <- function(a, b, mid, debug=F)
  {
    if (debug) cat("Split edge:", a, "- (", mid, ") -", b, fill=T)
    
    if (b != mid & a != mid) {
      old <- findHullEdge(a,b)
      # if old no longer exists find (possibly multiple) edges that have a/b as the source edge
      # then determine in which one "mid" is and swap that for a/b
      if (is.null(old)) {
        if (!is.null(findHullEdge(a,mid))) {
          if (!is.null(findHullEdge(mid,b))) {
            if (debug) cat("Segments", a, "-", mid, "and", mid, "-", b, "exist already", fill=T)
            return()
          }  
        }
        segmentSplits <- findHullSourceEdge(a,b)
        print(segmentSplits)
        segmentSplit <- which(sapply(hullEdges[segmentSplits], function(e) {
          pointSegmentCheck <- 
            intersect_2Segments(hullCoords[[mid]], hullCoords[[mid]], 
                                hullCoords[[e$vex1]], hullCoords[[e$vex2]])
          if (pointSegmentCheck$status != "intersect") return(FALSE)
          if (deltaEquals(distance(pointSegmentCheck$I0, hullCoords[[e$vex1]]), 0)) return(FALSE)
          if (deltaEquals(distance(pointSegmentCheck$I0, hullCoords[[e$vex2]]), 0)) return(FALSE)
          return(TRUE)
        }))
        if (length(segmentSplit) > 1) {
          print(hullEdges[segmentSplits])
          print(segmentSplit)
          stop(paste("There should be max one segment split, got", length(segmentSplit)))
        }
        if (length(segmentSplit) == 0) {
          if (debug) cat("Nothing to split - splits probably exists already", fill=T)
          return() # nothing to split
        }
        old <- segmentSplits[segmentSplit]
        if (debug) cat("Old segment",a,"-",b,"was split already so splitting",
            hullEdges[[old]]$vex1,"-",hullEdges[[old]]$vex2, "at", mid, fill=T)
        a <- hullEdges[[old]]$vex1
        b <- hullEdges[[old]]$vex2
      }
      if (!is.null(old)) {
        addHullEdge(a,mid,hullEdges[[old]]$F1,hullEdges[[old]]$F2,hullEdges[[old]]$srcVex1,hullEdges[[old]]$srcVex2,hullEdges[[old]]$onOriginal, debug=debug)
        addHullEdge(mid,b,hullEdges[[old]]$F1,hullEdges[[old]]$F2,hullEdges[[old]]$srcVex1,hullEdges[[old]]$srcVex2,hullEdges[[old]]$onOriginal, debug=debug)
        hullEdges[[old]] <<- NULL
      }
    }
  }
  
  # All face intersections, exclude face pairs that are connected anyway (poly$edgeToFaces)
  facePairs <- Filter(Negate(is.null), 
                      apply(combn(length(poly$faces),2), 2, function(fp) { 
                        if (!any(apply(poly$edgeToFaces, 1, function(edge) { 
                          return (2 == length(intersect(fp,edge))) }))) {
                        return(fp) }}))

  # TODO we could (perhaps) further optimize to exclude face pairs with a too large distance to eachother
  # TODO perhaps even if they have 1 point in common??
  
  # Find intersection lines of the face pairs - these are not segments yet but unbounded lines
  intersectionLines <- rbindlist(lapply(facePairs, function(facepair) {
    i <- intersect3D_2Planes(planeToNormalForm(poly$faces[[facepair[1]]], poly$coords), 
                             planeToNormalForm(poly$faces[[facepair[2]]], poly$coords))
    if (i$status == "intersect") {
        return(data.table(F1=facepair[1], F2=facepair[2], 
                        P0_x=i$P0[1], P0_y=i$P0[2], P0_z=i$P0[3], 
                        P1_x=i$P1[1], P1_y=i$P1[2], P1_z=i$P1[3]))
    }
  }))

  # for debugging restrict to one face
  if (debug) {
    if (nrow(intersectionLines) > 0) {
      if (!is.na(debugPrimaryFace)) {
        intersectionLines <- intersectionLines[F1==debugPrimaryFace | F2==debugPrimaryFace]
      }
      for (i in unique(unlist(intersectionLines[,c("F1","F2")]))) {
        if (is.na(debugPrimaryFace) | i != debugPrimaryFace) {
          drawPolygon(poly$faces[[i]], poly$coords, alpha=0.3, col = "grey", label=fToStr(i), drawlines=T, drawvertices=T)
        }
      }
    }
  }
  
  # Turn lines into segments by truncating to the edges of the face
  # TODO: for now assuming either face returns the same result. This is not true for complex cases
  # so we should do both then intersect (overlap) the resulting segments.
  dummy <- lapply(safeseq(nrow(intersectionLines)), function(i) {
    # intersection i is between F1 and F2, now truncate the lines to the boundaries of F1 and F2
    # TODO should probably take intersect of the two intersections... not assuming they stretch to the edges...
    for (f in c(intersectionLines[i]$F1, intersectionLines[i]$F2)) {
      faceVertices <- poly$faces[[f]]
      faceCoords <- poly$coords[faceVertices, ]
      skipNext <- F
      skipLast <- F
      # NB we're assuming the intersection yields exactly 0 or 2 points. In general this may not be true.
      intersectionSegmentP0 <- NULL
      intersectionSegmentP1 <- NULL
      for (j in seq(nrow(faceCoords))) {
        if (skipLast && j==nrow(faceCoords)) break;
        if (skipNext) { 
          skipNext <- F
          next
        }
        i_seg <- intersect_2Segments(c(intersectionLines[i]$P0_x, intersectionLines[i]$P0_y, intersectionLines[i]$P0_z), 
                                     c(intersectionLines[i]$P1_x, intersectionLines[i]$P1_y, intersectionLines[i]$P1_z), 
                                     faceCoords[j,], 
                                     faceCoords[(j%%nrow(faceCoords))+1,], 
                                     firstIsLine = T)
        if (i_seg$status == "intersect") {
          intersectionVertexIndex <- getHullCoordIdx(i_seg$I0)
          if (debug) drawDots(i_seg$I0, label=intersectionVertexIndex, color="red", radius=0.05)
          if (is.null(intersectionSegmentP0)) {
            intersectionSegmentP0 <- intersectionVertexIndex
          } else if (is.null(intersectionSegmentP1)) {
            intersectionSegmentP1 <- intersectionVertexIndex
            if (debug) cat("Add edge:", intersectionSegmentP0, "-", intersectionSegmentP1, 
                           "(", fToStr(intersectionLines[i]$F1), "|", fToStr(intersectionLines[i]$F2), ")", fill=T)
            addHullEdge(intersectionSegmentP0, intersectionSegmentP1,
                        intersectionLines[i]$F1, intersectionLines[i]$F2,
                        onOriginalPoly=FALSE, debug=debug)
          } else {
            stop(paste("Expected no more than 2 intersections of intersection line between",
                       fToStr(intersectionLines[i]$F1), "and", fToStr(intersectionLines[i]$F2), 
                       "with edges of", fToStr(intersectionLines[i]$F1), 
                       "but got additional intersection at", intersectionVertexIndex))
          }
          
          # intersection = segment end, then skip next round as this would start with the same point
          if (deltaEquals(0, distance(i_seg$I0, faceCoords[(j%%nrow(faceCoords))+1,]))) {
            skipNext <- T  
          } else {
            if (deltaEquals(0, distance(i_seg$I0, faceCoords[j,]))) {
              # intersection = segment start, then skip last round as this would start with the same point
              if (j == 1) skipLast <- T
            } else {
              # somewhere in the segment
              splitHullEdge(faceVertices[j], faceVertices[(j%%nrow(faceCoords))+1], intersectionVertexIndex, debug=debug)
            }
          }
          # TODO: for efficiency we could do below check however not having this check 
          # helps validate the assumption that there are no more than 2 intersection points
          # if (!is.null(intersectionSegmentP0) & !is.null(intersectionSegmentP1)) {
          #   break # both intersections found - assuming there are only 2! maybe we should not and assert?!
          # }
        }
      }
      # they can border on just a single point so only > 2 should be checked
    }
  })

  # Get intersections of the new segments, face by face
  allFaceIntersectionsSegments <- rbindlist(hullEdges)[(!onOriginal)]
  for (currentFace in unique(allFaceIntersectionsSegments$F1)) {
    if (debug) cat("Intersections of face segments for face", currentFace, fill=T)
    faceIntersectionsSegments <- allFaceIntersectionsSegments[F1 == currentFace | F2 == currentFace]
  
    # show new points
    if (debug) {
      if (length(hullCoords) > nrow(poly$coords)) {
        for (i in (nrow(poly$coords)+1):length(hullCoords)) {
          drawDots(hullCoords[[i]], color="red", radius = 0.03)
          drawTexts(spacing*hullCoords[[i]], text=i, color="red")
        }
      }
      # show face-face intersection segments
      # print(rbindlist(hullEdges))
      if (nrow(faceIntersectionsSegments) > 0) {
        drawSegments(spacing*t(as.data.table(hullCoords[faceIntersectionsSegments$vex1])), 
                     spacing*t(as.data.table(hullCoords[faceIntersectionsSegments$vex2])), color="red")
      }
    }
    
    # Intersect all those new segments amongst eachother
    if (nrow(faceIntersectionsSegments) > 0) {
      segmentPairs <- combn(nrow(faceIntersectionsSegments),2)
      for (i in safeseq(ncol(segmentPairs))) {
        seg1_start <- hullCoords[[faceIntersectionsSegments[segmentPairs[1,i]]$vex1]]
        seg1_end <- hullCoords[[faceIntersectionsSegments[segmentPairs[1,i]]$vex2]]
        seg2_start <- hullCoords[[faceIntersectionsSegments[segmentPairs[2,i]]$vex1]]
        seg2_end <- hullCoords[[faceIntersectionsSegments[segmentPairs[2,i]]$vex2]]
        
        intersectionNewSegs <- intersect_2Segments(seg1_start, seg1_end, seg2_start, seg2_end)
        
        if (intersectionNewSegs$status == "intersect") {
          idx <- getHullCoordIdx(intersectionNewSegs$I0)
          splitsSeg1 <- !deltaEquals(0, distance(intersectionNewSegs$I0, seg1_start)) & !deltaEquals(0, distance(intersectionNewSegs$I0, seg1_end))
          splitsSeg2 <- !deltaEquals(0, distance(intersectionNewSegs$I0, seg2_start)) & !deltaEquals(0, distance(intersectionNewSegs$I0, seg2_end))
          if (splitsSeg1 & splitsSeg2) {
            if (debug) cat("New inter-segment point:",idx,fill=T)
            drawDots(intersectionNewSegs$I0, label = idx, color="purple", radius = 0.02)
          }
          if (splitsSeg1) {
            splitHullEdge(faceIntersectionsSegments[segmentPairs[1,i]]$vex1, 
                          faceIntersectionsSegments[segmentPairs[1,i]]$vex2, idx, debug=debug)
          }
          if (splitsSeg2) {
            splitHullEdge(faceIntersectionsSegments[segmentPairs[2,i]]$vex1, 
                          faceIntersectionsSegments[segmentPairs[2,i]]$vex2, idx, debug=debug)
          }
        }
      }
    }
  }
  
  return (list(coords = as.matrix(t(as.data.table(hullCoords))), edges = rbindlist(hullEdges)))
}

# Find segment end points for the intersection of a line with one polygon
hull <- function(poly)
{
  poly <- greatDodecahedron
  poly <- compose(tetrahedron, dual(tetrahedron))
  #poly <- greatIcosahedron
  
  segments <- segmentation(tetrahedron) # no intersections at all - returns the same poly
  segments <- segmentation(cube) # some parallel faces
  segments <- segmentation(icosahedron) # faces connected at single point
  segments <- segmentation(greatDodecahedron, debug=T)
  segments <- segmentation(greatIcosahedron) # multiple intersections of face segments
  segments <- segmentation(compose(tetrahedron, dual(tetrahedron)))
  segments <- segmentation(greatStellatedDodecahedron) # does not work yet - gets multiple intersections 
  
  # segments <- list(coords = as.matrix(t(as.data.table(hullCoords))), edges = rbindlist(hullEdges))
  clear3d()
  edgesOriginal <- segments$edges[onOriginal & (srcVex1==vex1 & srcVex2==vex2)]
  edgesOnOriginalEdges <- segments$edges[onOriginal & !(srcVex1==vex1 & srcVex2==vex2)]
  edgesFromFaceIntersections <- segments$edges[(!onOriginal)]
  
  verticesOriginal <- which(rownames(segments$coords) != "")
  verticesOnOriginalEdges <- setdiff(unique(c(edgesOnOriginalEdges$vex1, edgesOnOriginalEdges$vex2)), verticesOriginal)
  verticesFromFaceIntersections <- setdiff(seq(length(hullCoords)), c(verticesOriginal, verticesOnOriginalEdges))
  
  drawDots(segments$coords[verticesOriginal,], label=verticesOriginal, color="green", radius=0.02)
  drawDots(segments$coords[verticesOnOriginalEdges,], label=verticesOnOriginalEdges, color="blue", radius=0.02)
  drawDots(segments$coords[verticesFromFaceIntersections,], label=verticesFromFaceIntersections, color="purple", radius=0.02)
  drawSegments(segments$coords[edgesOriginal$vex1,], segments$coords[edgesOriginal$vex2,], color="green")
  drawSegments(segments$coords[edgesOnOriginalEdges$vex1,], segments$coords[edgesOnOriginalEdges$vex2,], color="blue")
  drawSegments(segments$coords[edgesFromFaceIntersections$vex1,], segments$coords[edgesFromFaceIntersections$vex2,], color="purple")
  

  
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

