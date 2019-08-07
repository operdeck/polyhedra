# Discover new polyhedra

source("polyhedra.R")
source("draw.R")

discover <- function(p, progress=T, debug=F)
{
  discoveries <- list()
  
  # find unique distance pairs
  vexPairs <- combn(x = seq(nrow(p$coords)), m = 2)
  vexPairs <- data.table(t(vexPairs), 
                         dist=sapply(seq(ncol(vexPairs)), function(col) { return(distance(p$coords[vexPairs[1,col],], p$coords[vexPairs[2,col],]))}))
  vexPairs <- vexPairs[order(vexPairs$dist),]
  vexPairs$isUnique <- T
  row.names(vexPairs) <- NULL
  for (i in 2:nrow(vexPairs)) {
    if (deltaEquals(vexPairs$dist[i-1], vexPairs$dist[i])) { vexPairs$isUnique[i] <- F }
  }
  vexPairs <- vexPairs[(vexPairs$isUnique),][, -"isUnique"] # this is now a short list of vertex index pairs
  
  searchgrid <- expand.grid(facedim = 3:5, vertexdim = 3:20, vexPair = 1:nrow(vexPairs))
  searchgrid$isValid <- F
  searchgrid$description <- ""
  searchgrid$vertex1 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$V1[searchgrid$vexPair[i]]})
  searchgrid$vertex2 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$V2[searchgrid$vexPair[i]]})
  for (i in seq(nrow(searchgrid))) {
    if(debug) clear3d()
    
    # NB it also creates polys with holes in it
    # some of the ones from dodecahedron are new!
    s <- buildRegularPoly(p$coords, searchgrid$facedim[i], 
                             searchgrid$vertexdim[i], 
                             c(searchgrid$vertex1[i], searchgrid$vertex2[i]), debug=debug)
    if (!is.null(s)) {
      searchgrid$isValid[i] <- T
      searchgrid$description[i] <- description(s)
      
      if (progress | debug) {
        print(searchgrid[i, ])
      }
      
      discoveries[[length(discoveries)+1]] <- s
      
      if (debug) {
        clear3d()
        drawSinglePoly(s, debug=T)
      }
    }
  }
  return(list(specs=searchgrid[(searchgrid$isValid),], polyhedra=discoveries))
}

testDiscover <- function()
{
  discs <- discover(cube)
  clear3d()
  drawPoly(discs$polyhedra)
  
  discs <- discover(icosahedron, debug=F)
  clear3d()
  drawPoly(discs$polyhedra)

  discs <- discover(dodecahedron, debug=F)
  clear3d()
  drawPoly(discs$polyhedra)
  
  discs <- discover(rhombic(dodecahedron))
  clear3d()
  drawPoly(discs$polyhedra) # why not 5 x tetrahedra part of this?
  
  discover(dodecahedron, debug=T) # problems 
}


