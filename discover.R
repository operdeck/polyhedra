# Discover new polyhedra

source("polyhedra.R")
source("draw.R")

discover <- function(p, debug=F)
{
  clear3d()
  
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
  vexPairs <- vexPairs[(vexPairs$isUnique),] # this is now a short list of vertex index pairs
  
  searchgrid <- expand.grid(facedim = 3:5, vertexdim = 3:20, vexPair = 1:nrow(vexPairs))
  searchgrid$isPolygon <- F
  searchgrid$vertex1 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$V1[searchgrid$vexPair[i]]})
  searchgrid$vertex2 <- sapply(seq(nrow(searchgrid)), function(i) {vexPairs$V2[searchgrid$vexPair[i]]})
  for (i in seq(nrow(searchgrid))) {
    # NB it also creates polys with holes in it
    # some of the ones from dodecahedron are new!
    s <- buildRegularPoly(p$coords, searchgrid$facedim[i], 
                             searchgrid$vertexdim[i], 
                             c(searchgrid$vertex1[i], searchgrid$vertex2[i]), debug=debug)
    if (length(poly$faces) > 0) {
      print("Success!")
      searchgrid$isPolygon[i] <- T
      print(searchgrid[i, ])
      
      #print(poly)
      drawSinglePoly(poly, debug=T)
      
      #descr <- description(poly)
      descr<-"issues w discover(dodecahedron) 3rd"
      print(descr)
      searchgrid$isPolygon[i] <- T
      
      drawSinglePoly(poly, x=2*sum(searchgrid$isPolygon), label=paste("Discovery",sum(searchgrid$isPolygon),"from",p$name), debug=F)
    }
  }
}

testDiscover <- function()
{
  discover(rhombic(dodecahedron), debug=T)
 
  discover(dodecahedron, debug=T) # problems 
}


