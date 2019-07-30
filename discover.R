# Discover new polyhedra

source("polyhedra.R")
source("draw.R")

discover <- function(p, debug=F)
{
  clear3d()
  
  # find unique distance pairs
  vexPairs <- combn(x = seq(nrow(p$vertices)), m = 2)
  vexPairs <- data.table(t(vexPairs), 
                         dist=sapply(seq(ncol(vexPairs)), function(col) { return(distance(p$vertices[vexPairs[1,col],], p$vertices[vexPairs[2,col],]))}))
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
    poly <- buildRegularPoly(p$vertices, searchgrid$facedim[i], 
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

# if (debug) print("All done!")
#discover(archi(dodecahedron), debug=F)

stop("issues below:")

discover(dodecahedron, debug=T) # problems
poly <- buildRegularPoly(dodecahedron$vertices, 5, 6, c(1, 4))
clear3d()
drawSinglePoly(poly, debug = T) # this shows the error drawing {5/2} # but now goes into a loop?
topo <- getTopology(poly$faces, debug=T) # this errors out on edge 4-1 / into a loop

# face is "below" origin
faces <- poly$faces
debug <- T
simplified <- list(vertices=poly$vertices, faces=list(poly$faces[[1]], poly$faces[[2]]), name="debug")
print(simplified$faces[[1]])
print(simplified$faces[[2]])
drawSinglePoly(simplified, debug = T) # shows drawing error and orientation issue when using [[1]],[[2]]
getTopology(simplified$faces) # loops

# with a fix we get into a loop in 
getTopology(simplified$faces, debug=T)

