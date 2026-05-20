context("Draw")

test_that( "Assign colors", {
  seqColor <- function(n) { return(seq_len(n)) }
  seqColor2 <- function(n) { return(-seq_len(n)) }
  
  # one body, one type of face
  expect_equal( sort(assignColors(cube, colorProvider = seqColor)), c(1,2,3,4,5,6) )
  
  # one body, multiple types of faces
  expect_equal( length(unique(assignColors(truncate(quasi(dodecahedron)), colorProvider = seqColor))), 3 )
  
  # two bodies, one type of face, single provider
  expect_equal( length(unique(assignColors(compose(tetrahedron, dual(tetrahedron)), colorProvider = seqColor))), 2 )
  
  # two bodies, one type of face, two providers
  expect_equal( length(unique(assignColors(compose(tetrahedron, dual(tetrahedron)), 
                                           colorProvider = list(seqColor, seqColor) ))), 4 )
  expect_equal( length(unique(assignColors(compose(tetrahedron, dual(tetrahedron)), 
                                           colorProvider = list(seqColor, seqColor2) ))), 8 )

  # two bodies, two types of faces  
  expect_equal( length(unique(assignColors(compose(dodecahedron, dual(dodecahedron)), 
                                           colorProvider = list(seqColor, seqColor2)))), 32 )
})

