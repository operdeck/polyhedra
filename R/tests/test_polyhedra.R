context("Polyhedra operations")

test_that( "Topology", {
  expect_equal( length(cube$faces), 6 ) # faces
  expect_equal( length(cube$coords), 8 ) # vertices
  expect_equal( nrow(cube$edgeToFaces), 12 ) # edges = faces + vertices - 2
  
  expect_equal( length(icosahedron$faces), 20 )
  expect_equal( length(icosahedron$coords), 12 )
  expect_equal( nrow(icosahedron$edgeToFaces), 30 )
  
})

test_that( "Description", {
  # platonic
  expect_equal(description(cube), "{4,3}")
  expect_equal(description(dual(dodecahedron)), "{3,5}")
  
  # stars
  expect_equal(description(greatIcosahedron), "{3,5/2}")
  expect_equal(description(dual(greatIcosahedron)), "{5/2,3}")
  
  # archimedean
  expect_equal(name(truncate(tetrahedron)), "truncated Tetrahedron")
  expect_equal(description(truncate(tetrahedron)), "{3,6,6}")
  expect_equal(description(truncate(icosahedron)), "{5,6,6}")
  expect_equal(description(quasi(cube)), "{3,4,3,4}")
  expect_equal(description(quasi(dodecahedron)), description(quasi(icosahedron)))
  expect_equal(description(rhombic(cube)), "{3,4,4,4}")
  expect_equal(description(rhombic(dodecahedron)), "{3,4,5,4}")
  
  # compounds
  expect_equal(description(compose(cube, dual(cube))), "1x{3,4} + 1x{4,3}")
  compound5tetrahedra <- buildRegularPoly(dodecahedron$coords,
                                          polygonsize = 3,
                                          vertexsize = 3,
                                          exampleEdge = c(3, 8),
                                          name = "5 Tetrahedra")
  expect_equal(description(compound5tetrahedra), "{3,3}X5")
  expect_equal(description(quasi(compound5tetrahedra)), "{3,4}X5")
  expect_equal(description(truncate(compound5tetrahedra)), "{3,6,6}X5") # beautiful one - hard to layout
})
