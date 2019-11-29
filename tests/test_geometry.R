context("Geometry")

test_that("segment intersection", {
  
  # two overlapping segments
  seg1 <- matrix( c( c(1, 0, 0), c(3, 0, 0) ), ncol=3, byrow = T )
  seg2 <- matrix( c( c(2, 0, 0), c(4, 0, 0) ), ncol=3, byrow = T )
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "overlap")
  expect_equal(isect$I0, c(2, 0, 0))
  expect_equal(isect$I1, c(3, 0, 0))
  
  # two segments bordering eachother
  seg1 <- matrix( c( c(0, 2, 1), c(0, 4, 2) ), ncol=3, byrow = T )
  seg2 <- matrix( c( c(0, 5, 2.5), c(0, 4, 2) ), ncol=3, byrow = T )
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(0, 4, 2))
  
  # two segments intersecting eachother
  seg1 <- matrix( c( c(-1, 0, -1), c(-6, 0, -6) ), ncol=3, byrow = T )
  seg2 <- matrix( c( c(-3, 0, 0), c(-3, 0, -8) ), ncol=3, byrow = T )
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(-3, 0, -3))
  
  # two segments that would intersect eachother if the first one was extended
  seg1 <- matrix( c( c(-1, -1, 0), c(-2, -2, 0) ), ncol=3, byrow = T )
  seg2 <- matrix( c( c(-3, 0, 0), c(-3, -8, 0) ), ncol=3, byrow = T )
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,], firstIsLine = T)
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(-3, -3, 0))
  
  # two parallel segments that would never intersect eachother
  seg1 <- matrix( c( c(0, -1, -1), c(0, -6, -6) ), ncol=3, byrow = T )
  seg2 <- matrix( c( c(0, 2, 7), c(0, -3, 2) ), ncol=3, byrow = T )
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,], firstIsLine = T)
  expect_equal(isect$status, "disjoint")
  
  # a line segment in a face of a great dodecahedron vs some of the segments
  line_p0 <- c(-0.5257311, -0.2008114 , 0.0000000 ) # obtained by F4 intersect with F1
  line_p1 <- c(-0.07851752, -0.47720462, -0.72360680 )
  p9 <- c(-0.8506508,  0.0000000, -0.5257311)
  p11 <- c(-0.8506508,  0.0000000,  0.5257311 )
  p3 <- c(0.0000000, -0.5257311,  0.8506508 )
  p6 <- c(0.5257311, -0.8506508,  0.0000000)
  p2 <- c( 0.0000000,  0.5257311, -0.8506508  )
  isect <- intersect_2Segments(line_p0, line_p1, p9, p11, firstIsLine = T)
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, p11, tolerance=1e-6)
  isect <- intersect_2Segments(line_p0, line_p1, p9, p3, firstIsLine = T)
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(-5.257311e-01, -2.008114e-01, 0), tolerance=1e-6)
  isect <- intersect_2Segments(line_p0, line_p1, p6, p3, firstIsLine = T)
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments(line_p0, line_p1, p11, p3, firstIsLine = T)
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, p11, tolerance=1e-6)
  isect <- intersect_2Segments(line_p0, line_p1, p2, p9, firstIsLine = T)
  expect_equal(isect$status, "intersect") # TODO this should FAIL, totally different planes
  # need to check both expressions before returning anything
})