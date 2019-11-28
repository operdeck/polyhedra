context("Geometry")

test_that("segment intersection", {
  
  # two overlapping segments
  seg1 <- matrix( c( c(1, 0), c(3, 0) ), ncol=2, byrow = T )
  seg2 <- matrix( c( c(2, 0), c(4, 0) ), ncol=2, byrow = T )
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "overlap")
  expect_equal(isect$I0, c(2, 0))
  expect_equal(isect$I1, c(3, 0))
  
  # two segments bordering eachother
  seg1 <- matrix( c( c(2, 1), c(4, 2) ), ncol=2, byrow = T )
  seg2 <- matrix( c( c(5, 2.5), c(4, 2) ), ncol=2, byrow = T )
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(4, 2))
  
  # two segments intersecting eachother
  seg1 <- matrix( c( c(-1, -1), c(-6, -6) ), ncol=2, byrow = T )
  seg2 <- matrix( c( c(-3, 0), c(-3, -8) ), ncol=2, byrow = T )
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(-3, -3))
  
  # two segments that would intersect eachother if the first one was extended
  seg1 <- matrix( c( c(-1, -1), c(-2, -2) ), ncol=2, byrow = T )
  seg2 <- matrix( c( c(-3, 0), c(-3, -8) ), ncol=2, byrow = T )
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "disjoint")
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,], firstIsLine = T)
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(-3, -3))
  
  # two parallel segments that would never intersect eachother
  seg1 <- matrix( c( c(-1, -1), c(-6, -6) ), ncol=2, byrow = T )
  seg2 <- matrix( c( c(2, 7), c(-3, 2) ), ncol=2, byrow = T )
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,])
  expect_equal(isect$status, "disjoint")
  isect <- intersect2D_2Segments(seg1[1,], seg1[2,], seg2[1,], seg2[2,], firstIsLine = T)
  expect_equal(isect$status, "disjoint")
})