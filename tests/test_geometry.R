context("Geometry")

test_that("segment intersection", {
  
  p <- matrix(
    #  S1 p0    p1       S2 p0    01
    c( 1, 1, 1, 3, 3, 3, 0, 0, 2, 4, 4, 2,
       1, 1, 6, 3, 3, 8, 4, 4, 4, 5, 5, 5),
    ncol=3, byrow=T
  )
  
  # parallel lines/segments
  isect <- intersect_2Segments( p[1,], p[2,], p[5,], p[6,] )
  expect_equal(isect$status, "disjoint")
  
  # both are points
  isect <- intersect_2Segments( p[1,], p[1,], p[2,], p[2,] )
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments( p[5,], p[5,], p[5,], p[5,] )
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(1,1,6))
  
  # one is a point
  isect <- intersect_2Segments( p[1,], p[2,], p[7,], p[7,] )
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments( p[1,], p[2,], p[7,], p[7,], firstIsLine = T )
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(4,4,4))
  isect <- intersect_2Segments( p[1,], p[2,], p[6,], p[6,], firstIsLine = T )
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments( p[7,], p[7,], p[1,], p[2,] )
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments( p[7,], p[7,], p[1,], p[2,], firstIsLine = T )
  expect_equal(isect$status, "disjoint")
  
  # co-linear
  isect <- intersect_2Segments( p[1,], p[2,], p[7,], p[8,])
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments( c(1,0,0), c(6,0,0), c(3,0,0), c(4,0,0) )
  expect_equal(isect$status, "overlap")
  expect_equal(isect$I0, c(3,0,0))
  expect_equal(isect$I1, c(4,0,0))
  isect <- intersect_2Segments( c(1,0,0), c(3,0,0), c(3,0,0), c(4,0,0) )
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(3,0,0))
  
  isect <- intersect_2Segments( c(1,0,0), c(3,0,0), c(3,0,0), c(4,0,0), firstIsLine = T )
  expect_equal(isect$status, "overlap")
  expect_equal(isect$I0, c(3,0,0))
  expect_equal(isect$I1, c(4,0,0))
  
  isect <- intersect_2Segments( p[1,], p[2,], p[7,], p[8,], firstIsLine = T )
  expect_equal(isect$status, "overlap")
  expect_equal(isect$I0, c(4,4,4))
  expect_equal(isect$I1, c(5,5,5))
  
  isect <- intersect_2Segments( p[1,], p[2,], p[3,], p[4,] )
  expect_equal(isect$status, "intersect")
  expect_equal(isect$I0, c(2, 2, 2))
  
  # two overlapping segments
  # TODO can easiy add a test when not in the same plane by changing z=0 for one of them
  # overlap now can become intersect - not that trivial? Or is that addressed by choice of dim?
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
  expect_equal(isect$substatus, "parallel")
  
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
  
  # below some tests when one of the segments is just a single point
  mid <- c(1.672091e-17, 5.257311e-01, 2.008114e-01)
  a <- c(-1.672091e-17, 5.257311e-01, -2.008114e-01)
  b <- c(0.0000000, 0.5257311, -0.8506508)
  isect <- intersect_2Segments(mid, mid, a, b)
  expect_equal(isect$status, "disjoint")
  
  a <- c(1,0,0)
  b <- c(3,0,0)
  mid1 <- c(2,0,0)
  mid2 <- c(4,0,0)
  isect <- intersect_2Segments(mid1, mid1, a, b)
  expect_equal(isect$status, "intersect")
  isect <- intersect_2Segments(mid2, mid2, a, b)
  expect_equal(isect$status, "disjoint")
  a <- c(0,1,0)
  b <- c(0,3,0)
  mid1 <- c(0,2,0)
  mid2 <- c(0,4,0)
  isect <- intersect_2Segments(mid1, mid1, a, b)
  expect_equal(isect$status, "intersect")
  isect <- intersect_2Segments(mid2, mid2, a, b)
  expect_equal(isect$status, "disjoint")
  a <- c(0,0,1)
  b <- c(0,0,3)
  mid1 <- c(0,0,2)
  mid2 <- c(0,0,4)
  isect <- intersect_2Segments(mid1, mid1, a, b)
  expect_equal(isect$status, "intersect")
  isect <- intersect_2Segments(mid2, mid2, a, b)
  expect_equal(isect$status, "disjoint")
  
  
  p5 <- c(-0.5257311, -0.8506508, 0.0000000)
  p14 <- c(-5.257311e-01, -2.008114e-01, -1.110223e-16)
  p17 <- c(-5.257311e-01, 2.008114e-01, -1.110223e-16)
  p7 <- c(-0.5257311, 0.8506508, 0.0000000)
  isect <- intersect_2Segments(p17, p17, p5, p14)
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments(p17, p17, p14, p7)
  expect_equal(isect$status, "intersect")
  isect <- intersect_2Segments(p5, p14, p17, p17)
  expect_equal(isect$status, "disjoint")
  isect <- intersect_2Segments(p14, p7, p17, p17)
  expect_equal(isect$status, "intersect")
  # drawInit(T)
  # clear3d()
  # drawAxes()
  # drawDots(p5, "p5", color="red", radius=0.02)
  # drawDots(p14, "p14", color="red", radius=0.02)
  # drawDots(p17, "p17", color="black", radius=0.02)
  # drawDots(p7, "p7", color="red", radius=0.02)
  # drawSegments(p5, p14, color="blue")
  # drawSegments(p14, p7, color="green")
})