context("2D Layout routines")

test_that( "segmentation", {
  seg <- segmentation(compose(tetrahedron, dual(tetrahedron)), debug=F)
  expect_equal(36, nrow(seg$edges))
  expect_equal(24, sum(seg$edges$onOriginal)) # 8 pyramids of 3 sides
  expect_equal(8, sum(rownames(seg$coords)!="")) # 8 original vertices
  expect_equal(6, sum(rownames(seg$coords)=="")) # 6 new vertices
  
  seg <- segmentation(quasi(cube)) # identical to original
  expect_equal(24, nrow(seg$edges))
  expect_equal(24, sum(seg$edges$onOriginal)) 
  expect_equal(12, sum(rownames(seg$coords)!=""))
  expect_equal(0, sum(rownames(seg$coords)==""))
})

