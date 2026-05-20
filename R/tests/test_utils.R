context("R utils")

test_that( "Shift rotate", {
  expect_equal(shiftrotate(seq(3)), c(2,3,1))
  expect_equal(shiftrotate(seq(3),-1), c(3,1,2))
  expect_equal(shiftrotate(seq(3),2), c(3,1,2))
})

test_that( "Safe sequence", {
  expect_equal(length(safeseq(3)), 3)
  expect_equal(length(seq(0)), 2) # standard seq
  expect_equal(length(safeseq(0)), 0)
})

