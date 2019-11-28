library(testthat)

source("utils.R")
source("geometry.R")

source("tests/test_utils.R")
source("tests/test_geometry.R")

test_results <- test_dir("tests", reporter="summary")

#print(test_results)