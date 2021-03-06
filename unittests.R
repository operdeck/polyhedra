library(testthat)

source("utils.R")
source("geometry.R")
source("polyhedra.R")
source("layout.R")

source("tests/test_utils.R")
source("tests/test_geometry.R")
source("tests/test_polyhedra.R")
source("tests/test_draw.R")
source("tests/test_layout.R")

test_results <- test_dir("tests", reporter="summary")

#print(test_results)