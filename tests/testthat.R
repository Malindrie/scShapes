Sys.setenv("R_TESTS" = "")

library(testthat)
library(scShapes)

test_check("scShapes")
