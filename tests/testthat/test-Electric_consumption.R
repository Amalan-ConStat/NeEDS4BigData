context_start_file("Checking the Electric consumption data")
test_that("Comparing the number of rows",{
  expect_equal(nrow(Electric_consumption),2049280)
})
test_that("Comparing the number of columns",{
  expect_equal(ncol(Electric_consumption),4)
})
test_that("The class of the data",{
  expect_type(Electric_consumption,"list")
})
test_that("Mean of response variable",{
  expect_equal(round(mean(Electric_consumption$Intensity),5),1.10033)
})
