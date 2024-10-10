context_start_file("Checking the Bike sharing data")
test_that("Comparing the number of rows",{
  expect_equal(nrow(Bike_sharing),17379)
})
test_that("Comparing the number of columns",{
  expect_equal(ncol(Bike_sharing),4)
})
test_that("The class of the data",{
  expect_type(Bike_sharing,"list")
})
test_that("Number of response variables",{
  expect_equal(length(unique(Bike_sharing$Rented_Bikes)),869)
})
test_that("column names of bike sharing data",{
  expect_equal(colnames(Bike_sharing),c("Rented_Bikes","Temperature","Humidity","Windspeed"))
})
