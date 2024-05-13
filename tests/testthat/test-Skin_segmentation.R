context_start_file("Checking the Skin segmentation data")
test_that("Comparing the number of rows",{
  expect_equal(nrow(Skin_segmentation),245057)
})
test_that("Comparing the number of columns",{
  expect_equal(ncol(Skin_segmentation),4)
})
test_that("The class of the data",{
  expect_type(Skin_segmentation,"list")
})
test_that("Mean of response variable",{
  expect_equal(round(mean(Skin_segmentation$Skin_presence),1),0.8)
})

