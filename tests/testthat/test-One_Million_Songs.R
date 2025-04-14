context_start_file("Checking the One Million Songs data")
test_that("Comparing the number of rows",{
  expect_equal(nrow(One_Million_Songs),205032)
})
test_that("Comparing the number of columns",{
  expect_equal(ncol(One_Million_Songs),6)
})
test_that("The class of the data",{
  expect_type(One_Million_Songs,"list")
})
test_that("Mean of response variable",{
  expect_equal(round(mean(One_Million_Songs$Counts),5),2.87932)
})
test_that("column names of One Million Songs data",{
  expect_equal(colnames(One_Million_Songs),c("Counts","Duration","Loudness",
                                              "Tempo","Artist_Hotness","Song_Hotness"))
})
