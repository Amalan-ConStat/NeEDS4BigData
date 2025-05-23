context_start_file("Checking the One million songs data")
test_that("Comparing the number of rows",{
  expect_equal(nrow(One_million_songs),205032)
})
test_that("Comparing the number of columns",{
  expect_equal(ncol(One_million_songs),6)
})
test_that("The class of the data",{
  expect_type(One_million_songs,"list")
})
test_that("Mean of response variable",{
  expect_equal(round(mean(One_million_songs$Counts),5),2.87932)
})
test_that("column names of One million songs data",{
  expect_equal(colnames(One_million_songs),c("Counts","Duration","Loudness",
                                              "Tempo","Artist_Hotness","Song_Hotness"))
})
