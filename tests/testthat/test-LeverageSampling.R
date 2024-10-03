Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r<-100*c(6,10); Original_Data<-Full_Data$Complete_Data;

LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                 S_alpha = 0.95,
                 family = "linear") |> suppressMessages() ->Results

context_start_file("Checking the LeverageSampling function under linear regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),5)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r)*3,length(Beta)+2))
})

test_that("dimension of the Var Epsilon Estimates",{
  expect_equal(dim(Results$Variance_Epsilon_Estimates),c(length(r)*3,3))
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,2))
})

test_that("dimension of the Basic Leverage sample",{
  expect_equal(length(Results$Sample_Basic_Leverage),length(r))
})

test_that("dimension of the Shrinkage Leverage sample",{
  expect_equal(length(Results$Sample_Shrinkage_Leverage),length(r))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r<-100*c(6,10); Original_Data<-Full_Data$Complete_Data;

LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                 S_alpha = 0.95,
                 family = "logistic") |> suppressMessages()->Results

context_start_file("Checking the LeverageSampling function under logistic regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),4)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r)*3,length(Beta)+2))
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,2))
})

test_that("dimension of the Basic Leverage sample",{
  expect_equal(length(Results$Sample_Basic_Leverage),length(r))
})

test_that("dimension of the Shrinkage Leverage sample",{
  expect_equal(length(Results$Sample_Shrinkage_Leverage),length(r))
})

context_start_file("Checking the LeverageSampling for error messages")
test_that("Error on input for r",{
  expect_error(LeverageSampling(r = c(r,NA), Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                S_alpha = 0.95,
                                family = "logistic"),"NA or Infinite or NAN values in the r,S_alpha,N or family")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                S_alpha = 0.95,
                                family = "logistic"),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on size of X and Y",{
  expect_error(LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                                S_alpha = 0.95,
                                family = "logistic"),"The big data size N is not the same as of the size of X or Y")
})
test_that("Error on alphva value input",{
  expect_error(LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                S_alpha = 1.95,
                                family = "logistic"),"S_alpha value for shrinkage leverage scores are not in the range of zero and one or the length is more than one")
})
test_that("Error on character value of family",{
  expect_error(LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                S_alpha = 0.95,
                                family = "exp"),"Only the regression types 'linear','logistic' or 'poisson' are allowed")
})
