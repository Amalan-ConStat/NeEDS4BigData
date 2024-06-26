Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r1<-300; r2<-rep(c(600,1200),10); Original_Data<-Full_Data$Complete_Data;

ALoptimalGLMSub(r1 = r1, r2 = r2,
                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                family = "linear") |> suppressMessages()->Results

context_start_file("Checking the ALoptimalGLMSub function under linear regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),5)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2)*2,length(Beta)+2))
})

test_that("dimension of the Var Epsilon Estimates",{
  expect_equal(dim(Results$Variance_Epsilon_Estimates),c(length(r2)*2,3))
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,2))
})

test_that("dimension of the A-optimality subsample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(r2)+1))
})

test_that("dimension of the L-optimality subsample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(r2)+1))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r1<-300; r2<-rep(c(600,1200),10); Original_Data<-Full_Data$Complete_Data;

ALoptimalGLMSub(r1 = r1, r2 = r2,
                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                family = "logistic") |> suppressMessages()->Results

context_start_file("Checking the ALoptimalGLMSub function under logistic regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),5)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2)*2,length(Beta)+2))
})

test_that("dimension of the Utility Estimates",{
  expect_equal(dim(Results$Utility_Estimates),c(length(r2)*2,4))
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,2))
})

test_that("dimension of the A-optimality subsample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(r2)+1))
})

test_that("dimension of the L-optimality subsample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(r2)+1))
})

context_start_file("Checking the ALoptimalGLMSub for error messages")
test_that("Error on input for r2",{
  expect_error(ALoptimalGLMSub(r1 = r1, r2 = c(r2,NA),
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"NA or Infinite or NAN values in the r1,r2,N or family")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(ALoptimalGLMSub(r1 = r1, r2 = r2,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r1 and r2 input",{
  expect_error(ALoptimalGLMSub(r1 = r1+1000, r2 = r2,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(ALoptimalGLMSub(r1 = r1, r2 = r2,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                               family = "logistic"),"The big data size N is not the same as of the size of X or Y")
})

test_that("Error on character value of family",{
  expect_error(ALoptimalGLMSub(r1 = r1, r2 = r2,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "exp"),"Only the regression types 'linear','logistic' or 'poisson' are allowed")
})

