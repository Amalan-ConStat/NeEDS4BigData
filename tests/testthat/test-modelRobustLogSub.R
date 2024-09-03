All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
                 Assumed_Model_1=c("X0","X1","X2"),
                 Assumed_Model_2=c("X0","X1","X2","X2^2"),
                 Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family = "logistic"
Full_Data<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family) |> suppressWarnings()
r1<-300; r2<-rep(100*c(6,9,12),25); Original_Data<-Full_Data$Complete_Data;

modelRobustLogSub(r1 = r1, r2 = r2,
                  Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                  Alpha = rep(1/length(All_Models),length(All_Models)),
                  All_Combinations = All_Models,
                  All_Covariates = colnames(Original_Data)[-1]) |> suppressWarnings() |>
  suppressMessages()->Results

context_start_file("Checking the modelRobustLogSub function")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),7)
})

test_that("length of the Beta Estimates",{
  expect_equal(length(Results$Beta_Estimates),length(All_Models))
})

test_that("dim of the Model 1 Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates$Model_1),c(length(r2)*4,2+length(All_Models$Real_Model)))
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,length(All_Models)*2+2))
})

test_that("dimension of the Model 1 A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`$Model_1),length(r2)+1)
})

test_that("dimension of the Model 1 L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`$Model_1),length(r2)+1)
})

test_that("dimension of the Model 1 MR A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality_MR`$Model_1),length(r2)+1)
})

test_that("dimension of the Model 1 MR L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality_MR`$Model_1),length(r2)+1)
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family = "logistic"
Full_Data<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family) |> suppressWarnings()
r1<-300; r2<-rep(100*c(6,9,12),25); Original_Data<-Full_Data$Complete_Data;

context_start_file("Checking the modelRobustLogSub for error messages")
test_that("Error on input for r2",{
  expect_error(modelRobustLogSub(r1 = r1, r2 = NA,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Alpha = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "NA or Infinite or NAN values in the r1,r2,N,Alpha or All_Covariates")
})
Original_Data[100,]<-rep(NA,6)
test_that("Error on X input",{
  expect_error(modelRobustLogSub(r1 = r1, r2 = r2,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Alpha = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r1 and r2 input",{
  expect_error(modelRobustLogSub(r1 = r1+1000, r2 = r2,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Alpha = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(modelRobustLogSub(r1 = r1, r2 = r2,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                                 Alpha = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "The big data size N is not the same as of the size of X or Y")
})

test_that("Error on length on Alpha the a priori probability",{
  expect_error(modelRobustLogSub(r1 = r1, r2 = r2,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Alpha = c(rep(1/length(All_Models),length(All_Models)),0.1),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "No of models for averaging is not equal to the a priori alpha values")
})

