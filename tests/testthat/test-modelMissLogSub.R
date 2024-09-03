No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "logistic"
Full_Data<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)
r1<-300; r2<-100*c(6,9); Original_Data<-Full_Data$Full_Data;

Results<-modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                         X = as.matrix(Original_Data[,-1]),
                         N = Full_Data$N,Alpha = 10,
                         Beta_Estimate_Full = Full_Data$Beta$Estimate,
                         F_Estimate_Full = Full_Data$f$Real_GAM) |> suppressWarnings() |> suppressMessages()

context_start_file("Checking the modelMissLogSub function")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),8)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2)*5,length(Beta)-1+2))
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,5))
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(r2)+1))
})

test_that("dimension of the L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(r2)+1))
})

test_that("dimension of the RLmAMSE sample",{
  expect_equal(length(Results$Sample_RLmAMSE),c(length(r2)+1))
})

test_that("dimension of the RLmAMSE Log Odds sample",{
  expect_equal(length(Results$Sample_RLmAMSE_Log_Odds[[1]]),c(length(r2)+1))
})

test_that("dimension of the RLmAMSE Power sample",{
  expect_equal(length(Results$Sample_RLmAMSE_Power[[1]]),c(length(r2)+1))
})

No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "logistic"
Full_Data<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)
r1<-300; r2<-100*c(6,9); Original_Data<-Full_Data$Full_Data;

context_start_file("Checking the modelMissLogSub for error messages")
test_that("Error on input for r2",{
  expect_error(modelMissLogSub(r1 = r1, r2 = c(r2,NA),Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = Full_Data$N,Alpha = 10,
                               Beta_Estimate_Full = Full_Data$Beta$Estimate,
                               F_Estimate_Full = Full_Data$f$Real_GAM),
               "NA or Infinite or NAN values in the r1,r2,N or Alpha")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = Full_Data$N,Alpha = 10,
                               Beta_Estimate_Full = Full_Data$Beta$Estimate,
                               F_Estimate_Full = Full_Data$f$Real_GAM),
               "NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Full_Data
test_that("Error on r1 and r2 input",{
  expect_error(modelMissLogSub(r1 = r1+1000, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = Full_Data$N,Alpha = 10,
                               Beta_Estimate_Full = Full_Data$Beta$Estimate,
                               F_Estimate_Full = Full_Data$f$Real_GAM),
               "2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = Full_Data$N+1,Alpha = 10,
                               Beta_Estimate_Full = Full_Data$Beta$Estimate,
                               F_Estimate_Full = Full_Data$f$Real_GAM),
               "The big data size N is not the same as of the size of X or Y")
})

test_that("Error on Alpha the scaling factor",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = Full_Data$N,Alpha = 0.1,
                               Beta_Estimate_Full = Full_Data$Beta$Estimate,
                               F_Estimate_Full = Full_Data$f$Real_GAM),
               "Scaling factor alpha is not greater than one or the length is more than one")
})
