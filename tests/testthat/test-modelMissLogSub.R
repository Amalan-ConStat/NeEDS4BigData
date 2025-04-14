Beta<-c(-1,0.75,0.75,1); family <- "logistic"; N<-10000
X_1 <- replicate(2,stats::runif(n=N,min = -1,max = 1))

Temp<-Rfast::rowprods(X_1)
Misspecification <- (Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2)
X_Data <- cbind(X0=1,X_1);

Full_Data<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family)

r1<-300; r2<-rep(100*c(6,9),10); Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];

Results<-modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                         X = as.matrix(Original_Data[,-1]),
                         N = N,Alpha = 10, proportion=0.3) |> suppressWarnings() |> suppressMessages()

context_start_file("Checking the modelMissLogSub function")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),10)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2)*6,length(Beta)-1+2))
})

test_that("class of the Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates,"data.frame")
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,6))
})

test_that("class of the subsampling probability",{
  expect_s3_class(Results$Subsampling_Probability,"data.frame")
})

test_that("value in subsampling probability gte",{
  expect_gte(min(Results$Subsampling_Probability),0)
})

test_that("value in subsampling probability lte",{
  expect_lte(max(Results$Subsampling_Probability),1)
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(r2)+1))
})

test_that("dimension of the L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(r2)+1))
})

test_that("dimension of the L1-optimality sample",{
  expect_equal(length(Results$`Sample_L1-Optimality`),c(length(r2)+1))
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

context_start_file("Checking the modelMissLogSub for error messages")
test_that("Error on input for r2",{
  expect_error(modelMissLogSub(r1 = r1, r2 = c(r2,NA),Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 0.3),
               "NA or Infinite or NAN values in the r1,r2,N,Alpha or proportion")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 0.3),
               "NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)]
test_that("Error on r1 and r2 input",{
  expect_error(modelMissLogSub(r1 = r1+1000, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 0.3),
               "2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N+1,Alpha = 10, proportion = 0.3),
               "The big data size N is not the same as of the size of X or Y")
})

test_that("Error on Alpha the scaling factor",{
  expect_error(modelMissLogSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 0.1, proportion = 0.3),
               "Scaling factor alpha is not greater than one or the length is more than one")
})

test_that("Error on proportion value",{
  expect_error(modelMissLogSub(r1,r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 1.3),
               "Proportion should be a value higher than zero and less than or equal one")
})
test_that("Error on length of r1 or N",{
  expect_error(modelMissLogSub(r1,r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = c(N,2),Alpha = 10, proportion = 0.3),
               "proportion, r1 or N has a value greater than length one")
})
test_that("Error on model formula",{
  expect_error(modelMissLogSub(r1,r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 0.3, model=NULL),
               "The model formula for GAM is NA or NAN or NULL")
})
test_that("Wrong input on model formula",{
  expect_error(modelMissLogSub(r1,r2,Y = as.matrix(Original_Data[,1]),
                               X = as.matrix(Original_Data[,-1]),
                               N = N,Alpha = 10, proportion = 0.3, model="New"),
               "The model for GAM is not a valid formula or input")
})
