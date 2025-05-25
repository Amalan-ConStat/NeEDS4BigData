Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r0<-300; rf<-rep(c(600,1200),10); Original_Data<-Full_Data$Complete_Data;

ALoptimalGLMSub(r0 = r0, rf = rf,
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
  expect_equal(dim(Results$Beta_Estimates),c(length(rf)*2,length(Beta)+2))
})

test_that("class of the Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates,"data.frame")
})

test_that("dimension of the Var Epsilon Estimates",{
  expect_equal(dim(Results$Variance_Epsilon_Estimates),c(length(rf)*2,3))
})

test_that("class of the Var Epsilon Estimates",{
  expect_s3_class(Results$Variance_Epsilon_Estimates,"data.frame")
})

test_that("value in variance epsilon estimates",{
  expect_gt(min(Results$Variance_Epsilon_Estimates$`Var Epsilon`),0)
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,2))
})

test_that("class of the subsampling probability",{
  expect_s3_class(Results$Subsampling_Probability,"data.frame")
})

test_that("value in subsampling probability gte A-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`A-Optimality`),0)
})

test_that("value in subsampling probability lte A-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`A-Optimality`),1)
})

test_that("value in subsampling probability gte L-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`L-Optimality`),0)
})

test_that("value in subsampling probability lte L-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`L-Optimality`),1)
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(rf)+1))
})

test_that("dimension of the L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(rf)+1))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r0<-300; rf<-rep(c(600,1200),10); Original_Data<-Full_Data$Complete_Data;

ALoptimalGLMSub(r0 = r0, rf = rf,
                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                family = "logistic") |> suppressMessages()->Results

context_start_file("Checking the ALoptimalGLMSub function under logistic regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),4)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(rf)*2,length(Beta)+2))
})

test_that("class of the Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates,"data.frame")
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,2))
})

test_that("class of the subsampling probability",{
  expect_s3_class(Results$Subsampling_Probability,"data.frame")
})

test_that("value in subsampling probability gte A-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`A-Optimality`),0)
})

test_that("value in subsampling probability lte A-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`A-Optimality`),1)
})

test_that("value in subsampling probability gte L-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`L-Optimality`),0)
})

test_that("value in subsampling probability lte L-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`L-Optimality`),1)
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(rf)+1))
})

test_that("dimension of the L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(rf)+1))
})

context_start_file("Checking the ALoptimalGLMSub for error messages")
test_that("Error on input for r",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = c(rf,NA),
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"NA or Infinite or NAN values in the r0,rf,N or family")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r0 and r input",{
  expect_error(ALoptimalGLMSub(r0 = r0+1000, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "logistic"),"2*r0 cannot be greater than rf at any point")
})
test_that("Error on size of X and Y",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                               family = "logistic"),"The big data size N is not the same as of the size of X or Y")
})

test_that("Error on character value of family",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "exp"),"Only the regression types 'linear','logistic' or 'poisson' are allowed")
})
test_that("Error on length of r0, N or family",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = c("linear","logistic") ),"r0, N or family has a value greater than length one")
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"poisson"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r0<-300; rf<-rep(c(600,1200),10); Original_Data<-Full_Data$Complete_Data;

ALoptimalGLMSub(r0 = r0, rf = rf,
                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                family = "poisson") |> suppressMessages()->Results

context_start_file("Checking the ALoptimalGLMSub function under poisson regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),4)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(rf)*2,length(Beta)+2))
})

test_that("class of the Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates,"data.frame")
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(N,2))
})

test_that("class of the subsampling probability",{
  expect_s3_class(Results$Subsampling_Probability,"data.frame")
})

test_that("value in subsampling probability gte A-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`A-Optimality`),0)
})

test_that("value in subsampling probability lte A-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`A-Optimality`),1)
})

test_that("value in subsampling probability gte L-Opt",{
  expect_gte(min(Results$Subsampling_Probability$`L-Optimality`),0)
})

test_that("value in subsampling probability lte L-Opt",{
  expect_lte(max(Results$Subsampling_Probability$`L-Optimality`),1)
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(rf)+1))
})

test_that("dimension of the L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`),c(length(rf)+1))
})

context_start_file("Checking the ALoptimalGLMSub for error messages")
test_that("Error on input for r",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = c(rf,NA),
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "poisson"),"NA or Infinite or NAN values in the r0,rf,N or family")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "poisson"),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r0 and r input",{
  expect_error(ALoptimalGLMSub(r0 = r0+1000, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "poisson"),"2*r0 cannot be greater than rf at any point")
})
test_that("Error on size of X and Y",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                               family = "poisson"),"The big data size N is not the same as of the size of X or Y")
})

test_that("Error on character value of family",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = "exp"),"Only the regression types 'linear','logistic' or 'poisson' are allowed")
})
test_that("Error on length of r0, N or family",{
  expect_error(ALoptimalGLMSub(r0 = r0, rf = rf,
                               Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                               X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                               family = c("linear","poisson") ),"r0, N or family has a value greater than length one")
})
