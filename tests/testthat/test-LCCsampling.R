Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"logistic"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r1<-300; r2<-rep(100*c(6,9,12),50); Original_Data<-Full_Data$Complete_Data;

LCCsampling(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
            X = as.matrix(Original_Data[,-1]),
            N = nrow(Original_Data)) |> suppressMessages()->Results

context_start_file("Checking the LCCsampling function under logistic regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),4)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2),length(Beta)+2))
})

test_that("class of the Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates,"data.frame")
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,1))
})

test_that("class of the sampling probability",{
  expect_s3_class(Results$Sampling_Probability,"data.frame")
})

test_that("value in sampling probability gte",{
  expect_gte(min(Results$Sampling_Probability),0)
})

test_that("value in sampling probability lte",{
  expect_lte(max(Results$Sampling_Probability),1)
})

test_that("dimension of the LCC sample",{
  expect_equal(length(Results$Sample_LCC_Sampling),c(length(r2)+1))
})

context_start_file("Checking the LCCsampling for error messages")
test_that("Error on input for r2",{
  expect_error(LCCsampling(r1 = r1, r2 = c(r2,NA), Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = nrow(Original_Data)),"NA or Infinite or NAN values in the r1,r2 or N")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(LCCsampling(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = nrow(Original_Data)),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r1 and r2 input",{
  expect_error(LCCsampling(r1 = r1+1000, r2 = r2, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = nrow(Original_Data)),"2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(LCCsampling(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = nrow(Original_Data)+1),"The big data size N is not the same as of the size of X or Y")
})
test_that("Error on length of r1, N or family",{
  expect_error(LCCsampling(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                           X = as.matrix(Original_Data[,-1]),
                           N = c(nrow(Original_Data),2)),"r1 or N has a value greater than length one")
})
