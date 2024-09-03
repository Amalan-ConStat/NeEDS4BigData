Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"linear"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r1<-300; r2<-rep(100*c(6,9,12),50); Original_Data<-Full_Data$Complete_Data;

AoptimalGauLMSub(r1 = r1, r2 = r2,
                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                 X = as.matrix(Original_Data[,-1]),
                 N = nrow(Original_Data)) |> suppressMessages()->Results

context_start_file("Checking the AoptimalGauLMSub function under linear regression")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),4)
})

test_that("dimension of the Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates),c(length(r2),length(Beta)+2))
})

test_that("dimension of the Var Epsilon Estimates",{
  expect_equal(dim(Results$Variance_Epsilon_Estimates),c(length(r2),3))
})

test_that("dimension of the sampling probability",{
  expect_equal(dim(Results$Sampling_Probability),c(N,1))
})

test_that("dimension of the A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`),c(length(r2)+1))
})

context_start_file("Checking the AoptimalGauLMSub for error messages")
test_that("Error on input for r2",{
  expect_error(AoptimalGauLMSub(r1 = r1, r2 = c(r2,NA),
                                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),
                                N = nrow(Original_Data)),"NA or Infinite or NAN values in the r1,r2 or N")
})
Original_Data[100,]<-rep(NA,4)
test_that("Error on X input",{
  expect_error(AoptimalGauLMSub(r1 = r1, r2 = r2,
                                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),
                                N = nrow(Original_Data)),"NA or Infinite or NAN values in the Y or X")
})
Original_Data<-Full_Data$Complete_Data
test_that("Error on r1 and r2 input",{
  expect_error(AoptimalGauLMSub(r1 = r1+1000, r2 = r2,
                                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),
                                N = nrow(Original_Data)),"2*r1 cannot be greater than r2 at any point")
})
test_that("Error on size of X and Y",{
  expect_error(AoptimalGauLMSub(r1 = r1, r2 = r2,
                                Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                X = as.matrix(Original_Data[,-1]),
                                N = nrow(Original_Data)+1),"The big data size N is not the same as of the size of X or Y")
})
