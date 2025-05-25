indexes<-1:ceiling(nrow(One_million_songs)*0.5)
Original_Data<-One_million_songs[indexes,]
colnames(Original_Data)<-c("Y",paste0("X",1:ncol(Original_Data[,-1])))

# Scaling the covariate data
for (j in 2:4) {
  Original_Data[,j]<-scale(Original_Data[,j])
}

No_of_Variables<-ncol(Original_Data[,-1])
Squared_Terms<-paste0("X",1:No_of_Variables,"^2")
term_no <- 2
All_Models <- list(paste0("X",1:No_of_Variables))

Original_Data<-cbind(Original_Data,Original_Data[,-1]^2)
colnames(Original_Data)<-c("Y",paste0("X",1:No_of_Variables),
                            paste0("X",1:No_of_Variables,"^2"))

for (i in 1:No_of_Variables)
{
  x <- as.vector(combn(Squared_Terms,i,simplify = FALSE))
  for(j in 1:length(x))
  {
    All_Models[[term_no]] <- c(paste0("X",1:No_of_Variables),x[[j]])
    term_no <- term_no+1
  }
}
All_Models<-All_Models[c(1,27:32)]
names(All_Models)<-paste0("Model_",1:length(All_Models))
r0<-300; rf<-rep(100*c(6,9,12),5)

modelRobustPoiSub(r0 = r0, rf = rf,
                  Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                  Apriori_probs = rep(1/length(All_Models),length(All_Models)),
                  All_Combinations = All_Models,
                  All_Covariates = colnames(Original_Data)[-1]) |> suppressWarnings() |>
  suppressMessages()->Results

context_start_file("Checking the modelRobustPoiSub function")
test_that("type of the Results",{
  expect_type(Results,"list")
})

test_that("length of the Results output",{
  expect_equal(length(Results),6)
})

test_that("length of the Beta Estimates",{
  expect_equal(length(Results$Beta_Estimates),length(All_Models))
})

test_that("dim of the Model 1 Beta Estimates",{
  expect_equal(dim(Results$Beta_Estimates$Model_1),c(length(rf)*4,2+length(All_Models$Model_1)))
})

test_that("class of Model 1 Beta Estimates",{
  expect_s3_class(Results$Beta_Estimates$Model_1,"data.frame")
})

test_that("dimension of the subsampling probability",{
  expect_equal(dim(Results$Subsampling_Probability),c(nrow(Original_Data),length(All_Models)*2+2))
})

test_that("class of subsampling probability",{
  expect_s3_class(Results$Subsampling_Probability,"data.frame")
})

test_that("dimension of the Model 1 A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality`$Model_1),length(rf)+1)
})

test_that("dimension of the Model 1 L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality`$Model_1),length(rf)+1)
})

test_that("dimension of the Model 1 MR A-optimality sample",{
  expect_equal(length(Results$`Sample_A-Optimality_MR`$Model_1),length(rf)+1)
})

test_that("dimension of the Model 1 MR L-optimality sample",{
  expect_equal(length(Results$`Sample_L-Optimality_MR`$Model_1),length(rf)+1)
})

context_start_file("Checking the modelRobustPoiSub for error messages")
test_that("Error on input for r",{
  expect_error(modelRobustPoiSub(r0 = r0, rf = NA,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Apriori_probs = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "NA or Infinite or NAN values in the r0,rf,N,Apriori_probs or All_Covariates")
})
Temp_Data<-Original_Data
Temp_Data[100,]<-rep(NA,11)
test_that("Error on X input",{
  expect_error(modelRobustPoiSub(r0 = r0, rf = rf,
                                 Y = as.matrix(Temp_Data[,colnames(Temp_Data) %in% c("Y")]),
                                 X = as.matrix(Temp_Data[,-1]),N = nrow(Temp_Data),
                                 Apriori_probs = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Temp_Data)[-1]),
               "NA or Infinite or NAN values in the Y or X")
})

test_that("Error on r0 and r input",{
  expect_error(modelRobustPoiSub(r0 = r0+1000, rf = rf,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Apriori_probs = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "2*r0 cannot be greater than rf at any point")
})
test_that("Error on size of X and Y",{
  expect_error(modelRobustPoiSub(r0 = r0, rf = rf,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data)+1,
                                 Apriori_probs = rep(1/length(All_Models),length(All_Models)),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "The big data size N is not the same as of the size of X or Y")
})

test_that("Error on length on Apriori_probs the a priori probability",{
  expect_error(modelRobustPoiSub(r0 = r0, rf = rf,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Apriori_probs = c(rep(1/length(All_Models),length(All_Models)),0.1),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "No of models for averaging is not equal to the a priori probabilities")
})
test_that("Error on length of r0 or N",{
  expect_error(modelRobustPoiSub(r0 = c(r0,3), rf = rf,
                                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
                                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
                                 Apriori_probs = c(rep(1/length(All_Models),length(All_Models)),0.1),
                                 All_Combinations = All_Models,
                                 All_Covariates = colnames(Original_Data)[-1]),
               "r0 or N has a value greater than length one")
})
