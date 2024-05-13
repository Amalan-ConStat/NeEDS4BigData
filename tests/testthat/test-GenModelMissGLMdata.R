No_Of_Var<-2; Beta<-c(-1,2,2,1); Var_Epsilon<-0.5; N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "linear"
Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)

context_start_file("Checking the GenModelMissGLMdata function under linear regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),7)
})
test_that("dim of Results$Real_Full_Data output",{
  expect_equal(dim(Results$Real_Full_Data),c(N,length(Beta)+1))
})
test_that("dimension of Results$Full_Data",{
  expect_equal(dim(Results$Full_Data),c(N,length(Beta)))
})
test_that("names of Results",{
  expect_equal(names(Results),c("N","Beta","Variance_Epsilon","Xbeta",
                                "f","Real_Full_Data","Full_Data"))
})

No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "logistic"
Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)

context_start_file("Checking the GenModelMissGLMdata function under logistic regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),6)
})
test_that("dim of Results$Real_Full_Data output",{
  expect_equal(dim(Results$Real_Full_Data),c(N,length(Beta)+1))
})
test_that("dimension of Results$Full_Data",{
  expect_equal(dim(Results$Full_Data),c(N,length(Beta)))
})
test_that("names of Results",{
  expect_equal(names(Results),c("N","Beta","Xbeta","f",
                                "Real_Full_Data","Full_Data"))
})

No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "poisson"
Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)

context_start_file("Checking the GenModelMissGLMdata function under Poisson regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),6)
})
test_that("dim of Results$Real_Full_Data output",{
  expect_equal(dim(Results$Real_Full_Data),c(N,length(Beta)+1))
})
test_that("dimension of Results$Full_Data",{
  expect_equal(dim(Results$Full_Data),c(N,length(Beta)))
})
test_that("names of Results",{
  expect_equal(names(Results),c("N","Beta","Xbeta","f",
                                "Real_Full_Data","Full_Data"))
})

No_Of_Var<-2; Beta<-c(-1,2,2,1); Var_Epsilon<-0.5; N<-10000;
MisspecificationType <- "Type 2 Squared"; family <- "linear"
GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,
                    MisspecificationType,family)

context_start_file("Checking the GenModelMissGLMdata function for error")
test_that("Error on Results output for Beta",{
  expect_error(GenModelMissGLMdata(No_Of_Var,c(-1,2,NA,1),Var_Epsilon,N,MisspecificationType,family),
               "NA or Infinite or NAN values in the No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType or family")
})

family<-"hello"
test_that("Error on Results output for Family",{
  expect_error(GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family),
               "Only the regression types 'linear','logistic' or 'poisson' are allowed")
})

family<-"linear"; MisspecificationType<-"Type 4"
test_that("Error on Results output for misspecification type",{
  expect_error(GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family),
               "Only the misspecification types 'Type 1', 'Type 2 Squared', 'Type 2 Interaction', \n 'Type 3 Squared', 'Type 3 Interaction' are allowed")
})
