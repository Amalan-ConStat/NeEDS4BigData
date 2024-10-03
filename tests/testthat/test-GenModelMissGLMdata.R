Beta<-c(-1,0.75,0.75,1); Var_Epsilon<-0.5; family <- "linear"; N<-10000
X_1 <- replicate(2,stats::runif(n=N,min = -1,max = 1))

Temp<-Rfast::rowprods(X_1)
Misspecification <- (Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2)

X_Data <- cbind(X0=1,X_1);

Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family)

context_start_file("Checking the GenModelMissGLMdata function under linear regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),1)
})
test_that("dim of Results$Real_Full_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})

Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family="logistic")

context_start_file("Checking the GenModelMissGLMdata function under logistic regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),1)
})
test_that("dim of Results$Complete_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})

Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family="poisson")

context_start_file("Checking the GenModelMissGLMdata function under Poisson regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),1)
})
test_that("dim of Results$Complete_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})



context_start_file("Checking the GenModelMissGLMdata function for error")
test_that("Error on Results output for Beta",{
  expect_error(GenModelMissGLMdata(N,X_Data,Misspecification,Beta=c(-1,NA,0.75,1),Var_Epsilon,family),
               "NA or Infinite or NAN values in the Misspecification,Beta,Var_Epsilon or family")
})

family<-"hello"
test_that("Error on Results output for Family",{
  expect_error(GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family),
               "Only the regression types 'linear','logistic' or 'poisson' are allowed")
})
