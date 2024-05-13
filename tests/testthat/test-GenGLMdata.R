Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

context_start_file("Checking the GenGLMdata function under linear regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("length of Results$Basic output",{
  expect_equal(length(Results$Basic),7)
})
test_that("dimension of Results$Complete_Data",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Variance_Epsilon_Estimates",
                                      "Distribution","Distribution_Parameter","No_Of_Variables"))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

context_start_file("Checking the GenGLMdata function under logistic regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("length of Results$Basic output",{
  expect_equal(length(Results$Basic),6)
})
test_that("dimension of Results$Complete_Data",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Distribution",
                                      "Distribution_Parameter","No_Of_Variables"))
})

Dist<-"Normal";
No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"poisson"
Results<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)

context_start_file("Checking the GenGLMdata function under Poisson regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("length of Results$Basic output",{
  expect_equal(length(Results$Basic),6)
})
test_that("dimension of Results$Complete_Data",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+1))
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Distribution",
                                      "Distribution_Parameter","No_Of_Variables"))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,NA,1); N<-5000; Family<-"linear"

context_start_file("Checking the GenGLMdata function for error")
test_that("Error on Results output for Beta",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family),
               "NA or Infinite or NAN values in the Dist,Beta,No_Of_Var,N or family")
})

Family<-"hello"
test_that("Error on Results output for Family",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family),
               "NA or Infinite or NAN values in the Dist,Beta,No_Of_Var,N or family")
})

Family<-"linear"; Dist<-"Exp"; Beta<-c(-1,2,1)
test_that("Error on Results output for linear Dist",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family),
               "For linear regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
})

Family<-"logistic"; Dist<-"Exp"
test_that("Error on Results output for logistic Dist",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family),
               "For logistic regression select the distribution 'Exponential', \n 'Normal' or 'Uniform' to generate the covarate data")
})

Family<-"poisson"; Dist<-"Exp"
test_that("Error on Results output for poisson Dist",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family),
               "For poisson regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
})
