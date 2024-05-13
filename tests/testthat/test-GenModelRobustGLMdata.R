All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
                 Assumed_Model_1=c("X0","X1","X2"),
                 Assumed_Model_2=c("X0","X1","X2","X2^2"),
                 Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family<-"linear"
Results<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)

context_start_file("Checking the GenModelRobustGLMdata function under linear regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("dim of Results$Complete_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+2))
})
test_that("length of Results$Basic",{
  expect_equal(length(Results$Basic),8)
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Variance_Epsilon_Estimates",
                                      "Distribution","Distribution_Parameter","No_Of_Variables","All_Models"))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family = "logistic"
Results<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family) |> suppressWarnings()

context_start_file("Checking the GenModelRobustGLMdata function under logistic regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("dim of Results$Complete_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+2))
})
test_that("length of Results$Basic",{
  expect_equal(length(Results$Basic),7)
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Distribution",
                                      "Distribution_Parameter","No_Of_Variables","All_Models"))
})

Dist<-"Normal";
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family = "poisson"
Results<-GenModelRobustGLMdata(Dist,Dist_Par=NULL,No_Of_Var,Beta,N,All_Models,family)

context_start_file("Checking the GenModelRobustGLMdata function under Poisson regression")
test_that("type of the Results output",{
  expect_type(Results,"list")
})
test_that("length of Results output",{
  expect_equal(length(Results),2)
})
test_that("dim of Results$Complete_Data output",{
  expect_equal(dim(Results$Complete_Data),c(N,length(Beta)+2))
})
test_that("length of Results$Basic",{
  expect_equal(length(Results$Basic),6)
})
test_that("names of Results$Basic",{
  expect_equal(names(Results$Basic),c("N","Beta","Beta_Estimates","Distribution",
                                      "No_Of_Variables","All_Models"))
})

Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000; family<-"linear"

context_start_file("Checking the GenModelRobustGLMdata function for error")
test_that("Error on Results output for Beta",{
  expect_error(GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta=c(-1,2,NA,1),N,All_Models,family),
               "NA or Infinite or NAN values in the Dist,Beta,No_Of_Var,N or family")
})

family<-"hello"
test_that("Error on Results output for Family",{
  expect_error(GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family),
               "Only the regression types 'linear','logistic' or 'poisson' are allowed")
})

family<-"linear"; Dist<-"Exp";
test_that("Error on Results output for linear Dist",{
  expect_error(GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,family),
               "For linear regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
})

family<-"logistic"; Dist<-"Exp"
test_that("Error on Results output for logistic Dist",{
  expect_error(GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family),
               "For logistic regression select the distribution 'Exponential', \n 'Normal' or 'Uniform' to generate the covarate data")
})

family<-"poisson"; Dist<-"Exp"
test_that("Error on Results output for poisson Dist",{
  expect_error(GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family),
               "For poisson regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
})
