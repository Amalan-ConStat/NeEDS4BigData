Beta<-c(-1,0.75,0.75,1); Var_Epsilon<-0.5; family <- "linear"; N<-10000
X_1 <- replicate(2,stats::runif(n=N,min = -1,max = 1))

Temp<-Rfast::rowprods(X_1)
Misspecification <- (Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2)
X_Data <- cbind(X0=1,X_1);

Full_Data<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family)

r1<-300; r2<-rep(100*c(6,9),10); Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];

Results<-modelMissLinSub(r1 = r1, r2 = r2,
                         Y = as.matrix(Original_Data[,1]),
                         X = as.matrix(Original_Data[,-1]),
                         N = N, Alpha = 10, proportion = 0.3) |> suppressWarnings() |> suppressMessages()

plot_Results<-plot_AMSE(Results)

test_that("Testing for class error",{
  expect_error(plot_AMSE(Full_Data),
               "object is not from any of the classes ModelMisspecified")
})

test_that("class of the output",{
  expect_s3_class(plot_Results,"gg")
})
