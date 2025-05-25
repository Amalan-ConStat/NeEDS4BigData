Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"logistic"
Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)

r0<-300; rf<-rep(100*c(6,9,12),10); Original_Data<-Full_Data$Complete_Data;

LCCsampling(r0 = r0, rf = rf, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
            X = as.matrix(Original_Data[,-1]),
            N = nrow(Original_Data)) |> suppressMessages()->Results

plot_Results<-plot_Beta(Results)

test_that("Testing for class error",{
  expect_error(plot_Beta(Full_Data),
               "object is not from any of the classes LocalCaseControl, Leverage,\nA_L_OptimalSubsampling, AoptimalSubsampling, A_OptimalSamplingMC, \nModelRobust or ModelMisspecified")
})

test_that("class of the output",{
  expect_s3_class(plot_Results,"gg")
})
