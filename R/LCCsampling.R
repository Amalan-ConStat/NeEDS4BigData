#' Local case control sampling for logistic regression
#'
#' Using this function sample from big data under logistic regression to describe the data.
#' Sampling probabilities are obtained based on local case control method.
#'
#' @usage
#' LCCsampling(r0,rf,Y,X,N)
#'
#' @param r0      sample size for initial random sample
#' @param rf      final sample size including initial(r0) and case control(r) samples
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#'
#' @details
#' Two stage sampling algorithm for big data under logistic regression.
#'
#' First obtain a random sample of size \eqn{r_0} and estimate the model parameters.
#' Using the estimated parameters sampling probabilities are evaluated for local case control.
#'
#' Through the estimated sampling probabilities an optimal sample of size \eqn{r \ge r_0} is obtained.
#' Finally, the optimal sample is used and the model parameters are estimated.
#'
#' \strong{NOTE} : If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r \ge r_0} is not satisfied then an error message will be produced.
#'
#' If the big data \eqn{X,Y} has any missing values then an error message will be produced.
#'
#' The big data size \eqn{N} is compared with the sizes of \eqn{X,Y} and if they are not aligned an error
#' message will be produced.
#'
#' @return
#' The output of \code{LCCsampling} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters in a data.frame after sampling
#'
#' \code{Sample_LCC_Sampling} list of indexes for the initial and optimal samples obtained based on local case control sampling
#'
#' \code{Sampling_Probability} vector of calculated sampling probabilities for local case control sampling
#'
#' @references
#' \insertRef{fithian2015local}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"logistic"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r0<-300; rf<-rep(100*c(6,9,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' LCCsampling(r0 = r0, rf = rf, Y = as.matrix(Original_Data[,1]),
#'             X = as.matrix(Original_Data[,-1]),
#'             N = nrow(Original_Data))->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @export
LCCsampling<-function(r0,rf,Y,X,N){
  if(any(is.na(c(r0,rf,N))) | any(is.nan(c(r0,rf,N)))){
    stop("NA or Infinite or NAN values in the r0,rf or N")
  }

  if((length(r0) + length(N)) != 2){
    stop("r0 or N has a value greater than length one")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r0) > rf)){
    stop("2*r0 cannot be greater than rf at any point")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  n1 <- sum(Y)
  n0 <- N - n1
  PI.prop <- rep(1/(2*n0), N)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:N, size = r0, replace = TRUE, prob = PI.prop)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop,]
  pinv.prop <- 1/PI.prop[idx.prop]

  fit.prop <- .getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop_start <- fit.prop$par
  if(anyNA(beta.prop_start)){
    stop("There are NA or NaN values in the model parameters")
  }

  Xbeta_Final<-X%*% beta.prop_start
  P.prop  <- 1 - 1 / (1 + exp(Xbeta_Final))

  beta.LCC<-matrix(nrow = length(rf),ncol = ncol(X)+1 )
  Utility.LCC<-matrix(nrow = length(rf),ncol = 4 )
  Sample.LCC<-list()

  Sample.LCC[[1]]<-idx.prop

  beta.LCC[,1]<-Utility.LCC[,1]<-rf
  if(all(X[,1] == 1)){
    colnames(beta.LCC)<-c("rf",paste0("Beta_",0:(ncol(X)-1)))
  } else {
    colnames(beta.LCC)<-c("rf",paste0("Beta_",1:(ncol(X))))
  }

  ## local case control sampling
  PI.LCC <- abs(Y - P.prop)
  PI.LCC <- PI.LCC / sum(PI.LCC)

  message("Step 1 of the algorithm completed.\n")

  for(i in 1:length(rf)){
    ## local case control sampling
    idx.LCC <- sample(1:N, size = rf[i], replace = TRUE, prob = PI.LCC)

    x.LCC <- X[idx.LCC,]
    y.LCC <- Y[idx.LCC]
    pinv.LCC <- 1
    fit.LCC <- .getMLE(x=x.LCC, y=y.LCC, w=pinv.LCC)
    beta.prop <- fit.LCC$par+beta.prop_start

    Sample.LCC[[i+1]]<-idx.LCC;
    beta.LCC[i,-1] <- beta.prop

    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }
  }

  Full_SP<-cbind.data.frame(PI.LCC)
  colnames(Full_SP)<-c("Local case control")

  Sampling_Methods<-factor(c("Local case control"))
  Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(rf)),beta.LCC)

  names(Sample.LCC)<-c(r0,rf)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Sample_LCC_Sampling"=Sample.LCC,
            "Sampling_Probability"=Full_SP)

  class(ans)<-c("LocalCaseControl","logistic")
  return(ans)
}
