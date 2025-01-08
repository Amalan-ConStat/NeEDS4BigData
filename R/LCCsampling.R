#' Local case control sampling for logistic regression
#'
#' Using this function sample from big data under logistic regression to describe the data.
#' Sampling probabilities are obtained based on local case control method.
#'
#' @usage
#' LCCsampling(r1,r2,Y,X,N)
#'
#' @param r1      sample size for initial random sampling
#' @param r2      sample size for local case control sampling
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#'
#' @details
#' Two stage sampling algorithm for big data under logistic regression.
#'
#' First obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters sampling probabilities are evaluated for local case control.
#'
#' Through the estimated sampling probabilities an optimal sample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, the optimal sample is used and the model parameters are estimated.
#'
#' \strong{NOTE} : If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r_2 \ge r_1} is not satisfied then an error message will be produced.
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
#' \code{Utility_Estimates} estimated log scaled Information and variance for the estimated model parameters
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
#' r1<-300; r2<-rep(100*c(6,9,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' LCCsampling(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,1]),
#'             X = as.matrix(Original_Data[,-1]),
#'             N = nrow(Original_Data))->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @export
LCCsampling<-function(r1,r2,Y,X,N){
  if(any(is.na(c(r1,r2,N))) | any(is.nan(c(r1,r2,N)))){
    stop("NA or Infinite or NAN values in the r1,r2 or N")
  }

  if((length(r1) + length(N)) != 2){
    stop("r1 or N has a value greater than length one")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r1) > r2)){
    stop("2*r1 cannot be greater than r2 at any point")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  n1 <- sum(Y)
  n0 <- N - n1
  PI.prop <- rep(1/(2*n0), N)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:N, size = r1, replace = TRUE, prob = PI.prop)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop,]
  pinv.prop <- 1/PI.prop[idx.prop]

  fit.prop <- .getMLE(x=x.prop, y=y.prop, w=pinv.prop)
  beta.prop_start <- fit.prop$par
  if(anyNA(beta.prop_start)){
    stop("There are NA or NaN values in the model parameters")
  }

  P.prop  <- 1 - 1 / (1 + exp(X%*% beta.prop_start))

  beta.LCC<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  Sample.LCC<-list()

  Sample.LCC[[1]]<-idx.prop

  beta.LCC[,1]<-r2
  colnames(beta.LCC)<-c("r2",paste0("Beta",0:(ncol(X)-1)))

  ## local case control sampling
  PI.LCC <- abs(Y - P.prop)
  PI.LCC <- PI.LCC / sum(PI.LCC)

  message("Step 1 of the algorithm completed.\n")

  for(i in 1:length(r2)){
    ## local case control sampling
    idx.LCC <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.LCC)

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
  Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),beta.LCC)

  names(Sample.LCC)<-c(r1,r2)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Sample_LCC_Sampling"=Sample.LCC,
            "Sampling_Probability"=Full_SP)

  class(ans)<-c("LocalCaseControl","logistic")
  return(ans)
}
