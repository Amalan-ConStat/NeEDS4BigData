#' A-optimality criterion based subsampling under Gaussian Linear Models
#'
#' Using this function subsample from big data under Gaussian linear regression models
#' to describe the data. Subsampling probabilities are obtained based on the A-optimality criterions.
#'
#' @usage
#' AoptimalGauLMSub(r1,r2,Y,X,N)
#'
#' @param r1      subsample size for initial random sampling
#' @param r2      subsample size for optimal sampling
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#'
#' @details
#' Two stage subsampling algorithm for big data under Gaussian Linear Model.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-optimality criterion.
#'
#' Through the estimated subsampling probabilities an optimal subsample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated.
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
#' The output of \code{AoptimalGauLMSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters in a data.frame after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon in a data.frame after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal subsamples obtained based on A-Optimality criterion
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A-optimality criterion
#'
#' @references
#' \insertRef{lee2021fast}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"linear"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(100*c(6,9,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' AoptimalGauLMSub(r1 = r1, r2 = r2,
#'                  Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                  X = as.matrix(Original_Data[,-1]),
#'                  N = nrow(Original_Data))->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @export
AoptimalGauLMSub <- function(r1,r2,Y,X,N){
  if(any(is.na(c(r1,r2,N))) | any(is.nan(c(r1,r2,N)))){
    stop("NA or Infinite or NAN values in the r1,r2 or N")
  }

  if(any(is.na(cbind(Y,X))) | any(is.nan(cbind(Y,X)))){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r1) > r2)){
    stop("2*r1 cannot be greater than r2 at any point")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  PI.prop <- rep(1/N, N)
  idx.prop <- sample(1:N, r1, T)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop,]

  pinv.prop <- N
  pinv.prop <- 1/PI.prop[idx.prop]

  beta.prop<-solve(a=t(x.prop)%*%x.prop,b=t(x.prop)%*%y.prop)
  Xbeta_Final<-as.vector(X%*%beta.prop)
  Var.prop<-sum((Y-Xbeta_Final)^2)/N
  Epsilon.prop<-Y-Xbeta_Final

  if(anyNA(beta.prop)){
    stop("There are NA or NaN values in the model parameters")
  }

  Second <- (Epsilon.prop^2 - Var.prop)^2/(4 * N^2 * Var.prop)
  ML_Inv <- solve(t(X)%*%X)

  beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  Var_Epsilon<-matrix(nrow = length(r2),ncol = 2)
  Sample.mMSE<-list()

  Sample.mMSE[[1]]<-idx.prop

  beta.mMSE[,1]<-Var_Epsilon[,1]<-r2
  colnames(beta.mMSE)<-c("r2",paste0("Beta",0:(ncol(X)-1)))
  colnames(Var_Epsilon)<-c("r2","A-Optimality")

  ## mMSE
  PI.mMSE <- sqrt(Epsilon.prop^2 * rowSums((X %*% ML_Inv)^2) + Second)
  PI.mMSE <- PI.mMSE/sum(PI.mMSE)

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    ## mMSE
    idx.mMSE <- sample(1:N, r2[i]-r1, T, PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)

    pi4_r<-sqrt(r2[i]*pinv.mMSE^(-1))
    X_r4<-x.mMSE/pi4_r
    Y_r4<-y.mMSE/pi4_r
    beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4)
    Xbeta_Final<-as.vector(X%*%beta.prop)
    Var.prop<-sum((Y-Xbeta_Final)^2)/N

    Sample.mMSE[[i+1]]<-idx.mMSE;
    beta.mMSE[i,-1] <- beta.prop
    Var_Epsilon[i,2]<-Var.prop

    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }
  }

  Full_SP<-cbind.data.frame(PI.mMSE)
  colnames(Full_SP)<-c("A-Optimality")

  Subsampling_Methods<-factor(c("A-Optimality"))

  Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),beta.mMSE)

  Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     "Subsample"=r2,"Var Epsilon"=Var_Epsilon[,"A-Optimality"])

  names(Sample.mMSE)<-c(r1,r2)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("AoptimalSubsampling","linear")
  return(ans)
}
