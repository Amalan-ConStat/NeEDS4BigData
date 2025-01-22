#' A-optimality criteria based subsampling under measurement constraints for Generalised Linear Models
#'
#' Using this function sample from big data under linear, logistic and Poisson regression
#' to describe the data when response \eqn{y} is partially unavailable. Subsampling probabilities are
#' obtained based on the A-optimality criteria.
#'
#' @usage
#' AoptimalMCGLMSub(r1,r2,Y,X,N,family)
#'
#' @param r1      sample size for initial random sampling
#' @param r2      sample size for optimal sampling
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#' @param family  a character value for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Two stage subsampling algorithm for big data under Generalised Linear Models
#' (linear, logistic and Poisson regression) when the response is not available for
#' subsampling probability evaluation.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-optimality criteria.
#'
#' Through the estimated subsampling probabilities an optimal sample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, only the optimal sample is used and the model parameters are estimated.
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
#' A character value is provided for \code{family} and if it is not of the any three types an error message
#' will be produced.
#'
#' @return
#' The output of \code{AoptimalMCGLMSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters in a data.frame after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon in a data.frame after subsampling (valid only for linear regression)
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A-optimality criteria
#'
#' @references
#' \insertRef{zhang2021optimal}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"linear"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(100*c(6,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' AoptimalMCGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  family = "linear")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"logistic"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(100*c(6,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' AoptimalMCGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  family = "logistic")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal";
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-10000; Family<-"poisson"
#' Full_Data<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(100*c(6,12),50); Original_Data<-Full_Data$Complete_Data;
#'
#' AoptimalMCGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  family = "poisson")->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @importFrom matrixStats rowSums2
#' @export
AoptimalMCGLMSub <- function(r1,r2,Y,X,N,family){
  if(any(is.na(c(r1,r2,N,family))) | any(is.nan(c(r1,r2,N,family)))){
    stop("NA or Infinite or NAN values in the r1,r2,N or family")
  }

  if((length(r1) + length(N) + length(family)) != 3){
    stop("r1, N or family has a value greater than length one")
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

  if(!any(family == c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if(family %in% c("linear")){
    PI.prop <- rep(1/N, N)
    idx.prop <- sample(1:N, size = r1, replace = TRUE)

    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop,]

    pinv.prop <- N
    pinv.prop <- 1/PI.prop[idx.prop]

    beta.prop<-solve(a=crossprod(x.prop),b= crossprod(x.prop,y.prop))
    Xbeta_Final<-X%*%beta.prop
    Var.prop<-sum((Y-Xbeta_Final)^2)/N
    Epsilon.prop<-Y-Xbeta_Final

    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }

    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Var_Epsilon<-matrix(nrow = length(r2),ncol = 2)
    Sample.mMSE<-list()

    Sample.mMSE[[1]]<-idx.prop

    beta.mMSE[,1]<-Var_Epsilon[,1]<-r2

    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Var_Epsilon)<-c("r2","A-Optimality")

    ## mMSE
    PI.mMSE<-sqrt(matrixStats::rowSums2((X %*% solve(crossprod(X)))^2))
    PI.mMSE<-PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[idx.mMSE,]
      y.mMSE <- Y[idx.mMSE]
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE])

      pi4_r<-sqrt(r2[i]*pinv.mMSE^(-1))
      X_r4<-x.mMSE/pi4_r
      Y_r4<-y.mMSE/pi4_r
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final1<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Sample.mMSE[[i+1]]<-idx.mMSE
      beta.mMSE[i,-1] <- beta.prop
      Var_Epsilon[i,2]<-Var.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.mMSE)
    colnames(Full_SP)<-c("A-Optimality")

    Sampling_Methods<-factor(c("A-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),beta.mMSE)

    Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),
                                       "Sample"=rep(r2,times=length(Sampling_Methods)),
                                       "Var Epsilon"=c(Var_Epsilon[,"A-Optimality"]))

    names(Sample.mMSE)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Subsampling_Probability"=Full_SP)
    class(ans)<-c("A_OptimalSamplingMC","linear")
    return(ans)
  }
  if(family %in% "logistic"){
    n1 <- sum(Y)
    n0 <- N - n1
    PI.prop <- rep(1/(2*n0), N)
    PI.prop[Y==1] <- 1/(2*n1)
    idx.prop <- sample(1:N, size = r1, replace = TRUE, prob = PI.prop)

    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop,]
    pinv.prop <- 1/PI.prop[idx.prop]

    fit.prop <- .getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }

    P.prop  <- 1 - 1 / (1 + exp(X%*% beta.prop))

    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Sample.mMSE<-list()

    Sample.mMSE[[1]]<-idx.prop

    beta.mMSE[,1]<-r2
    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }


    ## mMSE
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- r1*solve(crossprod(x.prop,(x.prop * w.prop * pinv.prop)))

    PI.mMSE<-sqrt(P.prop*(1-P.prop))*sqrt(matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[idx.mMSE,]
      y.mMSE <- Y[idx.mMSE]
      fit.mMSE <- .getMLE(x=x.mMSE, y=y.mMSE,w=c(1 / PI.mMSE[idx.mMSE]))

      Sample.mMSE[[i+1]]<-idx.mMSE
      beta.mMSE[i,-1] <- fit.mMSE$par

      if(anyNA(fit.mMSE$par)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.mMSE)
    colnames(Full_SP)<-c("A-Optimality")

    Sampling_Methods<-factor(c("A-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),
                                beta.mMSE)

    names(Sample.mMSE)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Subsampling_Probability"=Full_SP)

    class(ans)<-c("A_OptimalSamplingMC","logistic")
    return(ans)
  }
  if(family %in% "poisson"){
    PI.prop <- rep(1/N, N)
    idx.prop <- sample(1:N, size = r1, replace = TRUE)

    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop,]

    pinv.prop <- N
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- stats::glm(y.prop~x.prop-1,family = "quasipoisson")

    beta.prop <- fit.prop$coefficients
    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }

    P.prop  <- exp(X %*% beta.prop)

    beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Sample.mMSE<-list()

    Sample.mMSE[[1]]<-idx.prop

    beta.mMSE[,1]<-r2
    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }

    ## mMSE
    w.prop <- P.prop[idx.prop]
    W.prop <- r1*solve(crossprod(x.prop,x.prop * w.prop * pinv.prop))

    PI.mMSE<-sqrt(P.prop)*sqrt(matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[idx.mMSE,]
      y.mMSE <- Y[idx.mMSE]
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE])

      fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1, family = "quasipoisson",weights=pinv.mMSE)
      Sample.mMSE[[i+1]]<-idx.mMSE

      beta.mMSE[i,-1] <-fit.mMSE$coefficients

      if(anyNA(fit.mMSE$coefficients)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.mMSE)
    colnames(Full_SP)<-c("A-Optimality")

    Sampling_Methods<-factor(c("A-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),
                                beta.mMSE)

    names(Sample.mMSE)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Subsampling_Probability"=Full_SP)
    class(ans)<-c("A_OptimalSamplingMC","poisson")
    return(ans)
  }
}
