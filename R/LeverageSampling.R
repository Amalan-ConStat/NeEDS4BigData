#' Basic and shrinkage leverage sampling for Generalised Linear Models
#'
#' Using this function sample from big data under linear, logistic and Poisson regression to describe the data.
#' Sampling probabilities are obtained based on the basic and shrinkage leverage method.
#'
#' @usage
#' LeverageSampling(r,Y,X,N,S_alpha,family)
#'
#' @param r       sample size
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#' @param S_alpha shrinkage factor in between 0 and 1
#' @param family  a character vector for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Leverage sampling algorithm for big data under Generalised Linear Models (linear, logistic and Poisson regression).
#'
#' First is to obtain a random sample of size \eqn{min(r)/2} and estimate the model parameters. Using the estimated parameters
#' leverage scores are evaluated for leverage sampling.
#'
#' Through the estimated leverage scores a sample of size \eqn{r} was obtained. Finally,
#' the sample of size \eqn{r} is used and the model parameters are estimated.
#'
#' \strong{NOTE} : If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r} is not satisfied then an error message will be produced.
#'
#' If the big data \eqn{X,Y} has any missing values then an error message will be produced.
#'
#' The big data size \eqn{N} is compared with the sizes of \eqn{X,Y} and if they are not aligned an error
#' message will be produced.
#'
#' If \eqn{0 < \alpha_{S} < 1} is not satisfied an error message will be produced.
#'
#' A character vector is provided for \code{family} and if it is not of the any three types an error message
#' will be produced.
#'
#' @return
#' The output of \code{LeverageSampling} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters in a data.frame after sampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon in a data.frame after sampling (valid only for linear regression)
#'
#' \code{Utility_Estimates} estimated D-(log scaled), A- and L- optimality values for the obtained subsamples
#'
#' \code{Sample_Basic_Leverage} list of indexes for the optimal samples obtained based on basic leverage
#'
#' \code{Sample_Shrinkage_Leverage} list of indexes for the optimal samples obtained based on shrinkage leverage
#'
#' \code{Sampling_Probability} matrix of calculated sampling probabilities for basic and shrinkage leverage
#'
#' @references
#' \insertRef{ma2014statistical}{NeEDS4BigData}
#' \insertRef{ma2015leveraging}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r<-rep(100*c(6,10),50); Original_Data<-Full_Data$Complete_Data;
#'
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  S_alpha = 0.95,
#'                  family = "linear")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r<-rep(100*c(6,10),50); Original_Data<-Full_Data$Complete_Data;
#'
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  S_alpha = 0.95,
#'                  family = "logistic")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal";
#' No_Of_Var<-2; Beta<-c(-1,0.5,0.5); N<-5000; Family<-"poisson"
#' Full_Data<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)
#'
#' r<-rep(100*c(6,10),50); Original_Data<-Full_Data$Complete_Data;
#'
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,1]),
#'                  X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                  S_alpha = 0.95,
#'                  family = "poisson")->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @export
LeverageSampling<-function(r,Y,X,N,S_alpha,family){
  if(any(is.na(c(r,S_alpha,N,family))) | any(is.nan(c(r,S_alpha,N,family)))){
    stop("NA or Infinite or NAN values in the r,S_alpha,N or family")
  }

  if((length(N) + length(family)) != 2){
    stop("N or family has a value greater than length one")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if(S_alpha >= 1 | S_alpha <= 0 | length(S_alpha) > 1){
    stop("S_alpha value for shrinkage leverage scores are not in the range of zero \nand one or the length is more than one")
  }

  if(!any(family == c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  S_alpha<-ifelse(is.null(S_alpha),0.9,S_alpha)
  if(family %in% c("linear")){
    svdf <- svd(X)
    U <- svdf$u
    PP <- apply(U, 1, crossprod)

    PI.blev <- PP / ncol(X)
    PI.slev <- S_alpha * PI.blev + (1-S_alpha) * 1 / N

    beta.blev<-beta.uwlev<-beta.slev<-matrix(nrow = length(r),ncol = ncol(X)+1 )
    Utility.blev<-Utility.uwlev<-Utility.slev<-matrix(nrow = length(r),ncol = 4 )
    Var_Epsilon<-matrix(nrow = length(r),ncol = 4)
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-Var_Epsilon[,1]<-
      Utility.blev[,1]<-Utility.uwlev[,1]<-Utility.slev[,1]<-r

    if(all(X[,1] == 1)){
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Var_Epsilon)<-c("r","Basic Leverage","Unweighted Leverage","Shrinkage Leverage")
    colnames(Utility.blev)<-colnames(Utility.uwlev)<-
      colnames(Utility.slev)<-c("r","D-optimality","A-optimality","L-optimality")

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample(1:N, size = r[i], replace = TRUE, prob = PI.blev)

      x.blev <- X[idx.blev,]
      y.blev <- Y[idx.blev]
      wgt <- 1 / PI.blev[idx.blev]

      Temp_Data <- data.frame(y = y.blev,x.blev)
      lm.blev <- stats::lm(y ~ . - 1, weights = wgt, data = Temp_Data)

      beta.prop<-stats::coefficients(lm.blev)
      Xbeta_Final<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Temp<-Var.prop*crossprod(x.blev)
      Temp_Inv<-solve(Temp)
      x.blev_t<-t(x.blev)
      Temp_Int<-Temp_Inv%*%x.blev_t
      Temp1<-x.blev%*%Temp_Int

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop
      Var_Epsilon[i,2]<-Var.prop
      Utility.blev[i,-1]<-c(log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # unweighted leverage sampling
      lm.uwlev <- stats::lm(y ~ . - 1, data = Temp_Data)

      beta.prop<-stats::coefficients(lm.uwlev)
      Xbeta_Final<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Temp<-Var.prop*crossprod(x.blev)
      Temp_Inv<-solve(Temp)
      Temp_Int<-Temp_Inv%*%x.blev_t
      Temp1<-x.blev%*%Temp_Int

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop
      Var_Epsilon[i,3]<-Var.prop
      Utility.uwlev[i,-1]<-c(log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # shrinkage leverage sampling
      idx.slev <- sample(N, size = r[i], replace = TRUE, prob = PI.slev)

      x.slev <- X[idx.slev,]
      y.slev <- Y[idx.slev]
      wgt <- 1 / PI.slev[idx.slev]

      Temp_Data <- data.frame(y = y.slev,x.slev)
      lm.slev <- stats::lm(y ~ . - 1, weights = wgt, data = Temp_Data)

      beta.prop<-stats::coefficients(lm.slev)
      Xbeta_Final<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Temp<-Var.prop*crossprod(x.slev)
      Temp_Inv<-solve(Temp)
      x.slev_t<-t(x.slev)
      Temp_Int<-Temp_Inv%*%x.slev_t
      Temp1<-x.slev%*%Temp_Int

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop
      Var_Epsilon[i,4]<-Var.prop
      Utility.slev[i,-1]<-c(log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.blev,PI.slev)
    colnames(Full_SP)<-c("Basic Leverage","Shrinkage Leverage")

    Sampling_Methods<-factor(c("Basic Leverage","Unweighted Leverage","Shrinkage Leverage"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                rbind(beta.blev,beta.uwlev,beta.slev))

    Utility_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                   rbind(Utility.blev,Utility.uwlev,Utility.slev))

    Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                       "Sample"=rep(r,times=length(Sampling_Methods)),
                                       "Var Epsilon"=c(Var_Epsilon[,"Basic Leverage"],
                                                       Var_Epsilon[,"Unweighted Leverage"],
                                                       Var_Epsilon[,"Shrinkage Leverage"]))

    names(Sample.blev)<-names(Sample.uwlev)<-names(Sample.slev)<-r

    message("Sampling completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
              "Sample_Basic_Leverage"=Sample.blev,
              "Sample_Shrinkage_Leverage"=Sample.slev,
              "Sampling_Probability"=Full_SP)
    class(ans)<-c("Leverage","linear")
    return(ans)
  }
  if(family %in% c("logistic")){
    r1<-round(min(r/2))
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

    Xbeta_Final <- X%*% beta.prop
    p.prop<-1 - 1 / (1 + exp(Xbeta_Final))
    w.prop<-as.vector(sqrt(p.prop*(1-p.prop)))
    X_bar<- w.prop*X
    XX_Inv<-solve(crossprod(X_bar))

    First_term <- X_bar %*% XX_Inv
    PP <- rowSums((First_term) * X_bar)

    PI.blev <- PP / sum(PP)
    PI.slev <- (S_alpha * PI.blev) + (1-S_alpha) / N

    beta.blev<-beta.uwlev<-beta.slev<-matrix(nrow = length(r),ncol = ncol(X)+1 )
    Utility.blev<-Utility.uwlev<-Utility.slev<-matrix(nrow = length(r),ncol = 4 )
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-
      Utility.blev[,1]<-Utility.uwlev[,1]<-Utility.slev[,1]<-r

    if(all(X[,1] == 1)){
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Utility.blev)<-colnames(Utility.uwlev)<-
      colnames(Utility.slev)<-c("r","D-optimality","A-optimality","L-optimality")

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample(N, size = r[i], replace = TRUE, PI.blev)

      x.blev <- X[idx.blev,]
      y.blev <- Y[idx.blev]
      wgt <- 1 / PI.blev[idx.blev]

      fit.blev <- .getMLE(x=x.blev, y=as.vector(y.blev), w=wgt)
      beta.prop <- fit.blev$par

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.blev %*% beta.blev[i,-1]
      pi<-c(1-1/(1 + exp(LP_data)))
      W<-pi*(1-pi)
      Mx<-crossprod(x.blev,(x.blev * W))
      Mx_Inv<-solve(Mx)
      x.blev_t<-t(x.blev)
      V_Mx_Temp <- Mx_Inv %*% x.blev_t
      V_Final<- x.blev%*% V_Mx_Temp

      Utility.blev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      # unweighted leverage sampling
      fit.uwlev <- .getMLE(x=x.blev, y=as.vector(y.blev), w=rep(N,r[i]))
      beta.prop<-fit.uwlev$par

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.blev %*% beta.uwlev[i,-1]
      pi<-c(1-1/(1 + exp(LP_data)))
      W<-pi*(1-pi)
      Mx<-crossprod(x.blev,(x.blev * W))
      Mx_Inv<-solve(Mx)
      V_Mx_Temp <- Mx_Inv %*% x.blev_t
      V_Final<- x.blev%*% V_Mx_Temp

      Utility.uwlev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      # shrinkage leverage sampling
      idx.slev <- sample(1:N, size = r[i], replace = TRUE, PI.slev)

      x.slev <- X[idx.slev,]
      y.slev <- Y[idx.slev]
      wgt <- 1 / PI.slev[idx.slev]

      fit.slev <- .getMLE(x=x.slev, y=as.vector(y.slev), w=wgt)
      beta.prop <- fit.slev$par

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.slev %*% beta.slev[i,-1]
      pi<-c(1-1/(1 + exp(LP_data)))
      W<-pi*(1-pi)
      Mx<-crossprod(x.slev,(x.slev * W))
      Mx_Inv<-solve(Mx)
      x.slev_t<-t(x.slev)
      V_Mx_Temp <- Mx_Inv %*% x.slev_t
      V_Final<- x.slev%*% V_Mx_Temp

      Utility.slev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))
    }

    Full_SP<-cbind.data.frame(PI.blev,PI.slev)
    colnames(Full_SP)<-c("Basic Leverage","Shrinkage Leverage")

    Sampling_Methods<-factor(c("Basic Leverage","Unweighted Leverage","Shrinkage Leverage"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                rbind(beta.blev,beta.uwlev,beta.slev))

    names(Sample.blev)<-names(Sample.uwlev)<-names(Sample.slev)<-r

    message("Sampling completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Sample_Basic_Leverage"=Sample.blev,
              "Sample_Shrinkage_Leverage"=Sample.slev,
              "Sampling_Probability"=Full_SP)
    class(ans)<-c("Leverage","logistic")
    return(ans)
  }
  if(family %in% c("poisson")){
    r1<-round(min(r/2))
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
    Xbeta_Final <- X%*% beta.prop
    p.prop<-exp(Xbeta_Final)
    w.prop<-as.vector(sqrt(p.prop))
    X_bar<- w.prop*X
    XX_Inv<-solve(crossprod(X_bar))

    # Precompute terms outside the loop
    First_term <- X_bar %*% XX_Inv
    PP <- rowSums((First_term) * X_bar)

    PI.blev <- PP / sum(PP)
    PI.slev <- (S_alpha * PI.blev) + ((1-S_alpha) / N)

    beta.blev<-beta.uwlev<-beta.slev<-matrix(nrow = length(r),ncol = ncol(X)+1 )
    Utility.blev<-Utility.uwlev<-Utility.slev<-matrix(nrow = length(r),ncol = 4 )
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-
      Utility.blev[,1]<-Utility.uwlev[,1]<-Utility.slev[,1]<-r
    if(all(X[,1] == 1)){
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.blev)<-colnames(beta.uwlev)<-
        colnames(beta.slev)<-c("r",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Utility.blev)<-colnames(Utility.uwlev)<-colnames(Utility.slev)<-
      c("r","D-optimality","A-optimality","L-optimality")

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample(1:N, size = r[i], replace = TRUE, prob = PI.blev)

      x.blev <- X[idx.blev,]
      y.blev <- Y[idx.blev]
      wgt <- 1 / PI.blev[idx.blev]

      fit.blev <-stats::glm(y.blev~x.blev-1, family = "quasipoisson",weights=wgt)
      beta.prop <- fit.blev$coefficients

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.blev %*% beta.blev[i,-1]
      pi<-c(exp(LP_data))
      W<-pi
      Mx<-crossprod(x.blev,(x.blev * W))
      Mx_Inv<-solve(Mx)
      x.blev_t<-t(x.blev)
      V_Mx_Temp <- Mx_Inv %*% x.blev_t
      V_Final<- x.blev%*% V_Mx_Temp

      Utility.blev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      # unweighted leverage sampling
      fit.uwlev <- stats::glm(y.blev~x.blev-1, family = "quasipoisson")
      beta.prop<-fit.uwlev$coefficients

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.blev %*% beta.uwlev[i,-1]
      pi<-c(exp(LP_data))
      W<-pi
      Mx<-crossprod(x.blev,(x.blev * W))
      Mx_Inv<-solve(Mx)
      V_Mx_Temp <- Mx_Inv %*% x.blev_t
      V_Final<- x.blev%*% V_Mx_Temp

      Utility.uwlev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      # shrinkage leverage sampling
      idx.slev <- sample(1:N, size = r[i], replace = TRUE, prob = PI.slev)

      x.slev <- X[idx.slev,]
      y.slev <- Y[idx.slev]
      wgt <- 1 / PI.slev[idx.slev]

      fit.slev <- stats::glm(y.slev~x.slev-1, family = "quasipoisson",weights = wgt)
      beta.prop <- fit.slev$coefficients

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.slev %*% beta.slev[i,-1]
      pi<-c(exp(LP_data))
      W<-pi
      Mx<-crossprod(x.slev,(x.slev * W))
      Mx_Inv<-solve(Mx)
      x.slev_t<-t(x.slev)
      V_Mx_Temp <- Mx_Inv %*% x.slev_t
      V_Final<- x.slev%*% V_Mx_Temp

      Utility.slev[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))
    }

    Full_SP<-cbind.data.frame(PI.blev,PI.slev)
    colnames(Full_SP)<-c("Basic Leverage","Shrinkage Leverage")

    Sampling_Methods<-factor(c("Basic Leverage","Unweighted Leverage","Shrinkage Leverage"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                rbind(beta.blev,beta.uwlev,beta.slev))

    Utility_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                   rbind(Utility.blev,Utility.uwlev,Utility.slev))

    names(Sample.blev)<-names(Sample.uwlev)<-names(Sample.slev)<-r

    message("Sampling completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Sample_Basic_Leverage"=Sample.blev,
              "Sample_Shrinkage_Leverage"=Sample.slev,
              "Sampling_Probability"=Full_SP)
    class(ans)<-c("Leverage","poisson")
    return(ans)
  }
}
