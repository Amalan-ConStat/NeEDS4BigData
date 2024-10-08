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
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
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
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
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
#' LeverageSampling(r = r, Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
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

  if(any(is.na(cbind(Y,X))) | any(is.nan(cbind(Y,X)))){
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
    Var_Epsilon<-matrix(nrow = length(r),ncol = 4)
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-Var_Epsilon[,1]<-r
    colnames(beta.blev)<-colnames(beta.uwlev)<-colnames(beta.slev)<-c("r",paste0("Beta",0:(ncol(X)-1)))
    colnames(Var_Epsilon)<-c("r","Basic Leverage","Unweighted Leverage","Shrinkage Leverage")

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample.int(N, size = r[i], replace = TRUE, PI.blev)

      x.blev <- as.matrix(X[idx.blev,])
      y.blev <- as.matrix(Y[idx.blev])
      wgt <- 1 / PI.blev[idx.blev]

      Temp_Data <- data.frame(y = y.blev,x.blev, wgt = wgt)
      lm.blev <- stats::lm(y ~ . - 1, weights = Temp_Data$wgt, data = Temp_Data[,-ncol(Temp_Data)])

      beta.prop<-stats::coefficients(lm.blev)
      Xbeta_Final<-as.vector(X%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop
      Var_Epsilon[i,2]<-Var.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # unweighted leverage sampling
      lm.uwlev <- stats::lm(y ~ . - 1, data = Temp_Data[,-ncol(Temp_Data)])

      beta.prop<-stats::coefficients(lm.uwlev)
      Xbeta_Final<-as.vector(X%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop
      Var_Epsilon[i,3]<-Var.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # shrinkage leverage sampling
      idx.slev <- sample.int(N, size = r[i], replace = TRUE, PI.slev)

      x.slev <- as.matrix(X[idx.slev,])
      y.slev <- as.matrix(Y[idx.slev])
      wgt <- 1 / PI.slev[idx.slev]

      Temp_Data <- data.frame(y = y.slev,x.slev, wgt = wgt)
      lm.slev <- stats::lm(y ~ . - 1, weights = wgt, data = Temp_Data[,-ncol(Temp_Data)])

      beta.prop<-stats::coefficients(lm.slev)
      Xbeta_Final<-as.vector(X%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop
      Var_Epsilon[i,4]<-Var.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.blev,PI.slev)
    colnames(Full_SP)<-c("Basic Leverage","Shrinkage Leverage")

    Sampling_Methods<-factor(c("Basic Leverage","Unweighted Leverage","Shrinkage Leverage"))

    Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                rbind(beta.blev,beta.uwlev,beta.slev))

    Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                       "Sample"=rep(r,times=length(Sampling_Methods)),
                                       "Var Epsilon"=c(Var_Epsilon[,"Basic Leverage"],
                                                       Var_Epsilon[,"Unweighted Leverage"],
                                                       Var_Epsilon[,"Shrinkage Leverage"]))

    names(Sample.blev)<-names(Sample.uwlev)<-names(Sample.slev)<-r

    message("Sampling completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
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
    idx.prop <- sample(1:N, r1, T, PI.prop)

    x.prop <- X[idx.prop,]
    y.prop <- Y[idx.prop,]
    pinv.prop <- 1/PI.prop[idx.prop]

    fit.prop <- .getMLE(x=x.prop, y=y.prop, w=pinv.prop)
    beta.prop <- fit.prop$par
    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }
    p.prop<-1 - 1 / (1 + exp(X%*% beta.prop))
    w.prop<-as.vector(sqrt(p.prop*(1-p.prop)))
    X_bar<- w.prop*X
    XX_Inv<-solve(t(X_bar)%*%X_bar)
    X_Temp<-X
    PP <-apply(X_Temp,1,function(X_Temp){
      P.prop  <- 1 - 1 / (1 + exp(X_Temp%*% beta.prop))
      W.prop <- sqrt(P.prop*(1-P.prop))
      x_bar <- W.prop%*%X_Temp
      x_bar%*%XX_Inv%*%t(x_bar)
      })

    PI.blev <- PP / sum(PP)
    PI.slev <- S_alpha * PI.blev + (1-S_alpha) * 1 / N

    beta.blev<-beta.uwlev<-beta.slev<-matrix(nrow = length(r),ncol = ncol(X)+1 )
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-r
    colnames(beta.blev)<-colnames(beta.uwlev)<-colnames(beta.slev)<-c("r",paste0("Beta",0:(ncol(X)-1)))

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample.int(N, size = r[i], replace = TRUE, PI.blev)

      x.blev <- as.matrix(X[idx.blev,])
      y.blev <- as.matrix(Y[idx.blev])
      wgt <- 1 / PI.blev[idx.blev]

      fit.blev <- .getMLE(x=x.blev, y=as.vector(y.blev), w=wgt)
      beta.prop <- fit.blev$par

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # unweighted leverage sampling
      fit.uwlev <- .getMLE(x=x.blev, y=as.vector(y.blev), w=N)
      beta.prop<-fit.uwlev$par

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # shrinkage leverage sampling
      idx.slev <- sample.int(N, size = r[i], replace = TRUE, PI.slev)

      x.slev <- as.matrix(X[idx.slev,])
      y.slev <- as.matrix(Y[idx.slev])
      wgt <- 1 / PI.slev[idx.slev]

      fit.slev <- .getMLE(x=x.slev, y=as.vector(y.slev), w=wgt)
      beta.prop <- fit.slev$par

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
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
    class(ans)<-c("Leverage","linear")
    return(ans)
  }
  if(family %in% c("poisson")){
    r1<-round(min(r/2))
    PI.prop <- rep(1/N, N)
    idx.prop <- sample(1:N, r1, T)

    x.prop<-X[idx.prop,]
    y.prop <- Y[idx.prop,]

    pinv.prop <- N
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- stats::glm(y.prop~x.prop-1,family = "poisson")

    beta.prop <- fit.prop$coefficients
    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }
    p.prop<-as.vector(sqrt(exp(X%*% beta.prop)))
    X_bar<- p.prop*X
    XX_Inv<-solve(t(X_bar)%*%X_bar)
    X_Temp<-X;
    PP <-apply(X_Temp,1,function(X_Temp){
      P.prop  <- sqrt(exp(X_Temp%*%beta.prop))
      x_bar <- P.prop%*%X_Temp
      x_bar%*%XX_Inv%*%t(x_bar)
    })

    PI.blev <- PP / sum(PP)
    PI.slev <- S_alpha * PI.blev + (1-S_alpha) * 1 / N

    beta.blev<-beta.uwlev<-beta.slev<-matrix(nrow = length(r),ncol = ncol(X)+1 )
    Sample.blev<-Sample.uwlev<-Sample.slev<-list()

    beta.blev[,1]<-beta.uwlev[,1]<-beta.slev[,1]<-r
    colnames(beta.blev)<-colnames(beta.uwlev)<-colnames(beta.slev)<-c("r",paste0("Beta",0:(ncol(X)-1)))

    message("Basic and shrinkage leverage probabilities calculated.\n")

    for (i in 1:length(r)) {
      # basic leverage sampling
      idx.blev <- sample.int(N, size = r[i], replace = TRUE, PI.blev)

      x.blev <- as.matrix(X[idx.blev,])
      y.blev <- as.matrix(Y[idx.blev])
      wgt <- 1 / PI.blev[idx.blev]

      fit.blev <-stats::glm(y.blev~x.blev-1, family = "poisson",weights=wgt)
      beta.prop <- fit.blev$coefficients

      Sample.blev[[i]]<-idx.blev
      beta.blev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # unweighted leverage sampling
      fit.uwlev <- stats::glm(y.blev~x.blev-1, family = "poisson")
      beta.prop<-fit.uwlev$coefficients

      Sample.uwlev[[i]]<-idx.blev
      beta.uwlev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      # shrinkage leverage sampling
      idx.slev <- sample.int(N, size = r[i], replace = TRUE, PI.slev)

      x.slev <- as.matrix(X[idx.slev,])
      y.slev <- as.matrix(Y[idx.slev])
      wgt <- 1 / PI.slev[idx.slev]

      fit.slev <- stats::glm(y.slev~x.slev-1, family = "poisson",weights = wgt)
      beta.prop <- fit.slev$coefficients

      Sample.slev[[i]]<-idx.slev
      beta.slev[i,-1] <- beta.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
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
    class(ans)<-c("Leverage","poisson")
    return(ans)
  }
}
