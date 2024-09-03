#' A- and L-optimality criteria based subsampling under Generalised Linear Models
#'
#' Using this function sample from big data under linear, logistic and Poisson regression
#' to describe the data. Sampling probabilities are obtained based on the A- and L-
#' optimality criteria.
#'
#' @usage
#' ALoptimalGLMSub(r1,r2,Y,X,N,family)
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
#' (linear, logistic and Poisson regression).
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters sampling probabilities are evaluated for A- and L-optimality criteria.
#'
#' Through the estimated sampling probabilities an optimal sample of size \eqn{r_2 \ge r_1} is obtained.
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
#' A character value is provided for \code{family} and if it is not of the any three types an error message
#' will be produced.
#'
#' @return
#' The output of \code{ALoptimalGLMSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters in a data.frame after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon in a data.frame after subsampling (valid only for linear regression)
#'
#' \code{Utility_Estimates} estimated log scaled Information and variance for the estimated model parameters
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
#'
#' \code{Sampling_Probability} matrix of calculated sampling probabilities for A- and L- optimality criteria
#'
#' @references
#' \insertRef{wang2018optimal}{NeEDS4BigData}
#' \insertRef{ai2021optimal}{NeEDS4BigData}
#' \insertRef{yao2021review}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(600,50); Original_Data<-Full_Data$Complete_Data;
#'
#' ALoptimalGLMSub(r1 = r1, r2 = r2,
#'                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                 family = "linear")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
#' Full_Data<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(600,50); Original_Data<-Full_Data$Complete_Data;
#'
#' ALoptimalGLMSub(r1 = r1, r2 = r2,
#'                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                 family = "logistic")->Results
#'
#' plot_Beta(Results)
#'
#' Dist<-"Normal";
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"poisson"
#' Full_Data<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)
#'
#' r1<-300; r2<-rep(600,50); Original_Data<-Full_Data$Complete_Data;
#'
#' ALoptimalGLMSub(r1 = r1, r2 = r2,
#'                 Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                 X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                 family = "poisson")->Results
#'
#' plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @importFrom psych tr
#' @importFrom matrixStats rowSums2
#' @export
ALoptimalGLMSub <- function(r1,r2,Y,X,N,family){
  if(any(is.na(c(r1,r2,N,family))) | any(is.nan(c(r1,r2,N,family)))){
    stop("NA or Infinite or NAN values in the r1,r2,N or family")
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

  if(!any(family == c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if(family %in% c("linear")){
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

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Var_Epsilon<-matrix(nrow = length(r2),ncol = 3)
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Var_Epsilon[,1]<-r2
    colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta",0:(ncol(X)-1)))
    colnames(Var_Epsilon)<-c("r2","A-Optimality","L-Optimality")

    ## mVc
    PI.mVc<-sqrt(Epsilon.prop^2 * matrixStats::rowSums2(X^2))
    PI.mVc<-PI.mVc/sum(PI.mVc)

    ## mMSE
    PI.mMSE<-sqrt(Epsilon.prop^2 * matrixStats::rowSums2((X %*% solve(t(X)%*%X))^2))
    PI.mMSE<-PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, r2[i]-r1, T, PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc<-c(1 / PI.mVc[idx.mVc], pinv.prop)

      pi4_r<-sqrt(r2[i]*pinv.mVc^(-1))
      X_r4<-x.mVc/pi4_r
      Y_r4<-y.mVc/pi4_r
      beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4)
      Xbeta_Final<-as.vector(X%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Sample.mVc[[i+1]]<-idx.mVc;
      beta.mVc[i,-1] <- beta.prop
      Var_Epsilon[i,3]<-Var.prop

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

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

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),rbind(beta.mMSE,beta.mVc))

    Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                       "Sample"=rep(r2,times=length(Subsampling_Methods)),
                                       "Var Epsilon"=c(Var_Epsilon[,"A-Optimality"],
                                                       Var_Epsilon[,"L-Optimality"]))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Sampling_Probability"=Full_SP)
    class(ans)<-c("A_L_OptimalSubsampling","linear")
    return(ans)
  }
  if(family %in% "logistic"){
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

    P.prop  <- 1 - 1 / (1 + exp(X%*% beta.prop))

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Utility_mVc<-Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Utility_mVc[,1]<-Utility_mMSE[,1]<-r2

    colnames(beta.mVc)<-colnames(beta.mMSE)<-c("r2",paste0("Beta",0:(ncol(X)-1)))
    colnames(Utility_mVc)<-colnames(Utility_mMSE)<-c("r2","Variance","Information")

    ## mVc
    PI.mVc<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2(X^2))
    PI.mVc <- PI.mVc/sum(PI.mVc)

    ## mMSE
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))

    PI.mMSE<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, r2[i]-r1, T, PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop), ]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      fit.mVc <-.getMLE(x=x.mVc, y=y.mVc,w=c(1 / PI.mVc[idx.mVc], pinv.prop))

      Sample.mVc[[i+1]]<-idx.mVc;
      beta.mVc[i,-1] <- fit.mVc$par

      if(anyNA(fit.mVc$par)){
        stop("There are NA or NaN values in the model parameters")
      }

      pi<-1-1/(1 + exp(x.mVc %*% beta.mVc[i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI.mVc[idx.mVc], pinv.prop))
      Mx<-solve((t(x.mVc) %*% (x.mVc * W)))
      Middle<-((as.vector(y.mVc)-as.vector(pi))*as.vector(c(1 / PI.mVc[idx.mVc], pinv.prop)))^2
      V_Temp<-(t(x.mVc) %*% (x.mVc * Middle) )
      V_Final<- Mx %*% V_Temp %*% Mx

      Utility_mVc[i,-1]<-c(psych::tr(V_Final),det(solve(V_Final)))

      ## mMSE
      idx.mMSE <- sample(1:N, r2[i]-r1, T, PI.mMSE)

      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      fit.mMSE <- .getMLE(x=x.mMSE, y=y.mMSE,w=c(1 / PI.mMSE[idx.mMSE], pinv.prop))

      Sample.mMSE[[i+1]]<-idx.mMSE;
      beta.mMSE[i,-1] <- fit.mMSE$par

      if(anyNA(fit.mMSE$par)){
        stop("There are NA or NaN values in the model parameters")
      }

      pi<-1-1/(1 + exp(x.mMSE %*% beta.mMSE[i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI.mMSE[idx.mMSE], pinv.prop))
      Mx<-solve((t(x.mMSE) %*% (x.mMSE * W)))
      Middle<-((as.vector(y.mMSE)-as.vector(pi))*as.vector(c(1 / PI.mMSE[idx.mMSE], pinv.prop)))^2
      V_Temp<-(t(x.mMSE) %*% (x.mMSE * Middle))
      V_Final<-Mx %*% V_Temp %*% Mx

      Utility_mMSE[i,-1]<-c(psych::tr(V_Final),det(solve(V_Final)))
    }

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                rbind(beta.mMSE,beta.mVc))

    Utility_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                   rbind(Utility_mMSE,Utility_mVc))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Sampling_Probability"=Full_SP)

    class(ans)<-c("A_L_OptimalSubsampling","logistic")
    return(ans)
  }
  if(family %in% "poisson"){
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

    P.prop  <- exp(X %*% beta.prop)

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Utility_mVc<-Utility_mMSE<-matrix(nrow = length(r2),ncol = 3 )
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Utility_mVc[,1]<-Utility_mMSE[,1]<-r2

    colnames(beta.mVc)<-colnames(beta.mMSE)<-c("r2",paste0("Beta",0:(ncol(X)-1)))
    colnames(Utility_mVc)<-colnames(Utility_mMSE)<-c("r2","Variance","Information")

    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * matrixStats::rowSums2(X^2))
    PI.mVc <- PI.mVc/sum(PI.mVc)

    ## mMSE
    w.prop <- P.prop[idx.prop]
    W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))

    PI.mMSE<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, r2[i]-r1, T, PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc<-c(1 / PI.mVc[idx.mVc], pinv.prop)

      fit.mVc <-stats::glm(y.mVc~x.mVc-1, family = "poisson",weights=pinv.mVc)
      Sample.mVc[[i+1]]<-idx.mVc;
      beta.mVc[i,-1] <- fit.mVc$coefficients

      if(anyNA(fit.mVc$coefficients)){
        stop("There are NA or NaN values in the model parameters")
      }

      pi<-c(exp(x.mVc %*% beta.mVc[i,-1]))
      pinv.mVc<-c(1 / PI.mVc[idx.mVc], pinv.prop)
      Mx<-solve(t(x.mVc) %*% (x.mVc * pi * pinv.mVc))
      V_Temp<-t(x.mVc) %*% (x.mVc*((as.vector(y.mVc)-pi)*pinv.mVc)^2)
      V_Final<- Mx %*% V_Temp %*% Mx

      Utility_mVc[i,-1]<-c(psych::tr(V_Final),det(solve(V_Final)))

      ## mMSE
      idx.mMSE <- sample(1:N, r2[i]-r1, T, PI.mMSE)

      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)

      fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1, family = "poisson",weights=pinv.mMSE)
      Sample.mMSE[[i+1]]<-idx.mMSE;
      beta.mMSE[i,-1] <-fit.mMSE$coefficients

      if(anyNA(fit.mMSE$coefficients)){
        stop("There are NA or NaN values in the model parameters")
      }

      pi<-c(exp(x.mMSE %*% beta.mMSE[i,-1]))
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      Mx<-solve(t(x.mMSE) %*% (x.mMSE*pi*pinv.mMSE))
      V_Temp<-t(x.mMSE) %*% (x.mMSE*((as.vector(y.mMSE)-pi)*pinv.mMSE)^2)
      V_Final<-Mx %*% V_Temp %*% Mx

      Utility_mMSE[i,-1]<-c(psych::tr(V_Final),det(solve(V_Final)))
    }

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                rbind(beta.mMSE,beta.mVc))

    Utility_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                   rbind(Utility_mMSE,Utility_mVc))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Sampling_Probability"=Full_SP)
    class(ans)<-c("A_L_OptimalSubsampling","poisson")
    return(ans)
  }
}

#' Generate data for Generalised Linear Models
#'
#' Function to simulate big data under linear, logistic and Poisson regression for sampling.
#' Covariate data X is through Normal or Uniform distribution for linear regression.
#' Covariate data X is through Exponential or Normal or Uniform distribution for logistic regression.
#' Covariate data X is through Normal or Uniform distribution for Poisson regression.
#'
#' @usage
#' GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,family)
#'
#' @param Dist        a character value for the distribution "Normal" or "Uniform or "Exponential"
#' @param Dist_Par    a list of parameters for the distribution that would generate data for covariate X
#' @param No_Of_Var   number of variables
#' @param Beta        a vector for the model parameters, including the intercept
#' @param N           the big data size
#' @param family      a character vector for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Big data for the Generalised Linear Models are generated by the "linear", "logistic" and "poisson"
#' regression types.
#'
#' We have limited the covariate data generation for
#' linear regression through normal and uniform distribution,
#' logistic regression through exponential, normal and uniform and
#' Poisson regression through normal and uniform distribution.
#'
#' @return
#' The output of \code{GenGLMData} gives a list of
#'
#' \code{Basic} a list of outputs based on the inputs and Beta Estimates for all models
#'
#' \code{Complete_Data} a matrix for Y and X
#'
#' @references
#' \insertRef{lee1996hierarchical}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"linear"
#' Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"logistic"
#' Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' Dist<-"Normal";
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000; Family<-"poisson"
#' Results<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)
#'
#' @import stats
#' @export
GenGLMdata<-function(Dist,Dist_Par,No_Of_Var,Beta,N,family){
  if(any(is.na(c(Dist,Beta,No_Of_Var,N,family))) | any(is.nan(c(Dist,No_Of_Var,Beta,N,family)))){
    stop("NA or Infinite or NAN values in the Dist,Beta,No_Of_Var,N or family")
  }

  if(any(is.na(unlist(Dist_Par))) | any(is.nan(unlist(Dist_Par)))){
    stop("NA or Infinite or NAN values in the Dist_Par")
  }

  if(!any(family == c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if(family == "linear"){
    if(!(Dist == "Normal" | Dist == "Uniform")){
      stop("For linear regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
    }
  }

  if(family == "logistic"){
    if(!(Dist == "Exponential" | Dist == "Normal" | Dist == "Uniform")){
      stop("For logistic regression select the distribution 'Exponential', \n 'Normal' or 'Uniform' to generate the covarate data")
    }
  }

  if(family == "poisson"){
    if(!(Dist == "Normal" | Dist == "Uniform")){
      stop("For poisson regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
    }
  }

  if(family %in% "linear"){
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }
    Complete_Data<-cbind(1,X);
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))
    Residual<-stats::rnorm(n=N,mean=0,sd=sqrt(Dist_Par$Error_Variance))
    Y <- Complete_Data%*%Beta + Residual

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    stats::lm(Y~.-1,data=as.data.frame(Complete_Data))->Results
    Beta_Estimate<-stats::coefficients(Results)
    Var_Epsilon_Estimate<-mean((Y-stats::fitted.values(Results))^2)

    Outputs<-list("Basic"=list("N"=N,"Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimate,
                               "Variance_Epsilon_Estimates"=Var_Epsilon_Estimate,
                               "Distribution"=Dist,
                               "Distribution_Parameter"=Dist_Par,
                               "No_Of_Variables"=No_Of_Var),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("A_L_OptimalSubsampling","linear")
    return(Outputs)
  }
  if(family %in% "logistic"){
    if(Dist %in% "Exponential"){
      X<-replicate(No_Of_Var,stats::rexp(n = N, rate = Dist_Par$Rate))
    }
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }

    Complete_Data<-cbind(1,X)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))

    Linear_Predictor_Data <- Complete_Data%*%Beta
    Pi_Data <- 1-1/(1+exp(Linear_Predictor_Data))
    Y <- stats::rbinom(n=N,size=1,prob = Pi_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    stats::glm(Y~.-1,data=as.data.frame(Complete_Data),family = "binomial")->Results
    Beta_Estimate<-stats::coefficients(Results)

    Outputs<-list("Basic"=list("N"=N,"Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimate,
                               "Distribution"=Dist,
                               "Distribution_Parameter"=Dist_Par,
                               "No_Of_Variables"=No_Of_Var),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("A_L_OptimalSubsampling","logistic")
    return(Outputs)
  }
  if(family %in% "poisson"){
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = 0, sd = 1))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = 0, max = 1))
    }

    Complete_Data<-cbind(1,X)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))

    Linear_Predictor_Data <- Complete_Data%*%Beta
    Lambda_Data <- exp(Linear_Predictor_Data)
    Y <- stats::rpois(n=N,lambda = Lambda_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    stats::glm(Y~.-1,data=as.data.frame(Complete_Data),family = "poisson")->Results
    Beta_Estimate<-stats::coefficients(Results)

    Outputs<-list("Basic"=list("N"=N,"Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimate,
                               "Distribution"=Dist,
                               "Distribution_Parameter"=Dist_Par,
                               "No_Of_Variables"=No_Of_Var),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("A_L_OptimalSubsampling","poisson")
    return(Outputs)
  }
}
