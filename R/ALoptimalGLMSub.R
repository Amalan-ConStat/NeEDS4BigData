#' A- and L-optimality criteria based subsampling under Generalised Linear Models
#'
#' Using this function sample from big data under linear, logistic and Poisson regression
#' to describe the data. Subsampling probabilities are obtained based on the A- and L-
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
#' Using the estimated parameters subsampling probabilities are evaluated for A- and L-optimality criteria.
#'
#' Through the estimated subsampling probabilities an optimal sample of size \eqn{r_2 \ge r_1} is obtained.
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
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon in a data.frame after subsampling
#'
#' \code{Utility_Estimates} estimated D-(log scaled), A- and L- optimality values for the obtained subsamples
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A- and L- optimality criteria
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
#' ALoptimalGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
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
#' ALoptimalGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
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
#' ALoptimalGLMSub(r1 = r1, r2 = r2,Y = as.matrix(Original_Data[,1]),
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

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Utility.mVc<-Utility.mMSE<-matrix(nrow = length(r2),ncol = 4 )
    Var_Epsilon<-matrix(nrow = length(r2),ncol = 3)
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Utility.mVc[,1]<-Utility.mMSE[,1]<-Var_Epsilon[,1]<-r2

    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Var_Epsilon)<-c("r2","A-Optimality","L-Optimality")
    colnames(Utility.mVc)<-colnames(Utility.mMSE)<-c("r2","D-optimality","A-optimality","L-optimality")

    ## mVc
    PI.mVc<-sqrt(Epsilon.prop^2 * matrixStats::rowSums2(X^2))
    PI.mVc<-PI.mVc/sum(PI.mVc)

    ## mMSE
    PI.mMSE<-sqrt(Epsilon.prop^2 * matrixStats::rowSums2((X %*% solve(crossprod(X)))^2))
    PI.mMSE<-PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc<-c(1 / PI.mVc[idx.mVc], pinv.prop)

      pi4_r<-sqrt(r2[i]*pinv.mVc^(-1))
      X_r4<-x.mVc/pi4_r
      Y_r4<-y.mVc/pi4_r
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Temp<-Var.prop*crossprod(x.mVc)
      Temp_Inv<-solve(Temp)
      x.mVc_t<-t(x.mVc)
      Temp_Int<-Temp_Inv%*%x.mVc_t
      Temp1<-x.mVc%*%Temp_Int

      Sample.mVc[[i+1]]<-idx.mVc
      beta.mVc[i,-1] <- beta.prop
      Var_Epsilon[i,3]<-Var.prop
      Utility.mVc[i,-1]<-c(log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }

      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)

      pi4_r<-sqrt(r2[i]*pinv.mMSE^(-1))
      X_r4<-x.mMSE/pi4_r
      Y_r4<-y.mMSE/pi4_r
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N

      Temp<-Var.prop*crossprod(x.mMSE)
      Temp_Inv<-solve(Temp)
      x.mMSE_t<-t(x.mMSE)
      Temp_Int<-Temp_Inv%*%x.mMSE_t
      Temp1<-x.mMSE%*%Temp_Int

      Sample.mMSE[[i+1]]<-idx.mMSE;
      beta.mMSE[i,-1] <- beta.prop
      Var_Epsilon[i,2]<-Var.prop
      Utility.mMSE[i,-1]<-c(log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),rbind(beta.mMSE,beta.mVc))

    Utility_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                   rbind(Utility.mMSE,Utility.mVc))

    Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                       "Sample"=rep(r2,times=length(Subsampling_Methods)),
                                       "Var Epsilon"=c(Var_Epsilon[,"A-Optimality"],
                                                       Var_Epsilon[,"L-Optimality"]))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Subsampling_Probability"=Full_SP)
    class(ans)<-c("A_L_OptimalSubsampling","linear")
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

    Xbeta_Final<-X%*%beta.prop
    P.prop  <- 1 - 1 / (1 + exp(Xbeta_Final))

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Utility.mVc<-Utility.mMSE<-matrix(nrow = length(r2),ncol = 4 )
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Utility.mVc[,1]<-Utility.mMSE[,1]<-r2

    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Utility.mVc)<-colnames(Utility.mMSE)<-c("r2","D-optimality","A-optimality","L-optimality")

    ## mVc
    PI.mVc<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2(X^2))
    PI.mVc <- PI.mVc/sum(PI.mVc)

    ## mMSE
    p.prop <- P.prop[idx.prop]
    w.prop <- p.prop * (1 - p.prop)
    W.prop <- solve(crossprod(x.prop,x.prop * w.prop * pinv.prop))

    PI.mMSE<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop), ]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      SS_Prob<-c(1 / PI.mVc[idx.mVc], pinv.prop)
      fit.mVc <-.getMLE(x=x.mVc, y=y.mVc,w=SS_Prob)

      Sample.mVc[[i+1]]<-idx.mVc;
      beta.mVc[i,-1] <- fit.mVc$par

      if(anyNA(fit.mVc$par)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.mVc %*% beta.mVc[i,-1]
      pi<-c(1-1/(1 + exp(LP_data)))
      W<-pi*(1-pi)
      Mx<-crossprod(x.mVc,(x.mVc * W))
      Mx_Inv<-solve(Mx)
      x.mVc_t<-t(x.mVc)
      V_Mx_Temp <- Mx_Inv %*% x.mVc_t
      V_Final<- x.mVc%*% V_Mx_Temp

      Utility.mVc[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      SS_Prob<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)
      fit.mMSE <- .getMLE(x=x.mMSE, y=y.mMSE,w=c(1 / PI.mMSE[idx.mMSE], pinv.prop))

      Sample.mMSE[[i+1]]<-idx.mMSE;
      beta.mMSE[i,-1] <- fit.mMSE$par

      if(anyNA(fit.mMSE$par)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.mMSE %*% beta.mMSE[i,-1]
      pi<-c(1-1/(1 + exp(LP_data)))
      W<-pi*(1-pi)
      Mx<-crossprod(x.mMSE,(x.mMSE * W))
      Mx_Inv<-solve(Mx)
      x.mMSE_t<-t(x.mMSE)
      V_Mx_Temp <- Mx_Inv %*% x.mMSE_t
      V_Final<- x.mMSE%*% V_Mx_Temp

      Utility.mMSE[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))
    }

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                rbind(beta.mMSE,beta.mVc))

    Utility_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                   rbind(Utility.mMSE,Utility.mVc))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Subsampling_Probability"=Full_SP)

    class(ans)<-c("A_L_OptimalSubsampling","logistic")
    return(ans)
  }
  if(family %in% "poisson"){
    PI.prop <- rep(1/N, N)
    idx.prop <- sample(1:N, size = r1, replace = TRUE)

    x.prop<-X[idx.prop,]
    y.prop <- Y[idx.prop,]

    pinv.prop <- N
    pinv.prop <- 1/PI.prop[idx.prop]
    fit.prop <- stats::glm(y.prop~x.prop-1,family = "quasipoisson")

    beta.prop <- fit.prop$coefficients
    if(anyNA(beta.prop)){
      stop("There are NA or NaN values in the model parameters")
    }

    Xbeta_Final<-X %*% beta.prop
    P.prop  <- exp(Xbeta_Final)

    beta.mVc<-beta.mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
    Utility.mVc<-Utility.mMSE<-matrix(nrow = length(r2),ncol = 4 )
    Sample.mMSE<-Sample.mVc<-list()

    Sample.mMSE[[1]]<-Sample.mVc[[1]]<-idx.prop

    beta.mVc[,1]<-beta.mMSE[,1]<-Utility.mVc[,1]<-Utility.mMSE[,1]<-r2

    if(all(X[,1] == 1)){
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
    } else {
      colnames(beta.mMSE)<-colnames(beta.mVc)<-c("r2",paste0("Beta_",1:(ncol(X))))
    }
    colnames(Utility.mVc)<-colnames(Utility.mMSE)<-c("r2","D-optimality","A-optimality","L-optimality")

    ## mVc
    PI.mVc <- sqrt((Y - P.prop)^2 * matrixStats::rowSums2(X^2))
    PI.mVc <- PI.mVc/sum(PI.mVc)

    ## mMSE
    w.prop <- P.prop[idx.prop]
    W.prop <- solve(crossprod(x.prop,x.prop * w.prop * pinv.prop))

    PI.mMSE<-sqrt((Y - P.prop)^2 * matrixStats::rowSums2((X%*%W.prop)^2))
    PI.mMSE <- PI.mMSE/sum(PI.mMSE)

    message("Step 1 of the algorithm completed.\n")

    for (i in 1:length(r2))
    {
      ## mVc
      idx.mVc <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mVc)

      x.mVc <- X[c(idx.mVc, idx.prop),]
      y.mVc <- Y[c(idx.mVc, idx.prop)]
      pinv.mVc<-c(1 / PI.mVc[idx.mVc], pinv.prop)

      fit.mVc <-stats::glm(y.mVc~x.mVc-1, family = "quasipoisson",weights=pinv.mVc)
      Sample.mVc[[i+1]]<-idx.mVc;
      beta.mVc[i,-1] <- fit.mVc$coefficients

      if(anyNA(fit.mVc$coefficients)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.mVc %*% beta.mVc[i,-1]
      pi<-c(exp(LP_data))
      W<-pi
      Mx<-crossprod(x.mVc,(x.mVc * W))
      Mx_Inv<-solve(Mx)
      x.mVc_t<-t(x.mVc)
      V_Mx_Temp <- Mx_Inv %*% x.mVc_t
      V_Final<- x.mVc %*% V_Mx_Temp

      Utility.mVc[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))

      ## mMSE
      idx.mMSE <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mMSE)

      x.mMSE <- X[c(idx.mMSE, idx.prop),]
      y.mMSE <- Y[c(idx.mMSE, idx.prop)]
      pinv.mMSE<-c(1 / PI.mMSE[idx.mMSE], pinv.prop)

      fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1, family = "quasipoisson",weights=pinv.mMSE)
      Sample.mMSE[[i+1]]<-idx.mMSE;
      beta.mMSE[i,-1] <-fit.mMSE$coefficients

      if(anyNA(fit.mMSE$coefficients)){
        stop("There are NA or NaN values in the model parameters")
      }

      LP_data<-x.mMSE %*% beta.mMSE[i,-1]
      pi<-c(exp(LP_data))
      W<-pi
      Mx<-crossprod(x.mMSE,(x.mMSE * W))
      Mx_Inv<-solve(Mx)
      x.mMSE_t<-t(x.mMSE)
      V_Mx_Temp <- Mx_Inv %*% x.mMSE_t
      V_Final<- x.mMSE %*% V_Mx_Temp

      Utility.mMSE[i,-1]<-c(log(det(Mx)),psych::tr(Mx_Inv),psych::tr(V_Final))
    }

    Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc)
    colnames(Full_SP)<-c("A-Optimality","L-Optimality")

    Subsampling_Methods<-factor(c("A-Optimality","L-Optimality"))

    Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                rbind(beta.mMSE,beta.mVc))

    Utility_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                   rbind(Utility.mMSE,Utility.mVc))

    names(Sample.mMSE)<-names(Sample.mVc)<-c(r1,r2)

    message("Step 2 of the algorithm completed.")

    ans<-list("Beta_Estimates"=Beta_Data,
              "Utility_Estimates"=Utility_Data,
              "Sample_A-Optimality"=Sample.mMSE,
              "Sample_L-Optimality"=Sample.mVc,
              "Subsampling_Probability"=Full_SP)
    class(ans)<-c("A_L_OptimalSubsampling","poisson")
    return(ans)
  }
}
