#' Model robust optimal subsampling for A- and L- optimality criterion under linear regression
#'
#' Using this function subsample from big data under linear regression when there are more than
#' one model to describe the data. Subsampling probabilities are obtained based on the A- and L-
#' optimality criterions.
#'
#' @usage
#' modelRobustLinSub(r1,r2,Y,X,N,Alpha,All_Combinations,All_Covariates)
#'
#' @param r1      subsample size for initial random sampling
#' @param r2      subsample size for optimal sampling
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#' @param Alpha   vector of alpha values that are used to obtain the model robust subsampling probabilities
#' @param All_Combinations list of possible models that can describe the data
#' @param All_Covariates all the covariates in the models
#'
#' @details
#' Two stage subsampling algorithm for big data under linear regression for multiple models that can
#' describe the big data.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters for all models.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-optimality criterion and
#' model averaging A-, L-optimality subsampling methods.
#'
#' Through the estimated subsampling probabilities an subsample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for all the models.
#'
#' \strong{NOTE} :  If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r_2 \ge r_1} is not satisfied then an error message will be produced.
#'
#' If the big data \eqn{X,Y} has any missing values then an error message will be produced.
#'
#' The big data size \eqn{N} is compared with the sizes of \eqn{X,Y} and
#' if they are not aligned an error message will be produced.
#'
#' If \eqn{0 < \alpha < 1} for the a priori probabilities are not satisfied an error message will be produced.
#'
#' @return
#' The output of \code{modelRobustLinSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters for each model in a list after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon for each model after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal subsamples obtained based on A-Optimality criterion
#'
#' \code{Sample_A-Optimality_MR} list of indexes for the initial and model robust optimal subsamples obtained based on A-Optimality criterion
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal subsamples obtained based on L-Optimality criterion
#'
#' \code{Sample_L-Optimality_MR} list of indexes for the initial and model robust optimal subsamples obtained based on L-Optimality criterion
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A- and L- optimality criterion
#'
#' @references
#' \insertRef{mahendran2023model}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000
#' All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
#'                  Assumed_Model_1=c("X0","X1","X2"),
#'                  Assumed_Model_2=c("X0","X1","X2","X2^2"),
#'                  Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))
#' family = "linear"
#'
#' Full_Data<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)
#'
#' r1<-300; r2<-rep(100*c(6,9,12),25); Original_Data<-Full_Data$Complete_Data;
#'
#' modelRobustLinSub(r1 = r1, r2 = r2,
#'                   Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                   X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                   Alpha = rep(1/length(All_Models),length(All_Models)),
#'                   All_Combinations = All_Models,
#'                   All_Covariates = colnames(Original_Data)[-1])->Results
#'
#' Beta_Plots<-plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @importFrom matrixStats rowSums2
#' @export
modelRobustLinSub <- function(r1,r2,Y,X,N,Alpha,All_Combinations,All_Covariates){
  if(any(is.na(c(r1,r2,N,Alpha,All_Covariates))) | any(is.nan(c(r1,r2,N,Alpha,All_Covariates)))){
    stop("NA or Infinite or NAN values in the r1,r2,N,Alpha or All_Covariates")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if(any(is.na(cbind(Y,X))) | any(is.nan(cbind(Y,X)))){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r1) > r2)){
    stop("2*r1 cannot be greater than r2 at any point")
  }

  if(length(Alpha) != length(All_Combinations)){
    stop("No of models for averaging is not equal to the a priori alpha values")
  }

  if(any(Alpha > 1) | any(Alpha < 0) | sum(Alpha) != 1){
    stop("A priori alpha is not inbetween zero and one or the sum is one")
  }

  PI.prop <- rep(1/N, N)
  idx.prop <- sample(1:N, r1, T)

  x.prop<-lapply(1:length(All_Combinations),function(j){
    X[idx.prop,All_Covariates %in% All_Combinations[[j]]]
  })

  y.prop <- Y[idx.prop,]

  pinv.prop <- N
  pinv.prop <- 1/PI.prop[idx.prop]

  fit.prop <- lapply(1:length(All_Combinations), function(j){
    beta.prop<-solve(a=t(x.prop[[j]])%*%x.prop[[j]],b=t(x.prop[[j]])%*%y.prop)
    Xbeta_Final<-as.vector(X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop)
    Var.prop<-sum((Y-Xbeta_Final)^2)/N
    Epsilon.prop<-Y-Xbeta_Final
    return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop,"Epsilon.prop"=Epsilon.prop))
  })

  beta.prop<-list()
  for (j in 1:length(All_Combinations))
  {
    beta.prop[[j]] <- fit.prop[[j]]$beta
    if(anyNA(beta.prop[[j]])){
      stop("There are NA or NaN values in the model parameters")
    }
  }

  # Single Model Results
  beta.mVc_Single<-beta.mMSE_Single<-list()
  Var_Epsilon_mVc<-Var_Epsilon_mMSE<-matrix(nrow = length(r2),ncol = length(All_Combinations) + 1)
  Sample.mMSE_Single<-Sample.mVc_Single<-list()

  # Model Robust Results
  beta.mVc_MR<-beta.mMSE_MR<-list()
  Var_Epsilon_mVc_MR<-Var_Epsilon_mMSE_MR<-matrix(nrow = length(r2),ncol = length(All_Combinations) + 1)
  Sample.mMSE_MR<-Sample.mVc_MR<-list()

  Var_Epsilon_mVc[,1]<-Var_Epsilon_mMSE[,1]<-Var_Epsilon_mVc_MR[,1]<-Var_Epsilon_mMSE_MR[,1]<-r2

  # For the models, Single and Model Robust
  for (a in 1:length(All_Combinations))
  {
    beta.mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Single Model Results
    beta.mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Sample.mMSE_Single[[a]]<-Sample.mVc_Single[[a]]<-list()

    beta.mVc_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Model Robust Results
    beta.mMSE_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Sample.mMSE_MR[[a]]<-Sample.mVc_MR[[a]]<-list()

    Sample.mMSE_Single[[a]][[1]]<-Sample.mVc_Single[[a]][[1]]<-
      Sample.mMSE_MR[[a]][[1]]<-Sample.mVc_MR[[a]][[1]]<-idx.prop

    colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
      colnames(beta.mMSE_MR[[a]])<-c("r2",paste0("beta_",0:(length(All_Combinations[[a]])-1)))
  }

  ## mVc
  PI_Single.mVc <- lapply(1:length(All_Combinations), function(j){ # Single Model Results
    PI.mVc<-sqrt(fit.prop[[j]]$Epsilon.prop^2 * matrixStats::rowSums2(X[,All_Covariates %in% All_Combinations[[j]] ]^2))
    return(PI.mVc/sum(PI.mVc))
  })
  PI_MR.mVc<-matrixStats::rowSums2(do.call(cbind,PI_Single.mVc)%*%diag(Alpha)) # Model Robust Results

  ## mMSE
  PI_Single.mMSE <- lapply(1:length(All_Combinations),function(j){ # Single Model Results
    PI.mMSE<-sqrt(fit.prop[[j]]$Epsilon.prop^2 *
                  matrixStats::rowSums2((X[,All_Covariates %in% All_Combinations[[j]] ] %*%
                                           solve(t(X[,All_Covariates %in% All_Combinations[[j]] ])%*%
                                                   X[,All_Covariates %in% All_Combinations[[j]] ]))^2))
    return(PI.mMSE/sum(PI.mMSE))
  })
  PI_MR.mMSE<-matrixStats::rowSums2(do.call(cbind,PI_Single.mMSE)%*%diag(Alpha))  # Model Robust Results

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    ## mVc
    idx_Single.mVc <- lapply(1:length(All_Combinations), function(j){
      sample(1:N, r2[i]-r1, T, PI_Single.mVc[[j]]) # Single Model Results
    })
    idx_MR.mVc <- sample(1:N, r2[i]-r1, T, PI_MR.mVc) # Model Robust Results

    x_Single.mVc <-lapply(1:length(All_Combinations),function(j){ # Single Model Results
      X[c(idx_Single.mVc[[j]], idx.prop),All_Covariates %in% All_Combinations[[j]] ]
    })
    y_Single.mVc <- lapply(1:length(All_Combinations),function(j){
      Y[c(idx_Single.mVc[[j]], idx.prop)]  # Single Model Results
    })

    x_MR.mVc <-  lapply(1:length(All_Combinations),function(j){
      X[c(idx_MR.mVc, idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Model Robust Results
    })
    y_MR.mVc <- Y[c(idx_MR.mVc, idx.prop)]
    pinv_MR.mVc<-c(1 / PI_MR.mVc[idx_MR.mVc], pinv.prop)

    fit_Single.mVc <-lapply(1:length(All_Combinations), function(j){ # Single Model Results
      pinv_Single.mVc<-c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)
      pi4_r<-sqrt(r2[i]*pinv_Single.mVc^(-1))
      X_r4<-x_Single.mVc[[j]]/pi4_r
      Y_r4<-y_Single.mVc[[j]]/pi4_r
      beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4)
      Xbeta_Final<-as.vector(X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    fit_MR.mVc <- lapply(1:length(All_Combinations),function(j){
      pi4_r<-sqrt(r2[i]*pinv_MR.mVc^(-1))
      X_r4<-x_MR.mVc[[j]]/pi4_r
      Y_r4<-y_MR.mVc/pi4_r
      beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4) # Model Robust Results
      Xbeta_Final<-as.vector(X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    for (j in 1:length(All_Combinations))
    {
      Sample.mVc_Single[[j]][[i+1]]<-idx_Single.mVc[[j]]
      Sample.mVc_MR[[j]][[i+1]]<-idx_MR.mVc

      beta.mVc_Single[[j]][i,] <- c(r2[i],fit_Single.mVc[[j]]$beta.prop)
      beta.mVc_MR[[j]][i,] <- c(r2[i],fit_MR.mVc[[j]]$beta.prop)

      Var_Epsilon_mVc[i,j+1]<-fit_Single.mVc[[j]]$Var.prop
      Var_Epsilon_mVc_MR[i,j+1]<-fit_MR.mVc[[j]]$Var.prop

      if(anyNA(fit_Single.mVc[[j]]$beta.prop) || anyNA(fit_MR.mVc[[j]]$beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    ## mMSE
    idx_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      sample(1:N, r2[i]-r1, T, PI_Single.mMSE[[j]]) # Single Model Results
    })
    idx_MR.mMSE <- sample(1:N, r2[i]-r1, T, PI_MR.mMSE) # Model Robust Results

    x_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      X[c(idx_Single.mMSE[[j]], idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Single Model Results
    })
    y_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      Y[c(idx_Single.mMSE[[j]], idx.prop)] # Single Model Results
    })

    x_MR.mMSE <- lapply(1:length(All_Combinations),function(j){
      X[c(idx_MR.mMSE, idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Model Robust Results
    })
    y_MR.mMSE <- Y[c(idx_MR.mMSE, idx.prop)] # Model Robust Results
    pinv_MR.mMSE<-c(1 / PI_MR.mMSE[idx_MR.mMSE], pinv.prop)

    fit_Single.mMSE <-lapply(1:length(All_Combinations), function(j){ # Single Model Results
      pinv_Single.mMSE<-c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)
      pi4_r<-sqrt(r2[i]*pinv_Single.mMSE^(-1))
      X_r4<-x_Single.mMSE[[j]]/pi4_r
      Y_r4<-y_Single.mMSE[[j]]/pi4_r
      beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4)
      Xbeta_Final<-as.vector(X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    fit_MR.mMSE <- lapply(1:length(All_Combinations),function(j){
      pi4_r<-sqrt(r2[i]*pinv_MR.mMSE^(-1))
      X_r4<-x_MR.mMSE[[j]]/pi4_r
      Y_r4<-y_MR.mMSE/pi4_r
      beta.prop<-solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4) # Model Robust Results
      Xbeta_Final<-as.vector(X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop)
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    for (j in 1:length(All_Combinations))
    {
      Sample.mMSE_Single[[j]][[i+1]]<-idx_Single.mMSE[[j]]
      Sample.mMSE_MR[[j]][[i+1]]<-idx_MR.mMSE

      beta.mMSE_Single[[j]][i,] <- c(r2[i],fit_Single.mMSE[[j]]$beta.prop)
      beta.mMSE_MR[[j]][i,] <- c(r2[i],fit_MR.mMSE[[j]]$beta.prop)

      Var_Epsilon_mMSE[i,j+1]<-fit_Single.mMSE[[j]]$Var.prop
      Var_Epsilon_mMSE_MR[i,j+1]<-fit_MR.mMSE[[j]]$Var.prop

      if(anyNA(fit_Single.mMSE[[j]]$beta.prop) || anyNA(fit_MR.mMSE[[j]]$beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }
  }

  if(anyNA(beta.mMSE_Single) || anyNA(beta.mVc_Single) || anyNA(beta.mMSE_MR) || anyNA(beta.mVc_MR)){
    stop("There are NA or NaN values")
  }

  Full_SP<-cbind.data.frame(do.call(cbind,PI_Single.mMSE),do.call(cbind,PI_Single.mVc),
                            PI_MR.mMSE,PI_MR.mVc)
  colnames(Full_SP)<-c(paste0(paste0("A-Optimality M",1:length(All_Combinations))),
                       paste0("L-Optimality M",1:length(All_Combinations)),
                       "MR A-Optimality","MR L-Optimality")

  Subsampling_Methods<-factor(c("A-Optimality","L-Optimality","MR A-Optimality","MR L-Optimality"))
  Beta_Data<-list()

  for (j in 1:length(All_Combinations))
  {
    Beta_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     rbind(beta.mMSE_Single[[j]],beta.mVc_Single[[j]],
                                           beta.mMSE_MR[[j]],beta.mVc_MR[[j]]))

    names(Sample.mVc_Single[[j]])<-names(Sample.mMSE_Single[[j]])<-
      names(Sample.mVc_MR[[j]])<-names(Sample.mMSE_MR[[j]])<-c(r1,r2)
  }

  Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     rbind(Var_Epsilon_mMSE,Var_Epsilon_mVc,
                                           Var_Epsilon_mMSE_MR,Var_Epsilon_mVc_MR))

  names(Beta_Data)<-paste0("Model_",1:length(All_Combinations))
  colnames(Var_Epsilon_Data)[-1]<-c("r2",paste0("Model_",1:length(All_Combinations)))

  names(Sample.mVc_Single)<-names(Sample.mVc_MR)<-paste0("Model_",1:length(All_Combinations))
  names(Sample.mMSE_Single)<-names(Sample.mMSE_MR)<-paste0("Model_",1:length(All_Combinations))

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
            "Sample_A-Optimality"=Sample.mMSE_Single,
            "Sample_A-Optimality_MR"=Sample.mMSE_MR,
            "Sample_L-Optimality"=Sample.mVc_Single,
            "Sample_L-Optimality_MR"=Sample.mVc_MR,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelRobust","linear")
  return(ans)
}

#' Generate data for Generalised Linear Models under model robust scenario
#'
#' Function to simulate big data under linear, logistic and Poisson regression for the model robust scenario
#' through a set of models.
#' Covariate data X is through Normal or Uniform distribution for linear regression.
#' Covariate data X is through Exponential or Normal or Uniform distribution for logistic regression.
#' Covariate data X is through Normal or Uniform distribution for Poisson regression.
#'
#' @usage
#' GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)
#'
#' @param Dist        a character value for the distribution "Normal" or "Uniform
#' @param Dist_Par    a list of parameters for the distribution that would generate data for covariate X
#' @param No_Of_Var   number of variables
#' @param Beta        a vector for the model parameters, including the intercept
#' @param N           the big data size
#' @param All_Models  a list that contains the possible models
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
#' For a given real model data is generated and then this data is modelled by All_Models.
#'
#' @return
#' The output of \code{GenModelRobustGLMdata} gives a list of
#'
#' \code{Basic} a list of outputs based on the inputs and Beta Estimates for all models
#'
#' \code{Complete_Data} a matrix for Y,X and X^2
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000
#' All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
#'                  Assumed_Model_1=c("X0","X1","X2"),
#'                  Assumed_Model_2=c("X0","X1","X2","X2^2"),
#'                  Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))
#' family<-"linear"
#'
#' Results<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000
#' All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
#'                  Assumed_Model_1=c("X0","X1","X2"),
#'                  Assumed_Model_2=c("X0","X1","X2","X2^2"),
#'                  Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))
#' family = "logistic"
#'
#' Results<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)
#'
#' Dist<-"Normal";
#' No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000
#' All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
#'                  Assumed_Model_1=c("X0","X1","X2"),
#'                  Assumed_Model_2=c("X0","X1","X2","X2^2"),
#'                  Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))
#' family = "poisson"
#'
#' Results<-GenModelRobustGLMdata(Dist,Dist_Par=NULL,No_Of_Var,Beta,N,All_Models,family)
#'
#' @import stats
#' @export
GenModelRobustGLMdata<-function(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family){
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

  if(family == "linear"){
    if(Dist == "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }
    if(Dist == "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }

    Complete_Data<-cbind(1,X,X^2);
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))
    Residual<-stats::rnorm(n=N,mean=0,sd=sqrt(Dist_Par$Error_Variance))
    Y <- Complete_Data[,colnames(Complete_Data) %in% All_Models$Real_Model]%*%Beta + Residual

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))

    Beta_Estimates<-list()
    Var_Epsilon_Estimates<-NULL
    for (j in 1:length(All_Models))
    {
      Temp_Data<-Complete_Data[,colnames(Complete_Data) %in% c("Y",All_Models[[j]])]

      stats::lm(Y~.-1,data=as.data.frame(Temp_Data))->Results
      Beta_Estimates[[j]]<-stats::coefficients(Results)
      Var_Epsilon_Estimates[j]<-mean((Y-stats::fitted.values(Results))^2)
    }
    names(Var_Epsilon_Estimates)<-paste0("Model ",1:length(All_Models))

    Outputs<-list("Basic"=list("N"=N,
                               "Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimates,
                               "Variance_Epsilon_Estimates"=Var_Epsilon_Estimates,
                               "Distribution"=Dist,
                               "Distribution_Parameter"=Dist_Par,
                               "No_Of_Variables"=No_Of_Var,
                               "All_Models"=All_Models),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("ModelRobust","linear")
    return(Outputs)
  }
  if(family == "logistic"){
    if(Dist %in% "Exponential"){
      X<-replicate(No_Of_Var,stats::rexp(n = N, rate = Dist_Par$Rate))
    }
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }

    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }

    Complete_Data<-cbind(1,X,X^2)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))

    Linear_Predictor_Data <- Complete_Data[,colnames(Complete_Data) %in% All_Models$Real_Model]%*%Beta
    Pi_Data <- 1-1/(1+exp(Linear_Predictor_Data))
    Y <- stats::rbinom(n=N,size=1,prob = Pi_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))

    Beta_Estimates<-list()
    for (j in 1:length(All_Models))
    {
      Temp_Data<-Complete_Data[,colnames(Complete_Data) %in% c("Y",All_Models[[j]])]

      stats::glm(Y~.-1,data=as.data.frame(Temp_Data),family = "binomial")->Results
      Beta_Estimates[[j]]<-stats::coefficients(Results)
    }

    Outputs<-list("Basic"=list("N"=N,
                               "Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimates,
                               "Distribution"=Dist,
                               "Distribution_Parameter"=Dist_Par,
                               "No_Of_Variables"=No_Of_Var,
                               "All_Models"=All_Models),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("ModelRobust","logistic")
    return(Outputs)
  }
  if(family=="poisson"){
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = 0, sd = 1))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = 0, max = 1))
    }

    Complete_Data<-cbind(1,X,X^2)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))

    Linear_Predictor_Data <- Complete_Data[,colnames(Complete_Data) %in% All_Models$Real_Model]%*%Beta
    Lambda_Data <- 1-1/(1+exp(Linear_Predictor_Data))
    Y <- stats::rpois(n=N,lambda = Lambda_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)),paste0("X",1:ncol(X),"^2"))

    Beta_Estimates<-list()
    for (j in 1:length(All_Models))
    {
      Temp_Data<-Complete_Data[,colnames(Complete_Data) %in% c("Y",All_Models[[j]])]

      stats::glm(Y~.-1,data=as.data.frame(Temp_Data),family = "poisson")->Results
      Beta_Estimates[[j]]<-stats::coefficients(Results)
    }

    Outputs<-list("Basic"=list("N"=N,
                               "Beta"=Beta,
                               "Beta_Estimates"=Beta_Estimates,
                               "Distribution"=Dist,
                               "No_Of_Variables"=No_Of_Var,
                               "All_Models"=All_Models),
                  "Complete_Data"=Complete_Data)

    class(Outputs)<-c("ModelRobust","poisson")
    return(Outputs)
  }
}
