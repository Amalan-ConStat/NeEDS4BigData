#' Model robust optimal subsampling for A- and L- optimality criterion under logistic regression
#'
#' Using this function subsample from big data under logistic regression when there are more than
#' one model to describe the data. Subsampling probabilities are obtained based on the A- and L-
#' optimality criterions.
#'
#' @usage
#' modelRobustLogSub(r1,r2,Y,X,N,Alpha,All_Combinations,All_Covariates)
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
#'Two stage subsampling algorithm for big data under logistic regression for multiple models that can
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
#' \code{Beta_Data} estimated model parameters for each model in a list after subsampling
#'
#' \code{Utility_Data} estimated Variance and Information of the model parameters after subsampling
#'
#' \code{Sample_L-optimality} list of indexes for the initial and optimal subsamples obtained based on L-optimality criterion
#'
#' \code{Sample_L-optimality_MR} list of indexes for the initial and model robust optimal subsamples obtained based on L-optimality criterion
#'
#' \code{Sample_A-optimality} list of indexes for the initial and optimal subsamples obtained based on A-optimality criterion
#'
#' \code{Sample_A-optimality_MR} list of indexes for the initial and model robust optimal subsamples obtained based on A-optimality criterion
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A- and L- optimality criterion
#'
#' @references
#' \insertRef{mahendran2023model}{NeEDS4BigData}
#'
#' @examples
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1)
#' No_Of_Var<-2; Beta<-c(-1,2,1,2); N<-10000
#' All_Models<-list(Real_Model=c("X0","X1","X2","X1^2"),
#'                  Assumed_Model_1=c("X0","X1","X2"),
#'                  Assumed_Model_2=c("X0","X1","X2","X2^2"),
#'                  Assumed_Model_3=c("X0","X1","X2","X1^2","X2^2"))
#' family = "logistic"
#'
#' Full_Data<-GenModelRobustGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,All_Models,family)
#'
#' r1<-300; r2<-rep(100*c(6,9,12),25); Original_Data<-Full_Data$Complete_Data;
#'
#' modelRobustLogSub(r1 = r1, r2 = r2,
#'                   Y = as.matrix(Original_Data[,colnames(Original_Data) %in% c("Y")]),
#'                   X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                   Alpha = rep(1/length(All_Models),length(All_Models)),
#'                   All_Combinations = All_Models,
#'                   All_Covariates = colnames(Original_Data)[-1])->Results
#'
#' Beta_Plots<-plot_Beta(Results)
#' Utility_Plots<-plot_Utility(Results)
#'
#' @importFrom Rdpack reprompt
#' @importFrom psych tr
#' @importFrom matrixStats rowSums2
#' @export
modelRobustLogSub <- function(r1,r2,Y,X,N,Alpha,All_Combinations,All_Covariates){
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

  n1 <- sum(Y)
  n0 <- N - n1
  PI.prop <- rep(1/(2*n0), N)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:N, r1, T, PI.prop)

  x.prop<-lapply(1:length(All_Combinations),function(j){
    X[idx.prop,All_Covariates %in% All_Combinations[[j]] ]
  })
  y.prop <- Y[idx.prop,]

  pinv.prop <- 1/PI.prop[idx.prop]

  fit.prop <- lapply(1:length(All_Combinations), function(j){
    .getMLE(x=x.prop[[j]], y=y.prop, w=pinv.prop)
  })

  beta.prop<-list()
  for (j in 1:length(All_Combinations))
  {
    beta.prop[[j]] <- fit.prop[[j]]$par
    if(anyNA(beta.prop[[j]])){
      stop("There are NA or NaN values in the model parameters")
    }
  }

  P.prop  <- lapply(1:length(All_Combinations),function(j){
    1 - 1 / (1 + exp(X[,All_Covariates %in% All_Combinations[[j]] ] %*% beta.prop[[j]]))
  })

  # Single Model Results
  beta.mVc_Single<-Utility_mVc_Single<-beta.mMSE_Single<-Utility_mMSE_Single<-list()
  Sample.mMSE_Single<-Sample.mVc_Single<-list()

  # Model Robust Results
  beta.mVc_MR<-Utility_mVc_MR<-beta.mMSE_MR<-Utility_mMSE_MR<-list()
  Sample.mMSE_MR<-Sample.mVc_MR<-list()

  # For the models, Single and Model Robust
  for (a in 1:length(All_Combinations))
  {
    beta.mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Single Model Results
    Utility_mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = 3 )
    beta.mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Utility_mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = 3 )
    Sample.mMSE_Single[[a]]<-Sample.mVc_Single[[a]]<-list()

    beta.mVc_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Model Robust Results
    Utility_mVc_MR[[a]]<-matrix(nrow = length(r2),ncol = 3 )
    beta.mMSE_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Utility_mMSE_MR[[a]]<-matrix(nrow = length(r2),ncol = 3 )
    Sample.mMSE_MR[[a]]<-Sample.mVc_MR[[a]]<-list()

    Sample.mMSE_Single[[a]][[1]]<-Sample.mVc_Single[[a]][[1]]<-
      Sample.mMSE_MR[[a]][[1]]<-Sample.mVc_MR[[a]][[1]]<-idx.prop

    colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
      colnames(beta.mMSE_MR[[a]])<-c("r2",paste0("beta_",0:(length(All_Combinations[[a]])-1)))

    colnames(Utility_mVc_Single[[a]])<-colnames(Utility_mMSE_Single[[a]])<-
      colnames(Utility_mVc_MR[[a]])<-colnames(Utility_mMSE_MR[[a]])<-c("r2","Variance","Information")
  }

  ## mVc
  PI_Single.mVc <- lapply(1:length(All_Combinations), function(j){ # Single Model Results
    PI.mVc<-sqrt((Y - P.prop[[j]])^2 * matrixStats::rowSums2(X[,All_Covariates %in% All_Combinations[[j]] ]^2))
    return(PI.mVc/sum(PI.mVc))
  })
  PI_MR.mVc<-matrixStats::rowSums2(do.call(cbind,PI_Single.mVc)%*%diag(Alpha)) # Model Robust Results

  ## mMSE
  p_Single.prop <- lapply(1:length(All_Combinations),function(j){
    P.prop[[j]][idx.prop] # Single Model Results
  })
  w_Single.prop <- lapply(1:length(All_Combinations),function(j){
    p_Single.prop[[j]] * (1 - p_Single.prop[[j]]) # Single Model Results
  })
  W_Single.prop <- lapply(1:length(All_Combinations),function(j){
    solve(t(x.prop[[j]]) %*% (x.prop[[j]] * w_Single.prop[[j]] * pinv.prop)) # Single Model Results
  })
  PI_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
    PI.mMSE<-sqrt((Y - P.prop[[j]])^2 * matrixStats::rowSums2((X[,All_Covariates %in% All_Combinations[[j]] ]%*%
                                                                 W_Single.prop[[j]])^2)) # Single Model Results
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

    x_Single.mVc <-lapply(1:length(All_Combinations),function(j){
      X[c(idx_Single.mVc[[j]], idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Single Model results
    })
    y_Single.mVc <- lapply(1:length(All_Combinations),function(j){
      Y[c(idx_Single.mVc[[j]], idx.prop)]  # Single Model Results
    })

    x_MR.mVc <-  lapply(1:length(All_Combinations),function(j){
      X[c(idx_MR.mVc, idx.prop),All_Covariates %in% All_Combinations[[j]] ] # Model Robust Results
    })
    y_MR.mVc <- Y[c(idx_MR.mVc, idx.prop)]

    fit_Single.mVc <-lapply(1:length(All_Combinations), function(j){
      .getMLE(x=x_Single.mVc[[j]], y=y_Single.mVc[[j]],
             w=c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)) # Single Model Results
    })

    fit_MR.mVc <- lapply(1:length(All_Combinations),function(j){
      .getMLE(x=x_MR.mVc[[j]], y=y_MR.mVc,w=c(1 / PI_MR.mVc[idx_MR.mVc], pinv.prop)) # Model Robust Results
    })

    for (j in 1:length(All_Combinations))
    {
      Sample.mVc_Single[[j]][[i+1]]<-idx_Single.mVc[[j]]
      Sample.mVc_MR[[j]][[i+1]]<-idx_MR.mVc

      beta.mVc_Single[[j]][i,] <- c(r2[i],fit_Single.mVc[[j]]$par)
      beta.mVc_MR[[j]][i,] <- c(r2[i],fit_MR.mVc[[j]]$par)

      if(anyNA(fit_Single.mVc[[j]]$par) || anyNA(fit_MR.mVc[[j]]$par)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    # Single Model Results
    V_Final<-lapply(1:length(All_Combinations),function(j){
      pi<-1-1/(1 + exp(x_Single.mVc[[j]] %*% beta.mVc_Single[[j]][i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop))
      Mx<-solve((t(x_Single.mVc[[j]]) %*% (x_Single.mVc[[j]] * W)))

      Middle<-((as.vector(y_Single.mVc[[j]])-as.vector(pi))*
                 as.vector(c(1 / PI_Single.mVc[[j]][idx_Single.mVc[[j]]], pinv.prop)))^2
      V_Temp<-(t(x_Single.mVc[[j]]) %*% (x_Single.mVc[[j]] * Middle) )

      Mx %*% V_Temp %*% Mx
    })

    for (j in 1:length(All_Combinations))
    {
      Utility_mVc_Single[[j]][i,]<-c(r2[i],psych::tr(V_Final[[j]]),det(solve(V_Final[[j]])))
    }

    # Model Robust results
    V_Final<-lapply(1:length(All_Combinations),function(j){
      pi<-1-1/(1+exp(x_MR.mVc[[j]] %*% beta.mVc_MR[[j]][i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI_MR.mVc[idx_MR.mVc], pinv.prop))
      Mx<-solve((t(x_MR.mVc[[j]]) %*% (x_MR.mVc[[j]] * W) ))

      Middle<-((as.vector(y_MR.mVc)-as.vector(pi))*as.vector(c(1 / PI_MR.mVc[idx_MR.mVc], pinv.prop)))^2
      V_Temp<-(t(x_MR.mVc[[j]]) %*% (x_MR.mVc[[j]] * Middle))

      Mx %*% V_Temp %*% Mx
    })

    for (j in 1:length(All_Combinations))
    {
      Utility_mVc_MR[[j]][i,]<-c(r2[i],psych::tr(V_Final[[j]]),det(solve(V_Final[[j]])))
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

    fit_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      .getMLE(x=x_Single.mMSE[[j]], y=y_Single.mMSE[[j]],
             w=c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)) # Single Model Results
    })

    fit_MR.mMSE <- lapply(1:length(All_Combinations), function(j){
      .getMLE(x=x_MR.mMSE[[j]], y=y_MR.mMSE,w=c(1 / PI_MR.mMSE[idx_MR.mMSE], pinv.prop)) # Model Robust Results
    })

    for (j in 1:length(All_Combinations))
    {
      # Single Model results
      Sample.mMSE_Single[[j]][[i+1]]<-idx_Single.mMSE[[j]]

      Sample.mMSE_MR[[j]][[i+1]]<-idx_MR.mMSE # Model Robust Results

      beta.mMSE_Single[[j]][i,] <- c(r2[i],fit_Single.mMSE[[j]]$par) # Single Model Results
      beta.mMSE_MR[[j]][i,] <- c(r2[i],fit_MR.mMSE[[j]]$par) # Model Robust Results

      if(anyNA(fit_Single.mMSE[[j]]$par) || anyNA(fit_MR.mMSE[[j]]$par)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    # Single Model
    V_Final<-lapply(1:length(All_Combinations),function(j){
      pi<-1-1/(1 + exp(x_Single.mMSE[[j]] %*% beta.mMSE_Single[[j]][i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop))
      Mx<-solve((t(x_Single.mMSE[[j]]) %*% (x_Single.mMSE[[j]] * W)))

      Middle<-((as.vector(y_Single.mMSE[[j]])-as.vector(pi))*
                 as.vector(c(1 / PI_Single.mMSE[[j]][idx_Single.mMSE[[j]]], pinv.prop)))^2
      V_Temp<-(t(x_Single.mMSE[[j]]) %*% (x_Single.mMSE[[j]] * Middle))

      Mx %*% V_Temp %*% Mx
    })

    for (j in 1:length(All_Combinations))
    {
      Utility_mMSE_Single[[j]][i,]<-c(r2[i],psych::tr(V_Final[[j]]),det(solve(V_Final[[j]])))
    }

    # Model Robust Results
    V_Final<-lapply(1:length(All_Combinations),function(j){
      pi<-1-1/(1+exp(x_MR.mMSE[[j]] %*% beta.mMSE_MR[[j]][i,-1]))
      W<-as.vector(pi*(1-pi)*c(1 / PI_MR.mMSE[idx_MR.mMSE], pinv.prop))
      Mx<-solve((t(x_MR.mMSE[[j]]) %*% (x_MR.mMSE[[j]] * W)))

      Middle<-((as.vector(y_MR.mMSE)-as.vector(pi))*as.vector(c(1 / PI_MR.mMSE[idx_MR.mMSE], pinv.prop)))^2
      V_Temp<-(t(x_MR.mMSE[[j]]) %*% (x_MR.mMSE[[j]] * Middle))

      Mx %*% V_Temp %*% Mx
    })

    for (j in 1:length(All_Combinations))
    {
      Utility_mMSE_MR[[j]][i,]<-c(r2[i],psych::tr(V_Final[[j]]),det(solve(V_Final[[j]])))
    }
  }

  Full_SP<-cbind.data.frame(do.call(cbind,PI_Single.mMSE),do.call(cbind,PI_Single.mVc),
                            PI_MR.mMSE,PI_MR.mVc)
  colnames(Full_SP)<-c(paste0(paste0("A-Optimality M",1:length(All_Combinations))),
                       paste0("L-Optimality M",1:length(All_Combinations)),
                       "MR A-Optimality","MR L-Optimality")

  Subsampling_Methods<-factor(c("A-Optimality","L-Optimality","MR A-Optimality","MR L-Optimality"))
  Beta_Data<-Utility_Data<-list()

  for (j in 1:length(All_Combinations))
  {
    Beta_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     rbind(beta.mMSE_Single[[j]],beta.mVc_Single[[j]],
                                           beta.mMSE_MR[[j]],beta.mVc_MR[[j]]))

    Utility_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                        rbind(Utility_mMSE_Single[[j]],Utility_mVc_Single[[j]],
                                              Utility_mMSE_MR[[j]],Utility_mVc_MR[[j]]))

    names(Sample.mVc_Single[[j]])<-names(Sample.mMSE_Single[[j]])<-
      names(Sample.mVc_MR[[j]])<-names(Sample.mMSE_MR[[j]])<-c(r1,r2)
  }

  names(Beta_Data)<-names(Utility_Data)<-paste0("Model_",1:length(All_Combinations))

  names(Sample.mVc_Single)<-names(Sample.mVc_MR)<-paste0("Model_",1:length(All_Combinations))
  names(Sample.mMSE_Single)<-names(Sample.mMSE_MR)<-paste0("Model_",1:length(All_Combinations))

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Utility_Estimates"=Utility_Data,
            "Sample_A-Optimality"=Sample.mMSE_Single,
            "Sample_A-Optimality_MR"=Sample.mMSE_MR,
            "Sample_L-Optimality"=Sample.mVc_Single,
            "Sample_L-Optimality_MR"=Sample.mVc_MR,
            "Subsampling_Probability"=Full_SP)

  class(ans)<-c("ModelRobust","logistic")
  return(ans)
}

.getMLE <- function(x, y, w) {
  beta <- rep(0, ncol(x))
  loop  <- 1
  Loop  <- 200
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) },
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop)
}
