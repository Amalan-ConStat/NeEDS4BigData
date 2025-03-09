#' Model robust optimal subsampling for A- and L- optimality criteria under linear regression
#'
#' Using this function sample from big data under linear regression when there are more than
#' one model to describe the data. Subsampling probabilities are obtained based on the A- and L-
#' optimality criteria.
#'
#' @usage
#' modelRobustLinSub(r1,r2,Y,X,N,Apriori_probs,All_Combinations,All_Covariates)
#'
#' @param r1      sample size for initial random sampling
#' @param r2      sample size for optimal sampling
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#' @param Apriori_probs   vector of a priori model probabilities that are used to obtain the model robust subsampling probabilities
#' @param All_Combinations list of possible models that can describe the data
#' @param All_Covariates all the covariates in the models
#'
#' @details
#' Two stage subsampling algorithm for big data under linear regression for multiple models that can
#' describe the big data.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters for all models.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-optimality criteria and
#' model averaging A-, L-optimality subsampling methods.
#'
#' Through the estimated subsampling probabilities a sample of size \eqn{r_2 \ge r_1} is obtained.
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
#' If \eqn{0 < \alpha_{q} < 1} for the a priori model probabilities are not satisfied an error message will be produced,
#' where \eqn{q=1,\ldots,Q} and \eqn{Q} is the number of models in the model set.
#'
#' @return
#' The output of \code{modelRobustLinSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters for each model in a list after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon for each model after subsampling
#'
#' \code{Utility_Estimates} estimated D-(log scaled), A- and L- optimality values for the obtained subsamples
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_A-Optimality_MR} list of indexes for the initial and model robust optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
#'
#' \code{Sample_L-Optimality_MR} list of indexes for the initial and model robust optimal samples obtained based on L-Optimality criteria
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A- and L- optimality criteria
#'
#' @references
#' \insertRef{mahendran2023model}{NeEDS4BigData}
#'
#' @examples
#' indexes<-1:ceiling(nrow(Electric_consumption)*0.005)
#' Original_Data<-cbind(Electric_consumption[indexes,1],1,
#'                      Electric_consumption[indexes,-1])
#' colnames(Original_Data)<-c("Y",paste0("X",0:ncol(Original_Data[,-c(1,2)])))
#' for (j in 3:5) {
#'   Original_Data[,j]<-scale(Original_Data[,j])
#' }
#'
#' No_of_Variables<-ncol(Original_Data[,-c(1,2)])
#' Squared_Terms<-paste0("X",1:No_of_Variables,"^2")
#' term_no <- 2
#' All_Models <- list(c("X0",paste0("X",1:No_of_Variables)))
#'
#' Original_Data<-cbind(Original_Data,Original_Data[,-c(1,2)]^2)
#' colnames(Original_Data)<-c("Y","X0",paste0("X",1:No_of_Variables),
#'                            paste0("X",1:No_of_Variables,"^2"))
#'
#' for (i in 1:No_of_Variables){
#'   x <- as.vector(combn(Squared_Terms,i,simplify = FALSE))
#'     for(j in 1:length(x)){
#'        All_Models[[term_no]] <- c("X0",paste0("X",1:No_of_Variables),x[[j]])
#'        term_no <- term_no+1
#'      }
#'    }
#'
#' All_Models<-All_Models[-c(5:7)]
#' names(All_Models)<-paste0("Model_",1:length(All_Models))
#'
#' r1<-300; r2<-rep(100*c(6,12),25);
#'
#' modelRobustLinSub(r1 = r1, r2 = r2, Y = as.matrix(Original_Data[,1]),
#'                   X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                   Apriori_probs = rep(1/length(All_Models),length(All_Models)),
#'                   All_Combinations = All_Models,
#'                   All_Covariates = colnames(Original_Data)[-1])->Results
#'
#' Beta_Plots<-plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @importFrom matrixStats rowSums2
#' @export
modelRobustLinSub <- function(r1,r2,Y,X,N,Apriori_probs,All_Combinations,All_Covariates){
  if(any(is.na(c(r1,r2,N,Apriori_probs,All_Covariates))) | any(is.nan(c(r1,r2,N,Apriori_probs,All_Covariates)))){
    stop("NA or Infinite or NAN values in the r1,r2,N,Apriori_probs or All_Covariates")
  }

  if((length(r1) + length(N)) != 2){
    stop("r1 or N has a value greater than length one")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r1) > r2)){
    stop("2*r1 cannot be greater than r2 at any point")
  }

  if(length(Apriori_probs) != length(All_Combinations)){
    stop("No of models for averaging is not equal to the a priori probabilities")
  }

  if(any(Apriori_probs > 1) | any(Apriori_probs < 0) | sum(Apriori_probs) != 1){
    stop("A priori probabilities not inbetween zero and one or the sum is one")
  }

  PI.prop <- rep(1/N, N)
  idx.prop <- sample(1:N, size = r1, replace = TRUE)

  x.prop<-lapply(1:length(All_Combinations),function(j){
    X[idx.prop,All_Covariates %in% All_Combinations[[j]]]
  })

  y.prop <- Y[idx.prop,]

  pinv.prop <- N
  pinv.prop <- 1/PI.prop[idx.prop]

  fit.prop <- lapply(1:length(All_Combinations), function(j){
    beta.prop<-solve(a=crossprod(x.prop[[j]]),b=crossprod(x.prop[[j]],y.prop))
    Xbeta_Final<-X[,All_Covariates %in% All_Combinations[[j]]]%*%beta.prop
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
  Utility.mVc_Single<-Utility.mMSE_Single<-list()
  Var_Epsilon_mVc<-Var_Epsilon_mMSE<-matrix(nrow = length(r2),ncol = length(All_Combinations) + 1)
  Sample.mMSE_Single<-Sample.mVc_Single<-list()

  # Model Robust Results
  beta.mVc_MR<-beta.mMSE_MR<-list()
  Utility.mVc_MR<-Utility.mMSE_MR<-list()
  Var_Epsilon_mVc_MR<-Var_Epsilon_mMSE_MR<-matrix(nrow = length(r2),ncol = length(All_Combinations) + 1)
  Sample.mMSE_MR<-Sample.mVc_MR<-list()

  Var_Epsilon_mVc[,1]<-Var_Epsilon_mMSE[,1]<-Var_Epsilon_mVc_MR[,1]<-Var_Epsilon_mMSE_MR[,1]<-r2

  # For the models, Single and Model Robust
  for (a in 1:length(All_Combinations))
  {
    beta.mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Single Model Results
    beta.mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Utility.mVc_Single[[a]]<-matrix(nrow = length(r2),ncol = 4) # Single Model Results
    Utility.mMSE_Single[[a]]<-matrix(nrow = length(r2),ncol = 4)
    Sample.mMSE_Single[[a]]<-Sample.mVc_Single[[a]]<-list()

    beta.mVc_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 ) # Model Robust Results
    beta.mMSE_MR[[a]]<-matrix(nrow = length(r2),ncol = length(All_Combinations[[a]])+1 )
    Utility.mVc_MR[[a]]<-matrix(nrow = length(r2),ncol = 4) # Model Robust Results
    Utility.mMSE_MR[[a]]<-matrix(nrow = length(r2),ncol = 4)
    Sample.mMSE_MR[[a]]<-Sample.mVc_MR[[a]]<-list()

    Sample.mMSE_Single[[a]][[1]]<-Sample.mVc_Single[[a]][[1]]<-
      Sample.mMSE_MR[[a]][[1]]<-Sample.mVc_MR[[a]][[1]]<-idx.prop

    if(all(x.prop[[a]][,1] == 1)){
      colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
        colnames(beta.mMSE_MR[[a]])<-c("r2",paste0("Beta_",0:(length(All_Combinations[[a]])-1)))
    } else {
      colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
        colnames(beta.mMSE_MR[[a]])<-c("r2",paste0("Beta_",1:(length(All_Combinations[[a]]))))
    }

    colnames(Utility.mVc_Single[[a]])<-colnames(Utility.mMSE_Single[[a]])<-
      colnames(Utility.mVc_MR[[a]])<-
      colnames(Utility.mMSE_MR[[a]])<-c("r2","D-optimality","A-optimality","L-optimality")

  }

  ## mVc
  PI_Single.mVc <- lapply(1:length(All_Combinations), function(j){ # Single Model Results
    PI.mVc<-sqrt(fit.prop[[j]]$Epsilon.prop^2 * matrixStats::rowSums2(X[,All_Covariates %in% All_Combinations[[j]] ]^2))
    return(PI.mVc/sum(PI.mVc))
  })
  PI_MR.mVc<-matrixStats::rowSums2(do.call(cbind,PI_Single.mVc)%*%diag(Apriori_probs)) # Model Robust Results

  ## mMSE
  # For efficient row-wise operations
  PI_Single.mMSE <- lapply(seq_along(All_Combinations), function(j) {
    X_sub <- X[, All_Covariates %in% All_Combinations[[j]], drop = FALSE]
    XtX_inv <- solve(crossprod(X_sub))
    row_sums_squared <- matrixStats::rowSums2((X_sub %*% XtX_inv)^2)
    PI.mMSE <- sqrt(fit.prop[[j]]$Epsilon.prop^2 * row_sums_squared)
    PI.mMSE / sum(PI.mMSE)
  })
  PI_MR.mMSE<-matrixStats::rowSums2(do.call(cbind,PI_Single.mMSE)%*%diag(Apriori_probs))  # Model Robust Results

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    ## mVc
    idx_Single.mVc <- lapply(1:length(All_Combinations), function(j){
      sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI_Single.mVc[[j]]) # Single Model Results
    })
    idx_MR.mVc <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI_MR.mVc) # Model Robust Results

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
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    fit_MR.mVc <- lapply(1:length(All_Combinations),function(j){
      pi4_r<-sqrt(r2[i]*pinv_MR.mVc^(-1))
      X_r4<-x_MR.mVc[[j]]/pi4_r
      Y_r4<-y_MR.mVc/pi4_r
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop
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

      Temp<-Var_Epsilon_mVc[i,j+1]*crossprod(x_Single.mVc[[j]])
      Temp_Inv<-solve(Temp)
      x.mVc_t<-t(x_Single.mVc[[j]])
      Temp_Int<-Temp_Inv%*%x.mVc_t
      Temp1<-x_Single.mVc[[j]]%*%Temp_Int
      Utility.mVc_Single[[j]][i,]<-c(r2[i],log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      Temp<-Var_Epsilon_mVc_MR[i,j+1]*crossprod(x_MR.mVc[[j]])
      Temp_Inv<-solve(Temp)
      x.mVc_t<-t(x_MR.mVc[[j]])
      Temp_Int<-Temp_Inv%*%x.mVc_t
      Temp1<-x_MR.mVc[[j]]%*%Temp_Int
      Utility.mVc_MR[[j]][i,]<-c(r2[i],log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      if(anyNA(fit_Single.mVc[[j]]$beta.prop) || anyNA(fit_MR.mVc[[j]]$beta.prop)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    ## mMSE
    idx_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      sample(1:N, size = r2[i]-r1, replace = T, prob = PI_Single.mMSE[[j]]) # Single Model Results
    })
    idx_MR.mMSE <- sample(1:N, size = r2[i]-r1, replace = T, prob = PI_MR.mMSE) # Model Robust Results

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
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop
      Var.prop<-sum((Y-Xbeta_Final)^2)/N
      return(list("beta.prop"=beta.prop,"Var.prop"=Var.prop))
    })

    fit_MR.mMSE <- lapply(1:length(All_Combinations),function(j){
      pi4_r<-sqrt(r2[i]*pinv_MR.mMSE^(-1))
      X_r4<-x_MR.mMSE[[j]]/pi4_r
      Y_r4<-y_MR.mMSE/pi4_r
      beta.prop<-solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4))
      Xbeta_Final<-X[,All_Covariates %in% All_Combinations[[j]] ]%*%beta.prop
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

      Temp<-Var_Epsilon_mMSE[i,j+1]*crossprod(x_Single.mMSE[[j]])
      Temp_Inv<-solve(Temp)
      x.mMSE_t<-t(x_Single.mMSE[[j]])
      Temp_Int<-Temp_Inv%*%x.mMSE_t
      Temp1<-x_Single.mMSE[[j]]%*%Temp_Int
      Utility.mMSE_Single[[j]][i,]<-c(r2[i],log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

      Temp<-Var_Epsilon_mMSE_MR[i,j+1]*crossprod(x_MR.mMSE[[j]])
      Temp_Inv<-solve(Temp)
      x.mMSE_t<-t(x_MR.mMSE[[j]])
      Temp_Int<-Temp_Inv%*%x.mMSE_t
      Temp1<-x_MR.mMSE[[j]]%*%Temp_Int
      Utility.mMSE_MR[[j]][i,]<-c(r2[i],log(det(Temp)),psych::tr(Temp_Inv),psych::tr(Temp1))

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
  Beta_Data<-Utility_Data<-list()

  for (j in 1:length(All_Combinations))
  {
    Beta_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     rbind(beta.mMSE_Single[[j]],beta.mVc_Single[[j]],
                                           beta.mMSE_MR[[j]],beta.mVc_MR[[j]]))

    Utility_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                        rbind(Utility.mMSE_Single[[j]],Utility.mVc_Single[[j]],
                                              Utility.mMSE_MR[[j]],Utility.mVc_MR[[j]]))

    names(Sample.mVc_Single[[j]])<-names(Sample.mMSE_Single[[j]])<-
      names(Sample.mVc_MR[[j]])<-names(Sample.mMSE_MR[[j]])<-c(r1,r2)
  }

  Var_Epsilon_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                     rbind(Var_Epsilon_mMSE,Var_Epsilon_mVc,
                                           Var_Epsilon_mMSE_MR,Var_Epsilon_mVc_MR))

  names(Beta_Data)<-names(Utility_Data)<-paste0("Model_",1:length(All_Combinations))
  colnames(Var_Epsilon_Data)[-1]<-c("r2",paste0("Model_",1:length(All_Combinations)))

  names(Sample.mVc_Single)<-names(Sample.mVc_MR)<-paste0("Model_",1:length(All_Combinations))
  names(Sample.mMSE_Single)<-names(Sample.mMSE_MR)<-paste0("Model_",1:length(All_Combinations))

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Utility_Estimates"=Utility_Data,
            "Variance_Epsilon_Estimates"=Var_Epsilon_Data,
            "Sample_A-Optimality"=Sample.mMSE_Single,
            "Sample_A-Optimality_MR"=Sample.mMSE_MR,
            "Sample_L-Optimality"=Sample.mVc_Single,
            "Sample_L-Optimality_MR"=Sample.mVc_MR,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelRobust","linear")
  return(ans)
}
