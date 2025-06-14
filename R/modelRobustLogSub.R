#' Model robust optimal subsampling for A- and L- optimality criteria under logistic regression
#'
#' Using this function sample from big data under logistic regression when there are more than
#' one model to describe the data. Subsampling probabilities are obtained based on the A- and L-
#' optimality criteria.
#'
#' @usage
#' modelRobustLogSub(r0,rf,Y,X,N,Apriori_probs,All_Combinations,All_Covariates)
#'
#' @param r0      sample size for initial random sample
#' @param rf      final sample size including initial(r0) and optimal(r) samples
#' @param Y       response data or Y
#' @param X       covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N       size of the big data
#' @param Apriori_probs    vector of a priori model probabilities that are used to obtain the model robust subsampling probabilities
#' @param All_Combinations list of possible models that can describe the data
#' @param All_Covariates   all the covariates in the models
#'
#' @details
#' Two stage subsampling algorithm for big data under logistic regression for multiple models that can
#' describe the big data.
#'
#' First stage is to obtain a random sample of size \eqn{r_0} and estimate the model parameters for all models.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-optimality criteria and
#' model averaging A-, L-optimality subsampling methods.
#'
#' Through the estimated subsampling probabilities a sample of size \eqn{r \ge r_0} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for all the models.
#'
#' \strong{NOTE} :  If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r \ge r_0} is not satisfied then an error message will be produced.
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
#' \code{Beta_Data} estimated model parameters for each model in a list after subsampling
#'
#' \code{Sample_L-optimality} list of indexes for the initial and optimal samples obtained based on L-optimality criteria
#'
#' \code{Sample_L-optimality_MR} list of indexes for the initial and model robust optimal samples obtained based on L-optimality criteria
#'
#' \code{Sample_A-optimality} list of indexes for the initial and optimal samples obtained based on A-optimality criteria
#'
#' \code{Sample_A-optimality_MR} list of indexes for the initial and model robust optimal samples obtained based on A-optimality criteria
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities for A- and L- optimality criteria
#'
#' @references
#' \insertRef{mahendran2023model}{NeEDS4BigData}
#'
#' @examples
#' indexes<-1:ceiling(nrow(Skin_segmentation)*0.25)
#' Original_Data<-cbind(Skin_segmentation[indexes,1],1,Skin_segmentation[indexes,-1])
#' colnames(Original_Data)<-c("Y",paste0("X",0:ncol(Original_Data[,-c(1,2)])))
#'
#' # Scaling the covariate data
#' for (j in 3:5) {
#'   Original_Data[,j]<-scale(Original_Data[,j])
#' }
#' No_of_Variables<-ncol(Original_Data[,-c(1,2)])
#'
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
#'   for(j in 1:length(x)){
#'     All_Models[[term_no]] <- c("X0",paste0("X",1:No_of_Variables),x[[j]])
#'     term_no <- term_no+1
#'   }
#' }
#' All_Models<-All_Models[-c(5:7)]
#' names(All_Models)<-paste0("Model_",1:length(All_Models))
#'
#' r0<-300; rf<-rep(100*c(6,9),50);
#'
#' modelRobustLogSub(r0 = r0, rf = rf, Y = as.matrix(Original_Data[,1]),
#'                   X = as.matrix(Original_Data[,-1]),N = nrow(Original_Data),
#'                   Apriori_probs = rep(1/length(All_Models),length(All_Models)),
#'                   All_Combinations = All_Models,
#'                   All_Covariates = colnames(Original_Data)[-1])->Results
#'
#' Beta_Plots<-plot_Beta(Results)
#'
#' @importFrom Rdpack reprompt
#' @importFrom psych tr
#' @importFrom matrixStats rowSums2
#' @export
modelRobustLogSub <- function(r0,rf,Y,X,N,Apriori_probs,All_Combinations,All_Covariates){
  if(any(is.na(c(r0,rf,N,Apriori_probs,All_Covariates))) | any(is.nan(c(r0,rf,N,Apriori_probs,All_Covariates)))){
    stop("NA or Infinite or NAN values in the r0,rf,N,Apriori_probs or All_Covariates")
  }

  if((length(r0) + length(N)) != 2){
    stop("r0 or N has a value greater than length one")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r0) > rf)){
    stop("2*r0 cannot be greater than rf at any point")
  }

  if(length(Apriori_probs) != length(All_Combinations)){
    stop("No of models for averaging is not equal to the a priori probabilities")
  }

  if(any(Apriori_probs > 1) | any(Apriori_probs < 0) | sum(Apriori_probs) != 1){
    stop("A priori probabilities not inbetween zero and one or the sum is one")
  }

  n1 <- sum(Y)
  n0 <- N - n1
  PI.prop <- rep(1/(2*n0), N)
  PI.prop[Y==1] <- 1/(2*n1)
  idx.prop <- sample(1:N, size = r0, replace = TRUE, prob = PI.prop)

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
  beta.mVc_Single<-beta.mMSE_Single<-list()
  Sample.mMSE_Single<-Sample.mVc_Single<-list()

  # Model Robust Results
  beta.mVc_MR<-beta.mMSE_MR<-list()
  Sample.mMSE_MR<-Sample.mVc_MR<-list()

  # For the models, Single and Model Robust
  for (a in 1:length(All_Combinations))
  {
    beta.mVc_Single[[a]]<-matrix(nrow = length(rf),ncol = length(All_Combinations[[a]])+1 ) # Single Model Results
    beta.mMSE_Single[[a]]<-matrix(nrow = length(rf),ncol = length(All_Combinations[[a]])+1 )
    Sample.mMSE_Single[[a]]<-Sample.mVc_Single[[a]]<-list()

    beta.mVc_MR[[a]]<-matrix(nrow = length(rf),ncol = length(All_Combinations[[a]])+1 ) # Model Robust Results
    beta.mMSE_MR[[a]]<-matrix(nrow = length(rf),ncol = length(All_Combinations[[a]])+1 )
    Sample.mMSE_MR[[a]]<-Sample.mVc_MR[[a]]<-list()

    Sample.mMSE_Single[[a]][[1]]<-Sample.mVc_Single[[a]][[1]]<-
      Sample.mMSE_MR[[a]][[1]]<-Sample.mVc_MR[[a]][[1]]<-idx.prop

    if(all(x.prop[[a]][,1] == 1)){
      colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
        colnames(beta.mMSE_MR[[a]])<-c("rf",paste0("Beta_",0:(length(All_Combinations[[a]])-1)))
    } else {
      colnames(beta.mVc_Single[[a]])<-colnames(beta.mMSE_Single[[a]])<-colnames(beta.mVc_MR[[a]])<-
        colnames(beta.mMSE_MR[[a]])<-c("rf",paste0("Beta_",1:(length(All_Combinations[[a]]))))
    }
  }

  ## mVc
  PI_Single.mVc <- lapply(1:length(All_Combinations), function(j){ # Single Model Results
    PI.mVc<-sqrt((Y - P.prop[[j]])^2 * matrixStats::rowSums2(X[,All_Covariates %in% All_Combinations[[j]] ]^2))
    return(PI.mVc/sum(PI.mVc))
  })
  PI_MR.mVc<-matrixStats::rowSums2(do.call(cbind,PI_Single.mVc)%*%diag(Apriori_probs)) # Model Robust Results

  ## mMSE
  p_Single.prop <- lapply(1:length(All_Combinations),function(j){
    P.prop[[j]][idx.prop] # Single Model Results
  })
  w_Single.prop <- lapply(1:length(All_Combinations),function(j){
    p_Single.prop[[j]] * (1 - p_Single.prop[[j]]) # Single Model Results
  })
  W_Single.prop <- lapply(1:length(All_Combinations),function(j){
    solve(crossprod(x.prop[[j]],x.prop[[j]] * w_Single.prop[[j]] * pinv.prop)) # Single Model Results
  })
  PI_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
    PI.mMSE<-sqrt((Y - P.prop[[j]])^2 * matrixStats::rowSums2((X[,All_Covariates %in% All_Combinations[[j]] ]%*%
                                                                 W_Single.prop[[j]])^2)) # Single Model Results
    return(PI.mMSE/sum(PI.mMSE))
  })
  PI_MR.mMSE<-matrixStats::rowSums2(do.call(cbind,PI_Single.mMSE)%*%diag(Apriori_probs))  # Model Robust Results

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(rf))
  {
    ## mVc
    idx_Single.mVc <- lapply(1:length(All_Combinations), function(j){
      sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI_Single.mVc[[j]]) # Single Model Results
    })
    idx_MR.mVc <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI_MR.mVc) # Model Robust Results

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

      beta.mVc_Single[[j]][i,] <- c(rf[i],fit_Single.mVc[[j]]$par)
      beta.mVc_MR[[j]][i,] <- c(rf[i],fit_MR.mVc[[j]]$par)

      if(anyNA(fit_Single.mVc[[j]]$par) || anyNA(fit_MR.mVc[[j]]$par)){
        stop("There are NA or NaN values in the model parameters")
      }
    }

    ## mMSE
    idx_Single.mMSE <- lapply(1:length(All_Combinations),function(j){
      sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI_Single.mMSE[[j]]) # Single Model Results
    })
    idx_MR.mMSE <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI_MR.mMSE) # Model Robust Results

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

      beta.mMSE_Single[[j]][i,] <- c(rf[i],fit_Single.mMSE[[j]]$par) # Single Model Results
      beta.mMSE_MR[[j]][i,] <- c(rf[i],fit_MR.mMSE[[j]]$par) # Model Robust Results

      if(anyNA(fit_Single.mMSE[[j]]$par) || anyNA(fit_MR.mMSE[[j]]$par)){
        stop("There are NA or NaN values in the model parameters")
      }
    }
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
    Beta_Data[[j]]<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(rf)),
                                     rbind(beta.mMSE_Single[[j]],beta.mVc_Single[[j]],
                                           beta.mMSE_MR[[j]],beta.mVc_MR[[j]]))

    names(Sample.mVc_Single[[j]])<-names(Sample.mMSE_Single[[j]])<-
      names(Sample.mVc_MR[[j]])<-names(Sample.mMSE_MR[[j]])<-c(r0,rf)
  }

  names(Beta_Data)<-paste0("Model_",1:length(All_Combinations))

  names(Sample.mVc_Single)<-names(Sample.mVc_MR)<-paste0("Model_",1:length(All_Combinations))
  names(Sample.mMSE_Single)<-names(Sample.mMSE_MR)<-paste0("Model_",1:length(All_Combinations))

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
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
    H <- crossprod(x,(pr * (1 - pr) * w * x))
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

.stratified_sample_indices <- function(response, n_total) {
  ind_0 <- which(response == 0)
  ind_1 <- which(response == 1)
  p_0 <- length(ind_0) / length(response)
  p_1 <- 1 - p_0  # or length(ind_1) / length(response)

  n_0 <- round(n_total * p_0)
  n_1 <- n_total - n_0  # to ensure n_total

  # Sampled indices
  samp_ind_0 <- sample(ind_0, size = n_0)
  samp_ind_1 <- sample(ind_1, size = n_1)

  # Return combined indices
  sort(c(samp_ind_0, samp_ind_1))
}
