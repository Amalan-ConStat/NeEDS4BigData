#' Subsampling under Poisson regression for a potentially misspecified model
#'
#' Using this function sample from big data under Poisson regression for a potentially misspecified model.
#' Subsampling probabilities are obtained based on the A- and L- optimality criteria
#' with the RLmAMSE (Reduction of Loss by minimizing the Average Mean Squared Error).
#'
#' @usage
#' modelMissPoiSub(r1,r2,Y,X,N,Alpha,proportion,model="Auto")
#'
#' @param r1          sample size for initial random sampling
#' @param r2          sample size for optimal sampling
#' @param Y           response data or Y
#' @param X           covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N           size of the big data
#' @param Alpha       scaling factor when using Log Odds or Power functions to magnify the probabilities
#' @param proportion  a proportion of the big data is used to help estimate AMSE values from the subsamples
#' @param model       formula for the model used in the GAM or the main effects ("s()"), squared term ("I()") or two-way interaction ("lo()") model
#'
#' @details
#' \strong{The article for this function is in preparation for publication. Please be patient.}
#'
#' Two stage subsampling algorithm for big data under Poisson regression for potential model misspecification.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-optimality criteria,
#' RLmAMSE and enhanced RLmAMSE (log-odds and power) subsampling methods.
#'
#' Through the estimated subsampling probabilities a sample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for A- and L-optimality,
#' while for RLmAMSE and enhanced RLmAMSE (log-odds and power) only the optimal sample is used.
#'
#' \strong{NOTE} :  If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r_2 \ge r_1} is not satisfied then an error message will be produced.
#'
#' If the big data \eqn{X,Y} has any missing values then an error message will be produced.
#'
#' The big data size \eqn{N} is compared with the sizes of \eqn{X,Y},F_estimate_Full and
#' if they are not aligned an error message will be produced.
#'
#' If \eqn{\alpha > 1} for the scaling vector is not satisfied an error message will be produced.
#'
#' If proportion is not in the region of \eqn{(0,1]} an error message will be produced.
#'
#' \code{model} is a formula input formed based on the covariates through the spline terms (s()),
#' squared term (I()), interaction terms (lo()) or automatically. If \code{model} is empty or NA
#' or NAN or not one of the defined inputs an error message is printed, as a default we have set
#' \code{model="Auto"}.
#'
#' @return
#' The output of \code{modelMissPoiSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{AMSE_Estimates} matrix of estimated AMSE values after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
#'
#' \code{Sample_RLmAMSE} list of indexes for the optimal samples obtained based on RLmAMSE
#'
#' \code{Sample_RLmAMSE_Log_Odds} list of indexes for the optimal samples obtained based on RLmAMSE with Log Odds function
#'
#' \code{Sample_RLmAMSE_Power} list of indexes for the optimal samples obtained based on RLmAMSE with Power function
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities
#'
#' @references
#' \insertRef{adewale2009robust}{NeEDS4BigData}
#' \insertRef{adewale2010robust}{NeEDS4BigData}
#'
#' @examples
#' Beta<-c(-1,0.75,0.75,1); family <- "poisson"; N<-10000
#' X_1 <- replicate(2,stats::runif(n=N,min = -1,max = 1))
#'
#' Temp<-Rfast::rowprods(X_1)
#' Misspecification <- (Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2)
#' X_Data <- cbind(X0=1,X_1);
#'
#' Full_Data<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family)
#'
#' r1<-300; r2<-rep(100*c(6,9),50);
#' Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#' \dontrun{
#' Results<-modelMissPoiSub(r1 = r1, r2 = r2,
#'                          Y = as.matrix(Original_Data[,1]),
#'                          X = as.matrix(Original_Data[,-1]),
#'                          N = N, Alpha = 10, proportion = 0.3)
#'
#' # parallel::stopCluster(cl)
#'
#' plot_Beta(Results)
#' plot_AMSE(Results)
#' }
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @import foreach
#' @importFrom gam s
#' @importFrom Rfast rowprods
#' @importFrom utils combn
#' @importFrom psych tr
#' @importFrom rlang is_formula
#' @export
modelMissPoiSub <- function(r1,r2,Y,X,N,Alpha,proportion,model="Auto"){
  if(any(is.na(c(r1,r2,N,Alpha,proportion))) | any(is.nan(c(r1,r2,N,Alpha,proportion)))){
    stop("NA or Infinite or NAN values in the r1,r2,N,Alpha or proportion")
  }

  if((length(r1) + length(N) + length(proportion)) != 3){
    stop("proportion, r1 or N has a value greater than length one")
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

  if(Alpha <= 1 | length(Alpha) > 1){
    stop("Scaling factor alpha is not greater than one or the length is more than one")
  }

  if(proportion >1 | proportion <=0){
    stop("Proportion should be a value higher than zero and less than or equal one")
  }
  if(proportion >= 0.5){
    message("50% or >=50% of the big data is used to help find AMSE for the subsamples, \nthis could take some time.")
  }

  if(anyNA(model) | any(is.nan(model)) | is.null(model) ){

    stop("The model formula for GAM is NA or NAN or NULL")

  } else if(model == "s()"){

    if(all(X[,1] == 1)){
      main_effects <- paste0("s(", colnames(X[,-1]),")")
      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
    } else {
      main_effects <- paste0("s(", colnames(X),")")
      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
    }

  } else if(model == "lo()"){

    if(all(X[,1] == 1)){
      main_effects <- paste0("s(", colnames(X[,-1]),")")
      two_way_interactions <- utils::combn(colnames(X[,-1]), 2,
                                           function(cols){paste0("lo(", paste(cols, collapse = "*"), ")")})

      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                          paste(two_way_interactions,collapse = " + ")))
    } else {
      main_effects <- paste0("s(",colnames(X),")")
      two_way_interactions <- utils::combn(colnames(X), 2,
                                           function(cols){paste0("lo(", paste(cols, collapse = "*"), ")")})

      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                          paste(two_way_interactions,collapse = " + ")))
    }

  } else if(model == "I()"){

    if(all(X[,1] == 1)){
      main_effects <- paste0("s(",colnames(X[,-1]),")")
      squared_terms <- paste0("s(I(",colnames(X[, -1]),"^2))")

      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                          paste(squared_terms,collapse = " + ")))
    } else {
      main_effects <- paste0("s(", colnames(X),")")
      squared_terms <- paste0("s(I(",colnames(X),"^2))")

      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                          paste(squared_terms,collapse = " + ")))
    }

  } else if(model == "Auto"){

    if(all(X[,1] == 1)){
      main_effects <- paste0("s(",colnames(X[, -1]),")")

      if(ncol(X[,-1]) == 2 ){
        two_way_interactions <- utils::combn(colnames(X[, -1]), 2,
                                             function(cols){paste0("lo(", paste(cols, collapse = "*"), ")")})

        my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                            paste(two_way_interactions,collapse = " + ")))
      } else{
        my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
      }
    } else{
      main_effects <- paste0("s(",colnames(X),")")

      if(ncol(X) == 2 ){
        two_way_interactions <- utils::combn(colnames(X), 2,
                                             function(cols){paste0("lo(", paste(cols, collapse = "*"), ")")})

        my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + "),"+",
                                            paste(two_way_interactions,collapse = " + ")))
      } else{
        my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
      }
    }

  } else if(rlang::is_formula(model)){
    my_formula <- model
  } else {
    stop("The model for GAM is not a valid formula or input")
  }

  PI.prop <- rep(1/N,N)
  idx.prop <- sample(1:N, size = r1, replace = TRUE, prob = PI.prop)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]

  pinv.prop <- rep(N,r1)

  fit.prop <- stats::glm(y.prop~x.prop-1,family="quasipoisson")
  beta.prop <- fit.prop$coefficients

  if (anyNA(beta.prop)){
    stop("There are NA or NaN values in the model parameters")
  }
  Xbeta_Final <- X %*% beta.prop
  Lambda.prop  <- exp(Xbeta_Final)

  #calculate f_hat
  Assumed_Data<-data.frame(Y=y.prop,x.prop)
  fit_GAM<-gam::gam(my_formula,data=Assumed_Data,family = "quasipoisson")
  Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X))
  Lambda_GAM<-exp(Xbeta_GAM)
  f_estimate<-Xbeta_GAM - Xbeta_Final

  if(proportion*N != r1){
    idx.proportion <- sample(1:N, size = ceiling(proportion*N), replace = TRUE, prob = PI.prop)

    Y_proportion <- Y[idx.proportion]
    X_proportion <- X[idx.proportion,]
    pinv.proportion <- rep(N,ceiling(proportion*N))

    fit.proportion <- stats::glm(Y_proportion~X_proportion-1,family="quasipoisson")
    beta.proportion <- fit.proportion$coefficients
    Xbeta_proportion <- X %*% beta.proportion

    Assumed_Data<-data.frame(Y=Y_proportion,X_proportion)
    fit_GAM_proportion<-gam::gam(my_formula,data=Assumed_Data,family = "quasipoisson")
    Xbeta_GAM_proportion<-gam::predict.Gam(fit_GAM_proportion,newdata = data.frame(X))

    F_Estimate_Full<-Xbeta_GAM_proportion - Xbeta_proportion
    Beta_Estimate_Full<-beta.proportion
  } else {
    Beta_Estimate_Full<- beta.prop ; F_Estimate_Full<-f_estimate
  }

  ## mVC
  PI.mVc <- sqrt((Y - Lambda.prop)^2 * rowSums(X^2))
  PI.mVc <- PI.mVc / sum(PI.mVc)

  # mMSE
  lambda.prop <- Lambda.prop[idx.prop]
  w.prop <- lambda.prop
  W.prop <- solve(crossprod(x.prop, x.prop * w.prop * pinv.prop))

  PI.mMSE <- sqrt((Y - Lambda.prop)^2 * rowSums((X%*%W.prop)^2))
  PI.mMSE <- PI.mMSE / sum(PI.mMSE)

  a<-NULL
  L_All<-foreach::foreach(a=1:N,.combine = rbind,.packages = "psych") %dopar% {
    X_r1<-X[c(idx.prop,a),]
    lambda_r1<-exp(X_r1%*%beta.prop)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(crossprod(X_r1,X_r1 * W_r1))
    Temp1<-(W_r1*X_r1)%*%H_r1

    f_r1<-f_estimate[c(idx.prop,a)]
    lambda_Tr1<-exp((X_r1 %*% beta.prop) + f_r1)
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <- crossprod(X_r1, X_r1 * W_Tr1)
    b_r1 <- crossprod(X_r1, lambda_Tr1-lambda_r1)

    Temp1_H <- Temp1 %*% H_Tr1
    diag_Temp <- rowSums(Temp1_H * Temp1)
    L1_r1 <- sum(diag_Temp)

    XH_r1 <- X_r1 %*% H_r1
    XH_b_r1 <- XH_r1 %*% b_r1
    diff <- XH_b_r1 - f_r1
    L2_r1 <- sum((W_r1 * diff)^2)

    L1_r1+L2_r1
  }

  lambda_r1<-exp(X[idx.prop,]%*%beta.prop)
  W_r1<-as.vector(lambda_r1)
  H_r1 <-solve(crossprod(X[idx.prop,], X[idx.prop,]*W_r1))
  Temp1<-(W_r1*X[idx.prop,])%*%H_r1

  lambda_Tr1<-exp((X[idx.prop,] %*% beta.prop) + f_estimate[idx.prop])
  W_Tr1<-as.vector(lambda_Tr1)
  H_Tr1 <- crossprod(X[idx.prop,], X[idx.prop,]*W_Tr1)
  b_r1 <- crossprod(X[idx.prop,], lambda_Tr1-lambda_r1)
  L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
  L2_r1 <- sum((W_r1*((X[idx.prop,]%*%H_r1%*%b_r1)-f_estimate[idx.prop]))^2)

  L_All_Temp<-L1_r1+L2_r1
  L_All_Final<-abs(L_All_Temp-L_All)

  # RLmAMSE
  PI.RLmAMSE <- (max(L_All_Final) - L_All_Final) / sum(max(L_All_Final) - L_All_Final)

  # RLmAMSE LO
  PI.RLmAMSE_LO<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    (1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1)/
      sum((1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1))},PI.RLmAMSE)

  # RLmAMSE Pow
  PI.RLmAMSE_Pow<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    PI.RLmAMSE^Alpha/sum(PI.RLmAMSE^Alpha)},PI.RLmAMSE)

  # For the Model with already available subsampling probabilities - mMSE and mVc
  beta_mVc<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  AMSE_Sample_mVc<-matrix(nrow = length(r2),ncol = 4); Sample.mVc<-list();

  beta_mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  AMSE_Sample_mMSE<-matrix(nrow = length(r2),ncol = 4); Sample.mMSE<-list();

  # For the Model with model robust subsampling probabilities RLmAMSE
  beta_RLmAMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  AMSE_Sample_RLmAMSE<-matrix(nrow = length(r2),ncol = 4); Sample.RLmAMSE<-list();

  # For the Model with model robust subsampling probabilities RLmAMSE LO
  beta_RLmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r2),ncol(X)+1)),simplify = FALSE)
  AMSE_Sample_RLmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r2),4)),simplify = FALSE)
  Sample.RLmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));

  # For the Model with model robust subsampling probabilities RLmAMSE Pow
  beta_RLmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r2),ncol(X)+1)),simplify = FALSE)
  AMSE_Sample_RLmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r2),4)),simplify = FALSE)
  Sample.RLmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));

  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.RLmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.RLmAMSE_LO[[j]][[1]]<-Sample.RLmAMSE_Pow[[j]][[1]]<-idx.prop
  }

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    # mVc
    idx.mVc <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mVc)

    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]

    fit.mVc <- stats::glm(y.mVc~x.mVc-1,family="quasipoisson",
                          weights = c(1 / PI.mVc[idx.mVc], pinv.prop))

    beta_mVc[i,] <- c(r2[i],fit.mVc$coefficients)

    idx.mVc->Sample.mVc[[i+1]]

    lambda_r1<-exp(x.mVc%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(crossprod(x.mVc, x.mVc * W_r1))
    Temp1<-(W_r1 * x.mVc)%*%H_r1

    lambda_Tr1<-exp((x.mVc %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.mVc, idx.prop)])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <- crossprod(x.mVc, x.mVc * W_Tr1)
    b_r1 <- crossprod(x.mVc, lambda_Tr1-lambda_r1)
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.mVc%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.mVc, idx.prop)]))^2)

    AMSE_Sample_mVc[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    # mMSE
    idx.mMSE <- sample(1:N, size = r2[i]-r1, replace = TRUE, prob = PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]

    fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1,family = "quasipoisson",
                           weights =c(1 / PI.mMSE[idx.mMSE], pinv.prop) )

    beta_mMSE[i,] <- c(r2[i],fit.mMSE$coefficients)

    idx.mMSE->Sample.mMSE[[i+1]]

    lambda_r1<-exp(x.mMSE%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(crossprod(x.mMSE, x.mMSE * W_r1))
    Temp1<-(W_r1*x.mMSE)%*%H_r1

    lambda_Tr1<-exp((x.mMSE %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.mMSE, idx.prop)])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <- crossprod(x.mMSE, x.mMSE * W_Tr1)
    b_r1 <- crossprod(x.mMSE, lambda_Tr1-lambda_r1)
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.mMSE%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.mMSE, idx.prop)]))^2)

    AMSE_Sample_mMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    # RLmAMSE
    idx.RLmAMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.RLmAMSE)

    x.RLmAMSE <- X[idx.RLmAMSE,]
    y.RLmAMSE <- Y[idx.RLmAMSE]

    fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                              weights = c(1 / PI.RLmAMSE[idx.RLmAMSE]))

    beta_RLmAMSE[i,] <- c(r2[i],fit.RLmAMSE$coefficients)

    idx.RLmAMSE->Sample.RLmAMSE[[i+1]]

    lambda_r1<-exp(x.RLmAMSE%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(crossprod(x.RLmAMSE, x.RLmAMSE * W_r1))
    Temp1<-(W_r1*x.RLmAMSE)%*%H_r1

    lambda_Tr1<-exp((x.RLmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[idx.RLmAMSE])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr1)
    b_r1 <- crossprod(x.RLmAMSE, lambda_Tr1-lambda_r1)
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.RLmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[idx.RLmAMSE]))^2)

    AMSE_Sample_RLmAMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    for (j in 1:length(Alpha))
    {
      # RLmAMSE Log Odds
      idx.RLmAMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.RLmAMSE_LO[,j])

      x.RLmAMSE <- X[idx.RLmAMSE,]
      y.RLmAMSE <- Y[idx.RLmAMSE]

      fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                                weights = c(1 / PI.RLmAMSE_LO[idx.RLmAMSE,j]))

      beta_RLmAMSE_LO[[j]][i,] <- c(r2[i],fit.RLmAMSE$coefficients)

      idx.RLmAMSE->Sample.RLmAMSE_LO[[j]][[i+1]]

      lambda_r1<-exp(x.RLmAMSE%*%Beta_Estimate_Full)
      W_r1<-as.vector(lambda_r1)
      H_r1 <-solve(crossprod(x.RLmAMSE, x.RLmAMSE * W_r1))
      Temp1<-(W_r1*x.RLmAMSE)%*%H_r1

      lambda_Tr1<-exp((x.RLmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[idx.RLmAMSE])
      W_Tr1<-as.vector(lambda_Tr1)
      H_Tr1 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr1)
      b_r1 <- crossprod(x.RLmAMSE, lambda_Tr1-lambda_r1)
      L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
      L2_r1 <- sum((W_r1*((x.RLmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[idx.RLmAMSE]))^2)

      AMSE_Sample_RLmAMSE_LO[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      # RLmAMSE Power
      idx.RLmAMSE <- sample(1:N, size = r2[i], replace = TRUE, prob = PI.RLmAMSE_Pow[,j])

      x.RLmAMSE <- X[idx.RLmAMSE,]
      y.RLmAMSE <- Y[idx.RLmAMSE]

      fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                                weights = c(1 / PI.RLmAMSE_Pow[idx.RLmAMSE,j]))

      beta_RLmAMSE_Pow[[j]][i,] <- c(r2[i],fit.RLmAMSE$coefficients)

      idx.RLmAMSE->Sample.RLmAMSE_Pow[[j]][[i+1]]

      lambda_r1<-exp(x.RLmAMSE%*%Beta_Estimate_Full)
      W_r1<-as.vector(lambda_r1)
      H_r1 <-solve(t(x.RLmAMSE) %*% (x.RLmAMSE * W_r1))
      Temp1<-(W_r1*x.RLmAMSE)%*%H_r1

      lambda_Tr1<-exp((x.RLmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[idx.RLmAMSE])
      W_Tr1<-as.vector(lambda_Tr1)
      H_Tr1 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr1)
      b_r1 <- crossprod(x.RLmAMSE, lambda_Tr1-lambda_r1)
      L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
      L2_r1 <- sum((W_r1*((x.RLmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[idx.RLmAMSE]))^2)

      AMSE_Sample_RLmAMSE_Pow[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)
    }
  }

  if(anyNA(beta_mVc) || anyNA(beta_mMSE) || anyNA(beta_RLmAMSE) || anyNA(beta_RLmAMSE_LO) || anyNA(beta_RLmAMSE_Pow))
  {
    stop("There are NA or NaN values")
  }

  Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc,PI.RLmAMSE,PI.RLmAMSE_LO,PI.RLmAMSE_Pow)
  colnames(Full_SP)<-c("A-Optimality","L-Optimality","RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),
                       paste0("RLmAMSE Power ",Alpha))

  Sampling_Methods<-factor(c("A-Optimality","L-Optimality","RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),
                                paste0("RLmAMSE Power ",Alpha)))

  # Beta Data
  Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),
                              rbind(beta_mMSE,beta_mVc,beta_RLmAMSE,
                                    do.call(rbind,beta_RLmAMSE_LO),
                                    do.call(rbind,beta_RLmAMSE_Pow)))

  if(all(X[,1] == 1)){
    colnames(Beta_Data)[-1]<-c("r2",paste0("Beta_",0:(ncol(X)-1)))
  } else {
    colnames(Beta_Data)[-1]<-c("r2",paste0("Beta_",1:(ncol(X))))
  }

  # AMSE Sample Data
  AMSE_Sample_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r2)),
                                        rbind(AMSE_Sample_mMSE,AMSE_Sample_mVc,
                                              AMSE_Sample_RLmAMSE,
                                              do.call(rbind,AMSE_Sample_RLmAMSE_LO),
                                              do.call(rbind,AMSE_Sample_RLmAMSE_Pow)))
  colnames(AMSE_Sample_Data)[-1]<-c("r2","Variance","Bias.2","AMSE")

  AMSE_Sample_Data[,-c(1,2)]<-AMSE_Sample_Data[,-c(1,2)]/r2

  all_r<-c(r1,r2)
  # Sample Data
  for(j in 1:length(Alpha)){
    names(Sample.RLmAMSE_LO[[j]])<-names(Sample.RLmAMSE_Pow[[j]])<-paste0("Alpha_",Alpha[j],"_",all_r)
  }

  names(Sample.mMSE)<-names(Sample.mVc)<-names(Sample.RLmAMSE)<-c(r1,r2)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "AMSE_Estimates"=AMSE_Sample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_RLmAMSE"=Sample.RLmAMSE,
            "Sample_RLmAMSE_Log_Odds"=Sample.RLmAMSE_LO,
            "Sample_RLmAMSE_Power"=Sample.RLmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","poisson")
  return(ans)
}
