#' Subsampling under Poisson regression for a potentially misspecified model
#'
#' Using this function sample from big data under Poisson regression for a potentially misspecified model.
#' Subsampling probabilities are obtained based on the A-, L- and L1- optimality criteria
#' with the RLmAMSE (Reduction of Loss by minimizing the Average Mean Squared Error).
#'
#' @usage
#' modelMissPoiSub(r0,r,Y,X,N,Alpha,proportion,model="Auto")
#'
#' @param r0          sample size for initial random sample
#' @param r           final sample size including initial(r0) and optimal(r1) samples
#' @param Y           response data or Y
#' @param X           covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N           size of the big data
#' @param Alpha       scaling factor when using Log Odds or Power functions to magnify the probabilities
#' @param proportion  a proportion of the big data is used to help estimate AMSE values from the subsamples
#' @param model       formula for the model used in the GAM or the default choice
#'
#' @details
#' \strong{The article for this function is in preparation for publication. Please be patient.}
#'
#' Two stage subsampling algorithm for big data under Poisson regression for potential model misspecification.
#'
#' First stage is to obtain a random sample of size \eqn{r_0} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-, L1-optimality criteria,
#' RLmAMSE and enhanced RLmAMSE (log-odds and power) subsampling methods.
#'
#' Through the estimated subsampling probabilities a sample of size \eqn{r_1 \ge r_0} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for A-, L-, L1-optimality,
#' RLmAMSE and enhanced RLmAMSE (log-odds and power).
#'
#' \strong{NOTE} :  If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r_1 \ge r_0} is not satisfied then an error message will be produced.
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
#' or NAN or not one of the defined inputs an error message is printed. As a default we have set
#' \code{model="Auto"}, which is the main effects model wit the spline terms.
#'
#' @return
#' The output of \code{modelMissPoiSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{Utility_Estimates} estimated A-, L- and L1- optimality values for the obtained subsamples
#'
#' \code{AMSE_Estimates} matrix of estimated AMSE values after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
#'
#' \code{Sample_L1-Optimality} list of indexes for the initial and optimal samples obtained based on L1-Optimality criteria
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
#' r0<-300; r<-rep(100*c(6,9),50);
#' Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#' \dontrun{
#' Results<-modelMissPoiSub(r0 = r0, r = r,
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
modelMissPoiSub <- function(r0,r,Y,X,N,Alpha,proportion,model="Auto"){
  if(any(is.na(c(r0,r,N,Alpha,proportion))) | any(is.nan(c(r0,r,N,Alpha,proportion)))){
    stop("NA or Infinite or NAN values in the r0,r,N,Alpha or proportion")
  }

  if((length(r0) + length(N) + length(proportion)) != 3){
    stop("proportion, r0 or N has a value greater than length one")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if(anyNA(Y) | anyNA(X) | any(is.nan(Y)) | any(is.nan(X)) ){
    stop("NA or Infinite or NAN values in the Y or X")
  }

  if(any((2*r0) > r)){
    stop("2*r0 cannot be greater than r at any point")
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

  } else if(model == "Auto"){

    if(all(X[,1] == 1)){
      main_effects <- paste0("s(",colnames(X[, -1]),")")
      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
    } else{
      main_effects <- paste0("s(",colnames(X),")")
      my_formula<-stats::as.formula(paste("Y ~ ",paste(main_effects,collapse = " + ")))
    }
  } else if(rlang::is_formula(model)){
    my_formula <- model
  } else {
    stop("The model for GAM is not a valid formula or input")
  }

  PI.prop <- rep(1/N,N)
  idx.prop <- sample(1:N, size = r0, replace = TRUE, prob = PI.prop)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]

  pinv.prop <- rep(N,r0)

  fit.prop <- stats::glm(y.prop~x.prop-1,family="quasipoisson")
  beta.prop <- fit.prop$coefficients

  if (anyNA(beta.prop)){
    stop("There are NA or NaN values in the model parameters")
  }
  Xbeta_Final <- X %*% beta.prop
  Lambda.prop  <- exp(Xbeta_Final)
  W_All <- as.vector(Lambda.prop)

  #calculate f_hat
  Assumed_Data<-data.frame(Y=y.prop,x.prop)
  fit_GAM<-gam::gam(my_formula,data=Assumed_Data,family = "quasipoisson")
  Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X))
  Lambda_GAM<-exp(Xbeta_GAM)
  f_estimate<-Xbeta_GAM - Xbeta_Final

  if(proportion*N != r0){
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

  ## idx.prop only
  f_estimate.prop<-f_estimate[idx.prop]
  H.prop1 <-crossprod(x.prop,(x.prop * w.prop))
  H.prop <-solve(H.prop1)

  lambda_Tr<-exp((Xbeta_Final) + f_estimate)
  W.Tr<-as.vector(lambda_Tr)
  lambda_Tr.prop<-lambda_Tr[idx.prop]
  w_Tr.prop<-W.Tr[idx.prop]
  H_Tr.prop <- crossprod(x.prop, (x.prop * w_Tr.prop))

  Temp.prop<-(w.prop*x.prop)%*%H.prop
  b_r.prop <- crossprod(x.prop, (lambda_Tr.prop-lambda.prop))
  L1_r0 <- tr(Temp.prop%*%H_Tr.prop%*%t(Temp.prop))
  Temp_diff<-((x.prop%*%H.prop%*%b_r.prop)-f_estimate.prop)
  L2_r0 <- sum((w.prop*Temp_diff)^2)
  L1_r0_1 <- tr(Temp.prop%*%H.prop1%*%t(Temp.prop))

  L_All_Temp<-cbind(L1_r0+L2_r0,L1_r0_1)
  r0_Temp<-(1+1/r0)

  a<-NULL
  L_All<-foreach::foreach(a=1:N,.combine = rbind,.packages = "psych") %dopar% {
    X_a<-matrix(X[a,],nrow=1)
    X_r0<-X[c(idx.prop,a),]
    f_r0<-c(f_estimate.prop,f_estimate[a])

    lambda_r0<-c(lambda.prop,Lambda.prop[a])
    W_r0<-c(w.prop,W_All[a])

    lambda_Tr0<-c(lambda_Tr.prop,lambda_Tr[a])
    W_Tr0<-c(w_Tr.prop,W.Tr[a])

    XtW_X<-H.prop1 + crossprod(X_a,(X_a * W_All[a]))
    H_r0 <-solve(XtW_X)
    Temp1<-(w.prop*x.prop)%*%H_r0

    H_Tr0<- H_Tr.prop + crossprod(X_a,(X_a*W.Tr[a]))

    Temp1_H <- Temp1 %*% H_Tr0
    diag_Temp <- rowSums(Temp1_H * Temp1)
    L1_r0 <- r0_Temp*sum(diag_Temp)

    XH_r0 <- x.prop %*% H_r0
    XH_b_r0 <- XH_r0 %*% b_r.prop
    diff <- r0_Temp*XH_b_r0 - f_estimate.prop
    L2_r0 <- sum((w.prop * diff)^2)

    Temp1_H <- Temp1 %*% XtW_X
    diag_Temp <- rowSums(Temp1_H * Temp1)
    L1_r0_1 <- r0_Temp*sum(diag_Temp)

    c(L1_r0+L2_r0,L1_r0_1)
  }

  # L1-optimality
  L_All_Final<-abs(L_All_Temp[,2]-L_All[,2])
  PI.L1_optimality <- (max(L_All_Final) - L_All_Final) / sum(max(L_All_Final) - L_All_Final)

  # RLmAMSE
  L_All_Final<-abs(L_All_Temp[,1]-L_All[,1])
  PI.RLmAMSE <- (max(L_All_Final) - L_All_Final) / sum(max(L_All_Final) - L_All_Final)

  # RLmAMSE LO
  PI.RLmAMSE_LO<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    (1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1)/
      sum((1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1))},PI.RLmAMSE)

  # RLmAMSE Pow
  PI.RLmAMSE_Pow<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    PI.RLmAMSE^Alpha/sum(PI.RLmAMSE^Alpha)},PI.RLmAMSE)

  # For the Model with already available subsampling probabilities mVc
  beta_mVc<-matrix(nrow = length(r),ncol = ncol(X)+1 )
  Utility_mVc<-matrix(nrow = length(r),ncol = 4 )
  AMSE_Sample_mVc<-matrix(nrow = length(r),ncol = 4); Sample.mVc<-list();

  # For the Model with already available subsampling probabilities mMSE
  beta_mMSE<-matrix(nrow = length(r),ncol = ncol(X)+1 )
  Utility_mMSE<-matrix(nrow = length(r),ncol = 4 )
  AMSE_Sample_mMSE<-matrix(nrow = length(r),ncol = 4); Sample.mMSE<-list();

  # For the Model with already available subsampling probabilities L1-optimality
  beta_L1_optimality<-matrix(nrow = length(r),ncol = ncol(X)+1 )
  Utility_L1_optimality<-matrix(nrow = length(r),ncol = 4 )
  AMSE_Sample_L1_optimality<-matrix(nrow = length(r),ncol = 4); Sample.L1_optimality<-list();

  # For the Model with model robust subsampling probabilities RLmAMSE
  beta_RLmAMSE<-matrix(nrow = length(r),ncol = ncol(X)+1 )
  Utility_RLmAMSE<-matrix(nrow = length(r),ncol = 4 )
  AMSE_Sample_RLmAMSE<-matrix(nrow = length(r),ncol = 4); Sample.RLmAMSE<-list();

  # For the Model with model robust subsampling probabilities RLmAMSE LO
  beta_RLmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r),ncol(X)+1)),simplify = FALSE)
  Utility_RLmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r),4)),simplify = FALSE)
  AMSE_Sample_RLmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r),4)),simplify = FALSE)
  Sample.RLmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(r)+1)));

  # For the Model with model robust subsampling probabilities RLmAMSE Pow
  beta_RLmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r),ncol(X)+1)),simplify = FALSE)
  Utility_RLmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r),4)),simplify = FALSE)
  AMSE_Sample_RLmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r),4)),simplify = FALSE)
  Sample.RLmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(r)+1)));

  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.L1_optimality[[1]]<-Sample.RLmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.RLmAMSE_LO[[j]][[1]]<-Sample.RLmAMSE_Pow[[j]][[1]]<-idx.prop
  }

  Utility_mVc[,1]<-Utility_mMSE[,1]<-Utility_L1_optimality[,1]<-
    Utility_RLmAMSE[,1]<-r

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r))
  {
    # mVc
    idx.mVc <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.mVc)

    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]

    fit.mVc <- stats::glm(y.mVc~x.mVc-1,family="quasipoisson",
                          weights = c(1 / PI.mVc[idx.mVc], pinv.prop))

    beta_mVc[i,] <- c(r[i],fit.mVc$coefficients)
    idx.mVc->Sample.mVc[[i+1]]

    Xbeta_r0<-x.mVc%*%Beta_Estimate_Full
    lambda_r0<-exp(Xbeta_r0)
    W_r0<-as.vector(lambda_r0)
    H_r0 <-solve(crossprod(x.mVc, x.mVc * W_r0))
    Temp1<-(W_r0 * x.mVc)%*%H_r0

    x.mVc_t<-t(x.mVc)
    Tempsy<-x.mVc%*%H_r0%*%x.mVc_t
    Tempsy1<-(W_r0 * Tempsy) * W_r0

    Utility_mVc[i,-1]<-c(psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

    f_r0<-F_Estimate_Full[c(idx.mVc, idx.prop)]
    lambda_Tr0<-exp(Xbeta_r0 + f_r0)
    W_Tr0<-as.vector(lambda_Tr0)
    H_Tr0 <- crossprod(x.mVc, x.mVc * W_Tr0)
    b_r0 <- crossprod(x.mVc, lambda_Tr0-lambda_r0)

    Temp1_t<-t(Temp1)
    L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
    Temp_diff<-((x.mVc%*%H_r0%*%b_r0) - f_r0)
    L2_r0 <- sum((W_r0*Temp_diff)^2)

    AMSE_Sample_mVc[i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)

    # mMSE
    idx.mMSE <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]

    fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1,family = "quasipoisson",
                           weights =c(1 / PI.mMSE[idx.mMSE], pinv.prop) )

    beta_mMSE[i,] <- c(r[i],fit.mMSE$coefficients)
    idx.mMSE->Sample.mMSE[[i+1]]

    Xbeta_r0<-x.mMSE%*%Beta_Estimate_Full
    lambda_r0<-exp(Xbeta_r0)
    W_r0<-as.vector(lambda_r0)
    H_r0 <-solve(crossprod(x.mMSE, x.mMSE * W_r0))
    Temp1<-(W_r0*x.mMSE)%*%H_r0

    x.mMSE_t<-t(x.mMSE)
    Tempsy<-x.mMSE%*%H_r0%*%x.mMSE_t
    Tempsy1<-(W_r0 * Tempsy) * W_r0

    Utility_mMSE[i,-1]<-c(psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

    f_r0<-F_Estimate_Full[c(idx.mMSE, idx.prop)]
    lambda_Tr0<-exp(Xbeta_r0 + f_r0)
    W_Tr0<-as.vector(lambda_Tr0)
    H_Tr0 <- crossprod(x.mMSE, x.mMSE * W_Tr0)
    b_r0 <- crossprod(x.mMSE, lambda_Tr0-lambda_r0)

    Temp1_t<-t(Temp1)
    L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
    Temp_diff<-((x.mMSE%*%H_r0%*%b_r0) - f_r0)
    L2_r0 <- sum((W_r0*Temp_diff)^2)

    AMSE_Sample_mMSE[i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)

    # L1-optimality
    idx.L1_optimality <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.L1_optimality)

    x.L1_optimality <- X[c(idx.L1_optimality, idx.prop),]
    y.L1_optimality <- Y[c(idx.L1_optimality, idx.prop)]

    fit.L1_optimality <- stats::glm(y.L1_optimality~x.L1_optimality-1,family = "quasipoisson",
                           weights =c(1 / PI.L1_optimality[idx.L1_optimality], pinv.prop) )

    beta_L1_optimality[i,] <- c(r[i],fit.L1_optimality$coefficients)
    idx.L1_optimality->Sample.L1_optimality[[i+1]]

    Xbeta_r0<-x.L1_optimality%*%Beta_Estimate_Full
    lambda_r0<-exp(Xbeta_r0)
    W_r0<-as.vector(lambda_r0)
    H_r0 <-solve(crossprod(x.L1_optimality, x.L1_optimality * W_r0))
    Temp1<-(W_r0*x.L1_optimality)%*%H_r0

    x.L1_optimality_t<-t(x.L1_optimality)
    Tempsy<-x.L1_optimality%*%H_r0%*%x.L1_optimality_t
    Tempsy1<-(W_r0 * Tempsy) * W_r0

    Utility_L1_optimality[i,-1]<-c(psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

    f_r0<-F_Estimate_Full[c(idx.L1_optimality, idx.prop)]
    lambda_Tr0<-exp(Xbeta_r0 + f_r0)
    W_Tr0<-as.vector(lambda_Tr0)
    H_Tr0 <- crossprod(x.L1_optimality, x.L1_optimality * W_Tr0)
    b_r0 <- crossprod(x.L1_optimality, lambda_Tr0-lambda_r0)

    Temp1_t<-t(Temp1)
    L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
    Temp_diff<-((x.L1_optimality%*%H_r0%*%b_r0) - f_r0)
    L2_r0 <- sum((W_r0*Temp_diff)^2)

    AMSE_Sample_L1_optimality[i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)

    # RLmAMSE
    idx.RLmAMSE <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.RLmAMSE)

    x.RLmAMSE <- X[c(idx.RLmAMSE,idx.prop),]
    y.RLmAMSE <- Y[c(idx.RLmAMSE,idx.prop)]

    fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                              weights = c(1 / PI.RLmAMSE[idx.RLmAMSE],pinv.prop))

    beta_RLmAMSE[i,] <- c(r[i],fit.RLmAMSE$coefficients)

    idx.RLmAMSE->Sample.RLmAMSE[[i+1]]

    Xbeta_r0<-x.RLmAMSE%*%Beta_Estimate_Full
    lambda_r0<-exp(Xbeta_r0)
    W_r0<-as.vector(lambda_r0)
    H_r0 <-solve(crossprod(x.RLmAMSE, x.RLmAMSE * W_r0))
    Temp1<-(W_r0*x.RLmAMSE)%*%H_r0

    x.RLmAMSE_t<-t(x.RLmAMSE)
    Tempsy<-x.RLmAMSE%*%H_r0%*%x.RLmAMSE_t
    Tempsy1<-(W_r0 * Tempsy) * W_r0

    Utility_RLmAMSE[i,-1]<-c(psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

    f_r0<-F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]
    lambda_Tr0<-exp((Xbeta_r0) + f_r0)
    W_Tr0<-as.vector(lambda_Tr0)
    H_Tr0 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr0)
    b_r0 <- crossprod(x.RLmAMSE, lambda_Tr0-lambda_r0)

    Temp1_t<-t(Temp1)
    L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
    Temp_diff<-((x.RLmAMSE%*%H_r0%*%b_r0) - f_r0)
    L2_r0 <- sum((W_r0*Temp_diff)^2)

    AMSE_Sample_RLmAMSE[i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)

    for (j in 1:length(Alpha))
    {
      # RLmAMSE Log Odds
      idx.RLmAMSE <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.RLmAMSE_LO[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE,idx.prop),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE,idx.prop)]

      fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                                weights = c(1 / PI.RLmAMSE_LO[idx.RLmAMSE,j],pinv.prop))

      beta_RLmAMSE_LO[[j]][i,] <- c(r[i],fit.RLmAMSE$coefficients)
      idx.RLmAMSE->Sample.RLmAMSE_LO[[j]][[i+1]]

      Xbeta_r0<-x.RLmAMSE%*%Beta_Estimate_Full
      lambda_r0<-exp(Xbeta_r0)
      W_r0<-as.vector(lambda_r0)
      H_r0 <-solve(crossprod(x.RLmAMSE, x.RLmAMSE * W_r0))
      Temp1<-(W_r0*x.RLmAMSE)%*%H_r0

      x.RLmAMSE_t<-t(x.RLmAMSE)
      Tempsy<-x.RLmAMSE%*%H_r0%*%x.RLmAMSE_t
      Tempsy1<-(W_r0 * Tempsy) * W_r0

      Utility_RLmAMSE_LO[[j]][i,]<-c(r[i],psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

      f_r0<-F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]
      lambda_Tr0<-exp((Xbeta_r0) + f_r0)
      W_Tr0<-as.vector(lambda_Tr0)
      H_Tr0 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr0)
      b_r0 <- crossprod(x.RLmAMSE, lambda_Tr0-lambda_r0)

      Temp1_t<-t(Temp1)
      L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
      Temp_diff<-((x.RLmAMSE%*%H_r0%*%b_r0) - f_r0)
      L2_r0 <- sum((W_r0*Temp_diff)^2)

      AMSE_Sample_RLmAMSE_LO[[j]][i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)

      # RLmAMSE Power
      idx.RLmAMSE <- sample(1:N, size = r[i]-r0, replace = TRUE, prob = PI.RLmAMSE_Pow[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE,idx.prop),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE,idx.prop)]

      fit.RLmAMSE <- stats::glm(y.RLmAMSE~x.RLmAMSE-1,family="quasipoisson",
                                weights = c(1 / PI.RLmAMSE_Pow[idx.RLmAMSE,j],pinv.prop))

      beta_RLmAMSE_Pow[[j]][i,] <- c(r[i],fit.RLmAMSE$coefficients)

      idx.RLmAMSE->Sample.RLmAMSE_Pow[[j]][[i+1]]

      Xbeta_r0<-x.RLmAMSE%*%Beta_Estimate_Full
      lambda_r0<-exp(Xbeta_r0)
      W_r0<-as.vector(lambda_r0)
      H_r0 <-solve(t(x.RLmAMSE) %*% (x.RLmAMSE * W_r0))
      Temp1<-(W_r0*x.RLmAMSE)%*%H_r0

      x.RLmAMSE_t<-t(x.RLmAMSE)
      Tempsy<-x.RLmAMSE%*%H_r0%*%x.RLmAMSE_t
      Tempsy1<-(W_r0 * Tempsy) * W_r0

      Utility_RLmAMSE_Pow[[j]][i,]<-c(r[i],psych::tr(H_r0),psych::tr(Tempsy),psych::tr(Tempsy1))

      f_r0<-F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]
      lambda_Tr0<-exp((Xbeta_r0) + f_r0)
      W_Tr0<-as.vector(lambda_Tr0)
      H_Tr0 <- crossprod(x.RLmAMSE, x.RLmAMSE * W_Tr0)
      b_r0 <- crossprod(x.RLmAMSE, lambda_Tr0-lambda_r0)

      Temp1_t<-t(Temp1)
      L1_r0 <- psych::tr(Temp1%*%H_Tr0%*%Temp1_t)
      Temp_diff<-((x.RLmAMSE%*%H_r0%*%b_r0) - f_r0)
      L2_r0 <- sum((W_r0*Temp_diff)^2)

      AMSE_Sample_RLmAMSE_Pow[[j]][i,]<-c(r[i],L1_r0,L2_r0,L1_r0+L2_r0)
    }
  }

  if(anyNA(beta_mVc) || anyNA(beta_mMSE) || anyNA(beta_L1_optimality) ||
     anyNA(beta_RLmAMSE) || anyNA(beta_RLmAMSE_LO) || anyNA(beta_RLmAMSE_Pow))
  {
    stop("There are NA or NaN values")
  }

  Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc,PI.L1_optimality,PI.RLmAMSE,PI.RLmAMSE_LO,PI.RLmAMSE_Pow)
  colnames(Full_SP)<-c("A-Optimality","L-Optimality","L1-Optimality",
                       "RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),
                       paste0("RLmAMSE Power ",Alpha))

  Sampling_Methods<-factor(c("A-Optimality","L-Optimality","L1-Optimality",
                             "RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),
                             paste0("RLmAMSE Power ",Alpha)))

  # Beta Data
  Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                              rbind(beta_mMSE,beta_mVc,beta_L1_optimality,beta_RLmAMSE,
                                    do.call(rbind,beta_RLmAMSE_LO),
                                    do.call(rbind,beta_RLmAMSE_Pow)))

  if(all(X[,1] == 1)){
    colnames(Beta_Data)[-1]<-c("r",paste0("Beta_",0:(ncol(X)-1)))
  } else {
    colnames(Beta_Data)[-1]<-c("r",paste0("Beta_",1:(ncol(X))))
  }

  # Utility Data
  Utility_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                 rbind(Utility_mMSE, Utility_mVc,Utility_L1_optimality,Utility_RLmAMSE,
                                       do.call(rbind,Utility_RLmAMSE_LO),
                                       do.call(rbind,Utility_RLmAMSE_Pow)))

  colnames(Utility_Data)<-c("Method","r","A-Optimality","L-Optimality","L1-Optimality")

  # AMSE Sample Data
  AMSE_Sample_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(r)),
                                        rbind(AMSE_Sample_mMSE,AMSE_Sample_mVc,
                                              AMSE_Sample_L1_optimality,AMSE_Sample_RLmAMSE,
                                              do.call(rbind,AMSE_Sample_RLmAMSE_LO),
                                              do.call(rbind,AMSE_Sample_RLmAMSE_Pow)))
  colnames(AMSE_Sample_Data)[-1]<-c("r","Variance","Bias.2","AMSE")

  AMSE_Sample_Data[,-c(1,2)]<-AMSE_Sample_Data[,-c(1,2)]/r

  all_r<-c(r0,r)
  # Sample Data
  for(j in 1:length(Alpha)){
    names(Sample.RLmAMSE_LO[[j]])<-names(Sample.RLmAMSE_Pow[[j]])<-paste0("Alpha_",Alpha[j],"_",all_r)
  }

  names(Sample.mMSE)<-names(Sample.mVc)<-names(Sample.L1_optimality)<-
    names(Sample.RLmAMSE)<-c(r0,r)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Utility_Estimates"=Utility_Data,
            "AMSE_Estimates"=AMSE_Sample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_L1-Optimality"=Sample.L1_optimality,
            "Sample_RLmAMSE"=Sample.RLmAMSE,
            "Sample_RLmAMSE_Log_Odds"=Sample.RLmAMSE_LO,
            "Sample_RLmAMSE_Power"=Sample.RLmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","poisson")
  return(ans)
}
