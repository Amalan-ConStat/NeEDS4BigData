#' Subsampling under linear regression for a potentially misspecified model
#'
#' Using this function sample from big data under linear regression for a potentially misspecified model.
#' Subsampling probabilities are obtained based on the A-, L- and L1- optimality criteria
#' with the RLmAMSE (Reduction of Loss by minimizing the Average Mean Squared Error).
#'
#' @usage
#' modelMissLinSub(r0,rf,Y,X,N,Alpha,proportion,model="Auto")
#'
#' @param r0          sample size for initial random sample
#' @param rf          final sample size including initial(r0) and optimal(r) samples
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
#' Two stage subsampling algorithm for big data under linear regression for potential model misspecification.
#'
#' First stage is to obtain a random sample of size \eqn{r_0} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-, L1-optimality criteria,
#' RLmAMSE and enhanced RLmAMSE (log-odds and power) subsampling methods.
#'
#' Through the estimated subsampling probabilities a sample of size \eqn{r \ge r_0} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for A-, L-, L1-optimality,
#' RLmAMSE and enhanced RLmAMSE (log-odds and power).
#'
#' \strong{NOTE} :  If input parameters are not in given domain conditions
#' necessary error messages will be provided to go further.
#'
#' If \eqn{r \ge r_0} is not satisfied then an error message will be produced.
#'
#' If the big data \eqn{X,Y} has any missing values then an error message will be produced.
#'
#' The big data size \eqn{N} is compared with the sizes of \eqn{X,Y},F_estimate_Full and
#' if they are not aligned an error message will be produced.
#'
#' If \eqn{\alpha > 1} for the scaling factor is not satisfied an error message will be produced.
#'
#' If proportion is not in the region of \eqn{(0,1]} an error message will be produced.
#'
#' \code{model} is a formula input formed based on the covariates through the spline terms (s()),
#' squared term (I()), interaction terms (lo()) or automatically. If \code{model} is empty or NA
#' or NAN or not one of the defined inputs an error message is printed. As a default we have set
#' \code{model="Auto"}, which is the main effects model wit the spline terms.
#'
#' @return
#' The output of \code{modelMissLinSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon after subsampling
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
#' \code{Sample_RLmAMSE} list of indexes for the optimal samples obtained based obtained based on RLmAMSE
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
#' Beta<-c(-1,0.75,0.75,1); Var_Epsilon<-0.5; family <- "linear"; N<-10000
#' X_1 <- replicate(2,stats::runif(n=N,min = -1,max = 1))
#'
#' Temp<-Rfast::rowprods(X_1)
#' Misspecification <- (Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2)
#' X_Data <- cbind(X0=1,X_1);
#' Full_Data<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family)
#'
#' r0<-300; rf<-rep(100*c(6,9),50);
#' Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#' \dontrun{
#' Results<-modelMissLinSub(r0 = r0, rf = rf,
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
#' @importFrom gam lo
#' @importFrom Rfast rowprods
#' @importFrom utils combn
#' @importFrom psych tr
#' @importFrom rlang is_formula
#' @export
modelMissLinSub <- function(r0,rf,Y,X,N,Alpha,proportion,model="Auto"){
  if(any(is.na(c(r0,rf,N,Alpha,proportion))) | any(is.nan(c(r0,rf,N,Alpha,proportion)))){
    stop("NA or Infinite or NAN values in the r0,rf,N,Alpha or proportion")
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

  if(any((2*r0) > rf)){
    stop("2*r0 cannot be greater than rf at any point")
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

  idx.prop <- sample(1:N, size = r0, replace = TRUE)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]

  beta.prop<-solve(a=crossprod(x.prop),b=crossprod(x.prop,y.prop))
  Xbeta_Final<-X%*%beta.prop
  Var.prop<-sum((Y-Xbeta_Final)^2)/N
  Epsilon.prop<-Y-Xbeta_Final

  #calculate f_hat
  Assumed_Data<-data.frame(Y=y.prop,x.prop)
  fit_GAM<-gam::gam(my_formula,data=Assumed_Data)
  Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X))
  f_estimate<-Xbeta_GAM-Xbeta_Final
  Var_GAM.prop<-sum((Y-Xbeta_GAM)^2)/N

  if(proportion*N != r0){
    idx.proportion <- sample(1:N, size = ceiling(proportion*N), replace = TRUE)
    Y_proportion<-Y[idx.proportion]
    X_proportion<-X[idx.proportion,]

    Proportion_Data<-data.frame(Y=Y_proportion,X_proportion)
    fit_GAM_Proportion<-gam::gam(my_formula,data=Proportion_Data)
    Xbeta_GAM_Proportion<-gam::predict.Gam(fit_GAM_Proportion,newdata = data.frame(X))

    beta_proportion<-solve(a=crossprod(X_proportion),b=crossprod(X_proportion,Y_proportion))
    Xbeta_proportion<-X%*%beta_proportion

    Var_GAM_Full<-sum((Y-Xbeta_GAM_Proportion)^2)/N
    Var_Full<-sum((Y-Xbeta_proportion)^2)/N
    F_Estimate_Full<-Xbeta_GAM_Proportion-Xbeta_proportion
  } else {
    Var_GAM_Full<-Var_GAM.prop ; Var_Full<-Var.prop ; F_Estimate_Full<-f_estimate
  }

  # mVc
  PI.mVc <- sqrt(Epsilon.prop^2 * rowSums((X)^2))
  PI.mVc <- PI.mVc / sum(PI.mVc)

  # mMSE
  PI.mMSE <- sqrt(Epsilon.prop^2 * rowSums((X %*% solve(crossprod(X)))^2))
  PI.mMSE <- PI.mMSE / sum(PI.mMSE)

  Tempsy_Var_Gam_Var<-Var_GAM.prop*Var.prop^(-2)
  Var_prop_inv <- Var.prop^(-1)
  x.prop_t<-t(x.prop)
  f_estimate.prop<-f_estimate[idx.prop]
  Xt_X.prop<-crossprod(x.prop)
  L1_r0_1_Temp<-(1+1/r0)*Tempsy_Var_Gam_Var
  r0_Temp<-(1+1/r0)

  a<-NULL
  L_All <- foreach::foreach(a = 1:N, .combine = rbind,.packages = "psych") %dopar% {
    X_r0<-matrix(X[a,],nrow=1)
    Xt_X<-Xt_X.prop + crossprod(X_r0)

    Temp_Solve<-solve(Xt_X)
    Temp_Xr0_Solve<-x.prop %*% Temp_Solve
    Temp1 <- tcrossprod(Temp_Xr0_Solve,x.prop)

    L1_r0 <- L1_r0_1_Temp*psych::tr(Temp1)
    Temp_f_r0 <- Temp1%*%f_estimate.prop
    Temp_diff <- (r0_Temp*Temp_f_r0 - f_estimate.prop)
    L2_r0 <- sum((Var_prop_inv*Temp_diff)^2)

    L1_r0_1 <- r0_Temp*Var_prop_inv*psych::tr(Temp1)

    c(L1_r0+L2_r0,L1_r0_1)
  }

  Temp1<-x.prop%*%solve(crossprod(x.prop))%*%x.prop_t
  L1_r0 <- Tempsy_Var_Gam_Var*psych::tr(Temp1)
  Temp_diff<-(Temp1%*%f_estimate.prop-f_estimate.prop)
  L2_r0 <- sum((Var_prop_inv*Temp_diff)^2)
  L1_r0_1 <- Var_prop_inv*psych::tr(Temp1)
  L_All_Temp<-cbind(L1_r0+L2_r0,L1_r0_1)

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

  Var_Epsilon<-matrix(nrow=length(rf),ncol=5)

  # For the Model with already available subsampling probabilities - mVc
  beta_mVc<-matrix(nrow = length(rf),ncol = ncol(X)+1);
  Utility_mVc<-matrix(nrow = length(rf),ncol = 4);
  AMSE_Sample_mVc<-matrix(nrow = length(rf),ncol = 4) ; Sample.mVc<-list();

  # For the Model with already available subsampling probabilities - mMSE
  beta_mMSE<-matrix(nrow = length(rf),ncol = ncol(X)+1);
  Utility_mMSE<-matrix(nrow = length(rf),ncol = 4);
  AMSE_Sample_mMSE<-matrix(nrow = length(rf),ncol = 4) ; Sample.mMSE<-list();

  # For the Model with already available subsampling probabilities - L1-optimality
  beta_L1_optimality<-matrix(nrow = length(rf),ncol = ncol(X)+1);
  Utility_L1_optimality<-matrix(nrow = length(rf),ncol = 4);
  AMSE_Sample_L1_optimality<-matrix(nrow = length(rf),ncol = 4) ; Sample.L1_optimality<-list();

  # For the Real and Model with model robust subsampling probabilities RLmAMSE
  beta_RLmAMSE<-matrix(nrow = length(rf),ncol = ncol(X)+1 );
  Utility_RLmAMSE<-matrix(nrow = length(rf),ncol = 4);
  AMSE_Sample_RLmAMSE<-matrix(nrow = length(rf),ncol = 4) ; Sample.RLmAMSE<-list()

  # For the Model with model robust subsampling probabilities RLmAMSE LO
  beta_RLmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(rf),ncol(X)+1)),simplify = FALSE) ;
  Utility_RLmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(rf),4)),simplify = FALSE) ;
  AMSE_Sample_RLmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(rf),4)),simplify = FALSE);
  Sample.RLmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(rf)+1)));
  Var_RLmAMSE_LO<-matrix(nrow = length(rf),ncol=length(Alpha))

  # For the Model with model robust subsampling probabilities RLmAMSE Pow
  beta_RLmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(rf),ncol(X)+1)),simplify = FALSE) ;
  Utility_RLmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(rf),4)),simplify = FALSE) ;
  AMSE_Sample_RLmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(rf),4)),simplify = FALSE);
  Sample.RLmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(rf)+1)));
  Var_RLmAMSE_Pow<-matrix(nrow = length(rf),ncol=length(Alpha))

  Tempsy_Var_Gam_Var_AMSE<-Var_GAM_Full*Var_Full^(-2)
  VarFull_Inv<-Var_Full^(-1)
  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.L1_optimality[[1]]<-Sample.RLmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.RLmAMSE_LO[[j]][[1]]<-Sample.RLmAMSE_Pow[[j]][[1]]<-idx.prop
  }

  Var_Epsilon[,1]<-rf
  Utility_mVc[,1]<-Utility_mMSE[,1]<-Utility_L1_optimality[,1]<-
    Utility_RLmAMSE[,1]<-rf

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(rf))
  {
    # mVc
    idx.mVc <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.mVc)

    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    w.mVc <- c(1/PI.mVc[idx.mVc], rep(N,r0))

    Temp_Inv<-solve(crossprod(x.mVc))
    x.mVc_t<-t(x.mVc)
    Temp1<-x.mVc%*%Temp_Inv%*%x.mVc_t

    Utility_mVc[i,-1]<-c(VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                          VarFull_Inv*psych::tr(Temp1))

    pi4_r<-sqrt(rf[i]*w.mVc^(-1))
    X_r4<-x.mVc/pi4_r
    Y_r4<-y.mVc/pi4_r
    beta_mVc[i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
    Xbeta_Final<-X%*%beta_mVc[i,-1]
    Var_Epsilon[i,2]<-sum((Y-Xbeta_Final)^2)/N

    L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.mVc,idx.prop)] -
                                 F_Estimate_Full[c(idx.mVc, idx.prop)]))^2)

    AMSE_Sample_mVc[i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

    idx.mVc->Sample.mVc[[i+1]]

    # mMSE
    idx.mMSE <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    w.mMSE <- c(1/PI.mMSE[idx.mMSE], rep(N,r0))

    Temp_Inv<-solve(crossprod(x.mMSE))
    x.mMSE_t<-t(x.mMSE)
    Temp1<-x.mMSE%*%Temp_Inv%*%x.mMSE_t

    Utility_mMSE[i,-1]<-c(VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                          VarFull_Inv*psych::tr(Temp1))

    pi4_r<-sqrt(rf[i]*w.mMSE^(-1))
    X_r4<-x.mMSE/pi4_r
    Y_r4<-y.mMSE/pi4_r
    beta_mMSE[i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
    Xbeta_Final<-X%*%beta_mMSE[i,-1]
    Var_Epsilon[i,3]<-sum((Y-Xbeta_Final)^2)/N

    L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.mMSE,idx.prop)]-
                                 F_Estimate_Full[c(idx.mMSE, idx.prop)]))^2)

    AMSE_Sample_mMSE[i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

    idx.mMSE->Sample.mMSE[[i+1]]

    # L1-optimality
    idx.L1_optimality <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.L1_optimality)

    x.L1_optimality <- X[c(idx.L1_optimality, idx.prop),]
    y.L1_optimality <- Y[c(idx.L1_optimality, idx.prop)]
    w.L1_optimality <- c(1/PI.L1_optimality[idx.L1_optimality], rep(N,r0))

    Temp_Inv<-solve(crossprod(x.L1_optimality))
    x.L1_optimality_t<-t(x.L1_optimality)
    Temp1<-x.L1_optimality%*%Temp_Inv%*%x.L1_optimality_t

    Utility_L1_optimality[i,-1]<-c(VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                                   VarFull_Inv*psych::tr(Temp1))

    pi4_r<-sqrt(rf[i]*w.L1_optimality^(-1))
    X_r4<-x.L1_optimality/pi4_r
    Y_r4<-y.L1_optimality/pi4_r
    beta_L1_optimality[i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
    Xbeta_Final<-X%*%beta_L1_optimality[i,-1]
    Var_Epsilon[i,4]<-sum((Y-Xbeta_Final)^2)/N

    L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.L1_optimality,idx.prop)]-
                                 F_Estimate_Full[c(idx.L1_optimality, idx.prop)]))^2)

    AMSE_Sample_L1_optimality[i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

    idx.L1_optimality->Sample.L1_optimality[[i+1]]

    # RLmAMSE
    idx.RLmAMSE <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.RLmAMSE)

    x.RLmAMSE <- X[c(idx.RLmAMSE, idx.prop),]
    y.RLmAMSE <- Y[c(idx.RLmAMSE, idx.prop)]
    w.RLmAMSE <- c(1/PI.RLmAMSE[idx.RLmAMSE],rep(N,r0))

    Temp_Inv<-solve(crossprod(x.RLmAMSE))
    x.RLmAMSE_t<-t(x.RLmAMSE)
    Temp1<-x.RLmAMSE%*%Temp_Inv%*%x.RLmAMSE_t

    Utility_RLmAMSE[i,-1]<-c(VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                                   VarFull_Inv*psych::tr(Temp1))

    pi4_r<-sqrt(rf[i]*w.RLmAMSE^(-1))
    X_r4<-x.RLmAMSE/pi4_r
    Y_r4<-y.RLmAMSE/pi4_r
    beta_RLmAMSE[i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
    Xbeta_Final<-X%*%beta_RLmAMSE[i,-1]
    Var_Epsilon[i,5]<-sum((Y-Xbeta_Final)^2)/N

    L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE,idx.prop)] -
                                 F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]))^2)

    AMSE_Sample_RLmAMSE[i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

    idx.RLmAMSE->Sample.RLmAMSE[[i+1]] # Model robust RLmAMSE

    for(j in 1:length(Alpha))
    {
      # RLmAMSE Log Odds
      idx.RLmAMSE <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.RLmAMSE_LO[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE,idx.prop),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE,idx.prop)]
      w.RLmAMSE <- c(1/PI.RLmAMSE_LO[idx.RLmAMSE,j],rep(N,r0))

      Temp_Inv<-solve(crossprod(x.RLmAMSE))
      x.RLmAMSE_t<-t(x.RLmAMSE)
      Temp1<-x.RLmAMSE%*%Temp_Inv%*%x.RLmAMSE_t
      Utility_RLmAMSE_LO[[j]][i,]<-c(rf[i],VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                                     VarFull_Inv*psych::tr(Temp1))

      pi4_r<-sqrt(rf[i]*w.RLmAMSE^(-1))
      X_r4<-x.RLmAMSE/pi4_r
      Y_r4<-y.RLmAMSE/pi4_r
      beta_RLmAMSE_LO[[j]][i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
      Xbeta_Final<-X%*%beta_RLmAMSE_LO[[j]][i,-1]
      Var_RLmAMSE_LO[i,j]<-sum((Y-Xbeta_Final)^2)/N

      L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
      L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]-
                                   F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]))^2)

      AMSE_Sample_RLmAMSE_LO[[j]][i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

      idx.RLmAMSE->Sample.RLmAMSE_LO[[j]][[i+1]] # Model robust RLmAMSE

      # RLmAMSE Power
      idx.RLmAMSE <- sample(1:N, size = rf[i]-r0, replace = TRUE, prob = PI.RLmAMSE_Pow[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE,idx.prop),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE,idx.prop)]
      w.RLmAMSE <- c(1/PI.RLmAMSE_Pow[idx.RLmAMSE,j],rep(N,r0))

      Temp_Inv<-solve(crossprod(x.RLmAMSE))
      x.RLmAMSE_t<-t(x.RLmAMSE)
      Temp1<-x.RLmAMSE%*%Temp_Inv%*%x.RLmAMSE_t
      Utility_RLmAMSE_Pow[[j]][i,]<-c(rf[i],VarFull_Inv*psych::tr(Temp_Inv),VarFull_Inv*psych::tr(Temp1),
                                      VarFull_Inv*psych::tr(Temp1))

      pi4_r<-sqrt(rf[i]*w.RLmAMSE^(-1))
      X_r4<-x.RLmAMSE/pi4_r
      Y_r4<-y.RLmAMSE/pi4_r
      beta_RLmAMSE_Pow[[j]][i,]<-c(rf[i],solve(a=crossprod(X_r4),b=crossprod(X_r4,Y_r4)))
      Xbeta_Final<-X%*%beta_RLmAMSE_Pow[[j]][i,-1]
      Var_RLmAMSE_Pow[i,j]<-sum((Y-Xbeta_Final)^2)/N

      L1_r0 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
      L2_r0 <- sum((VarFull_Inv*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE,idx.prop)] -
                                   F_Estimate_Full[c(idx.RLmAMSE,idx.prop)]))^2)

      AMSE_Sample_RLmAMSE_Pow[[j]][i,]<-c(rf[i],L1_r0,L2_r0,L1_r0+L2_r0)

      idx.RLmAMSE->Sample.RLmAMSE_Pow[[j]][[i+1]] # Model robust RLmAMSE
    }
  }

  if(anyNA(beta_mMSE) || anyNA(beta_mVc) || anyNA(beta_L1_optimality) ||
     anyNA(beta_RLmAMSE) || anyNA(beta_RLmAMSE_LO) || anyNA(beta_RLmAMSE_Pow) ){
    stop("There are NA or NaN values")
  }

  Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc,PI.L1_optimality,PI.RLmAMSE,PI.RLmAMSE_LO,PI.RLmAMSE_Pow)
  colnames(Full_SP)<-c("A-Optimality","L-Optimality","L1-Optimality",
                       "RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),paste0("RLmAMSE Power ",Alpha))

  Sampling_Methods<-factor(c("A-Optimality","L-Optimality","L1-Optimality",
                             "RLmAMSE",paste0("RLmAMSE Log Odds ",Alpha),
                             paste0("RLmAMSE Power ",Alpha)))

  # Beta Data
  Beta_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(rf)),
                              rbind(beta_mMSE,beta_mVc,beta_L1_optimality,beta_RLmAMSE,
                                    do.call(rbind,beta_RLmAMSE_LO),
                                    do.call(rbind,beta_RLmAMSE_Pow)))

  if(all(X[,1] == 1)){
    colnames(Beta_Data)[-1]<-c("rf",paste0("Beta_",0:(ncol(X)-1)))
  } else {
    colnames(Beta_Data)[-1]<-c("rf",paste0("Beta_",1:(ncol(X))))
  }

  # Utility Data
  Utility_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(rf)),
                                 rbind(Utility_mMSE, Utility_mVc,Utility_L1_optimality,Utility_RLmAMSE,
                                       do.call(rbind,Utility_RLmAMSE_LO),
                                       do.call(rbind,Utility_RLmAMSE_Pow)))

  colnames(Utility_Data)<-c("Method","rf","A-Optimality","L-Optimality","L1-Optimality")

  # Var Data
  Var_Data<-cbind.data.frame(Var_Epsilon,Var_RLmAMSE_LO,Var_RLmAMSE_Pow)

  colnames(Var_Data)<-c("rf","A-Optimality","L-Optimality","L1-Optimality","RLmAMSE",
                        paste0("RLmAMSE Log Odds ",Alpha),
                        paste0("RLmAMSE Power ",Alpha))

  # AMSE Sample Data
  AMSE_Sample_Data<-cbind.data.frame("Method"=rep(Sampling_Methods,each=length(rf)),
                                     rbind(AMSE_Sample_mMSE,AMSE_Sample_mVc,
                                           AMSE_Sample_L1_optimality,AMSE_Sample_RLmAMSE,
                                           do.call(rbind,AMSE_Sample_RLmAMSE_LO),
                                           do.call(rbind,AMSE_Sample_RLmAMSE_Pow)))
  colnames(AMSE_Sample_Data)[-1]<-c("rf","Variance","Bias.2","AMSE")

  AMSE_Sample_Data[,-c(1,2)]<-AMSE_Sample_Data[,-c(1,2)]/rf

  all_r<-c(r0,rf)
  # Sample Data
  for(j in 1:length(Alpha)){
    names(Sample.RLmAMSE_LO[[j]])<-names(Sample.RLmAMSE_Pow[[j]])<-paste0("Alpha_",Alpha[j],"_",all_r)
  }

  names(Sample.mMSE)<-names(Sample.mVc)<-names(Sample.L1_optimality)<-names(Sample.RLmAMSE)<-c(r0,rf)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Variance_Epsilon_Estimates"=Var_Data,
            "Utility_Estimates"=Utility_Data,
            "AMSE_Estimates"=AMSE_Sample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_L1-Optimality"=Sample.L1_optimality,
            "Sample_RLmAMSE"=Sample.RLmAMSE,
            "Sample_RLmAMSE_Log_Odds"=Sample.RLmAMSE_LO,
            "Sample_RLmAMSE_Power"=Sample.RLmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","linear")
  return(ans)
}
