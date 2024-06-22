#' Model misspecified subsampling under linear regression
#'
#' Using this function subsample from big data under linear regression when we assumed the model that
#' describes the data is misspecified. Subsampling probabilities are obtained based on the A- and L-
#' optimality criterions with the LmAMSE (Loss on mean of asymptotic mean squared error).
#'
#' @usage
#' modelMissLinSub(r1,r2,Y,X,N,Alpha,Var_GAM_Full,Var_Full,F_Estimate_Full)
#'
#' @param r1                subsample size for initial random sampling
#' @param r2                subsample size for optimal sampling
#' @param Y                 response data or Y
#' @param X                 covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N                 size of the big data
#' @param Alpha             scaling factor when using Log Odds or Power functions to magnify the probabilities
#' @param Var_GAM_Full      estimate of Var_Epsilon after fitting the GAM model
#' @param Var_Full          estimate of Var_Epsilon after fitting the linear regression model
#' @param F_Estimate_Full   estimate of f that is the difference of linear predictor on GAM and linear model
#'
#' @details
#' Two stage subsampling algorithm for big data under linear regression for potential model misspecification.
#'
#' First stage is to obtain a random sample of size \eqn{r_1} and estimate the model parameters.
#' Using the estimated parameters subsampling probabilities are evaluated for A-, L-optimality criterion,
#' LmAMSE and enhanced LmAMSE(log-odds and power) subsampling methods.
#'
#' Through the estimated subsampling probabilities an subsample of size \eqn{r_2 \ge r_1} is obtained.
#' Finally, the two samples are combined and the model parameters are estimated for A- and L-optimality,
#' while for LmAMSE and enhanced LmAMSE (log-odds and power) only the optimal subsample is used.
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
#' If \eqn{\alpha > 1} for the scaling factor is not satisfied an error message will be produced.
#'
#' @return
#' The output of \code{modelMissLinSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon after subsampling
#'
#' \code{Loss_Estimates} matrix of estimated LmAMSE values after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal subsamples obtained based on A-Optimality criterion
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal subsamples obtained based on L-Optimality criterion
#'
#' \code{Sample_LmAMSE} list of indexes for the optimal subsamples obtained based obtained based on LmAMSE
#'
#' \code{Sample_LmAMSE_Log_Odds} list of indexes for the optimal subsamples obtained based on LmAMSE with Log Odds function
#'
#' \code{Sample_LmAMSE_Power} list of indexes for the optimal subsamples obtained based on LmAMSE with Power function
#'
#' \code{Subsampling_Probability} matrix of calculated subsampling probabilities
#'
#' @references
#' \insertRef{adewale2009robust}{NeEDS4BigData}
#' \insertRef{adewale2010robust}{NeEDS4BigData}
#'
#' @examples
#' No_Of_Var<-2; Beta<-c(-1,2,2,1); Var_Epsilon<-0.5; N<-10000;
#' MisspecificationType <- "Type 2 Squared"; family <- "linear"
#'
#' Full_Data<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)
#'
#' r1<-300; r2<-100*c(6,9); Original_Data<-Full_Data$Full_Data;
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#' \dontrun{
#' Results<-modelMissLinSub(r1 = r1, r2 = r2,
#'                          Y = as.matrix(Original_Data[,1]),
#'                          X = as.matrix(Original_Data[,-1]),
#'                          N = Full_Data$N,
#'                          Alpha = 10 ,
#'                          Var_GAM_Full = Full_Data$Variance_Epsilon$Real_GAM,
#'                          Var_Full = Full_Data$Variance_Epsilon$Estimate,
#'                          F_Estimate_Full = Full_Data$f$Real_GAM)
#'
#' # parallel::stopCluster(cl)
#'
#' plot_Beta(Results)
#' plot_LmAMSE(Results)
#' }
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @import foreach
#' @importFrom gam s
#' @importFrom Rfast rowprods
#' @importFrom psych tr
#' @export
modelMissLinSub <- function(r1,r2,Y,X,N,Alpha,Var_GAM_Full,Var_Full,F_Estimate_Full){
  if(any(is.na(c(r1,r2,N,Alpha))) | any(is.nan(c(r1,r2,N,Alpha)))){
    stop("NA or Infinite or NAN values in the r1,r2,N or Alpha")
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

  if(Alpha <= 1 | length(Alpha) > 1){
    stop("Scaling factor alpha is not greater than one or the length is more than one")
  }

  idx.prop <- sample(1:N, r1, T)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]

  beta.prop<-solve(a=t(x.prop)%*%x.prop,b=t(x.prop)%*%y.prop)
  Xbeta_Final<-as.vector(X%*%beta.prop)
  Var.prop<-sum((Y-Xbeta_Final)^2)/N

  Xbeta.prop<-X%*%beta.prop
  Epsilon.prop<-Y-Xbeta.prop

  my_formula<-stats::as.formula(paste("Y ~ ",paste(paste0("s(X",1:ncol(x.prop[,-1]),")"),collapse = " + "),"+",
                                      paste(paste0("s(",paste0(colnames(x.prop[,-1]),collapse = "*"),")"),
                                            collapse = " + ")))

  #calculate f_hat
  Assumed_Data<-data.frame(Y=y.prop,x.prop)
  fit_GAM<-gam::gam(my_formula,data=Assumed_Data)
  Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X))
  f_estimate<-Xbeta_GAM-Xbeta.prop
  Var_GAM.prop<-sum((Y-Xbeta_GAM)^2)/N

  if(is.null(Var_GAM_Full) || is.null(Var_Full) || is.null(F_Estimate_Full) ||
     anyNA(Var_GAM_Full) || anyNA(Var_Full) || anyNA(F_Estimate_Full)){

    Var_GAM_Full<-Var_GAM.prop ; Var_Full<-Var.prop ; F_Estimate_Full<-f_estimate
    message("Var_GAM_Full, Var_Full and F_Estimate_Full from the initial sample is used.\n")

  }

  # mVc
  PI.mVc <- sqrt(Epsilon.prop^2 * rowSums((X)^2))
  PI.mVc <- PI.mVc / sum(PI.mVc)

  # mMSE
  PI.mMSE <- sqrt(Epsilon.prop^2 * rowSums((X %*% solve(t(X)%*%X))^2))
  PI.mMSE <- PI.mMSE / sum(PI.mMSE)

  Tempsy_Var_Gam_Var<-Var_GAM.prop*Var.prop^(-2)

  a<-NULL
  L_All <- foreach::foreach(a = 1:N, .combine = rbind,.packages = "psych") %dopar% {
    X_r1 <- rbind(x.prop, X[a, ])
    f_r1 <- f_estimate[c(idx.prop, a)]

    Temp_Solve<-solve(t(X_r1) %*% X_r1)
    Temp1 <- X_r1 %*% Temp_Solve %*% t(X_r1)
    L1_r1 <- Tempsy_Var_Gam_Var*psych::tr(Temp1)
    L2_r1 <- sum((Var.prop^(-1)*(Temp1%*%f_r1 - f_r1))^2)

    L1_r1 + L2_r1
  }

  Temp1<-X[idx.prop,]%*%solve(t(X[idx.prop,])%*%X[idx.prop,])%*%t(X[idx.prop,])
  L1_r1 <- Tempsy_Var_Gam_Var*psych::tr(Temp1)
  L2_r1 <- sum((Var.prop^(-1)*(Temp1%*%f_estimate[idx.prop]-f_estimate[idx.prop]))^2)

  L_All_Temp<-L1_r1+L2_r1
  L_All_Final<-abs(L_All_Temp-L_All)

  # LmAMSE
  PI.LmAMSE <- (max(L_All_Final) - L_All_Final) / sum(max(L_All_Final) - L_All_Final)

  # LmAMSE LO
  PI.LmAMSE_LO<-apply(as.matrix(Alpha),1,function(Alpha,PI.LmAMSE){
    (1+exp(-Alpha*log(PI.LmAMSE/(1-PI.LmAMSE))))^(-1)/
      sum((1+exp(-Alpha*log(PI.LmAMSE/(1-PI.LmAMSE))))^(-1))},PI.LmAMSE)

  # LmAMSE Pow
  PI.LmAMSE_Pow<-apply(as.matrix(Alpha),1,function(Alpha,PI.LmAMSE){
    PI.LmAMSE^Alpha/sum(PI.LmAMSE^Alpha)},PI.LmAMSE)

  Var_Epsilon<-matrix(nrow=length(r2),ncol=4)

  # For the Model with already available Sub-sampling probabilities - mVc
  beta_mVc<-matrix(nrow = length(r2),ncol = ncol(X)+1);
  Loss_Subsample_mVc<-matrix(nrow = length(r2),ncol = 4) ; Sample.mVc<-list();

  # For the Model with already available Sub-sampling probabilities - mMSE
  beta_mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1);
  Loss_Subsample_mMSE<-matrix(nrow = length(r2),ncol = 4) ; Sample.mMSE<-list();

  # For the Real and Model with model robust Sub-sampling probabilities LmAMSE
  beta_LmAMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 );
  Loss_Subsample_LmAMSE<-matrix(nrow = length(r2),ncol = 4) ; Sample.LmAMSE<-list()

  # For the Model with model robust Sub-sampling probabilities LmAMSE LO
  beta_LmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(r2),ncol(X)+1)),simplify = FALSE) ;
  Loss_Subsample_LmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(r2),4)),simplify = FALSE);
  Sample.LmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));
  Var_LmAMSE_LO<-matrix(nrow = length(r2),ncol=length(Alpha))

  # For the Model with model robust Sub-sampling probabilities LmAMSE Pow
  beta_LmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(r2),ncol(X)+1)),simplify = FALSE) ;
  Loss_Subsample_LmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(r2),4)),simplify = FALSE);
  Sample.LmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));
  Var_LmAMSE_Pow<-matrix(nrow = length(r2),ncol=length(Alpha))

  Tempsy_Var_Gam_Var_Loss<-Var_GAM_Full*Var_Full^(-2)
  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.LmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.LmAMSE_LO[[j]][[1]]<-Sample.LmAMSE_Pow[[j]][[1]]<-idx.prop
  }

  Var_Epsilon[,1]<-r2

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    # mVc
    idx.mVc <- sample(1:N, r2[i]-r1, T, PI.mVc)

    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]
    w.mVc <- c(1/PI.mVc[idx.mVc], rep(N,r1))

    pi4_r<-sqrt(r2[i]*w.mVc^(-1))
    X_r4<-x.mVc/pi4_r
    Y_r4<-y.mVc/pi4_r
    beta_mVc[i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
    Xbeta_Final<-as.vector(X%*%beta_mVc[i,-1])
    Var_Epsilon[i,2]<-sum((Y-Xbeta_Final)^2)/N

    Temp1<-x.mVc%*%solve(t(x.mVc)%*%x.mVc)%*%t(x.mVc)
    L1_r1 <- Tempsy_Var_Gam_Var_Loss*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.mVc,idx.prop)] -
                                   F_Estimate_Full[c(idx.mVc, idx.prop)]))^2)

    Loss_Subsample_mVc[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    idx.mVc->Sample.mVc[[i+1]]

    # mMSE
    idx.mMSE <- sample(1:N, r2[i]-r1, T, PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]
    w.mMSE <- c(1/PI.mMSE[idx.mMSE], rep(N,r1))

    pi4_r<-sqrt(r2[i]*w.mMSE^(-1))
    X_r4<-x.mMSE/pi4_r
    Y_r4<-y.mMSE/pi4_r
    beta_mMSE[i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
    Xbeta_Final<-as.vector(X%*%beta_mMSE[i,-1])
    Var_Epsilon[i,3]<-sum((Y-Xbeta_Final)^2)/N

    Temp1<-x.mMSE%*%solve(t(x.mMSE)%*%x.mMSE)%*%t(x.mMSE)
    L1_r1 <- Tempsy_Var_Gam_Var_Loss*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.mMSE,idx.prop)]-
                                   F_Estimate_Full[c(idx.mMSE, idx.prop)]))^2)

    Loss_Subsample_mMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    idx.mMSE->Sample.mMSE[[i+1]]

    # LmAMSE
    idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE)

    x.LmAMSE <- X[c(idx.LmAMSE),]
    y.LmAMSE <- Y[c(idx.LmAMSE)]
    w.LmAMSE <- c(1/PI.LmAMSE[idx.LmAMSE])

    pi4_r<-sqrt(r2[i]*w.LmAMSE^(-1))
    X_r4<-x.LmAMSE/pi4_r
    Y_r4<-y.LmAMSE/pi4_r
    beta_LmAMSE[i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
    Xbeta_Final<-as.vector(X%*%beta_LmAMSE[i,-1])
    Var_Epsilon[i,4]<-sum((Y-Xbeta_Final)^2)/N

    Temp1<-x.LmAMSE%*%solve(t(x.LmAMSE)%*%x.LmAMSE)%*%t(x.LmAMSE)
    L1_r1 <- Tempsy_Var_Gam_Var_Loss*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.LmAMSE)] -
                                   F_Estimate_Full[c(idx.LmAMSE)]))^2)

    Loss_Subsample_LmAMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    idx.LmAMSE->Sample.LmAMSE[[i+1]] # Model robust LmAMSE

    for(j in 1:length(Alpha))
    {
      # LmAMSE Log Odds
      idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE_LO[,j])

      x.LmAMSE <- X[c(idx.LmAMSE),]
      y.LmAMSE <- Y[c(idx.LmAMSE)]
      w.LmAMSE <- c(1/PI.LmAMSE_LO[idx.LmAMSE,j])

      pi4_r<-sqrt(r2[i]*w.LmAMSE^(-1))
      X_r4<-x.LmAMSE/pi4_r
      Y_r4<-y.LmAMSE/pi4_r
      beta_LmAMSE_LO[[j]][i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
      Xbeta_Final<-as.vector(X%*%beta_LmAMSE_LO[[j]][i,-1])
      Var_LmAMSE_LO[i,j]<-sum((Y-Xbeta_Final)^2)/N

      Temp1<-x.LmAMSE%*%solve(t(x.LmAMSE)%*%x.LmAMSE)%*%t(x.LmAMSE)
      L1_r1 <- Tempsy_Var_Gam_Var_Loss*psych::tr(Temp1)
      L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.LmAMSE)]-
                                     F_Estimate_Full[c(idx.LmAMSE)]))^2)

      Loss_Subsample_LmAMSE_LO[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      idx.LmAMSE->Sample.LmAMSE_LO[[j]][[i+1]] # Model robust LmAMSE

      # LmAMSE Power
      idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE_Pow[,j])

      x.LmAMSE <- X[c(idx.LmAMSE),]
      y.LmAMSE <- Y[c(idx.LmAMSE)]
      w.LmAMSE <- c(1/PI.LmAMSE_Pow[idx.LmAMSE,j])

      pi4_r<-sqrt(r2[i]*w.LmAMSE^(-1))
      X_r4<-x.LmAMSE/pi4_r
      Y_r4<-y.LmAMSE/pi4_r
      beta_LmAMSE_Pow[[j]][i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
      Xbeta_Final<-as.vector(X%*%beta_LmAMSE_Pow[[j]][i,-1])
      Var_LmAMSE_Pow[i,j]<-sum((Y-Xbeta_Final)^2)/N

      Temp1<-x.LmAMSE%*%solve(t(x.LmAMSE)%*%x.LmAMSE)%*%t(x.LmAMSE)
      L1_r1 <- Tempsy_Var_Gam_Var_Loss*psych::tr(Temp1)
      L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.LmAMSE)] -
                                     F_Estimate_Full[c(idx.LmAMSE)]))^2)

      Loss_Subsample_LmAMSE_Pow[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      idx.LmAMSE->Sample.LmAMSE_Pow[[j]][[i+1]] # Model robust LmAMSE
    }
  }

  if(anyNA(beta_mMSE) || anyNA(beta_mVc) || anyNA(beta_LmAMSE) ||
     anyNA(beta_LmAMSE_LO) || anyNA(beta_LmAMSE_Pow) ){
    stop("There are NA or NaN values")
  }

  Full_SP<-cbind.data.frame(PI.mMSE,PI.mVc,PI.LmAMSE,PI.LmAMSE_LO,PI.LmAMSE_Pow)
  colnames(Full_SP)<-c("A-Optimality","L-Optimality","LmAMSE",paste0("LmAMSE Log Odds ",Alpha),
                       paste0("LmAMSE Power ",Alpha))

  Subsampling_Methods<-factor(c("A-Optimality","L-Optimality","LmAMSE",paste0("LmAMSE Log Odds ",Alpha),
                                paste0("LmAMSE Power ",Alpha)))

  # Beta Data
  Beta_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                              rbind(beta_mMSE,beta_mVc,beta_LmAMSE,
                                    do.call(rbind,beta_LmAMSE_LO),
                                    do.call(rbind,beta_LmAMSE_Pow)))

  colnames(Beta_Data)[-1]<-c("r2",paste0("Beta",0:(ncol(X)-1)))

  # Var Data
  Var_Data<-cbind.data.frame(Var_Epsilon,Var_LmAMSE_LO,Var_LmAMSE_Pow)

  colnames(Var_Data)<-c("r2","A-Optimality","L-Optimality","LmAMSE",
                           paste0("LmAMSE Log Odds ",Alpha),
                           paste0("LmAMSE Power ",Alpha))

  # Loss Subsample Data
  Loss_Subsample_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                        rbind(Loss_Subsample_mMSE,Loss_Subsample_mVc,
                                              Loss_Subsample_LmAMSE,
                                              do.call(rbind,Loss_Subsample_LmAMSE_LO),
                                              do.call(rbind,Loss_Subsample_LmAMSE_Pow)))
  colnames(Loss_Subsample_Data)[-1]<-c("r2","Variance","Bias.2","Loss")

  Loss_Subsample_Data[,-c(1,2)]<-Loss_Subsample_Data[,-c(1,2)]/r2

  all_r<-c(r1,r2)
  # Sample Data
  for(j in 1:length(Alpha)){
    names(Sample.LmAMSE_LO[[j]])<-names(Sample.LmAMSE_Pow[[j]])<-paste0("Alpha_",Alpha[j],"_",all_r)
  }

  names(Sample.mMSE)<-names(Sample.mVc)<-names(Sample.LmAMSE)<-c(r1,r2)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Variance_Epsilon_Estimates"=Var_Data,
            "Loss_Estimates"=Loss_Subsample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_LmAMSE"=Sample.LmAMSE,
            "Sample_LmAMSE_Log_Odds"=Sample.LmAMSE_LO,
            "Sample_LmAMSE_Power"=Sample.LmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","linear")
  return(ans)
}

#' Generate data for Generalised Linear Models under model misspecification scenario
#'
#' Function to simulate big data under Generalised Linear Models for the model misspecification scenario through
#' any misspecification type.
#'
#' @usage
#' GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)
#'
#' @param No_Of_Var               number of variables
#' @param Beta                    a vector for the model parameters, including the intercept
#' @param Var_Epsilon             variance value for the residuals
#' @param N                       the big data size
#' @param MisspecificationType    a character vector referring to different types of misspecification
#' @param family                  a character vector for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Big data for the Generalised Linear Models are generated by the "linear", "logistic" and "poisson"
#' regression types under model misspecification.
#'
#' We have limited the covariate data generation through uniform distribution of limits \eqn{(-1,1)}.
#'
#' Different type of misspecifications are "Type 1", "Type 2 Squared", "Type 2 Interaction",
#' "Type 3 Squared" or "Type 3 Interaction".
#'
#' @return
#' The output of \code{GenModelMissGLMdata} gives a list of
#'
#' \code{N} the big data size
#' \code{Beta} a list of outputs(real and estimated) for the beta values
#' \code{Variance_Epsilon} a list of outputs(real and estimated) for the variance epsilon
#' \code{Xbeta} a list of outputs(real and estimated) for the linear predictor
#' \code{f} a list of outputs(real and estimated) misspecification
#' \code{Real_Full_Data} a matrix for Y,X and f(x)
#' \code{Full_Data} a matrix for Y and X
#'
#'
#' @examples
#' No_Of_Var<-2; Beta<-c(-1,2,2,1); Var_Epsilon<-0.5; N<-10000;
#' MisspecificationType <- "Type 2 Squared"; family <- "linear"
#'
#' Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)
#'
#' No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
#' MisspecificationType <- "Type 2 Squared"; family <- "logistic"
#'
#' Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)
#'
#' No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
#' MisspecificationType <- "Type 2 Squared"; family <- "poisson"
#'
#' Results<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)
#'
#' @import stats
#' @importFrom gam s
#' @export
GenModelMissGLMdata<-function(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)
{
  if(any(is.na(c(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family))) |
     any(is.nan(c(No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType,family)))){
    stop("NA or Infinite or NAN values in the No_Of_Var,Beta,Var_Epsilon,N,MisspecificationType or family")
  }

  if(!any(family %in% c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if(!any(MisspecificationType == c("Type 1","Type 2 Squared","Type 2 Interaction",
                                    "Type 3 Squared","Type 3 Interaction"))){
    stop("Only the misspecification types 'Type 1', 'Type 2 Squared', 'Type 2 Interaction', \n 'Type 3 Squared', 'Type 3 Interaction' are allowed")
  }

  X_1<-replicate(No_Of_Var,stats::runif(n=N,min = -1,max = 1))
  X_Data <- cbind(X0=1,X_1);
  colnames(X_Data)[-1]<-paste0("X",1:ncol(X_Data[,-1]))

  if(MisspecificationType == "Type 1"){
    X_Data_Real <- cbind(X_Data,"f(x)"=0)
  }
  if(MisspecificationType == "Type 2 Squared"){
    Temp<-rowSums(X_1^2)
    X_Data_Real <- cbind(X_Data,"f(x)"=(Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2))
  }
  if(MisspecificationType == "Type 2 Interaction"){
    Temp<-Rfast::rowprods(X_1)
    X_Data_Real <- cbind(X_Data,"f(x)"=(Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2))
  }
  if(MisspecificationType == "Type 3 Squared"){
    Temp<-rowSums(X_1^2)
    beta_dist<-stats::runif(1,min = 0.75,max = 1.25)
    X_Data_Real <- cbind(X_Data,"f(x)"=beta_dist*(Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2))
  }
  if(MisspecificationType == "Type 3 Interaction"){
    Temp<-Rfast::rowprods(X_1)
    beta_dist<-stats::runif(1,min = 0.75,max = 1.25)
    X_Data_Real <- cbind(X_Data,"f(x)"=beta_dist*(Temp-mean(Temp))/sqrt(mean(Temp^2)-mean(Temp)^2))
  }

  if(family == "linear"){
    Y_Data <- X_Data_Real%*%Beta + stats::rnorm(n = N,mean = 0,sd = sqrt(Var_Epsilon))

    if(MisspecificationType == "Type 1"){
      beta.prop_Real<-solve(a=t(X_Data_Real[,-ncol(X_Data_Real)])%*%X_Data_Real[,-ncol(X_Data_Real)],
                            b=t(X_Data_Real[,-ncol(X_Data_Real)])%*%Y_Data)
      Xbeta_Final_Real<-as.vector(X_Data_Real[,-ncol(X_Data_Real)]%*%beta.prop_Real)
      Var.prop_Real<-sum((Y_Data-Xbeta_Final_Real)^2)/N
    } else {
      beta.prop_Real<-solve(a=t(X_Data_Real)%*%X_Data_Real,b=t(X_Data_Real)%*%Y_Data)
      Xbeta_Final_Real<-as.vector(X_Data_Real%*%beta.prop_Real)
      Var.prop_Real<-sum((Y_Data-Xbeta_Final_Real)^2)/N
    }

    beta.prop<-solve(a=t(X_Data)%*%X_Data,b=t(X_Data)%*%Y_Data)
    Xbeta_Final<-as.vector(X_Data%*%beta.prop)
    Var.prop<-sum((Y_Data-Xbeta_Final)^2)/N
  }

  if(family == "logistic"){
    Pi_Data<-1-1/(1+exp(X_Data_Real%*%Beta))
    Y_Data <- stats::rbinom(N,1,Pi_Data)

    if(MisspecificationType == "Type 1"){
      parameter.propfits<-.getMLE(x=X_Data_Real[,-ncol(X_Data_Real)], y=Y_Data, w=rep(N,N))
      beta.prop_Real<-parameter.propfits$par
      Xbeta_Final_Real<-as.vector(X_Data_Real[,-ncol(X_Data_Real)]%*%beta.prop_Real)
    } else {
      parameter.propfits<-.getMLE(x=X_Data_Real, y=Y_Data, w=rep(N,N))
      beta.prop_Real<-parameter.propfits$par
      Xbeta_Final_Real<-as.vector(X_Data_Real%*%beta.prop_Real)
    }
    parameter.propfits<-.getMLE(x=X_Data, y=Y_Data, w=rep(N,N))
    beta.prop<-parameter.propfits$par
    Xbeta_Final<-as.vector(X_Data%*%beta.prop)
  }

  if(family == "poisson"){
    Lambda_Data<-exp(X_Data_Real%*%Beta)
    Y_Data <- stats::rpois(N,Lambda_Data)

    if(MisspecificationType == "Type 1"){
      Temp_Data<-cbind.data.frame(Y_Data,X_Data_Real[,-ncol(X_Data_Real)])
      colnames(Temp_Data)<-c("Y",paste0("X",0:(ncol(Temp_Data[,-1])-1)))
      parameter.propfits<-stats::glm(Y ~ . -1, data = Temp_Data, family = "poisson")
      beta.prop_Real<-parameter.propfits$coefficients
      Xbeta_Final_Real<-as.vector(X_Data_Real[,-ncol(X_Data_Real)]%*%beta.prop_Real)
    } else {
      Temp_Data<-cbind.data.frame(Y_Data,X_Data_Real)
      colnames(Temp_Data)<-c("Y",paste0("X",0:(ncol(Temp_Data[,-1])-1)))
      parameter.propfits<-stats::glm(Y ~ .-1, data = Temp_Data, family = "poisson")
      beta.prop_Real<-parameter.propfits$coefficients
      Xbeta_Final_Real<-as.vector(X_Data_Real%*%beta.prop_Real)
    }

    Temp_Data<-cbind.data.frame(Y_Data,X_Data)
    colnames(Temp_Data)<-c("Y",paste0("X",0:(ncol(Temp_Data[,-1])-1)))
    parameter.propfits<-stats::glm(Y ~ .-1, data = Temp_Data, family="poisson")
    beta.prop<-parameter.propfits$coefficients
    Xbeta_Final<-as.vector(X_Data%*%beta.prop)
  }

  my_formula<-stats::as.formula(paste("Y ~ ",paste(paste0("s(X",1:ncol(X_Data[,-1]),")"),collapse = " + "),"+",
                                      paste(paste0("s(",paste0(colnames(X_Data[,-1]),collapse = "*"),")"),
                                            collapse = " + ")))

  if(family == "linear"){
    # calculate f_hat
    Assumed_Data<-data.frame(Y=Y_Data,X_Data)
    fit_GAM<-gam::gam(formula = my_formula,data=Assumed_Data)
    Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X_Data))
    f_estimate<-Xbeta_GAM - Xbeta_Final
    Var_GAM.prop<-sum((Y_Data-Xbeta_GAM)^2)/N

    Real_Model_Data<-cbind(Y=Y_Data,X_Data_Real)
    Assumed_Model_Data<-cbind(Y=Y_Data,X_Data)

    colnames(Real_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])),"f(x)")
    colnames(Assumed_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])))

    Outputs<-list("N"=N,
                  "Beta"=list("Real"=Beta,"Real_Estimate"=t(beta.prop_Real),"Estimate"=t(beta.prop)),
                  "Variance_Epsilon"=list("Real"=Var_Epsilon,"Real_Estimate"=Var.prop_Real,
                                          "Real_GAM"=Var_GAM.prop,"Estimate"=Var.prop),
                  "Xbeta"=list("Real_Estimate"=Xbeta_Final_Real,"Real_GAM"=Xbeta_GAM,"Estimate"=Xbeta_Final),
                  "f"=list("Real"=X_Data_Real[,ncol(X_Data_Real)],"Real_GAM"=f_estimate),
                  "Real_Full_Data"=Real_Model_Data,
                  "Full_Data"=Assumed_Model_Data)
    class(Outputs)<-c("ModelMisspecified","linear")
  }
  if(family == "logistic"){
    # calculate f_hat
    Assumed_Data<-data.frame(Y=Y_Data,X_Data)
    fit_GAM<-gam::gam(formula = my_formula,data=Assumed_Data,family = "binomial")
    Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X_Data))
    f_estimate<-Xbeta_GAM - Xbeta_Final

    Real_Model_Data<-cbind(Y=Y_Data,X_Data_Real)
    Assumed_Model_Data<-cbind(Y=Y_Data,X_Data)

    colnames(Real_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])),"f(x)")
    colnames(Assumed_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])))

    Outputs<-list("N"=N,
                  "Beta"=list("Real"=Beta,"Real_Estimate"=beta.prop_Real,"Estimate"=beta.prop),
                  "Xbeta"=list("Real_Estimate"=Xbeta_Final_Real,"Real_GAM"=Xbeta_GAM,"Estimate"=Xbeta_Final),
                  "f"=list("Real"=X_Data_Real[,ncol(X_Data_Real)],"Real_GAM"=f_estimate),
                  "Real_Full_Data"=Real_Model_Data,
                  "Full_Data"=Assumed_Model_Data)
    class(Outputs)<-c("ModelMisspecified","logistic")
  }
  if(family == "poisson"){
    # calculate f_hat
    Assumed_Data<-data.frame(Y=Y_Data,X_Data)
    fit_GAM<-gam::gam(formula = my_formula,data=Assumed_Data,family = "poisson")
    Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X_Data))
    f_estimate<-Xbeta_GAM - Xbeta_Final

    Real_Model_Data<-cbind(Y=Y_Data,X_Data_Real)
    Assumed_Model_Data<-cbind(Y=Y_Data,X_Data)

    colnames(Real_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])),"f(x)")
    colnames(Assumed_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])))

    Outputs<-list("N"=N,
                  "Beta"=list("Real"=Beta,"Real_Estimate"=beta.prop_Real,"Estimate"=beta.prop),
                  "Xbeta"=list("Real_Estimate"=Xbeta_Final_Real,"Real_GAM"=Xbeta_GAM,"Estimate"=Xbeta_Final),
                  "f"=list("Real"=X_Data_Real[,ncol(X_Data_Real)],"Real_GAM"=f_estimate),
                  "Real_Full_Data"=Real_Model_Data,
                  "Full_Data"=Assumed_Model_Data)
    class(Outputs)<-c("ModelMisspecified","poisson")
  }
  return(Outputs)
}
