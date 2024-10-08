#' Subsampling under linear regression for a potentially misspecified model
#'
#' Using this function sample from big data under linear regression for a potentially misspecified model.
#' Subsampling probabilities are obtained based on the A- and L- optimality criteria
#' with the RLmAMSE (Reduction of Loss by minimizing the Average Mean Squared Error).
#'
#' @usage
#' modelMissLinSub(r1,r2,Y,X,N,Alpha,proportion)
#'
#' @param r1                sample size for initial random sampling
#' @param r2                sample size for optimal sampling
#' @param Y                 response data or Y
#' @param X                 covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N                 size of the big data
#' @param Alpha             scaling factor when using Log Odds or Power functions to magnify the probabilities
#' @param proportion        a proportion of the big data is used to help estimate AMSE values from the subsamples
#'
#' @details
#' Two stage subsampling algorithm for big data under linear regression for potential model misspecification.
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
#' If \eqn{\alpha > 1} for the scaling factor is not satisfied an error message will be produced.
#'
#' If proportion is not in the region of \eqn{(0,1]} an error message will be produced.
#'
#' @return
#' The output of \code{modelMissLinSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{Variance_Epsilon_Estimates} matrix of estimated variance for epsilon after subsampling
#'
#' \code{AMSE_Estimates} matrix of estimated AMSE values after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal samples obtained based on A-Optimality criteria
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal samples obtained based on L-Optimality criteria
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
#' r1<-300; r2<-rep(100*c(6,9),50);
#' Original_Data<-Full_Data$Complete_Data[,-ncol(Full_Data$Complete_Data)];
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#' \dontrun{
#' Results<-modelMissLinSub(r1 = r1, r2 = r2,
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
#' @importFrom psych tr
#' @export
modelMissLinSub <- function(r1,r2,Y,X,N,Alpha,proportion){
  if(any(is.na(c(r1,r2,N,Alpha,proportion))) | any(is.nan(c(r1,r2,N,Alpha,proportion)))){
    stop("NA or Infinite or NAN values in the r1,r2,N,Alpha or proportion")
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

  if(proportion >1 | proportion <=0){
    stop("Proportion should be a value higher than zero and less than or equal one")
  }
  if(proportion >= 0.5){
    message("50% or >=50% of the big data is used to help find AMSE for the subsamples, \nthis could take some time.")
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

  if(proportion*N != r1){
    idx.proportion <- sample(1:N, ceiling(proportion*N), T)
    Y_proportion<-Y[idx.proportion]
    X_proportion<-X[idx.proportion,]

    Proportion_Data<-data.frame(Y=Y_proportion,X_proportion)
    fit_GAM_Proportion<-gam::gam(my_formula,data=Proportion_Data)
    Xbeta_GAM_Proportion<-gam::predict.Gam(fit_GAM_Proportion,newdata = data.frame(X))

    beta_proportion<-solve(a=t(X_proportion)%*%X_proportion,b=t(X_proportion)%*%Y_proportion)
    Xbeta_proportion<-X%*%beta_proportion

    Var_GAM_Full<-sum((Y-Xbeta_GAM_Proportion)^2)/N
    Var_Full<-sum((Y-Xbeta_proportion)^2)/N
    F_Estimate_Full<-Xbeta_GAM_Proportion-Xbeta_proportion
  }
  else{
    Var_GAM_Full<-Var_GAM.prop ; Var_Full<-Var.prop ; F_Estimate_Full<-f_estimate
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

  # RLmAMSE
  PI.RLmAMSE <- (max(L_All_Final) - L_All_Final) / sum(max(L_All_Final) - L_All_Final)

  # RLmAMSE LO
  PI.RLmAMSE_LO<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    (1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1)/
      sum((1+exp(-Alpha*log(PI.RLmAMSE/(1-PI.RLmAMSE))))^(-1))},PI.RLmAMSE)

  # RLmAMSE Pow
  PI.RLmAMSE_Pow<-apply(as.matrix(Alpha),1,function(Alpha,PI.RLmAMSE){
    PI.RLmAMSE^Alpha/sum(PI.RLmAMSE^Alpha)},PI.RLmAMSE)

  Var_Epsilon<-matrix(nrow=length(r2),ncol=4)

  # For the Model with already available subsampling probabilities - mVc
  beta_mVc<-matrix(nrow = length(r2),ncol = ncol(X)+1);
  AMSE_Sample_mVc<-matrix(nrow = length(r2),ncol = 4) ; Sample.mVc<-list();

  # For the Model with already available subsampling probabilities - mMSE
  beta_mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1);
  AMSE_Sample_mMSE<-matrix(nrow = length(r2),ncol = 4) ; Sample.mMSE<-list();

  # For the Real and Model with model robust subsampling probabilities RLmAMSE
  beta_RLmAMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 );
  AMSE_Sample_RLmAMSE<-matrix(nrow = length(r2),ncol = 4) ; Sample.RLmAMSE<-list()

  # For the Model with model robust subsampling probabilities RLmAMSE LO
  beta_RLmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(r2),ncol(X)+1)),simplify = FALSE) ;
  AMSE_Sample_RLmAMSE_LO<-replicate(length(Alpha),array(dim = c(length(r2),4)),simplify = FALSE);
  Sample.RLmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));
  Var_RLmAMSE_LO<-matrix(nrow = length(r2),ncol=length(Alpha))

  # For the Model with model robust subsampling probabilities RLmAMSE Pow
  beta_RLmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(r2),ncol(X)+1)),simplify = FALSE) ;
  AMSE_Sample_RLmAMSE_Pow<-replicate(length(Alpha),array(dim = c(length(r2),4)),simplify = FALSE);
  Sample.RLmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));
  Var_RLmAMSE_Pow<-matrix(nrow = length(r2),ncol=length(Alpha))

  Tempsy_Var_Gam_Var_AMSE<-Var_GAM_Full*Var_Full^(-2)
  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.RLmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.RLmAMSE_LO[[j]][[1]]<-Sample.RLmAMSE_Pow[[j]][[1]]<-idx.prop
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
    L1_r1 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.mVc,idx.prop)] -
                                   F_Estimate_Full[c(idx.mVc, idx.prop)]))^2)

    AMSE_Sample_mVc[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

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
    L1_r1 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.mMSE,idx.prop)]-
                                   F_Estimate_Full[c(idx.mMSE, idx.prop)]))^2)

    AMSE_Sample_mMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    idx.mMSE->Sample.mMSE[[i+1]]

    # RLmAMSE
    idx.RLmAMSE <- sample(1:N, r2[i], T, PI.RLmAMSE)

    x.RLmAMSE <- X[c(idx.RLmAMSE),]
    y.RLmAMSE <- Y[c(idx.RLmAMSE)]
    w.RLmAMSE <- c(1/PI.RLmAMSE[idx.RLmAMSE])

    pi4_r<-sqrt(r2[i]*w.RLmAMSE^(-1))
    X_r4<-x.RLmAMSE/pi4_r
    Y_r4<-y.RLmAMSE/pi4_r
    beta_RLmAMSE[i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
    Xbeta_Final<-as.vector(X%*%beta_RLmAMSE[i,-1])
    Var_Epsilon[i,4]<-sum((Y-Xbeta_Final)^2)/N

    Temp1<-x.RLmAMSE%*%solve(t(x.RLmAMSE)%*%x.RLmAMSE)%*%t(x.RLmAMSE)
    L1_r1 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
    L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE)] -
                                   F_Estimate_Full[c(idx.RLmAMSE)]))^2)

    AMSE_Sample_RLmAMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    idx.RLmAMSE->Sample.RLmAMSE[[i+1]] # Model robust RLmAMSE

    for(j in 1:length(Alpha))
    {
      # RLmAMSE Log Odds
      idx.RLmAMSE <- sample(1:N, r2[i], T, PI.RLmAMSE_LO[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE)]
      w.RLmAMSE <- c(1/PI.RLmAMSE_LO[idx.RLmAMSE,j])

      pi4_r<-sqrt(r2[i]*w.RLmAMSE^(-1))
      X_r4<-x.RLmAMSE/pi4_r
      Y_r4<-y.RLmAMSE/pi4_r
      beta_RLmAMSE_LO[[j]][i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
      Xbeta_Final<-as.vector(X%*%beta_RLmAMSE_LO[[j]][i,-1])
      Var_RLmAMSE_LO[i,j]<-sum((Y-Xbeta_Final)^2)/N

      Temp1<-x.RLmAMSE%*%solve(t(x.RLmAMSE)%*%x.RLmAMSE)%*%t(x.RLmAMSE)
      L1_r1 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
      L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE)]-
                                     F_Estimate_Full[c(idx.RLmAMSE)]))^2)

      AMSE_Sample_RLmAMSE_LO[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      idx.RLmAMSE->Sample.RLmAMSE_LO[[j]][[i+1]] # Model robust RLmAMSE

      # RLmAMSE Power
      idx.RLmAMSE <- sample(1:N, r2[i], T, PI.RLmAMSE_Pow[,j])

      x.RLmAMSE <- X[c(idx.RLmAMSE),]
      y.RLmAMSE <- Y[c(idx.RLmAMSE)]
      w.RLmAMSE <- c(1/PI.RLmAMSE_Pow[idx.RLmAMSE,j])

      pi4_r<-sqrt(r2[i]*w.RLmAMSE^(-1))
      X_r4<-x.RLmAMSE/pi4_r
      Y_r4<-y.RLmAMSE/pi4_r
      beta_RLmAMSE_Pow[[j]][i,]<-c(r2[i],solve(a=t(X_r4)%*%X_r4,b=t(X_r4)%*%Y_r4))
      Xbeta_Final<-as.vector(X%*%beta_RLmAMSE_Pow[[j]][i,-1])
      Var_RLmAMSE_Pow[i,j]<-sum((Y-Xbeta_Final)^2)/N

      Temp1<-x.RLmAMSE%*%solve(t(x.RLmAMSE)%*%x.RLmAMSE)%*%t(x.RLmAMSE)
      L1_r1 <- Tempsy_Var_Gam_Var_AMSE*psych::tr(Temp1)
      L2_r1 <- sum((Var_Full^(-1)*(Temp1%*%F_Estimate_Full[c(idx.RLmAMSE)] -
                                     F_Estimate_Full[c(idx.RLmAMSE)]))^2)

      AMSE_Sample_RLmAMSE_Pow[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      idx.RLmAMSE->Sample.RLmAMSE_Pow[[j]][[i+1]] # Model robust RLmAMSE
    }
  }

  if(anyNA(beta_mMSE) || anyNA(beta_mVc) || anyNA(beta_RLmAMSE) ||
     anyNA(beta_RLmAMSE_LO) || anyNA(beta_RLmAMSE_Pow) ){
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

  colnames(Beta_Data)[-1]<-c("r2",paste0("Beta",0:(ncol(X)-1)))

  # Var Data
  Var_Data<-cbind.data.frame(Var_Epsilon,Var_RLmAMSE_LO,Var_RLmAMSE_Pow)

  colnames(Var_Data)<-c("r2","A-Optimality","L-Optimality","RLmAMSE",
                           paste0("RLmAMSE Log Odds ",Alpha),
                           paste0("RLmAMSE Power ",Alpha))

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
            "Variance_Epsilon_Estimates"=Var_Data,
            "AMSE_Estimates"=AMSE_Sample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_RLmAMSE"=Sample.RLmAMSE,
            "Sample_RLmAMSE_Log_Odds"=Sample.RLmAMSE_LO,
            "Sample_RLmAMSE_Power"=Sample.RLmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","linear")
  return(ans)
}
