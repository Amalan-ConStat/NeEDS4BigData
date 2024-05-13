#' Model misspecified subsampling under Poisson regression
#'
#' Using this function subsample from big data under Poisson regression when we assumed the model that
#' describes the data is misspecified. Subsampling probabilities are obtained based on the A- and L-
#' optimality criterions with the LmAMSE (Loss on mean of asymptotic mean squared error).
#'
#' @usage
#' modelMissPoiSub(r1,r2,Y,X,N,Alpha,Beta_Estimate_Full,F_Estimate_Full)
#'
#' @param r1                 subsample size for initial random sampling
#' @param r2                 subsample size for optimal sampling
#' @param Y                  response data or Y
#' @param X                  covariate data or X matrix that has all the covariates (first column is for the intercept)
#' @param N                  size of the big data
#' @param Alpha              scaling factor when using Log Odds or Power functions to magnify the probabilities
#' @param Beta_Estimate_Full estimate of Beta after fitting the Poisson model
#' @param F_Estimate_Full    estimate of f that is the difference of linear predictor on GAM and Poisson model
#'
#' @details
#' Two stage subsampling algorithm for big data under Poisson regression for potential model misspecification.
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
#' If \eqn{\alpha > 1} for the scaling vector is not satisfied an error message will be produced.
#'
#' @return
#' The output of \code{modelMissPoiSub} gives a list of
#'
#' \code{Beta_Estimates} estimated model parameters after subsampling
#'
#' \code{Loss_Estimates} matrix of estimated LmAMSE values after subsampling
#'
#' \code{Sample_A-Optimality} list of indexes for the initial and optimal subsamples obtained based on A-Optimality criterion
#'
#' \code{Sample_L-Optimality} list of indexes for the initial and optimal subsamples obtained based on L-Optimality criterion
#'
#' \code{Sample_LmAMSE} list of indexes for the optimal subsamples obtained based on LmAMSE
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
#' No_Of_Var<-2; Beta<-c(-1,2,2,1); N<-10000;
#' MisspecificationType <- "Type 2 Squared"; family <- "poisson"
#'
#' Full_Data<-GenModelMissGLMdata(No_Of_Var,Beta,Var_Epsilon=NULL,N,MisspecificationType,family)
#'
#' r1<-300; r2<-100*c(6,9); Original_Data<-Full_Data$Full_Data;
#'
#' # cl <- parallel::makeCluster(4)
#' # doParallel::registerDoParallel(cl)
#'
#' Results<-modelMissPoiSub(r1 = r1, r2 = r2,
#'                          Y = as.matrix(Original_Data[,1]),
#'                          X = as.matrix(Original_Data[,-1]),
#'                          N = Full_Data$N,
#'                          Alpha = 10,
#'                          Beta_Estimate_Full = Full_Data$Beta$Estimate,
#'                          F_Estimate_Full = Full_Data$f$Real_GAM)
#'
#' # parallel::stopCluster(cl)
#'
#' plot_Beta(Results)
#' plot_LmAMSE(Results)
#'
#' @importFrom Rdpack reprompt
#' @import stats
#' @import foreach
#' @importFrom gam s
#' @importFrom Rfast rowprods
#' @importFrom psych tr
#' @export
modelMissPoiSub <- function(r1,r2,Y,X,N,Alpha,Beta_Estimate_Full,F_Estimate_Full){
  if(any(is.na(c(r1,r2,N,Alpha,Beta_Estimate_Full))) | any(is.nan(c(r1,r2,N,Alpha,Beta_Estimate_Full)))){
    stop("NA or Infinite or NAN values in the r1,r2,N,Alpha or Beta_Estimate_Full")
  }

  if((N != nrow(X)) | (N != nrow(Y)) | nrow(X) != nrow(Y)){
    stop("The big data size N is not the same as of the size of X or Y")
  }

  if((length(F_Estimate_Full) != nrow(X)) | (length(F_Estimate_Full) != nrow(Y)) | N != length(F_Estimate_Full)){
    stop("The big data size N is not the same as of the size of F_Estimate_Full")
  }

  if(any(is.na(cbind(Y,X,F_Estimate_Full))) | any(is.nan(cbind(Y,X,F_Estimate_Full)))){
    stop("NA or Infinite or NAN values in the Y or X or F_Estimate_Full")
  }

  if(any((2*r1) > r2)){
    stop("2*r1 cannot be greater than r2 at any point")
  }

  if(Alpha <= 1 | length(Alpha) > 1){
    stop("Scaling factor alpha is not greater than one or the length is more than one")
  }

  PI.prop <- rep(1/N,N)
  idx.prop <- sample(1:N, r1, T,PI.prop)

  x.prop <- X[idx.prop,]
  y.prop <- Y[idx.prop]

  pinv.prop <- rep(N,r1)

  fit.prop <- stats::glm(y.prop~x.prop-1,family="poisson")
  beta.prop <- fit.prop$coefficients

  if (anyNA(beta.prop)){
    stop("There are NA or NaN values in the model parameters")
  }
  Lambda.prop  <- exp(X %*% beta.prop)

  my_formula<-stats::as.formula(paste("Y ~ ",paste(paste0("gam::s(X",1:ncol(x.prop[,-1]),")"),collapse = " + "),"+",
                                      paste(paste0("gam::s(",paste0(colnames(x.prop[,-1]),collapse = "*"),")"),
                                            collapse = " + ")))

  #calculate f_hat
  Assumed_Data<-data.frame(Y=y.prop,x.prop)
  fit_GAM<-gam::gam(my_formula,data=Assumed_Data,family = "poisson")
  Xbeta_GAM<-gam::predict.Gam(fit_GAM,newdata = data.frame(X))
  Lambda_GAM<-exp(Xbeta_GAM)
  f_estimate<-Xbeta_GAM - X %*% beta.prop

  if(is.null(Beta_Estimate_Full) || is.null(F_Estimate_Full)){
    Beta_Estimate_Full<- beta.prop ; F_Estimate_Full<-f_estimate
    message("Beta_Estimate_Full and F_Estimate_Full from the initial sample is used.")
  }

  ## mVC
  PI.mVc <- sqrt((Y - Lambda.prop)^2 * rowSums(X^2))
  PI.mVc <- PI.mVc / sum(PI.mVc)

  # mMSE
  lambda.prop <- Lambda.prop[idx.prop]
  w.prop <- lambda.prop
  W.prop <- solve(t(x.prop) %*% (x.prop * w.prop * pinv.prop))

  PI.mMSE <- sqrt((Y - Lambda.prop)^2 * rowSums((X%*%W.prop)^2))
  PI.mMSE <- PI.mMSE / sum(PI.mMSE)

  a<-NULL
  L_All<-foreach::foreach(a=1:N,.combine = rbind,.packages = "psych") %dopar% {
    X_r1<-rbind(X[idx.prop,],X[a,])
    lambda_r1<-exp(X_r1%*%beta.prop)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(t(X_r1) %*% (X_r1 * W_r1))
    Temp1<-(W_r1*X_r1)%*%H_r1

    f_r1<-c(f_estimate[idx.prop],f_estimate[a])
    lambda_Tr1<-exp((X_r1 %*% beta.prop) + f_r1)
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <-(t(X_r1)  %*% (X_r1 * W_Tr1))
    b_r1 <-(t(X_r1) %*% (lambda_Tr1-lambda_r1))
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((X_r1%*%H_r1%*%b_r1)-f_r1))^2)

    c(L1_r1+L2_r1)
  }

  lambda_r1<-exp(X[idx.prop,]%*%beta.prop)
  W_r1<-as.vector(lambda_r1)
  H_r1 <-solve(t(X[idx.prop,]) %*% (X[idx.prop,]*W_r1))
  Temp1<-(W_r1*X[idx.prop,])%*%H_r1

  lambda_Tr1<-exp((X[idx.prop,] %*% beta.prop) + f_estimate[idx.prop])
  W_Tr1<-as.vector(lambda_Tr1)
  H_Tr1 <-(t(X[idx.prop,]) %*% (X[idx.prop,]*W_Tr1))
  b_r1 <-(t(X[idx.prop,]) %*% (lambda_Tr1-lambda_r1))
  L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
  L2_r1 <- sum((W_r1*((X[idx.prop,]%*%H_r1%*%b_r1)-f_estimate[idx.prop]))^2)

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

  # For the Model with already available Sub-sampling probabilities - mMSE and mVc
  beta_mVc<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  Loss_Subsample_mVc<-matrix(nrow = length(r2),ncol = 4); Sample.mVc<-list();

  beta_mMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  Loss_Subsample_mMSE<-matrix(nrow = length(r2),ncol = 4); Sample.mMSE<-list();

  # For the Model with model robust Sub-sampling probabilities LmAMSE
  beta_LmAMSE<-matrix(nrow = length(r2),ncol = ncol(X)+1 )
  Loss_Subsample_LmAMSE<-matrix(nrow = length(r2),ncol = 4); Sample.LmAMSE<-list();

  # For the Model with model robust Sub-sampling probabilities LmAMSE LO
  beta_LmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r2),ncol(X)+1)),simplify = FALSE)
  Loss_Subsample_LmAMSE_LO<-replicate(length(Alpha),array(dim=c(length(r2),4)),simplify = FALSE)
  Sample.LmAMSE_LO<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));

  # For the Model with model robust Sub-sampling probabilities LmAMSE Pow
  beta_LmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r2),ncol(X)+1)),simplify = FALSE)
  Loss_Subsample_LmAMSE_Pow<-replicate(length(Alpha),array(dim=c(length(r2),4)),simplify = FALSE)
  Sample.LmAMSE_Pow<-replicate(length(Alpha),list(rep(list(NA),length(r2)+1)));

  Sample.mMSE[[1]]<-Sample.mVc[[1]]<-Sample.LmAMSE[[1]]<-idx.prop

  for (j in 1:length(Alpha)) {
    Sample.LmAMSE_LO[[j]][[1]]<-Sample.LmAMSE_Pow[[j]][[1]]<-idx.prop
  }

  message("Step 1 of the algorithm completed.\n")

  for (i in 1:length(r2))
  {
    # mVc
    idx.mVc <- sample(1:N, r2[i]-r1, T, PI.mVc)

    x.mVc <- X[c(idx.mVc, idx.prop),]
    y.mVc <- Y[c(idx.mVc, idx.prop)]

    fit.mVc <- stats::glm(y.mVc~x.mVc-1,family="poisson",weights = c(1 / PI.mVc[idx.mVc], pinv.prop))

    beta_mVc[i,] <- c(r2[i],fit.mVc$coefficients)

    idx.mVc->Sample.mVc[[i+1]]

    lambda_r1<-exp(x.mVc%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(t(x.mVc) %*% (x.mVc * W_r1))
    Temp1<-(W_r1 * x.mVc)%*%H_r1

    lambda_Tr1<-exp((x.mVc %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.mVc, idx.prop)])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <-(t(x.mVc) %*% (x.mVc * W_Tr1))
    b_r1 <-(t(x.mVc) %*% (lambda_Tr1-lambda_r1))
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.mVc%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.mVc, idx.prop)]))^2)

    Loss_Subsample_mVc[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    # mMSE
    idx.mMSE <- sample(1:N, r2[i]-r1, T, PI.mMSE)

    x.mMSE <- X[c(idx.mMSE, idx.prop),]
    y.mMSE <- Y[c(idx.mMSE, idx.prop)]

    fit.mMSE <- stats::glm(y.mMSE~x.mMSE-1,family = "poisson", weights =c(1 / PI.mMSE[idx.mMSE], pinv.prop) )

    beta_mMSE[i,] <- c(r2[i],fit.mMSE$coefficients)

    idx.mMSE->Sample.mMSE[[i+1]]

    lambda_r1<-exp(x.mMSE%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(t(x.mMSE) %*% (x.mMSE * W_r1))
    Temp1<-(W_r1*x.mMSE)%*%H_r1

    lambda_Tr1<-exp((x.mMSE %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.mMSE, idx.prop)])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <-(t(x.mMSE) %*% (x.mMSE * W_Tr1))
    b_r1 <-(t(x.mMSE) %*% (lambda_Tr1-lambda_r1))
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.mMSE%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.mMSE, idx.prop)]))^2)

    Loss_Subsample_mMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    # LmAMSE
    idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE)

    x.LmAMSE <- X[c(idx.LmAMSE),]
    y.LmAMSE <- Y[c(idx.LmAMSE)]

    fit.LmAMSE <- stats::glm(y.LmAMSE~x.LmAMSE-1,family="poisson",weights = c(1 / PI.LmAMSE[idx.LmAMSE]))

    beta_LmAMSE[i,] <- c(r2[i],fit.LmAMSE$coefficients)

    idx.LmAMSE->Sample.LmAMSE[[i+1]]

    lambda_r1<-exp(x.LmAMSE%*%Beta_Estimate_Full)
    W_r1<-as.vector(lambda_r1)
    H_r1 <-solve(t(x.LmAMSE) %*% (x.LmAMSE * W_r1))
    Temp1<-(W_r1*x.LmAMSE)%*%H_r1

    lambda_Tr1<-exp((x.LmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.LmAMSE)])
    W_Tr1<-as.vector(lambda_Tr1)
    H_Tr1 <-(t(x.LmAMSE) %*% (x.LmAMSE * W_Tr1))
    b_r1 <-(t(x.LmAMSE) %*% (lambda_Tr1-lambda_r1))
    L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
    L2_r1 <- sum((W_r1*((x.LmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.LmAMSE)]))^2)

    Loss_Subsample_LmAMSE[i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

    for (j in 1:length(Alpha))
    {
      # LmAMSE Log Odds
      idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE_LO[,j])

      x.LmAMSE <- X[c(idx.LmAMSE),]
      y.LmAMSE <- Y[c(idx.LmAMSE)]

      fit.LmAMSE <- stats::glm(y.LmAMSE~x.LmAMSE-1,family="poisson",weights = c(1 / PI.LmAMSE_LO[idx.LmAMSE,j]))

      beta_LmAMSE_LO[[j]][i,] <- c(r2[i],fit.LmAMSE$coefficients)

      idx.LmAMSE->Sample.LmAMSE_LO[[j]][[i+1]]

      lambda_r1<-exp(x.LmAMSE%*%Beta_Estimate_Full)
      W_r1<-as.vector(lambda_r1)
      H_r1 <-solve(t(x.LmAMSE) %*% (x.LmAMSE * W_r1))
      Temp1<-(W_r1*x.LmAMSE)%*%H_r1

      lambda_Tr1<-exp((x.LmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.LmAMSE)])
      W_Tr1<-as.vector(lambda_Tr1)
      H_Tr1 <-(t(x.LmAMSE) %*% (x.LmAMSE * W_Tr1))
      b_r1 <-(t(x.LmAMSE) %*% (lambda_Tr1-lambda_r1))
      L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
      L2_r1 <- sum((W_r1*((x.LmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.LmAMSE)]))^2)

      Loss_Subsample_LmAMSE_LO[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)

      # LmAMSE Power
      idx.LmAMSE <- sample(1:N, r2[i], T, PI.LmAMSE_Pow[,j])

      x.LmAMSE <- X[c(idx.LmAMSE),]
      y.LmAMSE <- Y[c(idx.LmAMSE)]

      fit.LmAMSE <- stats::glm(y.LmAMSE~x.LmAMSE-1,family="poisson",weights = c(1 / PI.LmAMSE_Pow[idx.LmAMSE,j]))

      beta_LmAMSE_Pow[[j]][i,] <- c(r2[i],fit.LmAMSE$coefficients)

      idx.LmAMSE->Sample.LmAMSE_Pow[[j]][[i+1]]

      lambda_r1<-exp(x.LmAMSE%*%Beta_Estimate_Full)
      W_r1<-as.vector(lambda_r1)
      H_r1 <-solve(t(x.LmAMSE) %*% (x.LmAMSE * W_r1))
      Temp1<-(W_r1*x.LmAMSE)%*%H_r1

      lambda_Tr1<-exp((x.LmAMSE %*% Beta_Estimate_Full) + F_Estimate_Full[c(idx.LmAMSE)])
      W_Tr1<-as.vector(lambda_Tr1)
      H_Tr1 <-(t(x.LmAMSE) %*% (x.LmAMSE * W_Tr1))
      b_r1 <-(t(x.LmAMSE) %*% (lambda_Tr1-lambda_r1))
      L1_r1 <- psych::tr(Temp1%*%H_Tr1%*%t(Temp1))
      L2_r1 <- sum((W_r1*((x.LmAMSE%*%H_r1%*%b_r1) - F_Estimate_Full[c(idx.LmAMSE)]))^2)

      Loss_Subsample_LmAMSE_Pow[[j]][i,]<-c(r2[i],L1_r1,L2_r1,L1_r1+L2_r1)
    }
  }

  if(anyNA(beta_mVc) || anyNA(beta_mMSE) || anyNA(beta_LmAMSE) || anyNA(beta_LmAMSE_LO) || anyNA(beta_LmAMSE_Pow))
  {
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

  # Loss Subsample Data
  Loss_Subsample_Data<-cbind.data.frame("Method"=rep(Subsampling_Methods,each=length(r2)),
                                        rbind(Loss_Subsample_mMSE,Loss_Subsample_mVc,
                                              Loss_Subsample_LmAMSE,
                                              do.call(rbind,Loss_Subsample_LmAMSE_LO),
                                              do.call(rbind,Loss_Subsample_LmAMSE_Pow)))
  colnames(Loss_Subsample_Data)[-1]<-c("r2","Variance","Bias.2","Loss")

  Loss_Subsample_Data[,-c(1,2)]<-Loss_Subsample_Data[,-c(1,2)]/r2

  # Sample Data
  for(j in 1:length(Alpha)){
    names(Sample.LmAMSE_LO[[j]])<-names(Sample.LmAMSE_Pow[[j]])<-c(length(r2)+1)
  }

  names(Sample.mMSE)<-names(Sample.mVc)<-names(Sample.LmAMSE)<-c(r1,r2)
  names(Sample.LmAMSE_LO)<-names(Sample.LmAMSE_Pow)<-paste0("Alpha_",Alpha)

  message("Step 2 of the algorithm completed.")

  ans<-list("Beta_Estimates"=Beta_Data,
            "Loss_Estimates"=Loss_Subsample_Data,
            "Sample_A-Optimality"=Sample.mMSE,
            "Sample_L-Optimality"=Sample.mVc,
            "Sample_LmAMSE"=Sample.LmAMSE,
            "Sample_LmAMSE_Log_Odds"=Sample.LmAMSE_LO,
            "Sample_LmAMSE_Power"=Sample.LmAMSE_Pow,
            "Subsampling_Probability"=Full_SP)
  class(ans)<-c("ModelMisspecified","poisson")
  return(ans)
}
