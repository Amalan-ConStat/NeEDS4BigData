#' Electric consumption data
#'
#' Hebrail and Berard (2012) described data which contains 2,049,280 completed measurements for a house
#' located at Sceaux, France between December 2006 and November 2010.
#' The log scale minute-averaged current intensity is selected as the response and
#' the covariates are active electrical energy (watt-hour) in the kitchen, the laundry room,
#' and electric water-heater and an air-conditioner.
#'
#' @format A data frame with 4 columns and 2,049,280 rows.
#' \describe{
#' \item{\code{Intensity}}{Minute-averaged current intensity}
#' \item{\code{EE_Kitchen}}{Active electrical energy (watt-hour) in the kitchen}
#' \item{\code{EE_Laundry}}{Active electrical energy (watt-hour) in the laundry room}
#' \item{\code{EE_WH_AC}}{Active electrical energy (watt-hour) of electric water-heater and an air-conditioner}
#' }
#'
#' @examples
#' nrow(Electric_consumption)
#'
#' @source
#' Extracted from
#'
#' Hebrail G, Berard A (2012) Individual Household Electric Power Consumption.
#' UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C58K54}
#'
"Electric_consumption"

#' Skin segmentation data
#'
#' Rajen and Abhinav (2012) addressed the challenge of detecting skin-like
#' regions in images as a component of the intricate process of facial recognition.
#' To achieve this goal, they curated the “Skin segmentation” data set, comprising
#' RGB (R-red, G-green, B-blue) values of randomly selected pixels from
#' N = 245,057 facial images, including 50,859 skin samples and 194,198 nonskin
#' samples, spanning diverse age groups, racial backgrounds, and genders.
#'
#' @format A data frame with 4 columns and 245,057 rows.
#' \describe{
#' \item{\code{Skin_presence}}{Skin presence in the randomly selected pixels}
#' \item{\code{Red}}{Red values in the randomly selected pixels}
#' \item{\code{Green}}{Green values in the randomly selected pixels}
#' \item{\code{Blue}}{Blue values in the randomly selected pixels}
#' }
#'
#' @examples
#' nrow(Skin_segmentation)
#'
#' @source
#' Extracted from
#'
#' Rajen B, Abhinav D (2012) Skin segmentation. UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C5T30C}
#'
"Skin_segmentation"

#' Bike sharing data
#'
#' Fanaee-T (2013) collected data to understand the bike sharing demands under
#' the rental and return process. The data in total contains 17,379 observations
#' where we consider the covariates temperature, humidity and windspeed
#' to model the response, the number of bikes rented hourly.
#'
#' @format A data frame with 4 columns and 17,379 rows.
#' \describe{
#' \item{\code{Rented_Bikes}}{Number of bikes rented hourly}
#' \item{\code{Temperature}}{Hourly temperature}
#' \item{\code{Humidity}}{Hourly humidity}
#' \item{\code{Windspeed}}{Hourly windspeed}
#' }
#'
#' @examples
#' nrow(Bike_sharing)
#'
#' @source
#' Extracted from
#'
#' Fanaee-T H (2013) Bike Sharing. UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C5W894}
#'
"Bike_sharing"

#' One Million Songs data
#'
#' This data set contains 1,019,318 unique users' music play counts in the Echo Nest,
#' which is available at "http://millionsongdataset.com/tasteprofile/".
#' As a basic step, it is interesting to predict the play counts using the song
#' information collected in the Million Song Dataset (Bertin-Mahieux et al. (2011)).
#' After cleaning up and feature engineering the data in total contains 205,032
#' observations where we consider the covariates duration, loudness, tempo, artist
#' hotness, song hotness, and album hotness to model the response,
#' the number of song counts.
#'
#' @format A data frame with 4 columns and 309,685 rows.
#' \describe{
#' \item{\code{Counts}}{Number of playback counts for songs}
#' \item{\code{Duration}}{Duration of the song}
#' \item{\code{Loudness}}{Loudness of the song}
#' \item{\code{Tempo}}{Tempo of the song}
#' \item{\code{Artist_Hotness}}{A value between 0 and 1}
#' \item{\code{Song_Hotness}}{A value between 0 and 1}
#' \item{\code{Album_Hotness}}{A value between 0 and 1}
#' }
#'
#' @examples
#' nrow(One_Million_Songs)
#'
#' @references
#' \insertRef{mcfee2012million}{NeEDS4BigData}
#' \insertRef{ai2021optimal}{NeEDS4BigData}
#'
"One_Million_Songs"

#' Generate data for Generalised Linear Models
#'
#' Function to simulate big data under linear, logistic and Poisson regression for sampling.
#' Covariate data X is through Normal, Multivariate Normal or Uniform distribution for linear regression.
#' Covariate data X is through Exponential, Normal, Multivariate Normal or Uniform distribution for logistic regression.
#' Covariate data X is through Normal or Uniform distribution for Poisson regression.
#'
#' @usage
#' GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,family)
#'
#' @param Dist        a character value for the distribution "Normal", "MVNormal", "Uniform or "Exponential"
#' @param Dist_Par    a list of parameters for the distribution that would generate data for covariate X
#' @param No_Of_Var   number of variables
#' @param Beta        a vector for the model parameters, including the intercept
#' @param N           the big data size
#' @param family      a character vector for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Big data for the Generalised Linear Models are generated by the "linear", "logistic" and "poisson"
#' regression types.
#'
#' We have limited the covariate data generation for
#' linear regression through normal, multivariate normal and uniform distribution,
#' logistic regression through exponential, normal, multivariate normal and uniform distribution
#' Poisson regression through normal and uniform distribution.
#'
#' @return
#' The output of \code{GenGLMdata} gives a list of
#'
#' \code{Complete_Data} a matrix for Y and X
#'
#' @references
#' \insertRef{lee1996hierarchical}{NeEDS4BigData}
#'
#' @examples
#' No_Of_Var<-2; Beta<-c(-1,2,1); N<-5000;
#'
#' # Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1,Error_Variance=0.5)
#' Dist<-"MVNormal";
#' Dist_Par<-list(Mean=rep(0,No_Of_Var),Variance=diag(rep(2,No_Of_Var)),Error_Variance=0.5)
#' # Dist<-"Uniform"; Dist_Par<-list(Min=0,Max=1)
#'
#' Family<-"linear"
#' Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' Dist<-"Normal"; Dist_Par<-list(Mean=0,Variance=1);
#' # Dist<-"MVNormal"; Dist_Par<-list(Mean=rep(0,No_Of_Var),Variance=diag(rep(2,No_Of_Var)))
#' # Dist<-"Exponential"; Dist_Par<-list(Rate=3)
#' # Dist<-"Uniform"; Dist_Par<-list(Min=0,Max=1)
#'
#' Family<-"logistic"
#' Results<-GenGLMdata(Dist,Dist_Par,No_Of_Var,Beta,N,Family)
#'
#' # Dist<-"Normal";
#' Dist<-"Uniform"; Family<-"poisson"
#' Results<-GenGLMdata(Dist,NULL,No_Of_Var,Beta,N,Family)
#'
#' @import stats
#' @importFrom mvnfast rmvn
#' @export
GenGLMdata<-function(Dist,Dist_Par,No_Of_Var,Beta,N,family){
  if(any(is.na(c(Dist,Beta,No_Of_Var,N,family))) | any(is.nan(c(Dist,No_Of_Var,Beta,N,family)))){
    stop("NA or Infinite or NAN values in the Dist,Beta,No_Of_Var,N or family")
  }

  if((length(N) + length(family)) != 2){
    stop("N or family has a value greater than length one")
  }

  if(any(is.na(unlist(Dist_Par))) | any(is.nan(unlist(Dist_Par)))){
    stop("NA or Infinite or NAN values in the Dist_Par")
  }

  if(!any(family == c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if(family == "linear"){
    if(!(Dist == "Normal" | Dist == "MVNormal" | Dist == "Uniform")){
      stop("For linear regression select the distribution 'Normal', 'MVNormal' \n or 'Uniform' to generate the covarate data")
    }
  }

  if(family == "logistic"){
    if(!(Dist == "Exponential" | Dist == "Normal" | Dist == "MVNormal" | Dist == "Uniform")){
      stop("For logistic regression select the distribution 'Exponential', \n 'Normal', 'MVNormal' or 'Uniform' to generate the covarate data")
    }
  }

  if(family == "poisson"){
    if(!(Dist == "Normal" | Dist == "Uniform")){
      stop("For poisson regression select the distribution 'Normal' \n or 'Uniform' to generate the covarate data")
    }
  }

  if(family %in% "linear"){
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }
    if(Dist %in% "MVNormal"){
      X<-mvnfast::rmvn(n = N, mu = Dist_Par$Mean, sigma = sqrt(Dist_Par$Variance))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }

    Complete_Data<-cbind(1,X);
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))
    Residual<-stats::rnorm(n=N,mean=0,sd=sqrt(Dist_Par$Error_Variance))
    Linear_Predictor_data <- Complete_Data%*%Beta
    Y <- Linear_Predictor_data + Residual

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    Outputs<-list("Complete_Data"=Complete_Data)

    class(Outputs)<-c("linear")
    return(Outputs)
  }
  if(family %in% "logistic"){
    if(Dist %in% "Exponential"){
      X<-replicate(No_Of_Var,stats::rexp(n = N, rate = Dist_Par$Rate))
    }
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = Dist_Par$Mean, sd = sqrt(Dist_Par$Variance)))
    }
    if(Dist %in% "MVNormal"){
      X<-mvnfast::rmvn(n = N, mu = Dist_Par$Mean, sigma = sqrt(Dist_Par$Variance))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = Dist_Par$Min, max = Dist_Par$Max))
    }

    Complete_Data<-cbind(1,X)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))

    Linear_Predictor_Data <- Complete_Data%*%Beta
    Pi_Data <- 1-1/(1+exp(Linear_Predictor_Data))
    Y <- stats::rbinom(n=N,size=1,prob = Pi_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    Outputs<-list("Complete_Data"=Complete_Data)

    class(Outputs)<-c("logistic")
    return(Outputs)
  }
  if(family %in% "poisson"){
    if(Dist %in% "Normal"){
      X<-replicate(No_Of_Var,stats::rnorm(n = N, mean = 0, sd = 1))
    }
    if(Dist %in% "Uniform"){
      X<-replicate(No_Of_Var,stats::runif(n = N, min = 0, max = 1))
    }

    Complete_Data<-cbind(1,X)
    colnames(Complete_Data)<-c(paste0("X",0:ncol(X)))

    Linear_Predictor_Data <- Complete_Data%*%Beta
    Lambda_Data <- exp(Linear_Predictor_Data)
    Y <- stats::rpois(n=N,lambda = Lambda_Data)

    Complete_Data<-cbind(Y,Complete_Data)
    colnames(Complete_Data)<-c("Y",paste0("X",0:ncol(X)))

    Outputs<-list("Complete_Data"=Complete_Data)

    class(Outputs)<-c("poisson")
    return(Outputs)
  }
}

#' Generate data for Generalised Linear Models under model misspecification scenario
#'
#' Function to simulate big data under Generalised Linear Models for the model misspecification scenario through
#' any misspecification type.
#'
#' @usage
#' GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family)
#'
#' @param N                       the big data size
#' @param X_Data                  a matrix for the covariate data
#' @param Misspecification        a vector of values for the misspecification
#' @param Beta                    a vector for the model parameters, including the intercept and misspecification term
#' @param Var_Epsilon             variance value for the residuals
#' @param family                  a character vector for "linear", "logistic" and "poisson" regression from Generalised Linear Models
#'
#' @details
#' Big data for the Generalised Linear Models are generated by the "linear", "logistic" and "poisson"
#' regression types under model misspecification.
#'
#' @return
#' The output of \code{GenModelMissGLMdata} gives a list of
#'
#' \code{Complete_Data} a matrix for Y,X and f(x)
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
#'
#' Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon,family)
#'
#' Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family="logistic")
#'
#' Results<-GenModelMissGLMdata(N,X_Data,Misspecification,Beta,Var_Epsilon=NULL,family="poisson")
#'
#' @import stats
#' @export
GenModelMissGLMdata<-function(N,X_Data,Misspecification,Beta,Var_Epsilon,family)
{
  if(any(is.na(c(Misspecification,Beta,Var_Epsilon,family))) |
     any(is.nan(c(Misspecification,Beta,Var_Epsilon,family)))){
    stop("NA or Infinite or NAN values in the Misspecification,Beta,Var_Epsilon or family")
  }

  if((length(N) + length(family)) != 2){
    stop("N or family has a value greater than length one")
  }

  if(!any(family %in% c("linear","logistic","poisson"))){
    stop("Only the regression types 'linear','logistic' or 'poisson' are allowed")
  }

  if((length(Misspecification) != nrow(X_Data)) | (N != length(Misspecification)) |
     (N != nrow(X_Data)) ){
    stop("The covariate data size N is not the same as of the length of Misspecification")
  }

  X_Data_Real<-cbind(X_Data,Misspecification)

  if(family == "linear"){
    if(is.null(Var_Epsilon) == TRUE){Var_Epsilon<-0.5}
    Linear_Predictor_data <- X_Data_Real%*%Beta
    Y_Data <- Linear_Predictor_data + stats::rnorm(n = N,mean = 0,sd = sqrt(Var_Epsilon))
  }

  if(family == "logistic"){
    Linear_Predictor_data <- X_Data_Real%*%Beta
    Pi_Data<-1-1/(1+exp(Linear_Predictor_data))
    Y_Data <- stats::rbinom(N,1,Pi_Data)
  }

  if(family == "poisson"){
    Linear_Predictor_data <- X_Data_Real%*%Beta
    Lambda_Data<-exp(Linear_Predictor_data)
    Y_Data <- stats::rpois(N,Lambda_Data)
  }

  Real_Model_Data<-cbind(Y=Y_Data,X_Data_Real)
  colnames(Real_Model_Data)<-c("Y","X0",paste0("X",1:ncol(X_Data[,-1])),"f(x)")

  Outputs<-list("Complete_Data"=Real_Model_Data)
  if(family == "linear"){
    class(Outputs)<-c("linear")
  }
  if(family == "logistic"){
    class(Outputs)<-c("logistic")
  }
  if(family == "poisson"){
    class(Outputs)<-c("poisson")
  }

  return(Outputs)
}
