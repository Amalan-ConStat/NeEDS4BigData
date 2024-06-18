#' Plotting model parameter (Beta) outputs after subsampling
#'
#' After using the subsampling methods we mostly obtain the estimated model parameter
#' estimates. Here, they are summarised as histogram plots.
#'
#' @usage
#' plot_Beta(object)
#'
#' @param object Any object after subsampling from our subsampling functions
#'
#' @details
#' For local case control sampling the facets are for subsample sizes and beta values.
#'
#' For leverage sampling the facets are for subsample sizes and beta values.
#'
#' For A- and L-optimality criterion subsampling under Generalised Linear Models
#' the facets are for subsample sizes and beta values.
#'
#' For A-optimality criterion subsampling under Gaussian Linear Models
#' the facets are for subsample sizes and beta values.
#'
#' For A-optimality criterion subsampling under Generalised Linear Models
#' with response variable not inclusive the facets are for subsample sizes and beta values.
#'
#' For A- and L-optimality criterion subsampling under Generalised Linear Models
#' where multiple models can describe the data the facets are for subsample sizes and beta values.
#'
#' For A- and L-optimality criterion and LmAMSE subsampling under Generalised Linear Models
#' with potential model misspecification the facets are for subsample sizes and beta values.
#'
#' @return
#' The output is a faceted ggplot result
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom ggh4x facet_grid2
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom dplyr group_by summarise
#' @export
plot_Beta<-function(object){
  UseMethod("plot_Beta",object)
}

#' @method plot_Beta LocalCaseControl
#' @export
plot_Beta.LocalCaseControl<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("darkred")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' @method plot_Beta Leverage
#' @export
plot_Beta.Leverage<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("maroon","red","darkred")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' @method plot_Beta A_L_OptimalSubsampling
#' @export
plot_Beta.A_L_OptimalSubsampling<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("red","darkred")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' @method plot_Beta AoptimalSubsampling
#' @export
plot_Beta.AoptimalSubsampling<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("red")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' @method plot_Beta A_OptimalSubsamplingMC
#' @export
plot_Beta.A_OptimalSubsamplingMC<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("red")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' @method plot_Beta ModelRobust
#' @export
plot_Beta.ModelRobust<-function(object){
  plot_beta<-list()
    for (i in 1:length(object$Beta_Estimates)) {
      Temp_Data<-data.frame(object$Beta_Estimates[[i]])
      Temp_Data |>
        tidyr::pivot_longer(cols=tidyr::starts_with("beta"),names_to="Beta",values_to="Values")->Temp_Data
      label_values<-0:(length(unique(Temp_Data$Beta))-1)
      Temp_Data$Beta<-factor(Temp_Data$Beta,
                             levels = paste0("beta_",label_values),
                             labels = paste0("beta[",label_values,"]"))

      ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
        ggplot2::geom_histogram()+
        ggh4x::facet_grid2(r2~Beta,scales = "free",labeller = ggplot2::label_parsed)+
        ggplot2::xlab(expression(paste(beta," values")))+
        ggplot2::ylab("Frequency")+
        ggplot2::scale_fill_manual(values = c("red","pink","darkgreen","green"))+
        ggplot2::ggtitle(names(object$Beta_Estimates)[i])+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "bottom")->plot_beta[[i]]
    }
    names(plot_beta)<-names(object$Beta_Estimates)
    return(plot_beta)
}

#' @method plot_Beta ModelMisspecified
#' @export
plot_Beta.ModelMisspecified<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)
  method_labels<-levels(Temp_Data$Method)
  Log_Odds_Labels<-method_labels[startsWith(method_labels,"LmAMSE Log Odds")]
  Power_Labels<-method_labels[startsWith(method_labels,"LmAMSE Power")]

  method_colors<-c("red","pink","lightgreen",paste0("springgreen",1:length(Log_Odds_Labels)),
                   paste0("seagreen",1:length(Power_Labels)))
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method))+
    ggplot2::geom_histogram()+
    ggh4x::facet_grid2(r2~Beta,scales = "free",independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Frequency")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_beta

  return(plot_beta)
}

#' Plotting A- and D-optimality for model parameter (Beta) outputs of the subsamples
#'
#' After using the subsampling methods we mostly obtain the estimated model parameter
#' estimates. Of these estimates we can obtain their respective A- and D-optimality
#' values and they are summarised as plots here. It should be noted that
#' A- and D-optimality represents the model parameter's tr(Variance) and log(det(Information)),
#' respectively.
#'
#' @usage
#' plot_Utility(object)
#'
#' @param object Any object after subsampling from our subsampling functions
#'
#' @details
#' For A- and L-optimality criterion subsampling under Generalised Linear Models
#' the facets are for model parameter variance and information.
#'
#' For A- and L-optimality criterion subsampling under Generalised Linear Models
#' where multiple models can describe the data the facets are for model parameter variance and information.
#'
#' @return
#' The output is a faceted ggplot result
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom ggh4x facet_grid2
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom dplyr group_by summarise
#' @export
plot_Utility<-function(object){
  UseMethod("plot_Utility",object)
}

#' @method plot_Utility A_L_OptimalSubsampling
#' @export
plot_Utility.A_L_OptimalSubsampling<-function(object){
  Temp_Data<-as.data.frame(object$Utility_Estimates)
  Temp_Data$Information<-log(Temp_Data$Information)
  colnames(Temp_Data)[4]<-c("log(Information)")
  Temp_Data |>
    tidyr::pivot_longer(cols=c("Variance","log(Information)"),names_to="Metric",values_to="Values")|>
    dplyr::group_by(.data$Method,.data$r2,.data$Metric) |>
    dplyr::summarise(Mean=mean(.data$Values),.groups = "drop")->Temp_Data

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$r2,y=.data$Mean,color=.data$Method,
                                              group=.data$Method,linetype=.data$Method))+
    ggplot2::geom_line()+
    ggplot2::geom_point(size=2)+
    ggplot2::scale_color_manual(values = c("red","darkred"))+
    ggplot2::scale_linetype_manual(values = c("dashed","dotted"))+
    ggplot2::scale_x_continuous(labels=unique(Temp_Data$r2),breaks=unique(Temp_Data$r2))+
    ggh4x::facet_grid2(~Metric,scales = "free",independent = "y")+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_utility

  return(plot_utility)
}

#' @method plot_Utility ModelRobust
#' @export
plot_Utility.ModelRobust<-function(object){
  Temp_Data<-list()
  for (i in 1:length(object$Utility_Estimates)) {
    Temp_Data[[i]]<-as.data.frame(object$Utility_Estimates[[i]])
    Temp_Data[[i]]$Information<-log(Temp_Data[[i]]$Information)
    colnames(Temp_Data[[i]])[4]<-c("log(Information)")
    Temp_Data[[i]] |>
      tidyr::pivot_longer(cols=c("Variance","log(Information)"),names_to="Metric",values_to="Values")|>
      dplyr::group_by(.data$Method,.data$r2,.data$Metric) |>
      dplyr::summarise(Mean=mean(.data$Values),.groups = "drop")->Temp_Data[[i]]
    Temp_Data[[i]]$Model<-c(names(object$Utility_Estimates)[i])
  }
  Temp_Data<-do.call(rbind,Temp_Data)
  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$r2,y=.data$Mean,color=.data$Method,
                                              group=.data$Method,linetype=.data$Method))+
    ggplot2::geom_line()+
    ggplot2::geom_point(size=2)+
    ggplot2::scale_color_manual(values = c("red","pink","darkgreen","green"))+
    ggplot2::scale_linetype_manual(values = c("dashed","dashed","dotted","dotted"))+
    ggplot2::scale_x_continuous(labels=unique(Temp_Data$r2),breaks=unique(Temp_Data$r2))+
    ggh4x::facet_grid2(Metric~Model,scales = "free",independent = "y")+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_utility

  return(plot_utility)
}

#' Plotting loss function or mAMSE outputs for the subsamples under model misspecification
#'
#' After using the subsampling methods under potential model misspecification we obtain
#' their respective loss or mAMSE values. They are summarised as plots here.
#'
#' @usage
#' plot_LmAMSE(object)
#'
#' @param object Any object after subsampling from our subsampling function under potential model misspecification
#'
#' @details
#' For A- and L-optimality criterion and LmAMSE subsampling under Generalised Linear Models
#' with potential model misspecification the facets are for variance and bias^2 of mAMSE values.
#'
#' @return
#' The output is a faceted ggplot result
#'
#' @import ggplot2
#' @importFrom rlang .data
#' @importFrom ggh4x facet_grid2
#' @importFrom tidyr pivot_longer starts_with
#' @importFrom dplyr group_by summarise
#' @export
plot_LmAMSE<-function(object){
  UseMethod("plot_LmAMSE",object)
}

#' @method plot_LmAMSE ModelMisspecified
#' @export
plot_LmAMSE.ModelMisspecified<-function(object){
  Temp_Data<-data.frame(object$Loss_Estimates) |>
    tidyr::pivot_longer(cols = c("Variance","Bias.2","Loss"),names_to = "Metric",values_to = "Values") |>
    dplyr::group_by(.data$Method,.data$r2,.data$Metric) |>
    dplyr::summarise(Mean=mean(.data$Values),.groups = "drop")

  Temp_Data$Metric<-factor(Temp_Data$Metric,levels = c("Bias.2","Variance","Loss"),
                           labels = c("Bias^2","Variance","Loss"))

  method_labels<-levels(Temp_Data$Method)

  Log_Odds_Labels<-method_labels[startsWith(method_labels,"LmAMSE Log Odds")]
  Power_Labels<-method_labels[startsWith(method_labels,"LmAMSE Power")]
  method_colors<-c("red","pink","lightgreen",paste0("springgreen",1:length(Log_Odds_Labels)),
                   paste0("seagreen",1:length(Log_Odds_Labels)))
  method_linetypes<-c(rep("dashed",2),"solid",rep("dotted",length(Log_Odds_Labels)),
                      rep("dotdash",length(Power_Labels)))

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$r2,y=.data$Mean,color=.data$Method,
                                              linetype=.data$Method,group=.data$Method))+
    ggplot2::geom_point(size=2)+
    ggplot2::geom_line()+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::scale_linetype_manual(values = method_linetypes)+
    ggplot2::scale_x_continuous(labels=unique(Temp_Data$r2),breaks=unique(Temp_Data$r2))+
    ggplot2::facet_wrap(~Metric,scales = "free")+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom")->plot_lmAMSE

  return(plot_lmAMSE)
}
