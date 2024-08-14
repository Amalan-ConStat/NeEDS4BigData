#' Plotting model parameter outputs after subsampling
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
#' For A- and L-optimality criteria subsampling under Generalised Linear Models
#' the facets are for subsample sizes and beta values.
#'
#' For A-optimality criteria subsampling under Gaussian Linear Models
#' the facets are for subsample sizes and beta values.
#'
#' For A-optimality criteria subsampling under Generalised Linear Models
#' with response variable not inclusive the facets are for subsample sizes and beta values.
#'
#' For A- and L-optimality criteria subsampling under Generalised Linear Models
#' where multiple models can describe the data the facets are for subsample sizes and beta values.
#'
#' For A- and L-optimality criteria and LmAMSE subsampling under Generalised Linear Models
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

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))->plot_beta

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

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))->plot_beta

  return(plot_beta)
}

#' @method plot_Beta A_L_OptimalSubsampling
#' @export
plot_Beta.A_L_OptimalSubsampling<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("darkred","red")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))->plot_beta

  return(plot_beta)
}

#' @method plot_Beta AoptimalSubsampling
#' @export
plot_Beta.AoptimalSubsampling<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("darkred")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))->plot_beta

  return(plot_beta)
}

#' @method plot_Beta A_OptimalSubsamplingMC
#' @export
plot_Beta.A_OptimalSubsamplingMC<-function(object){
  Temp_Data<-data.frame(object$Beta_Estimates)
  Temp_Data |>
    tidyr::pivot_longer(cols=tidyr::starts_with("Beta"),names_to="Beta",values_to="Values")->Temp_Data
  label_values<-0:(length(unique(Temp_Data$Beta))-1)

  method_colors<-c("darkred")
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r2~Beta,scales = "free", independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))->plot_beta

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

      Mean_Data <- Temp_Data |>
        dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
        dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

      method_colors<-c("darkred","red","darkgreen","green")

      ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
        ggplot2::geom_density(alpha=0.4)+
        ggh4x::facet_grid2(r2~Beta,scales = "free",labeller = ggplot2::label_parsed)+
        ggplot2::xlab(expression(paste(beta," values")))+
        ggplot2::ylab("Density")+
        ggplot2::scale_fill_manual(values = method_colors)+
        ggplot2::ggtitle(names(object$Beta_Estimates)[i])+
        ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
        ggplot2::scale_color_manual(values = method_colors)+
        ggplot2::theme_bw()+
        ggplot2::theme(legend.position = "bottom",
                       axis.text.x = ggplot2::element_text(angle = 30))->plot_beta[[i]]
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
  Log_Odds_Labels<-method_labels[startsWith(method_labels,"RLmAMSE Log Odds")]
  Power_Labels<-method_labels[startsWith(method_labels,"RLmAMSE Power")]

  method_colors<-c("darkred","red","lightgreen",rep("green",length(Log_Odds_Labels)),
                   rep("darkgreen",length(Power_Labels)))
  Temp_Data$Beta<-factor(Temp_Data$Beta,
                         levels = paste0("Beta",label_values),
                         labels = paste0("beta[",label_values,"]"))

  Mean_Data <- Temp_Data |>
    dplyr::group_by(.data$Method,.data$r2,.data$Beta) |>
    dplyr::summarise(Mean = mean(.data$Values),.groups = "drop")

  ggplot2::ggplot(data=Temp_Data,ggplot2::aes(x=.data$Values,fill=.data$Method,color=.data$Method))+
    ggplot2::geom_density(alpha=0.4)+
    ggh4x::facet_grid2(r2~Beta,scales = "free",independent = "x",labeller = ggplot2::label_parsed)+
    ggplot2::xlab(expression(paste(beta," values")))+
    ggplot2::ylab("Density")+
    ggplot2::scale_fill_manual(values = method_colors)+
    ggplot2::geom_vline(data=Mean_Data,ggplot2::aes(xintercept = .data$Mean, color = .data$Method))+
    ggplot2::scale_color_manual(values = method_colors)+
    ggplot2::theme_bw()+
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))+
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))->plot_beta

  return(plot_beta)
}

#' Plotting AMSE outputs for the subsamples under model misspecification
#'
#' After using the subsampling methods under potential model misspecification we obtain
#' their respective AMSE values for the predictions. They are summarised as plots here.
#'
#' @usage
#' plot_AMSE(object)
#'
#' @param object Any object after subsampling from our subsampling function under potential model misspecification
#'
#' @details
#' For A- and L-optimality criteria and RLmAMSE subsampling under Generalised Linear Models
#' with potential model misspecification the facets are for variance and bias^2 of AMSE values.
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
plot_AMSE<-function(object){
  UseMethod("plot_AMSE",object)
}

#' @method plot_AMSE ModelMisspecified
#' @export
plot_AMSE.ModelMisspecified<-function(object){
  Temp_Data<-data.frame(object$AMSE_Estimates) |>
    tidyr::pivot_longer(cols = c("Variance","Bias.2","AMSE"),names_to = "Metric",values_to = "Values") |>
    dplyr::group_by(.data$Method,.data$r2,.data$Metric) |>
    dplyr::summarise(Mean=mean(.data$Values),.groups = "drop")

  Temp_Data$Metric<-factor(Temp_Data$Metric,levels = c("Bias.2","Variance","AMSE"),
                           labels = c("Bias^2","Variance","AMSE"))

  method_labels<-levels(Temp_Data$Method)

  Log_Odds_Labels<-method_labels[startsWith(method_labels,"RLmAMSE Log Odds")]
  Power_Labels<-method_labels[startsWith(method_labels,"RLmAMSE Power")]
  method_colors<-c("darkred","red","lightgreen",rep("green",length(Log_Odds_Labels)),
                   rep("darkgreen",length(Power_Labels)))
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
    ggplot2::theme(legend.position = "bottom",
                   axis.text.x = ggplot2::element_text(angle = 30))+
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))->plot_AMSE

  return(plot_AMSE)
}
