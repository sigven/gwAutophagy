library(dplyr)
library(ggplot2)
library(ggsci)
library(matrixStats)
library(reshape2)
library(ggrepel)

load_all_bf_data <- function(){

  df_BF_overall <- fst::read_fst(
    file.path(here::here(), "data","processed","df_BF_overall.fst"))
  df_BF_starv <- fst::read_fst(
    file.path(here::here(), "data","processed","df_BF_starv.fst"))
  df_BF_repl <- fst::read_fst(
    file.path(here::here(), "data","processed","df_BF_repl.fst"))
  df_BF_temporal <- fst::read_fst(
    file.path(here::here(), "data","processed","df_BF_temporal.fst"))
  gene_ids_bf <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_ids_bf.fst"))

  return(list(
    df_BF_overall = df_BF_overall,
    df_BF_starv = df_BF_starv,
    df_BF_repl = df_BF_repl,
    df_BF_temporal = df_BF_temporal,
    gene_ids_bf = gene_ids_bf
  ))

}
load_all_response_data <- function(){

  df_DNN_preds <- fst::read_fst(
    file.path(here::here(), "data","processed","df_DNN_preds.fst"))
  df_DS_curvefits <- fst::read_fst(
    file.path(here::here(), "data","processed","df_DS_curvefits.fst"))
  df_DS_parms <- fst::read_fst(
    file.path(here::here(), "data","processed","df_DS_parms.fst"))
  df_DS_parms_ctrs <- fst::read_fst(
    file.path(here::here(), "data","processed","df_DS_parms_ctrs.fst"))
  regr_df_all <- fst::read_fst(
    file.path(here::here(), "data","processed","regr_df_all.fst"))
  gene_ids_kinetic <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_ids_kinetic.fst"))

  df_DS_parms_comb <- dplyr::bind_rows(
    df_DS_parms,df_DS_parms_ctrs)

  return(list(
    df_DS_parms_comb = df_DS_parms_comb,
    df_DNN_preds = df_DNN_preds,
    df_DS_curvefits = df_DS_curvefits,
    df_DS_parms = df_DS_parms,
    df_DS_parms_ctrs = df_DS_parms_ctrs,
    regr_df_all = regr_df_all,
    gene_ids_kinetic = gene_ids_kinetic
  ))
}


load_kinetic_response_data <- function(){

  kinetic_response_plot_input <- readRDS(
    file.path(here::here(), "data","processed","gw_autoph_kinetic_plot_input.rds"))
  gene_ids_kinetic <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_ids_kinetic.fst"))

  return(list(
    per_ko = kinetic_response_plot_input,
    gene_ids = gene_ids_kinetic
  ))

}

load_autophagy_competence_data <- function(){

  autophagy_competence_plot_input <- readRDS(
    file.path(here::here(), "data","processed","gw_autoph_competence_plot_input.rds"))
  gene_ids_bf <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_ids_bf.fst"))

  return(list(
    per_ko = autophagy_competence_plot_input,
    gene_ids = gene_ids_bf
  ))

}

# load_gene_response_data <- function(id = "RTG3 YBL103C", gw_autoph_data = NULL){
#
#   y_pred <- gw_autoph_data$df_DS_curvefits |>
#     dplyr::mutate(
#       id_gene_orf = stringr::str_replace(
#         paste(.data$Gene, .data$ORF),"NA ","")) |>
#     dplyr::filter(.data$id_gene_orf == id)
#
#   ## if NROW(y_pred == 0)?
#
#   #If mutant in rec plate, only evaluate rec mutants
#   if(any(grepl("Rec",y_pred$Plate))){
#     y_pred <- y_pred |>
#       dplyr::filter(stringr::str_detect(.data$Plate,"Rec"))
#   }
#   #else what?
#
#   #Find representative response
#   y_pred <- y_pred |>
#     dplyr::group_by(Time) |>
#     dplyr::mutate(d=abs(P1_30_fit-median(P1_30_fit))) |>
#     dplyr::group_by(Plate, Position) |>
#     dplyr::mutate(d=mean(d))
#   y_pred <- y_pred[which(y_pred$d==min(y_pred$d)),]
#   y_raw <- gw_autoph_data$df_DNN_preds |>
#     dplyr::filter(
#       paste(Plate, Position) == paste(y_pred$Plate, y_pred$Position)[1])
#
#   slibrary <- ifelse(y_raw$Type== "KO" | !(is.na(y_raw$Plate_controls)),"KO","DAmP")[1]
#   plate_id <- unique(y_raw$Plate)
#   y_pred.ctr <- gw_autoph_data$df_DS_curvefits[
#     which(is.na(gw_autoph_data$df_DS_curvefits$Plate_controls)),] |>
#     #subset(Type == Library) |>
#     dplyr::filter(
#       .data$Plate == plate_id) |>
#     dplyr::group_by(.data$Time) |>
#     dplyr::summarise(
#       P1_30_fit=median(.data$P1_30_fit), .groups = "drop")
#
#   plate_position <- toupper(stringr::str_replace_all(
#     paste(y_pred$Plate, y_pred$Position)[1]," ",""))
#   double_sigmoidal_model <- gw_autoph_data$regr_df_all |>
#     dplyr::filter(startsWith(toupper(.data$id2), plate_position))
#
#   return(list(
#     y_pred = y_pred,
#     y_raw = y_raw,
#     y_pred.ctr = y_pred.ctr,
#     double_sigmoidal_model = double_sigmoidal_model,
#     slibrary = slibrary,
#     id = id
#   ))
#
# }
#

plot_autophagy_competence <- function(competence_data = NULL, Model = "30"){

  p <- NULL
  if(Model == "30"){
    p <- ggplot2::ggplot(competence_data$BF_response,
                         ggplot2::aes(log_BFt_VAM6.ATG1_30, log_BFt_WT.VAM6_30, col=TimeR))+
      ggplot2::geom_vline(xintercept = 0, lty=1, col="black", size=0.1)+
      ggplot2::geom_hline(yintercept = 0, lty=1, col="black", size=0.1)+
      scico::scale_color_scico(palette = 'lisbon')+
      ggplot2::scale_size(range = c(0, 1.5))+
      #geom_path(data=BF_response_ctr, aes(),col="gray",lty=2, alpha = 1) +
      ggplot2::geom_point(
        data=competence_data$BF_response_ctr,aes(size=log_BFt_WT.ATG1_30, pch="Control"),size=1.8, alpha = 1) +
      ggplot2::geom_path(ggplot2::aes(),col="black",lty=2, alpha = 1) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = log_BFt_VAM6.ATG1_30.shift, yend = log_BFt_WT.VAM6_30.shift, size=log_BFt_WT.ATG1_30),
        arrow = arrow(angle = 35, length = unit(0.11, "inches"), type = "closed"), alpha = 1, size=1.3) +
      ggplot2::labs(x="log BF (VAM6:ATG1) 30", y="log BF (WT:VAM6) 30", col="Time", size="BF (WT:ATG1)",
                    title=paste(competence_data$id, competence_data$Library), shape="")+
      ggplot2::theme_bw(base_size = 20, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=20, face="bold"),
        axis.title = ggplot2::element_text(size=18),
        axis.text = ggplot2::element_text(size=18),
        legend.text = ggplot2::element_text(size=18),
        strip.background = ggplot2::element_blank(),
        legend.position = "top")
  }else{
    p <- ggplot2::ggplot(competence_data$BF_response,
                         ggplot2::aes(log_BFt_VAM6.ATG1_22, log_BFt_WT.VAM6_22, col=TimeR))+
      ggplot2::geom_vline(xintercept = 0, lty=1, col="black", size=0.1)+
      ggplot2::geom_hline(yintercept = 0, lty=1, col="black", size=0.1)+
      scale_color_scico(palette = 'lisbon')+
      ggplot2::scale_size(range = c(0, 1.5))+
      #geom_path(data=BF_response_ctr,aes(),col="gray",lty=2, alpha = 1) +
      ggplot2::geom_point(
        data=competence_data$BF_response_ctr,aes(size=log_BFt_WT.ATG1_22, pch="Control"),size=1.8, alpha = 1) +
      ggplot2::geom_path(ggplot2::aes(),col="black",lty=2, alpha = 1) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = log_BFt_VAM6.ATG1_22.shift, yend = log_BFt_WT.VAM6_22.shift, size=log_BFt_WT.ATG1_22),
        arrow = arrow(angle = 35, length = unit(0.11, "inches"), type = "closed"), alpha = 1, size=3) +
      ggplot2::labs(x="log BF (VAM6:ATG1) 22", y="log BF (WT:VAM6) 22", col="Time", size="BF (WT:ATG1)",
                    title=paste(competence_data$id, competence_data$Library), shape="")+
      ggplot2::theme_bw(base_size = 20, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=20, face="bold"),
        axis.title = ggplot2::element_text(size=18),
        axis.text = ggplot2::element_text(size=18),
        legend.text = ggplot2::element_text(size=18),
        legend.position = "top")
  }
  return(p)


}

plot_response_kinetics <- function(response_data = NULL){

  y_pred <- NULL
  y_raw <- NULL
  y_pred.ctr <- NULL
  double_sigmoidal_model <- NULL
  slibrary <- NULL
  id = NULL
  if(!is.null(response_data)){
    if("y_pred" %in% names(response_data)){
      y_pred <- response_data$y_pred
    }
    if("y_raw" %in% names(response_data)){
      y_raw <- response_data$y_raw
    }
    if("y_pred.ctr" %in% names(response_data)){
      y_pred.ctr <- response_data$y_pred.ctr
    }
    if("double_sigmoidal_model" %in% names(response_data)){
      double_sigmoidal_model <- response_data$double_sigmoidal_model
    }
    if("slibrary" %in% names(response_data)){
      slibrary <- response_data$slibrary
    }
    if("id" %in% names(response_data)){
      id <- response_data$id
    }
  }

  alpha_ribbon <- 0.1
  alpha_segment <- 1
  size_segment <- 0.4
  size_line <- 0.4
  arrow_length <- 0.15

  p <- ggplot2::ggplot()+
    ggplot2::geom_vline(xintercept=0, lty=2)+
    ggplot2::geom_vline(xintercept=12, lty=2)+
    ggplot2::geom_hline(
      yintercept=0, lty=1, linewidth=size_line)+
    ggplot2::geom_hline(
      yintercept=double_sigmoidal_model[,"maximum_y"],
      lty=1, linewidth=size_line)+
    ggplot2::geom_hline(
      yintercept=double_sigmoidal_model[,"finalAsymptoteIntensity"],
      lty=1, linewidth=size_line)+
    ggplot2::geom_segment(
      ggplot2::aes(x=double_sigmoidal_model[,"startPoint_x"], y=0,
                   xend=double_sigmoidal_model[,"reachMaximum_x"],
                   yend=double_sigmoidal_model[,"reachMaximum_y"] ),
      lty=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=double_sigmoidal_model[,"startDeclinePoint_x"],
                   y=double_sigmoidal_model[,"startDeclinePoint_y"],
                   xend=double_sigmoidal_model[,"endDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"endDeclinePoint_y"] ),
      lty=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=double_sigmoidal_model[,"midPoint1_y"],
                   xend=double_sigmoidal_model[,"midPoint1_x"],
                   yend=double_sigmoidal_model[,"midPoint1_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"midPoint2_y"],
                   xend=double_sigmoidal_model[,"midPoint2_x"],
                   yend=double_sigmoidal_model[,"midPoint2_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=0,
                   xend=double_sigmoidal_model[,"startPoint_x"],
                   yend=0),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"maximum_y"],
                   xend=double_sigmoidal_model[,"startDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"maximum_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=double_sigmoidal_model[,"reachMaximum_y"],
                   xend=double_sigmoidal_model[,"reachMaximum_x"],
                   yend=double_sigmoidal_model[,"reachMaximum_y"]),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"endDeclinePoint_y"],
                   xend=double_sigmoidal_model[,"endDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"endDeclinePoint_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment)+
    geom_line(
      ggplot2::aes(x= y_pred.ctr$Time, y=y_pred.ctr$P1_30_fit, col="Control"))+
    geom_ribbon(
      ggplot2::aes(x=y_pred.ctr$Time,
                   y=colMeans(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                   ymin=colMins(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                   ymax=colMaxs(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit))),
      fill="yellow", alpha=0.1)+
    geom_ribbon(
      ggplot2::aes(x=c(double_sigmoidal_model[,"startPoint_x"],
                       double_sigmoidal_model[,"reachMaximum_x"]),
                   y=c(0, double_sigmoidal_model[,"reachMaximum_y"]),
                   ymin=c(0,0),
                   ymax=c(double_sigmoidal_model[,"reachMaximum_y"],
                          double_sigmoidal_model[,"reachMaximum_y"])),fill="black",
      alpha=0.1)+
    ggplot2::geom_ribbon(
      ggplot2::aes(x=c(double_sigmoidal_model[,"startDeclinePoint_x"],
                       double_sigmoidal_model[,"endDeclinePoint_x"]),
                   y=c(double_sigmoidal_model[,"startDeclinePoint_y"],
                       double_sigmoidal_model[,"endDeclinePoint_y"]),
                   ymin=c(double_sigmoidal_model[,"endDeclinePoint_y"],double_sigmoidal_model[,"endDeclinePoint_y"]),
                   ymax=c(double_sigmoidal_model[,"startDeclinePoint_y"],
                          double_sigmoidal_model[,"startDeclinePoint_y"])),fill="black",
      alpha=0.1)+
    ggplot2::geom_line(
      ggplot2::aes(x= y_pred$Time, y=y_pred$P1_30_fit, col="Model"))+
    ggplot2::geom_point(
      ggplot2::aes(x=ifelse(
        grepl("Starv",y_raw$Phase),
        y_raw$TimeR,y_raw$TimeR+12),y_raw$P1_30, col="Data"))+
    ggplot2::ylim(0,100)+
    ggplot2::xlim(-1,20)+
    ggplot2::labs(x="Time (hours)", y="Autophagy (%)", title=paste(id, slibrary))+
    ggplot2::scale_color_manual(
      values=ggsci::pal_jco("default", alpha = 1)(7)[c(7,1,4)])+
    ggplot2::theme_bw(base_size = 18, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size=20, face="bold"),
      axis.title = ggplot2::element_text(size=18),
      axis.text = ggplot2::element_text(size=18),
      legend.text = ggplot2::element_text(size=18),
      legend.title = ggplot2::element_blank()
    )

  return(p)
}

