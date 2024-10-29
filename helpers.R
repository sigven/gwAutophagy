library(dplyr)
library(ggplot2)
library(ggsci)
library(matrixStats)
library(reshape2)
library(ggrepel)


show_gene_info <- function(primary_id = NULL,
                           gene_info = NULL){

  res <- list()
  res[['name']] <- '<i>No information available</i>'
  res[['human_orthologs']] <- '<i>No information available</i>'
  res[['description']] <- ''
  res[['sgd_link']] <- '<i>No information available</i>'
  if(!is.null(primary_id) & !is.null(gene_info)){
    gene_info_id <- gene_info |>
      dplyr::filter(primary_identifier == primary_id)
    if(NROW(gene_info_id) > 0){
      res[['name']] <- gene_info_id$genename
      res[['description']] <- gene_info_id$sgd_description
      res[['human_orthologs']] <-
        gene_info_id$human_ortholog_links
      if(!is.na(gene_info_id$sgd_id)){
        res[['sgd_link']] <- paste0(
          "<a href='https://www.yeastgenome.org/locus/",
          stringr::str_replace(
            gene_info_id$sgd_id,"SGD:",""),
            "' target='_blank'>",res[['name']],
          "</a>")
      }else{
        res[['sgd_link']] <- res[['name']]
      }
    }
  }
  return(res)

}


load_kinetic_response_data <- function(){

  df_DS_curvefits <- fst::read_fst(
    file.path(here::here(), "data","processed","ds_curvefits.fst"))
  df_DS_parms <- fst::read_fst(
    file.path(here::here(), "data","processed","ds_parms.fst"))
  df_DS_parms_ctrs <- fst::read_fst(
    file.path(here::here(), "data","processed","ds_parms_ctrs.fst"))
  kinetic_response_plot_input <- readRDS(
    file.path(here::here(), "data","processed","gw_autoph_kinetic_plot_input.rds"))
  gene_info_kinetic <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_info_kinetic.fst"))
  gene_info_kinetic_multi <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_info_kinetic_multi.fst"))
  df_DNN_preds <- fst::read_fst(
    file.path(here::here(), "data","processed","dnn_preds.fst"))
  df_DS_parms_comb <- dplyr::bind_rows(
    df_DS_parms,df_DS_parms_ctrs)


  return(list(
    ds_curvefits = df_DS_curvefits,
    ds_parms = df_DS_parms,
    ds_parms_ctrs = df_DS_parms_ctrs,
    ds_parms_comb = df_DS_parms_comb,
    dnn_preds = df_DNN_preds,
    per_ko = kinetic_response_plot_input,
    gene_info_kinetic = gene_info_kinetic,
    gene_info_kinetic_multi = gene_info_kinetic_multi
  ))

}

load_autophagy_competence_data <- function(){

  autophagy_competence_plot_input <- readRDS(
    file.path(here::here(), "data","processed","gw_autoph_competence_plot_input.rds"))
  gene_info_bf <- fst::read_fst(
    file.path(here::here(), "data","processed","gene_info_bf.fst"))
  df_DNN_preds <- fst::read_fst(
    file.path(here::here(), "data","processed","dnn_preds.fst"))
  df_BF_overall <- fst::read_fst(
    file.path(here::here(), "data","processed","BF_overall.fst"))
  df_BF_temporal <- fst::read_fst(
    file.path(here::here(), "data","processed","BF_temporal.fst"))

  return(list(
    bf_overall = df_BF_overall,
    dnn_preds = df_DNN_preds,
    bf_temporal = df_BF_temporal,
    per_ko = autophagy_competence_plot_input,
    gene_info_bf = gene_info_bf

  ))

}

plot_autophagy_competence_multi <- function(
    competence_data = NULL,
    primary_identifiers = NULL,
    dnn_model = "30",
    user_x = "Autophagosome formation",
    user_y = "Autophagosome clearance"){

  #Select x and y
  variables <- colnames(
    competence_data[['bf_overall']])[5:16]

  X <- variables[1]
  Y <- variables[1]
  if(dnn_model == "30"){
    if(user_x == "Autophagosome formation"){
      X <- variables[3]
    }
    if(user_x == "Autophagosome clearance"){
      X <- variables[2]
    }
    if(user_y == "Autophagosome formation"){
      Y <- variables[3]
    }
    if(user_y == "Autophagosome clearance"){
      Y <- variables[2]
    }
  }

  if(dnn_model == "22"){
    if(user_x == "Autophagosome clearance"){
      X <- variables[5]
    }
    if(user_x == "Autophagosome formation"){
      X <- variables[6]
    }
    if(user_x == "Overall autophagy"){
      X <- variables[4]
    }
    if(user_y == "Autophagosome clearance"){
      Y <- variables[5]
    }
    if(user_y == "Autophagosome formation"){
      Y <- variables[6]
    }
    if(user_y == "Overall autophagy"){
      Y <- variables[4]
    }
  }

  lab_x <- ""
  lab_y <- ""
  if(grepl("VAM6.ATG1", X, fixed = T)){
    lab_x <- paste0("<br><b>Autophagosome formation</b><br><br>",
                    "<i>log BF (VAM6:ATG1) - DNN model '",
                    dnn_model,"'</i>")
  }else if(grepl("WT.ATG1", X, fixed = T)){
    lab_x <- paste0("<br><b>Overall autophagy</b><br><br>",
                    "<i>log BF (WT:ATG1) - DNN model '",
                    dnn_model,"'</i>")
  }else{
    lab_x <- paste0("<br><b>Autophagosome clearance</b><br><br>",
                    "<i>log BF (WT:VAM6) - DNN model '",
                    dnn_model,"'</i>")
  }
  if(grepl("VAM6.ATG1", Y, fixed = T)){
    lab_y <- paste0("<b>Autophagosome formation</b><br><br>",
                    "<i>log BF (VAM6:ATG1) - DNN model '",
                    dnn_model,"'</i><br>")

  }else if(grepl("WT.ATG1", Y, fixed = T)){
    lab_y <- paste0("<b>Overall autophagy</b><br><br>",
                    "<i>log BF (WT:ATG1) - DNN model '",
                    dnn_model,"'</i><br>")
  }else{
    lab_y <- paste0("<b>Autophagosome clearance</b><br><br>",
                    "<i>log BF (WT:VAM6) - DNN model '",
                    dnn_model,"'</i><br>")
  }

  Positions <- c()
  for(i in primary_identifiers){
    BF_response <- competence_data[['bf_temporal']] |>
      dplyr::filter(primary_identifier == i)
      #subset(gsub("NA ", "", paste(Gene, ORF)) == i)
    #If mutant in rec plate, only evaluate rec mutants
    if(any(grepl("Rec",BF_response$Plate))){
      BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
    }
    #Find representative response
    BF_response <- BF_response |>
      dplyr::group_by(TimeR) |>
      dplyr::mutate(d=abs(log_BFt_WT.ATG1_30-median(log_BFt_WT.ATG1_30)) +
               abs(log_BFt_WT.VAM6_30-median(log_BFt_WT.VAM6_30)) +
               abs(log_BFt_VAM6.ATG1_30-median(log_BFt_VAM6.ATG1_30)) +
               abs(log_BFt_WT.ATG1_22-median(log_BFt_WT.ATG1_22)) +
               abs(log_BFt_WT.VAM6_22-median(log_BFt_WT.VAM6_22)) +
               abs(log_BFt_VAM6.ATG1_22-median(log_BFt_VAM6.ATG1_22))) |>
      dplyr::group_by(Plate, Position) |>
      dplyr::mutate(d=mean(d))
    BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
    Positions <- c(Positions,paste(BF_response$Plate, BF_response$Position)[1])
  }

  mat <- competence_data[['bf_overall']] |>
    as.data.frame()
  mat_select <- competence_data[['bf_overall']] |>
    subset(paste(Plate, Position) %in% Positions) |>
    as.data.frame()

  mat$X <- mat[,X]
  mat$Y <- mat[,Y]
  mat_select$X <- mat_select[,X]
  mat_select$Y <- mat_select[,Y]

  mat$Type <- competence_data[['dnn_preds']]$Type[match(
    paste(mat$Gene, mat$ORF),paste(
      competence_data[['dnn_preds']]$Gene,
      competence_data[['dnn_preds']]$ORF))]


  p <- ggplot2::ggplot(mat[is.na(mat$Plate_controls),], ggplot2::aes(X, Y)) +
    ggplot2::geom_point(col="lightgray", size=0.9, pch=16, alpha=0.8) +
    # ggplot2::stat_density_2d(data=mat[is.na(mat$Plate_controls),],
    #                 aes(lty=Type), col="grey30", alpha=1) +
    ggplot2::stat_density_2d(data=mat[which(mat$Plate_controls=="+"),] |>
                               dplyr::mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
                    ggplot2::aes(fill=Gene, group=Gene, alpha = ..level..),
                    geom = "polygon", col=NA) +
    ggplot2::geom_point(
      data=mat[which(!is.na(mat$Reference_sets) & !grepl("ORF",mat$Reference_sets)),],
               ggplot2::aes(col=Reference_sets), size=1.3) +
    ggplot2::geom_point(data=mat_select, ggplot2::aes(), pch=21) +
    ggrepel::geom_text_repel(data=mat_select, ggplot2::aes(label=Gene),
                    force=2, size=5,
                    max.overlaps = Inf,
                    min.segment.length = 0,
                    segment.size = 0.3) +
    ggplot2::labs(
      x = lab_x,
      y = lab_y,
      color="Reference sets", linetype="Library")+
    ggsci::scale_fill_jama() +
    ggsci::scale_color_d3() +
    ggplot2::scale_linetype_manual(values=c(2,1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "top",
          plot.margin = ggplot2::margin(2, 1, 2, 1, "cm"),
          legend.text = ggplot2::element_text(size=18),
          legend.title = ggplot2::element_text(size=18),
          axis.title.y = ggtext::element_markdown(size=18),
          axis.title.x = ggtext::element_markdown(size=18),
          axis.text = ggplot2::element_text(size=18)) +
    ggplot2::guides(alpha = FALSE)

  return(p)

  # ggplot2::ggplot(mat[is.na(mat$Plate_controls),], ggplot2::aes(X, Y)) +
  #   ggplot2::geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8) +
  #   ggplot2::stat_density_2d(data=mat[is.na(mat$Plate_controls),],
  #                   ggplot2::aes(), col="grey30", alpha=1) +
  #   ggplot2::stat_density_2d(
  #     data=mat[which(mat$Plate_controls=="+"),] |>
  #       dplyr::mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
  #     ggplot2::aes(fill=Gene, group=Gene, alpha = ..level..), geom = "polygon", col=NA) +
  #   ggplot2::geom_point(
  #     data=mat[which(!is.na(mat$Reference_sets) &
  #                      !grepl("ORF",mat$Reference_sets)),],
  #     ggplot2::aes(col=Reference_sets), size=1) +
  #   ggplot2::geom_point(data=mat_select, ggplot2::aes(), pch=21) +
  #   geom_text_repel(data=mat_select,aes(label=Gene),
  #                   force=2, size=4.2,
  #                   max.overlaps = Inf,
  #                   min.segment.length = 0,
  #                   segment.size = 0.3) +
  #   ggplot2::labs(x=gsub("\\.",":",
  #               gsub("_"," ",X)),
  #        y=gsub("\\.",":",gsub("_"," ",Y)),
  #        color="Reference sets",
  #        linetype="Library") +
  #   ggsci::scale_fill_jama() +
  #   ggsci::scale_color_d3() +
  #   ggplot2::scale_linetype_manual(values=c(2,1)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(axis.title = ggplot2::element_text(size=16),
  #                  axis.text = ggplot2::element_text(size=16),
  #                  panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank()) +
  #   guides(alpha = FALSE)
  #
  # ggplot(mat[is.na(mat$Plate_controls),], aes(X, Y)) +
  #   ggplot2::geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8) +
  #   ggplot2::stat_density_2d(data=mat[is.na(mat$Plate_controls),],
  #                   ggplot2::aes(), col="grey30", alpha=1) +
  #   ggplot2::stat_density_2d(data=mat[!is.na(mat$Plate_controls) &
  #                                       !is.na(mat$Gene) & mat$Gene != "WT",],
  #                   ggplot2::aes(fill=Gene, group=Gene, alpha = ..level..),
  #                   geom = "polygon", col=NA) +
  #   ggplot2::geom_point(data=mat[which(!is.na(mat$Reference_sets) &
  #                                        !grepl("ORF",mat$Reference_sets)),],
  #              ggplot2::aes(col=Reference_sets), size=1) +
  #   ggplot2::geom_point(data=mat_select, ggplot2::aes(), pch=21) +
  #   ggrepel::geom_text_repel(data=mat_select, ggplot2::aes(label=Gene),
  #                   force=2, size=4.2,
  #                   max.overlaps = Inf,
  #                   min.segment.length = 0,
  #                   segment.size = 0.3) +
  #   ggplot2::labs(
  #     x=gsub("\\.",":",gsub("_"," ",X)),
  #     y=gsub("\\.",":",gsub("_"," ",Y)),
  #     color="Reference sets") +
  #   ggsci::scale_fill_jama() +
  #   ggsci::scale_color_d3() +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(panel.grid.major = element_blank(),
  #                  panel.grid.minor = element_blank()) +
  #   guides(alpha = FALSE)
  #
  # ggplot2::ggplot(mat[is.na(mat$Plate_controls),], ggplot2::aes(X, Y)) +
  #   ggplot2::geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8) +
  #   ggplot2::stat_density_2d(data=mat[is.na(mat$Plate_controls),],
  #                   ggplot2::aes(lty=Type), col="grey30", alpha=1) +
  #   ggplot2::stat_density_2d(data=mat[which(mat$Plate_controls=="+"),] |>
  #                              dplyr::mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
  #                   ggplot2::aes(fill=Gene, group=Gene, alpha = ..level..),
  #                   geom = "polygon", col=NA) +
  #   ggplot2::geom_point(data=mat[which(!is.na(mat$Reference_sets) &
  #                                        !grepl("ORF",mat$Reference_sets)),],
  #              ggplot2::aes(col=Reference_sets), size=1) +
  #   ggplot2::geom_point(data=mat_select, ggplot2::aes(), pch=21) +
  #   ggrepel::geom_text_repel(data=mat_select,aes(label=Gene),
  #                   force=2, size=4.2,
  #                   max.overlaps = Inf,
  #                   min.segment.length = 0,
  #                   segment.size = 0.3) +
  #   ggplot2::labs(
  #     x=gsub("\\.",":",gsub("_"," ",X)),
  #     y=gsub("\\.",":",gsub("_"," ",Y)),
  #     color="Reference sets", linetype="Library") +
  #   ggsci::scale_fill_jama() +
  #   ggsci::scale_color_d3() +
  #   ggplot2::scale_linetype_manual(values=c(2,1)) +
  #   ggplot2::theme_bw() +
  #   ggplot2::theme(
  #     axis.title = ggplot2::element_text(size=16),
  #     axis.text = ggplot2::element_text(size=16),
  #     legend.text = ggplot2::element_text(size=14),
  #     legend.title = ggplot2::element_text(size=16),
  #     panel.grid.major = element_blank(),
  #     panel.grid.minor = element_blank()) +
  #   ggplot2::guides(alpha = FALSE)



}


plot_autophagy_competence <- function(competence_data = NULL,
                                      dnn_model = "30"){

  p <- NULL

  x_lab <-
    paste0("<br><b>Autophagosome formation</b><br><br>",
           "<i>log BFt (VAM6:ATG1), DNN-model '",
           dnn_model,"'</i>")
  y_lab <-
    paste0("<b>Autophagosome clearance</b><br><br>",
           "<i>log BFt (WT:VAM6), DNN-model '",
           dnn_model,"'</i><br>")

  if(dnn_model == "30"){
    p <- ggplot2::ggplot(
      competence_data$BF_response,
      ggplot2::aes(
        log_BFt_VAM6.ATG1_30,
        log_BFt_WT.VAM6_30, col=TimeR)) +
      ggplot2::geom_vline(xintercept = 0, lty=1, col="black", size=0.1) +
      ggplot2::geom_hline(yintercept = 0, lty=1, col="black", size=0.1) +
      scico::scale_color_scico(palette = 'lisbon') +
      ggplot2::scale_size(range = c(0, 1.5)) +
      #geom_path(data=BF_response_ctr, aes(),col="gray",lty=2, alpha = 1) +
      ggplot2::geom_point(
        data=competence_data$BF_response_ctr,
        ggplot2::aes(
          size=log_BFt_WT.ATG1_30, pch="Control"),size=1.8, alpha = 1) +
      ggplot2::geom_path(ggplot2::aes(),col="black",lty=2, alpha = 1) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = log_BFt_VAM6.ATG1_30.shift,
        yend = log_BFt_WT.VAM6_30.shift,
        size=log_BFt_WT.ATG1_30),
        arrow = arrow(
          angle = 35,
          length = unit(0.11, "inches"),
          type = "closed"),
        alpha = 1,
        size=1.3) +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        col="Time", size="BF (WT:ATG1)",
        title=paste(competence_data$id, competence_data$Library), shape="") +
      ggplot2::theme_bw(base_size = 20, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=20, face="bold"),
        axis.text = ggplot2::element_text(size=18),
        plot.margin = ggplot2::margin(2, 1, 2, 1, "cm"),
        axis.title.y = ggtext::element_markdown(size = 18),
        axis.title.x = ggtext::element_markdown(size = 18),
        legend.text = ggplot2::element_text(size=18, family = "Helvetica"),
        strip.background = ggplot2::element_blank(),
        legend.position = "top")
  }else{
    p <- ggplot2::ggplot(
      competence_data$BF_response,
      ggplot2::aes(log_BFt_VAM6.ATG1_22,
                   log_BFt_WT.VAM6_22,
                   col=TimeR)) +
      ggplot2::geom_vline(
        xintercept = 0, lty=1, col="black", size=0.1) +
      ggplot2::geom_hline(
        yintercept = 0, lty=1, col="black", size=0.1) +
      scico::scale_color_scico(palette = 'lisbon') +
      ggplot2::scale_size(range = c(0, 1.5)) +
      ggplot2::geom_point(
        data=competence_data$BF_response_ctr,
        ggplot2::aes(size=log_BFt_WT.ATG1_22, pch="Control"),size=1.8, alpha = 1) +
      ggplot2::geom_path(ggplot2::aes(),col="black",lty=2, alpha = 1) +
      ggplot2::geom_segment(ggplot2::aes(
        xend = log_BFt_VAM6.ATG1_22.shift,
        yend = log_BFt_WT.VAM6_22.shift,
        size=log_BFt_WT.ATG1_22),
        arrow = arrow(
          angle = 35,
          length = unit(0.11, "inches"),
          type = "closed"),
        alpha = 1, size=1.3) +
      ggplot2::labs(
        x = x_lab,
        y = y_lab,
        col="Time",
        size="BF (WT:ATG1)",
        title=paste(competence_data$id, competence_data$Library), shape="") +
      ggplot2::theme_bw(base_size = 20, base_family = "Helvetica") +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(size=20, face="bold"),
        axis.text = ggplot2::element_text(size=18),
        plot.margin = ggplot2::margin(2, 1, 2, 1, "cm"),
        axis.title.y = ggtext::element_markdown(size = 18),
        axis.title.x = ggtext::element_markdown(size = 18),
        legend.text = ggplot2::element_text(size=18),
        legend.position = "top")
  }
  return(p)


}

plot_response_kinetics_multi <- function(
    response_data = NULL,
    primary_identifiers = "RTG1 / YOL067C",
    Filter = TRUE,
    user_x = "T50 +N",
    user_y = "T50 -N",
    show_library_type_contour = FALSE,
    Value = "Perturbation"){

  # Kinetic parameters
  #Select one or several IDs for text repel
  #id=IDs_unique[grep("RTG",IDs_unique)]

  #Select x and y
  #X <- unique(response_data[['ds_parms']]$Parameter)[12]
  #Y <- unique(response_data[['ds_parms']]$Parameter)[11]
  X <- user_x
  Y <- user_y

  Positions <- c()
  for(i in primary_identifiers){
    y_pred <- response_data[['ds_curvefits']] |>
      dplyr::filter(primary_identifier == i)
      #subset(gsub("NA ", "", paste(Gene, ORF)) == i)
    #If mutant in rec plate, only evaluate rec mutants
    if(any(grepl("Rec",y_pred$Plate))){
      y_pred <- y_pred[which(grepl("Rec",y_pred$Plate)),]
    }
    #Find representative response
    y_pred <- y_pred %>%
      dplyr::group_by(Time) %>%
      dplyr::mutate(d=abs(P1_30_fit-median(P1_30_fit))) %>%
      dplyr::group_by(Plate, Position) %>%
      dplyr::mutate(d=mean(d))
    y_pred <- y_pred[which(y_pred$d==min(y_pred$d)),]
    Positions <- c(
      Positions,paste(y_pred$Plate, y_pred$Position)[1])
  }

  mat <- dcast(response_data[['ds_parms']][which(
    response_data[['ds_parms']]$Parameter %in% c(X,Y)),],
               Plate+Position+ORF+Gene+ Reference_sets~Parameter,
               value.var = Value)

  mat_select <- dcast(response_data[['ds_parms_comb']][which(
    response_data[['ds_parms_comb']]$Parameter %in% c(X,Y)),] |>
      subset(paste(Plate, Position) %in% Positions),
    Plate+Position+ORF+Gene + Reference_sets+Plate_controls~Parameter,
    value.var = Value)

  mat$X <- mat[,X]
  mat$Y <- mat[,Y]
  mat_select$X <- mat_select[,X]
  mat_select$Y <- mat_select[,Y]

  #Filter, replace manual control of scale?
  if(Filter){
    x_lower_bound <- as.numeric(quantile(mat$X, 1-0.99, na.rm = T)) -
      5 * IQR(mat$X, na.rm = T)
    x_upper_bound <- as.numeric(quantile(mat$X, 0.99, na.rm = T)) +
      5 * IQR(mat$X, na.rm = T)
    y_lower_bound <- as.numeric(quantile(mat$Y, 1-0.99, na.rm = T)) -
      5 * IQR(mat$Y, na.rm = T)
    y_upper_bound <- as.numeric(quantile(mat$Y, 0.99, na.rm = T)) +
      5 * IQR(mat$Y, na.rm = T)
  }else{
    x_lower_bound <- min(mat$X, na.rm = T)
    x_upper_bound <- max(mat$X, na.rm = T)
    y_lower_bound <- min(mat$Y, na.rm = T)
    y_upper_bound <- max(mat$Y, na.rm = T)
  }

  mat$Type <- response_data[['dnn_preds']]$Type[match(paste(
    mat$Gene, mat$ORF),paste(
      response_data[['dnn_preds']]$Gene,
      response_data[['dnn_preds']]$ORF))]

  p <- ggplot2::ggplot(mat, ggplot2::aes(X, Y))+
    ggplot2::geom_point(col="lightgray", size=1.4, pch=19, alpha=0.8)+
    ggplot2::geom_point(
      data=mat[which(!is.na(mat$Reference_sets)),],
      ggplot2::aes(col=Reference_sets), size=1.5)+
    ggplot2::geom_point(
      data=mat_select, ggplot2::aes(), pch=19)+
    ggrepel::geom_text_repel(
      data=mat_select,
      ggplot2::aes(label=Gene),
      force=2, size=5.5,
      nudge_x = 0.3, nudge_y = 0.3,
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.size = 0.9)+
    ggplot2::xlim(x_lower_bound,x_upper_bound)+
    ggplot2::ylim(y_lower_bound,y_upper_bound)+
    ggplot2::labs(x=X, y=Y) +
    ggsci::scale_fill_jama() +
    ggsci::scale_color_d3() +
    ggplot2::theme_bw(base_size = 18, base_family = "Helvetica") +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size=18),
      axis.text = ggplot2::element_text(size=18),
      #plot.margin = ggplot2::margin(1, 1, 2, 2, "cm"),
      legend.text = ggplot2::element_text(size=18),
      plot.margin = ggplot2::margin(1, 1, 1, 1, "cm"),
      legend.position = "top",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())

  if(show_library_type_contour == T){
    p <- p + ggplot2::stat_density_2d(
      data=mat[,],
      ggplot2::aes(lty=Type), col="grey30", alpha=1)
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

  p <- ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept=0, lty=2) +
    ggplot2::geom_vline(xintercept=12, lty=2) +
    ggplot2::geom_hline(
      yintercept=0, lty=1, linewidth=size_line) +
    ggplot2::geom_hline(
      yintercept=double_sigmoidal_model[,"maximum_y"],
      lty=1, linewidth=size_line) +
    ggplot2::geom_hline(
      yintercept=double_sigmoidal_model[,"finalAsymptoteIntensity"],
      lty=1, linewidth=size_line) +
    ggplot2::geom_segment(
      ggplot2::aes(x=double_sigmoidal_model[,"startPoint_x"], y=0,
                   xend=double_sigmoidal_model[,"reachMaximum_x"],
                   yend=double_sigmoidal_model[,"reachMaximum_y"] ),
      lty=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=double_sigmoidal_model[,"startDeclinePoint_x"],
                   y=double_sigmoidal_model[,"startDeclinePoint_y"],
                   xend=double_sigmoidal_model[,"endDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"endDeclinePoint_y"] ),
      lty=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=double_sigmoidal_model[,"midPoint1_y"],
                   xend=double_sigmoidal_model[,"midPoint1_x"],
                   yend=double_sigmoidal_model[,"midPoint1_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"midPoint2_y"],
                   xend=double_sigmoidal_model[,"midPoint2_x"],
                   yend=double_sigmoidal_model[,"midPoint2_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=0,
                   xend=double_sigmoidal_model[,"startPoint_x"],
                   yend=0),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"maximum_y"],
                   xend=double_sigmoidal_model[,"startDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"maximum_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=0, y=double_sigmoidal_model[,"reachMaximum_y"],
                   xend=double_sigmoidal_model[,"reachMaximum_x"],
                   yend=double_sigmoidal_model[,"reachMaximum_y"]),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    ggplot2::geom_segment(
      ggplot2::aes(x=12, y=double_sigmoidal_model[,"endDeclinePoint_y"],
                   xend=double_sigmoidal_model[,"endDeclinePoint_x"],
                   yend=double_sigmoidal_model[,"endDeclinePoint_y"] ),
      lineend="round",arrow=arrow(length=unit(arrow_length, "inches"), ends = "both"),
      alpha=1, size=size_segment) +
    geom_line(
      ggplot2::aes(x= y_pred.ctr$Time, y=y_pred.ctr$P1_30_fit, col="Control")) +
    geom_ribbon(
      ggplot2::aes(x=y_pred.ctr$Time,
                   y=colMeans(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                   ymin=colMins(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                   ymax=colMaxs(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit))),
      fill="yellow", alpha=0.1) +
    geom_ribbon(
      ggplot2::aes(x=c(double_sigmoidal_model[,"startPoint_x"],
                       double_sigmoidal_model[,"reachMaximum_x"]),
                   y=c(0, double_sigmoidal_model[,"reachMaximum_y"]),
                   ymin=c(0,0),
                   ymax=c(double_sigmoidal_model[,"reachMaximum_y"],
                          double_sigmoidal_model[,"reachMaximum_y"])),fill="black",
      alpha=0.1) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x=c(double_sigmoidal_model[,"startDeclinePoint_x"],
                       double_sigmoidal_model[,"endDeclinePoint_x"]),
                   y=c(double_sigmoidal_model[,"startDeclinePoint_y"],
                       double_sigmoidal_model[,"endDeclinePoint_y"]),
                   ymin=c(double_sigmoidal_model[,"endDeclinePoint_y"],double_sigmoidal_model[,"endDeclinePoint_y"]),
                   ymax=c(double_sigmoidal_model[,"startDeclinePoint_y"],
                          double_sigmoidal_model[,"startDeclinePoint_y"])),fill="black",
      alpha=0.1) +
    ggplot2::geom_line(
      ggplot2::aes(x= y_pred$Time, y=y_pred$P1_30_fit, col="Model")) +
    ggplot2::geom_point(
      ggplot2::aes(x=ifelse(
        grepl("Starv",y_raw$Phase),
        y_raw$TimeR,y_raw$TimeR+12),y_raw$P1_30, col="Data")) +
    ggplot2::ylim(0,100) +
    ggplot2::xlim(-1,20) +
    ggplot2::labs(x="Time (hours)", y="Autophagy (%)", title=paste(id, slibrary)) +
    ggplot2::scale_color_manual(
      values=ggsci::pal_jco("default", alpha = 1)(7)[c(7,1,4)]) +
    ggplot2::theme_bw(base_size = 18, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size=20, face="bold"),
      axis.title = ggplot2::element_text(size=18),
      axis.text = ggplot2::element_text(size=18),
      plot.margin = ggplot2::margin(1, 1, 2, 2, "cm"),
      legend.text = ggplot2::element_text(size=18),
      legend.title = ggplot2::element_blank()
    )

  return(p)
}

