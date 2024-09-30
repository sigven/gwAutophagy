## Script to prepare data, utilize https://www.fstpackage.org/ for faster read/write

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - gwAutophagyApp - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <-
  log4r::logger(
    threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("GW_AUTOPHAGY_LOGGER" = log4r_logger)

autophagy_competence_dataprep_per_gene <- function(id, gw_autoph_data = NULL){

  BF_response <- gw_autoph_data$df_BF_temporal |>
    subset(gsub("NA ", "", paste(Gene, ORF)) == id)

    #If mutant in rec plate, only evaluate rec mutants
  if(any(grepl("Rec",BF_response$Plate))){
    BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
  }
  #Find representative response
  BF_response <- BF_response |>
    group_by(TimeR) |>
    mutate(d=abs(log_BFt_WT.ATG1_30-median(log_BFt_WT.ATG1_30))+
             abs(log_BFt_WT.VAM6_30-median(log_BFt_WT.VAM6_30))+
             abs(log_BFt_VAM6.ATG1_30-median(log_BFt_VAM6.ATG1_30))+
             abs(log_BFt_WT.ATG1_22-median(log_BFt_WT.ATG1_22))+
             abs(log_BFt_WT.VAM6_22-median(log_BFt_WT.VAM6_22))+
             abs(log_BFt_VAM6.ATG1_22-median(log_BFt_VAM6.ATG1_22))) |>
    dplyr::group_by(Plate, Position) |>
    dplyr::mutate(d = mean(d))
  BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
  BF_response_ctr <- df_BF_temporal[df_BF_temporal$Plate == BF_response$Plate[1],] |>
    dplyr::group_by(TimeR) |>
    dplyr::summarise_if(is.numeric, mean)


  Time.shift <- unique(BF_response$TimeR) + 1
  names(Time.shift) <- unique(BF_response$TimeR)
  Time.shift[which.max(Time.shift)] <- 0

  BF_response$shift <- Time.shift[as.character(BF_response$TimeR)]
  BF_response <- data.frame(BF_response)
  BF_response <- BF_response |>
    mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
           log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
           log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
           log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
           log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
           log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)
  BF_response_ctr$shift <- Time.shift[as.character(BF_response_ctr$TimeR)]
  BF_response_ctr <- data.frame(BF_response_ctr)
  BF_response_ctr <- BF_response_ctr |>
    mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
           log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
           log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
           log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
           log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
           log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)

  Library <- df_DNN_preds$Type[gsub("NA ", "",
                                    paste(gw_autoph_data$df_DNN_preds$Gene, gw_autoph_data$df_DNN_preds$ORF)) == id][1]

  return(list(
    BF_response = BF_response,
    BF_response_ctr = BF_response_ctr,
    Library = Library,
    id = id
  ))

}

#autophagy_execution_dataprep_per_gene
#kinetic_response_dataprep_per_gene
kinetic_response_dataprep_per_gene <- function(id = "RTG3 YBL103C", gw_autoph_data = NULL){
  y_pred <- gw_autoph_data$df_DS_curvefits |>
    dplyr::mutate(
      id_gene_orf = stringr::str_replace(
        paste(.data$Gene, .data$ORF),"NA ","")) |>
    dplyr::filter(.data$id_gene_orf == id)

  ## if NROW(y_pred == 0)?

  #If mutant in rec plate, only evaluate rec mutants
  if(any(grepl("Rec",y_pred$Plate))){
    y_pred <- y_pred |>
      dplyr::filter(stringr::str_detect(.data$Plate,"Rec"))
  }
  #else what?

  #Find representative response
  y_pred <- y_pred |>
    dplyr::group_by(Time) |>
    dplyr::mutate(d=abs(P1_30_fit-median(P1_30_fit))) |>
    dplyr::group_by(Plate, Position) |>
    dplyr::mutate(d=mean(d))
  y_pred <- y_pred[which(y_pred$d==min(y_pred$d)),]
  y_raw <- gw_autoph_data$df_DNN_preds |>
    dplyr::filter(
      paste(Plate, Position) == paste(y_pred$Plate, y_pred$Position)[1])

  slibrary <- ifelse(y_raw$Type== "KO" | !(is.na(y_raw$Plate_controls)),"KO","DAmP")[1]
  plate_id <- unique(y_raw$Plate)
  y_pred.ctr <- gw_autoph_data$df_DS_curvefits[
    which(is.na(gw_autoph_data$df_DS_curvefits$Plate_controls)),] |>
    #subset(Type == Library) |>
    dplyr::filter(
      .data$Plate == plate_id) |>
    dplyr::group_by(.data$Time) |>
    dplyr::summarise(
      P1_30_fit=median(.data$P1_30_fit), .groups = "drop")

  plate_position <- toupper(stringr::str_replace_all(
    paste(y_pred$Plate, y_pred$Position)[1]," ",""))
  double_sigmoidal_model <- gw_autoph_data$regr_df_all |>
    dplyr::filter(startsWith(toupper(.data$id2), plate_position))

  return(list(
    y_pred = y_pred,
    y_raw = y_raw,
    y_pred.ctr = y_pred.ctr,
    double_sigmoidal_model = double_sigmoidal_model,
    slibrary = slibrary,
    id = id
  ))
}

## Read autophagy predictions
excel_file_autophagy_predictions <- file.path(
  here::here(),"data-raw","Chica_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx")
load(file.path(here::here(), "data-raw","Double sigmoidal curve fitting_v2.RData"))

df_DNN_preds <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 2)
df_DS_curvefits <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 3)
df_DS_parms <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 5)
df_DS_parms_ctrs <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 6)
df_DS_parms_comb <- dplyr::bind_rows(
  df_DS_parms, df_DS_parms_ctrs)

IDs_unique <- gsub("NA ", "", unique(paste(df_DS_parms_comb$Gene,df_DS_parms_comb$ORF)))
gene_ids_kinetic <- data.frame('id' = IDs_unique)


#Read Bayes factor data
excel_file_bfactors <- file.path(
  here::here(),"data-raw","Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx")

df_BF_overall <- read_excel(
  excel_file_bfactors, sheet = 2)
df_BF_starv <- read_excel(
  excel_file_bfactors, sheet = 3)
df_BF_repl <- read_excel(
  excel_file_bfactors, sheet = 4)
df_BF_temporal <- read_excel(
  excel_file_bfactors, sheet = 5)

#ID list
IDs_unique <- gsub("NA ", "", unique(paste(df_BF_repl$Gene,df_BF_repl$ORF)))
gene_ids_bf <- data.frame('id' = IDs_unique)


regr_df_all <- data.frame()
for(n in names(regression.results)){
  if(is.data.frame(regression.results[[n]])){
    regr_df <- as.data.frame(regression.results[[n]]) |>
      dplyr::mutate(id = n, id2 = stringr::str_replace_all(n," ","")) |>
      dplyr::select(c("id", "id2", dplyr::everything()))

    regr_df_all <- dplyr::bind_rows(regr_df_all, regr_df)
  }
}


gw_autoph_data <- list()
gw_autoph_data$df_DS_curvefits <- df_DS_curvefits
gw_autoph_data$df_DNN_preds <- df_DNN_preds
gw_autoph_data$regr_df_all <- regr_df_all
gw_autoph_data$df_BF_overall <- df_BF_overall
gw_autoph_data$df_BF_starv <- df_BF_starv
gw_autoph_data$df_BF_repl <- df_BF_repl
gw_autoph_data$df_BF_temporal <- df_BF_temporal



gw_autoph_kinetic_plot_input <- list()
gw_autoph_competence_plot_input <- list()

i <- 1
for(id in unique(gene_ids_kinetic$id)){
  gw_autoph_kinetic_plot_input[[id]] <-
    kinetic_response_dataprep_per_gene(id = id, gw_autoph_data = gw_autoph_data)
  if(i %% 100 == 0){
    log4r::info(log4r_logger, msg = paste0(
      "Preparing kinetic plot input - processed ",i, " gene identifiers.. "))
  }
  i <- i + 1
}

i <- 1
for(id in unique(gene_ids_bf$id)){
  gw_autoph_competence_plot_input[[id]] <-
    autophagy_competence_dataprep_per_gene(id = id, gw_autoph_data = gw_autoph_data)
  if(i %% 100 == 0){
    log4r::info(log4r_logger, msg = paste0(
      "Preparing autophagy competence plot input - processed ",i, " gene identifiers.. "))
  }
  i <- i + 1
}

saveRDS(gw_autoph_competence_plot_input, file.path(
  here::here(), "data","processed","gw_autoph_competence_plot_input.rds"))
saveRDS(gw_autoph_kinetic_plot_input, file.path(
  here::here(), "data","processed","gw_autoph_kinetic_plot_input.rds"))
fst::write_fst(
  df_DNN_preds, file.path(here::here(), "data","processed","df_DNN_preds.fst"))
fst::write_fst(
  df_DS_curvefits, file.path(here::here(), "data","processed","df_DS_curvefits.fst"))
fst::write_fst(
  df_DS_parms, file.path(here::here(), "data","processed","df_DS_parms.fst"))
fst::write_fst(
  df_DS_parms_ctrs, file.path(here::here(), "data","processed","df_DS_parms_ctrs.fst"))
fst::write_fst(
  regr_df_all, file.path(here::here(), "data","processed","regr_df_all.fst"))
fst::write_fst(
  gene_ids_kinetic, file.path(here::here(), "data","processed","gene_ids_kinetic.fst")
)

fst::write_fst(
  gene_ids_bf, file.path(here::here(), "data","processed","gene_ids_bf.fst")
)
fst::write_fst(
  df_BF_overall, file.path(here::here(), "data","processed","df_BF_overall.fst"))
fst::write_fst(
  df_BF_starv, file.path(here::here(), "data","processed","df_BF_starv.fst"))
fst::write_fst(
  df_BF_repl, file.path(here::here(), "data","processed","df_BF_repl.fst"))
fst::write_fst(
  df_BF_temporal, file.path(here::here(), "data","processed","df_BF_temporal.fst"))


