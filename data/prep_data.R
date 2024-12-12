## Script to prepare data for the AutoDRY Shiny app,
## utilize https://www.fstpackage.org/ for faster read/write

my_log4r_layout <- function(level, ...) {
  paste0(format(Sys.time()), " - gwAutophagyApp - ",
         level, " - ", ..., "\n", collapse = "")
}

log4r_logger <-
  log4r::logger(
    threshold = "INFO", appenders = log4r::console_appender(my_log4r_layout))

# this gets passed on to all the log4r_* functions inside the pkg
options("GW_AUTOPHAGY_LOGGER" = log4r_logger)


clean_internal_gene_orfs <- function(dataset = NULL){

  dataset <- dataset |>
    dplyr::mutate(record_idx = dplyr::row_number())

  ## Case 0: Only Gene is present (ORF is missing) - WT
  tmp0 <- dataset |>
    dplyr::filter(is.na(ORF) & !is.na(Gene))

  ## Case 1: Both Gene and ORF are present - confident entry
  tmp1 <- dataset |>
    dplyr::filter(!is.na(ORF) & !is.na(Gene)) |>
    dplyr::mutate(Gene2 = Gene, ORF2 = ORF)

  ## Case 2: Only ORF is present - Gene is missing)
  ## - here, we exclude records for which records with the same ORF
  ##  (and Gene) are part of tmp1 dataset; some 44 ORFs carry this property
  ##   Rationale: ORF is a unique identifier, users will likely get confused
  ##.           if they can choose between 'YCL023C' and 'YCL023C / ATG22'
  tmp2 <- dataset |>
    dplyr::filter(!is.na(ORF) & is.na(Gene)) |>
    dplyr::anti_join(tmp1, by = "ORF") |>
    dplyr::mutate(ORF2 = ORF, Gene2 = Gene) |>
    dplyr::distinct()

  if(NROW(tmp0) > 0){
    tmp0 <- tmp0 |>
      dplyr::mutate(Gene2 = NA, ORF2 = "WT")
  }

  dataset_final <- data.frame()
  if(NROW(tmp0) > 0){
    dataset_final <-
      dplyr::bind_rows(
        tmp0, tmp1, tmp2)
  }else{
    dataset_final <-
      dplyr::bind_rows(
        tmp1, tmp2)
  }
  dataset_final <-
   dataset_final |>
    dplyr::arrange(record_idx) |>
    dplyr::select(-c("record_idx")) |>
    dplyr::mutate(
      screen_id = gsub(
        "NA / ", "",
        paste(Gene2, ORF2, sep = " / "))) |>
    dplyr::arrange(screen_id)

  ## skip duplicated/merged mutants
  merged_duplicate_mutants <-
    readr::read_tsv("data-raw/duplicated_merged_mutants.txt",
                    show_col_types = F, col_names = F)
  colnames(merged_duplicate_mutants) <- c("screen_id")
  dataset_final <- dataset_final |>
    dplyr::anti_join(
      merged_duplicate_mutants, by = "screen_id")

  return(dataset_final)

}

autophagy_competence_per_ko <- function(id, gw_autoph_data = NULL){

  #Get data for the mutant
  BF_response <- gw_autoph_data$BF_temporal |>
    dplyr::filter(primary_identifier == id)

    #If mutant in rec plate, only evaluate rec mutants
  if(any(grepl("Rec",BF_response$Plate))){
    BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
  }
  #Find representative response
  BF_response <- BF_response |>
    dplyr::group_by(TimeR) |>
    dplyr::mutate(d=abs(log_BFt_WT.ATG1-median(log_BFt_WT.ATG1))+
             abs(log_BFt_WT.VAM6 - median(log_BFt_WT.VAM6))+
             abs(log_BFt_VAM6.ATG1 - median(log_BFt_VAM6.ATG1))) |>
    dplyr::group_by(Plate, Position) |>
    dplyr::mutate(d = mean(d), NCells = mean(NCells))
  BF_response <- BF_response[which(BF_response$d == min(BF_response$d)),]
  BF_response <- BF_response[which(BF_response$NCells == max(BF_response$NCells)),]

  BF_response_ctr <- gw_autoph_data$BF_temporal[
    gw_autoph_data$BF_temporal$Plate == BF_response$Plate[1],] |>
    dplyr::group_by(TimeR) |>
    dplyr::summarise_if(is.numeric, mean)

  if("entrezgene" %in% colnames(BF_response_ctr)){
    BF_response_ctr <- BF_response_ctr |>
      dplyr::select(-c("entrezgene"))
  }


  Time.shift <- unique(BF_response$TimeR) + 1
  names(Time.shift) <- unique(BF_response$TimeR)
  Time.shift[which.max(Time.shift)] <- 0

  BF_response$shift <- Time.shift[as.character(BF_response$TimeR)]
  BF_response <- data.frame(BF_response)
  BF_response <- BF_response |>
    dplyr::mutate(log_BFt_WT.ATG1.shift = (log_BFt_WT.ATG1[match(shift, TimeR)]-log_BFt_WT.ATG1)/3+log_BFt_WT.ATG1,
           log_BFt_WT.VAM6.shift =  (log_BFt_WT.VAM6[match(shift, TimeR)]-log_BFt_WT.VAM6)/3 + log_BFt_WT.VAM6,
           log_BFt_VAM6.ATG1.shift =  (log_BFt_VAM6.ATG1[match(shift, TimeR)]-log_BFt_VAM6.ATG1)/3 + log_BFt_VAM6.ATG1)
  BF_response_ctr$shift <- Time.shift[as.character(BF_response_ctr$TimeR)]
  BF_response_ctr <- data.frame(BF_response_ctr)
  BF_response_ctr <- BF_response_ctr |>
    dplyr::mutate(log_BFt_WT.ATG1.shift = (log_BFt_WT.ATG1[match(shift, TimeR)]-log_BFt_WT.ATG1)/3+log_BFt_WT.ATG1,
           log_BFt_WT.VAM6.shift =  (log_BFt_WT.VAM6[match(shift, TimeR)]-log_BFt_WT.VAM6)/3 + log_BFt_WT.VAM6,
           log_BFt_VAM6.ATG1.shift =  (log_BFt_VAM6.ATG1[match(shift, TimeR)]-log_BFt_VAM6.ATG1)/3 + log_BFt_VAM6.ATG1)
  Library <- gw_autoph_data$dnn_preds[gw_autoph_data$dnn_preds$primary_identifier == id,]$Type[1]

  return(list(
    BF_response = BF_response,
    BF_response_ctr = BF_response_ctr,
    Library = Library,
    id = id
  ))

}

## Prepare kinetic response data per gene mutant
kinetic_response_per_ko <- function(
    id = "RTG3 / YBL103C", gw_autoph_data = NULL){

  y_pred <- gw_autoph_data$ds_curvefits |>
    dplyr::filter(primary_identifier == id)
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
  y_raw <- gw_autoph_data$dnn_preds |>
    subset(paste(Plate, Position) %in% paste(y_pred$Plate, y_pred$Position)) |>
    dplyr::group_by(Plate, Position) |>
    dplyr::mutate(NCells = mean(NCells ))
  y_raw <- y_raw[which(y_raw$NCells==max(y_raw$NCells)),]
  y_pred <- y_pred |>
    subset(paste(Plate, Position) %in% paste(y_raw$Plate, y_raw$Position))

  slibrary <- ifelse(y_raw$Type== "KO" | !(is.na(y_raw$Plate_controls)),"KO","DAmP")[1]
  plate_id <- unique(y_raw$Plate)
  y_pred.ctr <- gw_autoph_data$ds_curvefits[
    which(is.na(gw_autoph_data$ds_curvefits$Plate_controls)),] |>
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


map_yeast_gene_identifiers <- function(
    screen_ids = NULL,
    update = T,
    datestamp_genome_alliance = "v7.4.0_2024-10-09"){

  datestamp <- as.character(Sys.Date())
  destfile_gene_info <- file.path(
    here::here(), "data-raw",
    paste0("Saccharomyces_cerevisiae.gene_info_",
           datestamp, ".gz"))
  destfile_gene_description <- file.path(
    here::here(), "data-raw",
    paste0("GENE-DESCRIPTION-TSV_SGD_",
           datestamp_genome_alliance,".tsv.gz"))
  if(update == T){
      remote_url1 = paste0(
        "https://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/",
        "Fungi/Saccharomyces_cerevisiae.gene_info.gz")
      if(!file.exists(destfile_gene_info)){
        download.file(remote_url1, destfile_gene_info)
      }
      remote_url2 = paste0(
        "https://fms.alliancegenome.org/",
        "download/GENE-DESCRIPTION-TSV_SGD.tsv.gz")
      if(!file.exists(destfile_gene_description)){
        download.file(remote_url2, destfile_gene_description)
      }

  }

  ncbi_gene_info <- readr::read_tsv(
    file=destfile_gene_info,
    show_col_types = F) |>
    janitor::clean_names() |>
    tidyr::separate_rows(synonyms, sep="\\|") |>
    dplyr::rename(entrezgene = gene_id) |>
    dplyr::mutate(id_gene_orf = paste(
      symbol, locus_tag, sep=" / ")) |>
    dplyr::mutate(gene_alias = dplyr::if_else(
      symbol == locus_tag, symbol,
      as.character(id_gene_orf))) |>
    dplyr::mutate(sgd_id = stringr::str_replace_all(
      db_xrefs,"\\|.+$","")) |>
    dplyr::select(
      entrezgene, symbol, locus_tag,
      description, synonyms,
      gene_alias, sgd_id) |>
    dplyr::rename(genename = description) |>
    dplyr::distinct()


  ## Get different gene designations/aliases
  ## from NCBI's gene info resources
  primary_gene_info <- ncbi_gene_info |>
    dplyr::select(gene_alias, sgd_id,
                  entrezgene, genename) |>
    dplyr::distinct() |>
    dplyr::mutate(alias_type = "symbol_locus_tag") |>
    dplyr::mutate(primary_identifier = gene_alias) |>
    dplyr::filter(
      !stringr::str_detect(gene_alias, "NEWENTRY"))

  gene_alias0 <- ncbi_gene_info |>
    dplyr::select(symbol, sgd_id, entrezgene) |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::rename(gene_alias = symbol) |>
    dplyr::mutate(alias_type = "symbol") |>
    dplyr::left_join(
      dplyr::select(
        primary_gene_info, primary_identifier,
        sgd_id, entrezgene, genename),
      by = c("sgd_id","entrezgene")
    ) |>
    dplyr::distinct()

  gene_alias1 <- ncbi_gene_info |>
    dplyr::select(locus_tag, sgd_id, entrezgene) |>
    dplyr::filter(locus_tag != "-") |>
    dplyr::rename(gene_alias = locus_tag) |>
    dplyr::mutate(gene_alias = paste(
      gene_alias, paste(gene_alias, gene_alias, sep=" / "),sep=";")) |>
    tidyr::separate_rows(gene_alias, sep=";") |>
    dplyr::mutate(alias_type = "locus_tag") |>
    dplyr::left_join(
      dplyr::select(
        primary_gene_info, primary_identifier,
        sgd_id, entrezgene, genename),
      by = c("sgd_id","entrezgene")
    ) |>
    dplyr::distinct()

  gene_alias2 <- ncbi_gene_info |>
    dplyr::filter(synonyms != "-") |>
    dplyr::select(synonyms, symbol, locus_tag, sgd_id, entrezgene) |>
    dplyr::mutate(gene_alias = dplyr::if_else(
      !startsWith(synonyms, stringr::str_sub(symbol,0,4)),
      paste(synonyms, locus_tag, sep=" / "),
      as.character(NA))) |>
    dplyr::filter(!is.na(gene_alias)) |>
    dplyr::select(gene_alias, sgd_id, entrezgene) |>
    dplyr::mutate(alias_type = "synonym_locus_tag") |>
    dplyr::left_join(
      dplyr::select(
        primary_gene_info, primary_identifier,
        sgd_id, entrezgene, genename),
      by = c("sgd_id","entrezgene")
    )

  gene_alias3 <- ncbi_gene_info |>
    dplyr::filter(synonyms != "-") |>
    dplyr::select(synonyms, sgd_id, entrezgene) |>
    dplyr::mutate(gene_alias = synonyms) |>
    dplyr::filter(!is.na(gene_alias)) |>
    dplyr::select(gene_alias, sgd_id, entrezgene) |>
    dplyr::mutate(alias_type = "synonym") |>
    dplyr::left_join(
      dplyr::select(
        primary_gene_info, primary_identifier,
        sgd_id, entrezgene, genename),
      by = c("sgd_id","entrezgene")
    )

  ## Get gene descriptors from Genome Alliance
  sgd_descriptors <- readr::read_tsv(
    file = destfile_gene_description,
    show_col_types = F,
    comment = "#", col_names = F) |>
    dplyr::distinct()
  colnames(sgd_descriptors) <-
    c("sgd_id","gene_alias","sgd_description")
  #sgd_descriptors$sgd_description <- NULL
  gene_alias4 <- sgd_descriptors |>
    dplyr::select(-c("sgd_description")) |>
    dplyr::mutate(alias_type = "sgd_designator") |>
    dplyr::left_join(
      dplyr::select(
        primary_gene_info, primary_identifier,
        sgd_id, entrezgene, genename),
      by = c("sgd_id")
    ) |>
    dplyr::distinct()

  all_gene_info <- dplyr::bind_rows(
    primary_gene_info, gene_alias0, gene_alias1,
    gene_alias2, gene_alias3, gene_alias4) |>
    dplyr::arrange(entrezgene) |>
    dplyr::group_by(dplyr::across(-c("alias_type"))) |>
    dplyr::summarise(alias_type = paste(unique(alias_type), collapse=";"),
                     .groups = "drop") |>
    dplyr::distinct()

  dup_aliases <-
    plyr::count(all_gene_info$gene_alias) |>
    dplyr::filter(freq > 1) |>
    dplyr::rename(gene_alias = x) |>
    dplyr::select(gene_alias)

  all_gene_info <- all_gene_info |>
    dplyr::anti_join(dup_aliases, by = "gene_alias") |>
    dplyr::distinct()

  id_map <- screen_ids |>
    dplyr::left_join(all_gene_info,
                     by = c("screen_id" = "gene_alias")) |>
    dplyr::distinct() |>
    dplyr::mutate(primary_identifier = dplyr::if_else(
      !is.na(sgd_id) & is.na(primary_identifier),
      screen_id,
      primary_identifier
    )) |>
    dplyr::left_join(
      dplyr::select(
        sgd_descriptors, sgd_id, sgd_description),
      by = "sgd_id"
    )


  not_found <- id_map |>
    dplyr::filter(
      is.na(sgd_id)) |>
    dplyr::distinct()

  rescued_set1 <- suppressWarnings(not_found |>
    dplyr::select(-c("entrezgene","primary_identifier",
                     "genename", "alias_type",
                     "sgd_description","sgd_id")) |>
    tidyr::separate(screen_id, c("id1","id2"),
                    sep =" / ", remove = F) |>
    dplyr::left_join(all_gene_info,
                     by = c("id1" = "gene_alias")) |>
    dplyr::filter(!is.na(sgd_id)) |>
    dplyr::select(-c("id1","id2")) |>
    dplyr::left_join(
      dplyr::select(
        sgd_descriptors, sgd_id, sgd_description),
      by = "sgd_id"
    ) |>
    dplyr::mutate(
      primary_identifier = dplyr::if_else(
        is.na(primary_identifier),
        screen_id,
        primary_identifier
      )
    ))

  not_found <- not_found |>
    dplyr::anti_join(rescued_set1, by = "screen_id")

  rescued_set2 <- suppressWarnings(
    not_found |>
      dplyr::select(
        -c("entrezgene","primary_identifier",
           "genename", "alias_type",
           "sgd_description","sgd_id")) |>
      tidyr::separate(screen_id, c("id1","id2"),
                      sep =" / ", remove = F) |>
      dplyr::left_join(all_gene_info,
                       by = c("id2" = "gene_alias")) |>
      dplyr::filter(!is.na(sgd_id)) |>
      dplyr::select(-c("id1","id2")) |>
      dplyr::left_join(
        dplyr::select(
          sgd_descriptors, sgd_id, sgd_description),
        by = "sgd_id"
      ) |>
      dplyr::mutate(
        primary_identifier = dplyr::if_else(
          is.na(primary_identifier),
          screen_id,
          primary_identifier
        )
      ))

  not_found <- not_found |>
    dplyr::anti_join(rescued_set2, by = "screen_id") |>
    dplyr::mutate(primary_identifier = screen_id) |>
    dplyr::mutate(sgd_description = "No description available")

  id_map_final <- dplyr::bind_rows(
    dplyr::filter(id_map, !is.na(sgd_id)),
    rescued_set1, rescued_set2, not_found) |>
    dplyr::distinct() |>
    dplyr::mutate(genename = dplyr::if_else(
      is.na(genename),
      "No genename available",
      genename
    )) |>
    dplyr::mutate(sgd_description = dplyr::if_else(
      is.na(sgd_description),
      paste0(
        "No description available [<i>Alliance ",
        "of Genome Resources/SGD - ", datestamp_genome_alliance,
        "</i>]"),
      paste0(sgd_description,
             " [<i>Alliance of Genome Resources/SGD - ",
             datestamp_genome_alliance, "</i>]")
    )) |>
    dplyr::distinct() |>
    dplyr::arrange(primary_identifier)

  return(id_map_final)


}

####--- Start preparing data ----####

datestamp_genome_alliance <- "v7.4.0_2024-10-09"

####--- 1. Read autophagy predictions ---####

excel_file_autophagy_predictions <- file.path(
  here::here(),"data-raw",
  "Chica_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx")
load(file.path(here::here(),
               "data-raw","Double sigmoidal curve fitting_v2.RData"))


autophagy_preds <- list()

autophagy_preds[['dnn']] <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 2)
autophagy_preds[['ds_curvefits']] <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 3)
autophagy_preds[['ds_parms']] <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 5)
autophagy_preds[['ds_parms_ctrs']] <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 6)
autophagy_preds[['ds_parms_comb']] <- dplyr::bind_rows(
  autophagy_preds[['ds_parms']], autophagy_preds[['ds_parms_ctrs']])
autophagy_preds[['profile']] <- readxl::read_excel(
  excel_file_autophagy_predictions, sheet = 7)

for(elem in c('dnn','ds_curvefits','ds_parms','profile',
              'ds_parms_ctrs','ds_parms_comb')){
  autophagy_preds[[elem]] <-
    clean_internal_gene_orfs(dataset = autophagy_preds[[elem]])
  autophagy_preds[[elem]] <-
    map_yeast_gene_identifiers(screen_ids = autophagy_preds[[elem]])
  if(elem != 'ds_parms_comb'){
    autophagy_preds[[elem]] <- autophagy_preds[[elem]] |>
      dplyr::select(
        -dplyr::any_of(c("sgd_id","entrezgene",
                       "genename","alias_type","sgd_description")))
  }
}
gene_info_kinetic <-
  data.frame(
    'orf_gene_id' = autophagy_preds[['ds_parms_comb']]$primary_identifier,
    'genename' = autophagy_preds[['ds_parms_comb']]$genename,
    'entrezgene' = autophagy_preds[['ds_parms_comb']]$entrezgene,
    'screen_id' = autophagy_preds[['ds_parms_comb']]$screen_id,
    'sgd_id' = autophagy_preds[['ds_parms_comb']]$sgd_id,
    'Gene' = autophagy_preds[['ds_parms_comb']]$Gene,
    'ORF' = autophagy_preds[['ds_parms_comb']]$ORF,
    'Gene2' = autophagy_preds[['ds_parms_comb']]$Gene2,
    'ORF2' = autophagy_preds[['ds_parms_comb']]$ORF2,
    'primary_identifier' = autophagy_preds[['ds_parms_comb']]$primary_identifier,
    'sgd_description' = autophagy_preds[['ds_parms_comb']]$sgd_description
  ) |>
  dplyr::distinct() |>
  dplyr::left_join(
    dplyr::select(autophagy_preds$profile,
                  primary_identifier, Profile),
    by = "primary_identifier") |>
  dplyr::rename(response_profile = Profile) |>
  dplyr::mutate(response_profile = dplyr::if_else(
    is.na(response_profile),
    "Not available",
    stringr::str_replace(response_profile, "[1-6]: ", "")
  )) |>
  dplyr::distinct()

####---- 2. Identify gene/ORFs with missing/limited data----#####
missing_multi <-
  autophagy_preds$ds_parms_comb |>
  dplyr::select(primary_identifier, Parameter) |>
  dplyr::distinct() |>
  dplyr::group_by(primary_identifier) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::filter(n < 18) |>
  dplyr::select(primary_identifier)

gene_info_kinetic_multi <-
  gene_info_kinetic |>
  dplyr::filter(!(
    primary_identifier %in% missing_multi$primary_identifier))

autophagy_preds[['ds_parms_comb']] <-
  autophagy_preds[['ds_parms_comb']] |>
  dplyr::select(
    -dplyr::any_of(
      c("sgd_id","entrezgene",
        "genename","alias_type",
        "sgd_description")))


## Read Bayes factor (BF) data - model specific (22 vs. 30)
excel_file_bfactors <- file.path(
   here::here(),"data-raw",
   "Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx")

bfactors_data <- list()
bfactors_data[['starv']] <- readxl::read_excel(
  excel_file_bfactors, sheet = 3)
bfactors_data[['repl']] <- readxl::read_excel(
   excel_file_bfactors, sheet = 4)

####---- 3. Read average BF model data----#####
load(file.path(
  here::here(), "data-raw",
  "Bayes_factor_data_for_web_portal_averaged_20241206.RData"))

bfactors_data[['overall']] <- as.data.frame(df_BF_average |>
  dplyr::mutate(Plate_controls = dplyr::if_else(
    nchar(Plate_controls) == 0,
    as.character(NA),
    Plate_controls
  )) |>
  dplyr::mutate(Gene = dplyr::if_else(
    nchar(Gene) == 0,
    as.character(NA),
    Gene
  )) |>
  dplyr::mutate(Reference_sets = dplyr::if_else(
    nchar(Reference_sets) == 0,
    as.character(NA),
    Reference_sets
  )) |>
  dplyr::select(-c("Red","Green","GxR","qRed","qGreen",
                   "Red.late","Green.late","GxR.late",
                   "P1_adj","P1_binary_adj"))
)

bfactors_data[['temporal']] <- as.data.frame(df_BF_TW_average |>
  dplyr::mutate(Plate_controls = dplyr::if_else(
    nchar(Plate_controls) == 0,
    as.character(NA),
    Plate_controls
  )) |>
  dplyr::mutate(Gene = dplyr::if_else(
    nchar(Gene) == 0,
    as.character(NA),
    Gene
  )) |>
  dplyr::mutate(ORF = dplyr::if_else(
    nchar(ORF) == 0,
    as.character(NA),
    ORF
  ))
)

rownames(bfactors_data$temporal) <- NULL
rownames(bfactors_data$overall) <- NULL

for(elem in c('overall','repl','temporal')){
  bfactors_data[[elem]] <- clean_internal_gene_orfs(
    dataset = bfactors_data[[elem]])
  bfactors_data[[elem]] <- map_yeast_gene_identifiers(
    screen_ids = bfactors_data[[elem]])
  if(elem != 'repl'){
    bfactors_data[[elem]] <- bfactors_data[[elem]] |>
      dplyr::select(
        -dplyr::any_of(
          c("sgd_id","entrezgene",
            "genename","alias_type","sgd_description")))

  }
}

gene_info_bf <-
  data.frame(
    'orf_gene_id' = bfactors_data[['repl']]$primary_identifier,
    'genename' = bfactors_data[['repl']]$genename,
    'entrezgene' = bfactors_data[['repl']]$entrezgene,
    'screen_id' = bfactors_data[['repl']]$screen_id,
    'sgd_id' = bfactors_data[['repl']]$sgd_id,
    'Gene' = bfactors_data[['repl']]$Gene,
    'ORF' = bfactors_data[['repl']]$ORF,
    'Gene2' = bfactors_data[['repl']]$Gene2,
    'ORF2' = bfactors_data[['repl']]$ORF2,
    'primary_identifier' = bfactors_data[['repl']]$primary_identifier,
    'sgd_description' = bfactors_data[['repl']]$sgd_description
  ) |>
  dplyr::distinct()

bfactors_data[['repl']] <- bfactors_data[['repl']] |>
  dplyr::select(
    -dplyr::any_of(c("sgd_id","entrezgene",
                     "genename","alias_type","sgd_description")))




regr_df_all <- data.frame()
for(n in names(regression.results)){
  if(is.data.frame(regression.results[[n]])){
    regr_df <- as.data.frame(regression.results[[n]]) |>
      dplyr::mutate(
        id = n, id2 = stringr::str_replace_all(n," ","")) |>
      dplyr::select(
        c("id", "id2", dplyr::everything()))

    regr_df_all <- dplyr::bind_rows(regr_df_all, regr_df)
  }
}

####---- 4. Prepare mutant-specific data for plots----#####
gw_autoph_data <- list()
gw_autoph_data$ds_curvefits <- autophagy_preds[['ds_curvefits']]
gw_autoph_data$dnn_preds <- autophagy_preds[['dnn']]
gw_autoph_data$regr_df_all <- regr_df_all
gw_autoph_data$BF_overall <- bfactors_data[['overall']]
gw_autoph_data$BF_temporal <- bfactors_data[['temporal']]

gw_autoph_kinetic_plot_input <- list()
gw_autoph_competence_plot_input <- list()

## Loop through all gene identifiers and prepare input data for plots
## A) Kinetic response
i <- 1
for(id in unique(gene_info_kinetic$orf_gene_id)){
  gw_autoph_kinetic_plot_input[[id]] <-
    kinetic_response_per_ko(id = id, gw_autoph_data = gw_autoph_data)
  if(i %% 100 == 0){
    log4r::info(log4r_logger, msg = paste0(
      "Preparing kinetic plot input - processed ",
      i, " gene identifiers.. "))
  }
  i <- i + 1
}

## B) Autophagy competence (BF)
i <- 1
for(id in unique(gene_info_bf$orf_gene_id)){
  gw_autoph_competence_plot_input[[id]] <-
    autophagy_competence_per_ko(
      id = id, gw_autoph_data = gw_autoph_data)
  if(i %% 100 == 0){
    log4r::info(log4r_logger, msg = paste0(
      "Preparing autophagy competence plot input - processed ",
      i, " gene identifiers.. "))
  }
  i <- i + 1
}

####---- 5. Gather gene information for plots----#####
ortholog_links <- readRDS(
  file = file.path(
    here::here(), "data-raw", "yeast_ortholog_map.rds")) |>
  tidyr::separate(
    value, c("entrez_human","symbol",
             "orthology_support","orthology_match",
             "primary_match"),sep="\\|", remove=F) |>
  dplyr::mutate(
    entrez_link = dplyr::if_else(
      primary_match == "Yes",
      paste0(
        "<b><a href='https://www.ncbi.nlm.nih.gov/gene/",
        entrez_human,"' target='_blank'>",
        symbol,"</a></b>"),
      paste0(
        "<a href='https://www.ncbi.nlm.nih.gov/gene/",
        entrez_human,"' target='_blank'>",
        symbol,"</a>")
      )
  ) |>
  dplyr::select(
    -dplyr::any_of(c("value","entrez_human","symbol",
                     "orthology_support","orthology_match",
                     "primary_match"))) |>
  dplyr::group_by(entrezgene) |>
  dplyr::summarise(
    human_ortholog_links = paste0(entrez_link, collapse = " | "),
    .groups = "drop") |>
  dplyr::mutate(human_ortholog_links = paste0(
    human_ortholog_links,
    " [<i>provided by Alliance of Genome Resources - ",
    datestamp_genome_alliance, "</i>]")) |>
  dplyr::distinct() |>
  dplyr::mutate(entrezgene = as.numeric(entrezgene))

gene_info_kinetic <- gene_info_kinetic |>
  dplyr::left_join(ortholog_links,
                   by=c("entrezgene"="entrezgene")) |>
  dplyr::mutate(
    human_ortholog_links = dplyr::if_else(
      is.na(human_ortholog_links),
      paste0(
        "<i>No data on human orthologs available ",
        "[provided by Alliance of Genome Resources - ",
        datestamp_genome_alliance,"]</i>."),
      human_ortholog_links)
  ) |>
  dplyr::mutate(sgd_description = stringr::str_replace_all(
    sgd_description,"\\. Orthologous to .+\\.$","."))


gene_info_bf <- gene_info_bf |>
  dplyr::left_join(ortholog_links,
                   by=c("entrezgene"="entrezgene")) |>
  dplyr::mutate(
    human_ortholog_links = dplyr::if_else(
      is.na(human_ortholog_links),
      paste0(
        "<i>No data on human orthologs available ",
        "[provided by Alliance of Genome Resources - ",
        datestamp_genome_alliance,"]</i>."),
      human_ortholog_links)
  ) |>
  dplyr::mutate(sgd_description = stringr::str_replace_all(
    sgd_description,"\\. Orthologous to .+\\.$","."))

## change name of parameters
for(elem in c('ds_parms','ds_parms_ctrs')){
  if(!is.null(autophagy_preds[[elem]])){
    autophagy_preds[[elem]]$Parameter <-
      gsub("A_","Autophagy ",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      ifelse(!grepl("Param",autophagy_preds[[elem]]$Parameter),
             gsub("slope","Slope (tangent) ",autophagy_preds[[elem]]$Parameter),
             gsub("slope","Slope (sigmoid) ",autophagy_preds[[elem]]$Parameter))
    autophagy_preds[[elem]]$Parameter <-
      gsub("Param","",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      gsub("1","-N",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      gsub("2","+N ",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      gsub("_"," ",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      gsub("starvation","-N",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      gsub("replenishment","+N",autophagy_preds[[elem]]$Parameter)
    autophagy_preds[[elem]]$Parameter <-
      stringr::str_trim(autophagy_preds[[elem]]$Parameter)

    autophagy_preds[[elem]] <- autophagy_preds[[elem]] |>
      dplyr::semi_join(
        gene_info_kinetic_multi,
        by = c("primary_identifier"))
  }
}

autophagy_preds[['ds_curvefits']] <-
  autophagy_preds[['ds_curvefits']] |>
  dplyr::semi_join(
    gene_info_kinetic_multi,
    by = c("primary_identifier"))


####---- 6. Save processed data ----####

## - These data files will be loaded initially in the
##   shiny app to improve response time of plots etc.
saveRDS(gw_autoph_competence_plot_input, file.path(
  here::here(), "data","processed",
  "gw_autoph_competence_plot_input.rds"))
saveRDS(gw_autoph_kinetic_plot_input, file.path(
  here::here(), "data","processed",
  "gw_autoph_kinetic_plot_input.rds"))
fst::write_fst(
  autophagy_preds$dnn, file.path(
    here::here(), "data","processed",
    "dnn_preds.fst"))
fst::write_fst(
  autophagy_preds$ds_curvefits, file.path(
    here::here(), "data","processed",
    "ds_curvefits.fst"))
fst::write_fst(
  autophagy_preds$ds_parms, file.path(
    here::here(), "data","processed",
    "ds_parms.fst"))
fst::write_fst(
  autophagy_preds$ds_parms_ctrs, file.path(
    here::here(), "data","processed",
    "ds_parms_ctrs.fst"))
fst::write_fst(
  gw_autoph_data$regr_df_all, file.path(
    here::here(), "data","processed",
    "regr_df_all.fst"))
fst::write_fst(
  gene_info_kinetic, file.path(
    here::here(), "data","processed",
    "gene_info_kinetic.fst"))
fst::write_fst(
  gene_info_kinetic_multi, file.path(
    here::here(), "data","processed",
    "gene_info_kinetic_multi.fst"))


fst::write_fst(
  gene_info_bf, file.path(
    here::here(), "data","processed","gene_info_bf.fst"))
fst::write_fst(
  gw_autoph_data$BF_overall, file.path(
    here::here(), "data","processed","BF_overall.fst"))
fst::write_fst(
  gw_autoph_data$BF_temporal, file.path(
    here::here(), "data","processed","BF_temporal.fst"))


