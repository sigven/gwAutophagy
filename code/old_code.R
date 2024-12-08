autophagy_competence_per_ko_old <- function(id, gw_autoph_data = NULL){

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
    dplyr::mutate(d=abs(log_BFt_WT.ATG1_30-median(log_BFt_WT.ATG1_30))+
                    abs(log_BFt_WT.VAM6_30-median(log_BFt_WT.VAM6_30))+
                    abs(log_BFt_VAM6.ATG1_30-median(log_BFt_VAM6.ATG1_30))+
                    abs(log_BFt_WT.ATG1_22-median(log_BFt_WT.ATG1_22))+
                    abs(log_BFt_WT.VAM6_22-median(log_BFt_WT.VAM6_22))+
                    abs(log_BFt_VAM6.ATG1_22-median(log_BFt_VAM6.ATG1_22))) |>
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
    dplyr::mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
                  log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
                  log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
                  log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
                  log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
                  log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)
  BF_response_ctr$shift <- Time.shift[as.character(BF_response_ctr$TimeR)]
  BF_response_ctr <- data.frame(BF_response_ctr)
  BF_response_ctr <- BF_response_ctr |>
    dplyr::mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
                  log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
                  log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
                  log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
                  log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
                  log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)

  Library <- gw_autoph_data$dnn_preds[gw_autoph_data$dnn_preds$primary_identifier == id,]$Type[1]

  return(list(
    BF_response = BF_response,
    BF_response_ctr = BF_response_ctr,
    Library = Library,
    id = id
  ))

}

