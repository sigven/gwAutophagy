library(readxl)
library(dplyr)
library(ggplot2)
library(ggsci)
library(matrixStats)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(scico)

#Read data
data_path <- file.path(here::here(),"data-raw")
load(file.path(data_path,"Bayes_factor_data_for_web_portal_averaged_20241206.RData"))
#load("C:/Users/aramn/OneDrive/Desktop/GW autophagy revision Bayes factors/Bayes factor data for web portal_averaged_241206.RData")

df_BF_overall <- df_BF_average |>
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

df_BF_temporal <- df_BF_TW_average |>
  dplyr::mutate(Plate_controls = dplyr::if_else(
    nchar(Plate_controls) == 0,
    as.character(NA),
    Plate_controls
  )) |>
  dplyr::mutate(Gene = dplyr::if_else(
    nchar(Gene) == 0,
    as.character(NA),
    Gene
  ))
rm(df_BF_average, df_BF_TW_average)

df_DNN_preds <- readxl::read_excel(
  path = file.path(data_path,"Chica_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx"),
                           sheet=2)
  #paste0(path,"/Chica_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx"),
  #                         sheet=2)


#ID list
library(org.Sc.sgd.db)
map <- as.data.frame(org.Sc.sgdCOMMON2ORF)
df_BF_overall$Gene <- map$gene_name[match(df_BF_overall$ORF, map$systematic_name)]
df_BF_temporal$Gene <- map$gene_name[match(df_BF_temporal$ORF, map$systematic_name)]
df_DNN_preds$Gene <- map$gene_name[match(df_DNN_preds$ORF, map$systematic_name)]

IDs_unique <- gsub("NA ", "", unique(paste(df_BF_overall$Gene,df_BF_overall$ORF)))


# Response kinetics----
#Select ID
id=IDs_unique[grep("TOR1",IDs_unique)]

BF_response <- df_BF_temporal %>% subset(gsub("NA ", "", paste(Gene, ORF)) == id)
#If mutant in rec plate, only evaluate rec mutants
if(any(grepl("Rec",BF_response$Plate))){
  BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
}
#Find representative response
BF_response <- BF_response %>%
  group_by(TimeR) %>%
  mutate(d=abs(log_BFt_WT.ATG1-median(log_BFt_WT.ATG1))+
           abs(log_BFt_WT.VAM6-median(log_BFt_WT.VAM6))+
           abs(log_BFt_VAM6.ATG1-median(log_BFt_VAM6.ATG1))) %>%
  group_by(Plate, Position) %>%
  mutate(d=mean(d),
         NCells = mean(NCells))
BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
BF_response <- BF_response[which(BF_response$NCells==max(BF_response$NCells)),]


BF_response_ctr <- df_BF_temporal[df_BF_temporal$Plate == BF_response$Plate[1],] %>%
  group_by(TimeR) %>% summarise_if(is.numeric, mean)


Time.shift <- unique(BF_response$TimeR) + 1
names(Time.shift) <- unique(BF_response$TimeR)
Time.shift[which.max(Time.shift)] <- 0

BF_response$shift <- Time.shift[as.character(BF_response$TimeR)]
BF_response <- data.frame(BF_response)
BF_response <- BF_response %>%
  mutate(log_BFt_WT.ATG1.shift = (log_BFt_WT.ATG1[match(shift, TimeR)]-log_BFt_WT.ATG1)/3+log_BFt_WT.ATG1,
         log_BFt_WT.VAM6.shift =  (log_BFt_WT.VAM6[match(shift, TimeR)]-log_BFt_WT.VAM6)/3 + log_BFt_WT.VAM6,
         log_BFt_VAM6.ATG1.shift =  (log_BFt_VAM6.ATG1[match(shift, TimeR)]-log_BFt_VAM6.ATG1)/3 + log_BFt_VAM6.ATG1)
BF_response_ctr$shift <- Time.shift[as.character(BF_response_ctr$TimeR)]
BF_response_ctr <- data.frame(BF_response_ctr)
BF_response_ctr <- BF_response_ctr %>%
  mutate(log_BFt_WT.ATG1.shift = (log_BFt_WT.ATG1[match(shift, TimeR)]-log_BFt_WT.ATG1)/3+log_BFt_WT.ATG1,
         log_BFt_WT.VAM6.shift =  (log_BFt_WT.VAM6[match(shift, TimeR)]-log_BFt_WT.VAM6)/3 + log_BFt_WT.VAM6,
         log_BFt_VAM6.ATG1.shift =  (log_BFt_VAM6.ATG1[match(shift, TimeR)]-log_BFt_VAM6.ATG1)/3 + log_BFt_VAM6.ATG1)

Library <- df_DNN_preds$Type[gsub("NA ", "", paste(df_DNN_preds$Gene, df_DNN_preds$ORF)) == id][1]

ggplot(BF_response,
       aes(log_BFt_VAM6.ATG1, log_BFt_WT.VAM6, col=TimeR))+
  geom_vline(xintercept = 0, lty=1, col="black", size=0.1)+
  geom_hline(yintercept = 0, lty=1, col="black", size=0.1)+
  scale_color_scico(palette = 'lisbon')+
  scale_size(range = c(0, 1.5))+
  #geom_path(data=BF_response_ctr, aes(),col="gray",lty=2, alpha = 1) +
  geom_point(data=BF_response_ctr,aes(size=log_BFt_WT.ATG1, pch="Control"),size=0.2, alpha = 1) +
  geom_path(aes(),col="black",lty=2, alpha = 1) +
  geom_segment(aes(xend = log_BFt_VAM6.ATG1.shift, yend = log_BFt_WT.VAM6.shift, size=log_BFt_WT.ATG1),
               arrow = arrow(angle = 15, length = unit(0.07, "inches"), type = "closed"), alpha = 1, size=1) +
  labs(x="log BF (VAM6:ATG1) 30\nAutophagosome formation", y="Autophagosome clearance\nlog BF (WT:VAM6) 30", col="Time", size="BF (WT:ATG1)", title=paste(id, Library), shape="")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(size=12),
        strip.background = element_blank(),
        legend.position = "top")


# BF overall plots----
df <- df_BF_overall
df$Plate_controls[df$Plate_controls != "+"] <- NA
df$Reference_sets[df$Reference_sets == ""] <- NA

#Select one or several IDs for text repel
#id=IDs_unique[grep("tor1|tor2|tco89|kog1|lst8",IDs_unique, ignore.case = T)]

id=IDs_unique[grep("tco89",IDs_unique, ignore.case = T)]

#Select x and y
variables <- colnames(df)[5:7]
X <- variables[3]
Y <- variables[2]

Library_adjust <- TRUE

Positions <- c()
for(i in id){
  BF_response <- df_BF_temporal %>% subset(gsub("NA ", "", paste(Gene, ORF)) == i)
  #If mutant in rec plate, only evaluate rec mutants
  if(any(grepl("Rec",BF_response$Plate))){
    BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
  }
  #Find representative response
  BF_response <- BF_response %>%
    group_by(TimeR) %>%
    mutate(d=abs(log_BFt_WT.ATG1-median(log_BFt_WT.ATG1))+
             abs(log_BFt_WT.VAM6-median(log_BFt_WT.VAM6))+
             abs(log_BFt_VAM6.ATG1-median(log_BFt_VAM6.ATG1))) %>%
    group_by(Plate, Position) %>%
    mutate(d=mean(d),
           NCells = mean(NCells))
  BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
  BF_response <- BF_response[which(BF_response$NCells==max(BF_response$NCells)),]
  Positions <- c(Positions,paste(BF_response$Plate, BF_response$Position)[1])
}

mat <- df %>% as.data.frame()
mat_select <- df %>% subset(paste(Plate, Position) %in% Positions) %>% as.data.frame()

mat$X <- mat[,X]
mat$Y <- mat[,Y]
mat_select$X <- mat_select[,X]
mat_select$Y <- mat_select[,Y]

mat$Type <- df_DNN_preds$Type[match(paste(mat$Gene, mat$ORF),paste(df_DNN_preds$Gene, df_DNN_preds$ORF))]

if(grepl("VAM6.ATG1", X, fixed = T)){
  lab_x <- "Autophagosome formation"
}else if(grepl("WT.ATG1", X, fixed = T)){
  lab_x <- "Overall autophagy"
}else{
  lab_x <- "Autophagosome clearance"
}
if(grepl("VAM6.ATG1", Y, fixed = T)){
  lab_y <- "Autophagosome formation"
}else if(grepl("WT.ATG1", Y, fixed = T)){
  lab_y <- "Overall autophagy"
}else{
  lab_y <- "Autophagosome clearance"
}

if(Library_adjust){
  mat_lib <- mat[is.na(mat$Plate_controls),] %>%
    mutate(mean_x = mean(X, na.rm=T), mean_y=mean(Y, na.rm=T),
           sd_x = sd(X, na.rm=T), sd_y=sd(Y, na.rm=T)) %>%
    group_by(Type) %>%
    summarise(mean_type_x = mean(X, na.rm=T), mean_type_y=mean(Y, na.rm=T),
              sd_type_x = sd(X, na.rm=T), sd_type_y=sd(Y, na.rm=T),
              mean_x = mean(mean_x, na.rm=T), mean_y=mean(mean_y, na.rm=T),
              sd_x = mean(sd_x, na.rm=T), sd_y=mean(sd_y, na.rm=T))

  mat$Type[which(mat$Plate_controls == "+")] <- "KO"
  mat$Type[which(mat$Plate_controls=="+" & is.na(mat$ORF))] <- "WT"

  mat_wt <- mat[which(mat$Plate_controls=="+" & is.na(mat$ORF)),] %>%
    mutate(mean_x = mean(X), mean_y=mean(Y),
           sd_x = sd(X), sd_y=sd(Y)) %>%
    group_by(Type) %>%
    summarise(mean_type_x = mean(X, na.rm=T), mean_type_y=mean(Y, na.rm=T),
              sd_type_x = sd(X, na.rm=T), sd_type_y=sd(Y, na.rm=T),
              mean_x = mean(mean_x, na.rm=T), mean_y=mean(mean_y, na.rm=T),
              sd_x = mean(sd_x, na.rm=T), sd_y=mean(sd_y, na.rm=T))
  mat_wt[,6:9] <- mat_lib[1,6:9]
  mat_wt[1,4:5] <- mat_wt[1,8:9]
  mat_lib <- rbind(mat_lib, mat_wt)

  mat$X <-  mat$X + (mat_lib$mean_x-mat_lib$mean_type_x)[match(mat$Type,mat_lib$Type)]
  mat$Y <- mat$Y + (mat_lib$mean_y-mat_lib$mean_type_y)[match(mat$Type,mat_lib$Type)]
  mat$X <- mat$X*((mat_lib$sd_x/mat_lib$sd_type_x)[match(mat$Type,mat_lib$Type)])
  mat$Y <- mat$Y*((mat_lib$sd_y/mat_lib$sd_type_y)[match(mat$Type,mat_lib$Type)])
}
mat$Type <- factor(mat$Type, levels=c("KO","DAmP"))

ggplot(mat[is.na(mat$Plate_controls),], aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8)+
  stat_density_2d(data=mat[is.na(mat$Plate_controls),],
                  aes(lty=Type), col="grey30", alpha=1)+
  stat_density_2d(data=mat[which(mat$Plate_controls=="+"),] %>% mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
                  aes(fill=Gene, group=Gene, alpha = ..level..), geom = "polygon", col=NA)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets) & !grepl("ORF",mat$Reference_sets)),],
             aes(col=Reference_sets), size=1)+
  geom_point(data=mat_select,aes(), pch=21)+
  geom_text_repel(data=mat_select,aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  labs(x=gsub("\\.",":",gsub("_"," ",X)), y=gsub("\\.",":",gsub("_"," ",Y)), color="Reference sets", linetype="Library")+
  scale_fill_jama()+
  scale_color_d3()+
 # scale_linetype_manual(values=c(2,1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(alpha = FALSE)

ggplot(mat[is.na(mat$Plate_controls),], aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8)+
  # stat_density_2d(data=mat[is.na(mat$Plate_controls),],
  #                 aes(lty=Type), col="grey30", alpha=1)+
  stat_density_2d(data=mat[which(mat$Plate_controls=="+"),] %>% mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
                  aes(fill=Gene, group=Gene, alpha = ..level..), geom = "polygon", col=NA)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets) & !grepl("ORF",mat$Reference_sets)),],
             aes(col=Reference_sets), size=1)+
  geom_point(data=mat_select,aes(), pch=21)+
  geom_text_repel(data=mat_select,aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  labs(x=paste0(gsub("\\.",":",gsub("_"," ",X)),"\n",lab_x), y=paste0(lab_y,"\n",gsub("\\.",":",gsub("_"," ",Y))), color="Reference sets", linetype="Library")+
  scale_fill_jama()+
  scale_color_d3()+
  scale_linetype_manual(values=c(2,1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(alpha = FALSE)

ggplot(mat[is.na(mat$Plate_controls),], aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8)+
  stat_density_2d(data=mat[is.na(mat$Plate_controls),],
                  aes(), col="grey30", alpha=1)+
  stat_density_2d(data=mat[which(mat$Plate_controls=="+"),] %>% mutate(Gene = ifelse(is.na(Gene), "WT", Gene)),
                  aes(fill=Gene, group=Gene, alpha = ..level..), geom = "polygon", col=NA)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets) & !grepl("ORF",mat$Reference_sets)),],
             aes(col=Reference_sets), size=1)+
  geom_point(data=mat_select,aes(), pch=21)+
  geom_text_repel(data=mat_select,aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  labs(x=gsub("\\.",":",gsub("_"," ",X)), y=gsub("\\.",":",gsub("_"," ",Y)), color="Reference sets", linetype="Library")+
  scale_fill_jama()+
  scale_color_d3()+
  scale_linetype_manual(values=c(2,1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(alpha = FALSE)

ggplot(mat[is.na(mat$Plate_controls),], aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8)+
  stat_density_2d(data=mat[is.na(mat$Plate_controls),],
                  aes(), col="grey30", alpha=1)+
  stat_density_2d(data=mat[!is.na(mat$Plate_controls) & !is.na(mat$Gene) & mat$Gene != "WT",],
                  aes(fill=Gene, group=Gene, alpha = ..level..), geom = "polygon", col=NA)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets) & !grepl("ORF",mat$Reference_sets)),],
             aes(col=Reference_sets), size=1)+
  geom_point(data=mat_select,aes(), pch=21)+
  geom_text_repel(data=mat_select,aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  labs(x=gsub("\\.",":",gsub("_"," ",X)), y=gsub("\\.",":",gsub("_"," ",Y)), color="Reference sets")+
  scale_fill_jama()+
  scale_color_d3()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(alpha = FALSE)






