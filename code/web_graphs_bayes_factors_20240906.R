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
path <- "data-raw"

df_BF_overall <- read_excel(paste0(path,"/Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx"),
                            sheet=2)
df_BF_starv <- read_excel(paste0(path,"/Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx"),
                          sheet=3)
df_BF_repl <- read_excel(paste0(path,"/Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx"),
                         sheet=4)
df_BF_temporal <- read_excel(paste0(path,"/Chica_et_al_Data_File_6_Autophagy_Bayes_factors.xlsx"),
                             sheet=5)
df_DNN_preds <- read_excel(paste0(path,"/Chica_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx"),
                           sheet=2)

#ID list
IDs_unique <- gsub("NA ", "", unique(paste(df_BF_repl$Gene,df_BF_repl$ORF)))

# Response kinetics----
#Select ID
id=IDs_unique[grep("RTG2",IDs_unique)]
#Select model
Model = c("22","30")[2]

BF_response <- df_BF_temporal %>% subset(gsub("NA ", "", paste(Gene, ORF)) == id)
#If mutant in rec plate, only evaluate rec mutants
if(any(grepl("Rec",BF_response$Plate))){
  BF_response <- BF_response[which(grepl("Rec",BF_response$Plate)),]
}
#Find representative response
BF_response <- BF_response %>%
  group_by(TimeR) %>%
  mutate(d=abs(log_BFt_WT.ATG1_30-median(log_BFt_WT.ATG1_30))+
           abs(log_BFt_WT.VAM6_30-median(log_BFt_WT.VAM6_30))+
           abs(log_BFt_VAM6.ATG1_30-median(log_BFt_VAM6.ATG1_30))+
           abs(log_BFt_WT.ATG1_22-median(log_BFt_WT.ATG1_22))+
           abs(log_BFt_WT.VAM6_22-median(log_BFt_WT.VAM6_22))+
           abs(log_BFt_VAM6.ATG1_22-median(log_BFt_VAM6.ATG1_22))) %>%
  group_by(Plate, Position) %>%
  mutate(d=mean(d))
BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
BF_response_ctr <- df_BF_temporal[df_BF_temporal$Plate == BF_response$Plate[1],] %>%
  group_by(TimeR) %>% summarise_if(is.numeric, mean)


Time.shift <- unique(BF_response$TimeR) + 1
names(Time.shift) <- unique(BF_response$TimeR)
Time.shift[which.max(Time.shift)] <- 0

BF_response$shift <- Time.shift[as.character(BF_response$TimeR)]
BF_response <- data.frame(BF_response)
BF_response <- BF_response %>%
  mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
         log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
         log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
         log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
         log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
         log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)
BF_response_ctr$shift <- Time.shift[as.character(BF_response_ctr$TimeR)]
BF_response_ctr <- data.frame(BF_response_ctr)
BF_response_ctr <- BF_response_ctr %>%
  mutate(log_BFt_WT.ATG1_30.shift = (log_BFt_WT.ATG1_30[match(shift, TimeR)]-log_BFt_WT.ATG1_30)/3+log_BFt_WT.ATG1_30,
         log_BFt_WT.VAM6_30.shift =  (log_BFt_WT.VAM6_30[match(shift, TimeR)]-log_BFt_WT.VAM6_30)/3 + log_BFt_WT.VAM6_30,
         log_BFt_VAM6.ATG1_30.shift =  (log_BFt_VAM6.ATG1_30[match(shift, TimeR)]-log_BFt_VAM6.ATG1_30)/3 + log_BFt_VAM6.ATG1_30,
         log_BFt_WT.ATG1_22.shift = (log_BFt_WT.ATG1_22[match(shift, TimeR)]-log_BFt_WT.ATG1_22)/3+log_BFt_WT.ATG1_22,
         log_BFt_WT.VAM6_22.shift =  (log_BFt_WT.VAM6_22[match(shift, TimeR)]-log_BFt_WT.VAM6_22)/3 + log_BFt_WT.VAM6_22,
         log_BFt_VAM6.ATG1_22.shift =  (log_BFt_VAM6.ATG1_22[match(shift, TimeR)]-log_BFt_VAM6.ATG1_22)/3 + log_BFt_VAM6.ATG1_22)

Library <- df_DNN_preds$Type[gsub("NA ", "", paste(df_DNN_preds$Gene, df_DNN_preds$ORF)) == id][1]

if(Model == "30"){
  ggplot(BF_response,
         aes(log_BFt_VAM6.ATG1_30, log_BFt_WT.VAM6_30, col=TimeR))+
    geom_vline(xintercept = 0, lty=1, col="black", size=0.1)+
    geom_hline(yintercept = 0, lty=1, col="black", size=0.1)+
    scale_color_scico(palette = 'lisbon')+
    scale_size(range = c(0, 1.5))+
    #geom_path(data=BF_response_ctr, aes(),col="gray",lty=2, alpha = 1) +
    geom_point(data=BF_response_ctr,aes(size=log_BFt_WT.ATG1_30, pch="Control"),size=0.2, alpha = 1) +
    geom_path(aes(),col="black",lty=2, alpha = 1) +
    geom_segment(aes(xend = log_BFt_VAM6.ATG1_30.shift, yend = log_BFt_WT.VAM6_30.shift, size=log_BFt_WT.ATG1_30),
                 arrow = arrow(angle = 15, length = unit(0.07, "inches"), type = "closed"), alpha = 1, size=1) +
    labs(x="log BF (VAM6:ATG1) 30", y="log BF (WT:VAM6) 30", col="Time", size="BF (WT:ATG1)", title=paste(id, Library), shape="")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12),
          strip.background = element_blank(),
          legend.position = "top")
}else{
  ggplot(BF_response,
         aes(log_BFt_VAM6.ATG1_22, log_BFt_WT.VAM6_22, col=TimeR))+
    geom_vline(xintercept = 0, lty=1, col="black", size=0.1)+
    geom_hline(yintercept = 0, lty=1, col="black", size=0.1)+
    scale_color_scico(palette = 'lisbon')+
    scale_size(range = c(0, 1.5))+
    #geom_path(data=BF_response_ctr,aes(),col="gray",lty=2, alpha = 1) +
    geom_point(data=BF_response_ctr,aes(size=log_BFt_WT.ATG1_22, pch="Control"),size=0.2, alpha = 1) +
    geom_path(aes(),col="black",lty=2, alpha = 1) +
    geom_segment(aes(xend = log_BFt_VAM6.ATG1_22.shift, yend = log_BFt_WT.VAM6_22.shift, size=log_BFt_WT.ATG1_22),
                 arrow = arrow(angle = 15, length = unit(0.07, "inches"), type = "closed"), alpha = 1, size=1) +
    labs(x="log BF (VAM6:ATG1) 22", y="log BF (WT:VAM6) 22", col="Time", size="BF (WT:ATG1)", shape="")+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          plot.title = element_text(size=12),
          strip.background = element_blank(),
          legend.position = "top")
}




# BF overall plots----
#Select one or several IDs for text repel
id=IDs_unique[grep("RTG",IDs_unique)]

#Select x and y
variables <- colnames(df_BF_overall)[5:16]
X <- variables[3]
Y <- variables[2]


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
    mutate(d=abs(log_BFt_WT.ATG1_30-median(log_BFt_WT.ATG1_30))+
             abs(log_BFt_WT.VAM6_30-median(log_BFt_WT.VAM6_30))+
             abs(log_BFt_VAM6.ATG1_30-median(log_BFt_VAM6.ATG1_30))+
             abs(log_BFt_WT.ATG1_22-median(log_BFt_WT.ATG1_22))+
             abs(log_BFt_WT.VAM6_22-median(log_BFt_WT.VAM6_22))+
             abs(log_BFt_VAM6.ATG1_22-median(log_BFt_VAM6.ATG1_22))) %>%
    group_by(Plate, Position) %>%
    mutate(d=mean(d))
  BF_response <- BF_response[which(BF_response$d==min(BF_response$d)),]
  Positions <- c(Positions,paste(BF_response$Plate, BF_response$Position)[1])
}

mat <- df_BF_overall %>% as.data.frame()
mat_select <- df_BF_overall %>% subset(paste(Plate, Position) %in% Positions) %>% as.data.frame()

mat$X <- mat[,X]
mat$Y <- mat[,Y]
mat_select$X <- mat_select[,X]
mat_select$Y <- mat_select[,Y]

mat$Type <- df_DNN_preds$Type[match(paste(mat$Gene, mat$ORF),paste(df_DNN_preds$Gene, df_DNN_preds$ORF))]

#Dot plots
lab_x <- ""
lab_y <- ""
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
  #labs(x=gsub("\\.",":",gsub("_"," ",X)), y=gsub("\\.",":",gsub("_"," ",Y)), color="Reference sets", linetype="Library")+
  labs(x=paste0(
    gsub("\\.",":",gsub("_"," ",X)),"\n",lab_x),
    y=paste0(lab_y,"\n",gsub("\\.",":",gsub("_"," ",Y))),
    color="Reference sets", linetype="Library")+
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
  scale_linetype_manual(values=c(2,1))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  guides(alpha = FALSE)





