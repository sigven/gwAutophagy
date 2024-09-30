library(readxl)
library(dplyr)
library(ggplot2)
library(ggsci)
library(matrixStats)
library(reshape2)
library(ggrepel)

#Read data
#path <- "C:/Users/aramn/OneDrive/Desktop/GW autophagy_scripts for web portal"
load(file.path(
  here::here(),"code",
  "code_data_ana", "Double sigmoidal curve fitting_v2.RData"))
data_excel <- file.path(here::here(), "code","code_data_ana",
                        "ChiCa_et_al_Data_File_S2_Autophagy_predictions_Statistics_Clustering.xlsx")

df_DNN_preds <- readxl::read_excel(
  data_excel, sheet=2)

df_DS_curvefits <- readxl::read_excel(
  data_excel, sheet=3)

df_DS_parms <- readxl::read_excel(
  data_excel, sheet=5)

df_DS_parms_ctrs <- readxl::read_excel(
  data_excel, sheet=6)

df_DS_parms_comb <- dplyr::bind_rows(
  df_DS_parms,df_DS_parms_ctrs)

#fst::write_fst(df_DNN_preds, path = "data/df_DNN_preds.fst")
#fst::write_fst(df_DS_curvefits, path = "data/df_DS_curvefits.fst")
#fst::write_fst(df_DS_parms, path = "data/df_DS_parms.fst")
#fst::write_fst(df_DS_parms_ctrs, path = "data/df_DS_parms_ctrs.fst")
#fst::write_fst(df_DS_parms_comb, path = "data/df_DS_parms_comb.fst")


curve_fitting_data <- list()
curve_fitting_data$DNN_preds <- df_DNN_preds



#ID list
IDs_unique <- gsub("NA ", "", unique(
  paste(df_DS_parms_comb$Gene,df_DS_parms_comb$ORF)))

# Response kinetics----
#Select ID
id=IDs_unique[grep("RTG3",IDs_unique)]

y_pred <- df_DS_curvefits |> subset(gsub("NA ", "", paste(Gene, ORF)) == id)
#If mutant in rec plate, only evaluate rec mutants
if(any(grepl("Rec",y_pred$Plate))){
  y_pred <- y_pred[which(grepl("Rec",y_pred$Plate)),]
}
#Find representative response
y_pred <- y_pred |>
  dplyr::group_by(Time) |>
  dplyr::mutate(d=abs(P1_30_fit-median(P1_30_fit))) |>
  dplyr::group_by(Plate, Position) |>
  dplyr::mutate(d=mean(d))
y_pred <- y_pred[which(y_pred$d==min(y_pred$d)),]
y_raw <- df_DNN_preds |>
  subset(paste(Plate, Position) == paste(y_pred$Plate, y_pred$Position)[1])

Library <- ifelse(y_raw$Type== "KO" | !(is.na(y_raw$Plate_controls)),"KO","DAmP")[1]
Plate_id <- unique(y_raw$Plate)
y_pred.ctr <- df_DS_curvefits[which(is.na(df_DS_curvefits$Plate_controls)),] |>
  #subset(Type == Library) |>
  subset(Plate == Plate_id) |>
  dplyr::group_by(Time) |>
  dplyr::summarise(P1_30_fit=median(P1_30_fit))

doubleSigmoidalModel <- regression.results[[which(grepl(gsub(" ", "", paste(y_pred$Plate, y_pred$Position)[1]),gsub(" ","",names(regression.results)), ignore.case = T))]]


ggplot()+
  geom_vline(xintercept=0, lty=2)+
  geom_vline(xintercept=12, lty=2)+
  geom_hline(yintercept=0, lty=1, size=0.2)+
  geom_hline(yintercept=doubleSigmoidalModel[,"maximum_y"], lty=1, size=0.2)+
  geom_hline(yintercept=doubleSigmoidalModel[,"finalAsymptoteIntensity"], lty=1, size=0.2)+
  geom_segment(aes(x=doubleSigmoidalModel[,"startPoint_x"], y=0,
                   xend=doubleSigmoidalModel[,"reachMaximum_x"],
                   yend=doubleSigmoidalModel[,"reachMaximum_y"] ), lty=1, size=0.2)+
  geom_segment(aes(x=doubleSigmoidalModel[,"startDeclinePoint_x"],
                   y=doubleSigmoidalModel[,"startDeclinePoint_y"],
                   xend=doubleSigmoidalModel[,"endDeclinePoint_x"],
                   yend=doubleSigmoidalModel[,"endDeclinePoint_y"] ), lty=1, size=0.2)+
  geom_segment(aes(x=0, y=doubleSigmoidalModel[,"midPoint1_y"],
                   xend=doubleSigmoidalModel[,"midPoint1_x"],
                   yend=doubleSigmoidalModel[,"midPoint1_y"] ),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_segment(aes(x=12, y=doubleSigmoidalModel[,"midPoint2_y"],
                   xend=doubleSigmoidalModel[,"midPoint2_x"],
                   yend=doubleSigmoidalModel[,"midPoint2_y"] ),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_segment(aes(x=0, y=0,
                   xend=doubleSigmoidalModel[,"startPoint_x"],
                   yend=0),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_segment(aes(x=12, y=doubleSigmoidalModel[,"maximum_y"],
                   xend=doubleSigmoidalModel[,"startDeclinePoint_x"],
                   yend=doubleSigmoidalModel[,"maximum_y"] ),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_segment(aes(x=0, y=doubleSigmoidalModel[,"reachMaximum_y"],
                   xend=doubleSigmoidalModel[,"reachMaximum_x"],
                   yend=doubleSigmoidalModel[,"reachMaximum_y"]),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_segment(aes(x=12, y=doubleSigmoidalModel[,"endDeclinePoint_y"],
                   xend=doubleSigmoidalModel[,"endDeclinePoint_x"],
                   yend=doubleSigmoidalModel[,"endDeclinePoint_y"] ),
               lineend="round",arrow=arrow(length=unit(0.1, "inches"), ends = "both"),alpha=1, size=0.2)+
  geom_line(aes(x= y_pred.ctr$Time, y=y_pred.ctr$P1_30_fit, col="Control"))+
  geom_ribbon(aes(x=y_pred.ctr$Time,
                  y=colMeans(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                  ymin=colMins(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit)),
                  ymax=colMaxs(rbind(y_pred.ctr$P1_30_fit,y_pred$P1_30_fit))),fill="yellow", alpha=0.1)+
  geom_ribbon(aes(x=c(doubleSigmoidalModel[,"startPoint_x"],
                      doubleSigmoidalModel[,"reachMaximum_x"]),
                  y=c(0, doubleSigmoidalModel[,"reachMaximum_y"]),
                  ymin=c(0,0),
                  ymax=c(doubleSigmoidalModel[,"reachMaximum_y"],
                         doubleSigmoidalModel[,"reachMaximum_y"])),fill="black", alpha=0.1)+
  geom_ribbon(aes(x=c(doubleSigmoidalModel[,"startDeclinePoint_x"],
                      doubleSigmoidalModel[,"endDeclinePoint_x"]),
                  y=c(doubleSigmoidalModel[,"startDeclinePoint_y"],
                      doubleSigmoidalModel[,"endDeclinePoint_y"]),
                  ymin=c(doubleSigmoidalModel[,"endDeclinePoint_y"],doubleSigmoidalModel[,"endDeclinePoint_y"]),
                  ymax=c(doubleSigmoidalModel[,"startDeclinePoint_y"],
                         doubleSigmoidalModel[,"startDeclinePoint_y"])),fill="black", alpha=0.1)+
  geom_line(aes(x= y_pred$Time, y=y_pred$P1_30_fit, col="Model"))+
  geom_point(aes(x=ifelse(grepl("Starv",y_raw$Phase),y_raw$TimeR,y_raw$TimeR+12),y_raw$P1_30, col="Data"))+
  ylim(0,100)+
  xlim(-1,20)+
  labs(x="Time (hours)", y="Autophagy (%)", title=paste(id, Library))+
  scale_color_manual(values=pal_jco("default", alpha = 1)(7)[c(7,1,4)])+
  theme_bw()


# Kinetic parameters----
#Select one or several IDs for text repel
id=IDs_unique[grep("RTG",IDs_unique)]

#Select x and y
X <- unique(df_DS_parms$Parameter)[12]
Y <- unique(df_DS_parms$Parameter)[11]

#Select value or perturbation
Value = c("Value", "Perturbation")[1]

#Filter
Filter=TRUE

Positions <- c()
for(i in id){
  y_pred <- df_DS_curvefits |> subset(gsub("NA ", "", paste(Gene, ORF)) == i)
  #If mutant in rec plate, only evaluate rec mutants
  if(any(grepl("Rec",y_pred$Plate))){
    y_pred <- y_pred[which(grepl("Rec",y_pred$Plate)),]
  }
  #Find representative response
  y_pred <- y_pred |>
    dplyr::group_by(Time) |>
    dplyr::mutate(d=abs(P1_30_fit-median(P1_30_fit))) |>
    dplyr::group_by(Plate, Position) |>
    dplyr::mutate(d=mean(d))
  y_pred <- y_pred[which(y_pred$d==min(y_pred$d)),]
  Positions <- c(Positions,paste(y_pred$Plate, y_pred$Position)[1])
}

mat <- dcast(df_DS_parms[which(df_DS_parms$Parameter %in% c(X,Y)),],
             Plate+Position+ORF+Gene+ Reference_sets~Parameter,
             value.var = Value)
mat_select <- dcast(df_DS_parms_comb[which(df_DS_parms_comb$Parameter %in% c(X,Y)),] |> subset(paste(Plate, Position) %in% Positions),
                    Plate+Position+ORF+Gene + Reference_sets+Plate_controls~Parameter,
                    value.var = Value)
mat$X <- mat[,X]
mat$Y <- mat[,Y]
mat_select$X <- mat_select[,X]
mat_select$Y <- mat_select[,Y]

#Filter, replace manual control of scale?
if(Filter){
  x_lower_bound <- quantile(mat$X, 1-0.99, na.rm = T) - 5*IQR(mat$X, na.rm = T)
  x_upper_bound <- quantile(mat$X, 0.99, na.rm = T) + 5*IQR(mat$X, na.rm = T)
  y_lower_bound <- quantile(mat$Y, 1-0.99, na.rm = T) - 5*IQR(mat$Y, na.rm = T)
  y_upper_bound <- quantile(mat$Y, 0.99, na.rm = T) + 5*IQR(mat$Y, na.rm = T)
}else{
  x_lower_bound <- min(mat$X)
  x_upper_bound <- max(mat$X)
  y_lower_bound <- min(mat$Y)
  y_upper_bound <- max(mat$Y)
}

ggplot(mat, aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16, alpha=0.8)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets)),],aes(col=Reference_sets), size=1)+
  geom_point(data=mat_select,aes(), pch=21)+
  geom_text_repel(data=mat_select,aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  xlim(x_lower_bound,x_upper_bound)+
  ylim(y_lower_bound,y_upper_bound)+
  labs(x=X, y=Y)+
  scale_fill_jama()+
  scale_color_d3()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggplot(mat, aes(X, Y))+
  geom_point(col="lightgray", size=0.7, pch=16)+
  geom_point(data=mat[which(!is.na(mat$Reference_sets)),],aes(col=Reference_sets), size=1)+
  geom_text_repel(data=mat[which(!is.na(mat$Reference_sets)),],aes(label=Gene),
                  force=2, size=2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.size = 0.3)+
  xlim(x_lower_bound,x_upper_bound)+
  ylim(y_lower_bound,y_upper_bound)+
  labs(x=X, y=Y)+
  scale_fill_jama()+
  scale_color_d3()+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
