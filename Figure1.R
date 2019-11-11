options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("plotly")
library(cowplot)


setwd("E:/DO/R-Analyses/Summer_2019/") 
load("P_L_Lipids.RData")

# SELECT TISSUE 
tis <- "liver" #For Figure 1: "plasma" #For Figure S1: "liver"

#Select normalized dataset
if(tis=="plasma") {
  dataset <- dataset.plasma.lipids$norm
  } else if (tis=="liver") {
  dataset <- dataset.liver.lipids$norm}

lipids <- lipids[which(lipids$tissue==tis),]
rownames(lipids) <- lipids$identifier
lipid_qtls <- lipid_qtls[which(lipid_qtls$tissue==tis),]

#removes int.St. and other hiccups
dataset <- as.data.frame(t(dataset[,which(colnames(dataset)%in%lipids$identifier)]))
dataset$identifier <- rownames(dataset)

#combine extended info from lipids table and quant. data from dataset
data <- merge(lipids,dataset,by = "identifier")
rownames(data) <- data$identifier

#Do calcs on quant info
start <- which(colnames(data)=="DO021")
end <- length(data[1,])

data$max <- apply(data[,start:end],1,max) #ADD COL FOR MAXIMUM VALUE
data$min <- apply(data[,start:end],1,min) #ADD COL FOR MINIMUM VALUE
data$mean <- apply(data[,start:end],1,mean) #ADD COL FOR AVERAGE VALUE
data$DR <- data$max - data$min #ADD COL FOR DYNAMIC RANGE
data$SD <- apply(data[,start:end],1,sd) #ADD COL FOR STANDARD DEVIATION

###############################################################################
#DYNAMIC RANGE#

ggplotly(ggplot(data, aes(identifier,SD,color=category)) + geom_point()) #helps you pick out candidates

ID <- subset(data, identifier=="GM3.NGNA_d18.1_22.0_9.799_1253.81067_plus")
#Plasma ID: "PC_41.4_11.136_852.65076_plus"
#Plasma UNK: "GM3.NGNA_d18.1_22.0_9.799_1253.81067_plus"
#Liver ID: "TG_16.0_18.1_20.0_1_20.103_906.84802_plus"
#Liver UNK: "UNK_5.124_465.30432_minus"

#UPDATE DATASET SELECTION HERE BASED ON TISSUE
mrg <- data.frame("Sex" = ifelse(dataset.plasma.lipids$covar[,1]==0,coon_grey,"grey"), "Lip" = as.numeric(ID[,24:407]))
mrg <- mrg[order(mrg$Lip),]
mrg <- data.frame("Rank" = c(1:384), mrg)

pdf(file=paste0(ID$analyte,"_DR.pdf"), pointsize = 8, width = 3.75, height = 2, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), cex.main = 1, font.main = 1, cex.lab = 1)
plot(mrg$Rank, mrg$Lip, col=mrg$Sex, pch=16, bty="l",
     ylim=c(6,16.75), 
     xlab="384 DO mice, rank-ordered",
     ylab="Log2(norm. peak int.)")
title(ID$analyte, line = -1, adj = 0.04)
legend("bottomright", legend=c("Male", "Female"),
       col=c("grey",coon_grey), pch=16, bty = "n", cex=1)
dev.off()

###############################################################################
theme_set(theme_classic())
pdf(file=paste0(tis,"_mzrt.pdf"), pointsize = 8, width = 4.5, height = 3.75, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), cex.main = 1, font.main = 1, cex.lab = 1)

ggplot() +

geom_point(data=subset(data,id.status!="identified"),
           aes(x=rt,y=mz,color="Unidentified",alpha=0.75)) + #making sure the hand-IDd ones are shown as UNK
geom_point(data=subset(data,(id.status=="identified"|category=="Fatty.Acyl")),
           aes(x=rt,y=mz,color=as.factor(category)),alpha=0.75) + 

scale_color_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
scale_alpha(guide = 'none') + #remove unuseful legend
scale_x_continuous(name="RT [min]") +
scale_y_continuous(name="m/z") +
theme(legend.position = "none")

dev.off()

