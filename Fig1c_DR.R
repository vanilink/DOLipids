options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("plotly")
library(cowplot)

setwd("E:/DO/R-Analyses/Feb_2019/") 
load("E:/DO/R-Analyses/Nov_2018/DO_Lipidomics.RData") #has Coon colors

#GOAL: make dynamic range plot for both tissues (for Fig 2)

# PREPARE DATASET FOR PLOTTING ################################################
plasma_lipids <- lipids[which(lipids$tissue=="plasma"),]
liver_lipids <- lipids[which(lipids$tissue=="liver"),]

dataset.key <- c("pheno"=7, "norm"=8, "raw"=9)
this.key <- dataset.key["norm"] #CAN ADD PHENO, norm OR RAW VALUES

this.dataset <- dataset.plasma.lipids[[this.key]]
forplot.plasma <- cbind(plasma_lipids, t(this.dataset[,which(colnames(this.dataset)%in%plasma_lipids$identifier)])) 

this.dataset <- dataset.liver.lipids[[this.key]]
forplot.liver <- cbind(liver_lipids, t(this.dataset[,which(colnames(this.dataset)%in%liver_lipids$identifier)])) 

forplot <- rbind(forplot.liver, forplot.plasma)
rm(dataset.key, this.key, this.dataset, forplot.plasma, forplot.liver, plasma_lipids, liver_lipids)

l <- length(forplot[1,])
forplot$max <- apply(forplot[,24:l],1,max) #ADD A COL FOR THE MAXIMUM VALUE
forplot$min <- apply(forplot[,24:l],1,min)
forplot$range <- forplot$max - forplot$min
forplot$mean <- apply(forplot[,24:l],1,mean) #ADD A COL FOR THE AVERAGE VALUE
forplot$stddev <- apply(forplot[,24:l],1,sd)
rownames(forplot) <- forplot$identifier
###############################################################################

ggplotly(ggplot(forplot, aes(identifier,stddev,color=category)) + 
           geom_point())

tis = "liver"

UNK <- subset(subset(forplot,id.status!="identified"),tissue==tis)
rownames(UNK) <- UNK$identifier
subset(UNK, stddev>1.5)$identifier

UNK <- subset(forplot, identifier=="UNK_9.863_766.0896_plus") #"GM2.NGNA_d18.1_22.0_9.756_1456.88892_plus") #"UNK_9.85_828.60968_plus" UNK_8.282_1197.74463_plus
#subset(forplot,stddev==max(subset(subset(forplot,id.status!="identified"),tissue==tis)$stddev))
  
ID <- subset(subset(forplot,id.status=="identified"),tissue==tis)
rownames(ID) <- ID$identifier
subset(ID, stddev>1.8)$identifier

ID <- subset(forplot, identifier=="TG_54.6_2_16.785_896.7688_plus")#"PC_16.1_22.6_2_8.21_804.55225_plus")
#subset(forplot,stddev==max(subset(subset(forplot,id.status=="identified"),tissue==tis)$stddev))

mrg <- rbind("Sex" = dataset.plasma.lipids$covar[,1], UNK[,24:407])
mrg <- mrg[,order(mrg[2,])]
mrg <- rbind("Rank" = c(1:384), mrg)
mrg <- as.data.frame(t(mrg))

ggplot(mrg, aes(x=Rank, y=mrg[,3], color=as.logical(Sex))) + geom_point() +
  scale_x_continuous(name="384 DO mice, rank-ordered") +
  scale_y_continuous(name="Log2 (normalized peak intensity)", limits=c(5,17)) +
  ggtitle(names(mrg)[3]) + theme(legend.position = "none")

plot(c(1:384),sort(mrg[,3]), #sort() if rank-ordered
     ylim=c(5,21), #xaxt="n",
     xlab="384 DO mice, rank-ordered",
     ylab="Log2 (normalized peak intensity)", 
     main="Plasma Lipidomics Dynamic Range", #change according to tissue
     col=mrg$Sex, pch=16,
     mgp=c(1.8,0.5,0))

points(c(1:384),sort(ID[24:407]), #sort() if rank-ordered
       col=coon_turq, pch=16) #change color according to lipid class

legend(-8, 22.5, 
       legend=c(UNK$identifier,
                ID$identifier),
       col=c("grey",coon_turq), pch=16, #change color according to lipid class
       bty = "n", cex=0.75)


############### ALL DR ####################
test <- cbind(forplot, "DR"=2^(forplot$max - forplot$min))

ggplot(subset(test,id.status=="identified"), aes(x=2^(max-min),fill=tissue)) + #geom_histogram(binwidth = 5) + 
  geom_density(alpha=0.8) + xlim(0,200) +
  geom_vline(aes(xintercept=mean(2^(max-min), na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1) +
  scale_fill_manual("Tissue", values=c(coon_yel,coon_red))
  #scale_fill_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, "grey"))
  
ggplot(test, aes(x=category,y=(max-min),fill=category,color=tissue)) + geom_violin() + ylim(0,11.5) +
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, "grey")) +
  scale_color_manual("Tissue", values=c("black","white"))
