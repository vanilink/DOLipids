options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("plotly")
library(cowplot)

setwd("E:/DO/R-Analyses/Feb_2019/") 
load("E:/DO/R-Analyses/Nov_2018/DO_Lipidomics.RData") #has Coon colors

#GOAL: make m/z vs RT plot for both tissues (for Fig 2)
# by color and size by intensity

# PREPARE DATASET FOR PLOTTING ################################################
plasma_lipids <- lipids[which(lipids$tissue=="plasma"),]
liver_lipids <- lipids[which(lipids$tissue=="liver"),]

dataset.key <- c("pheno"=7, "norm"=8, "raw"=9)
this.key <- dataset.key["raw"] #CAN ADD PHENO, norm OR RAW VALUES

this.dataset <- dataset.plasma.lipids[[this.key]]
forplot.plasma <- cbind(plasma_lipids, t(this.dataset[,which(colnames(this.dataset)%in%plasma_lipids$identifier)])) 

this.dataset <- dataset.liver.lipids[[this.key]]
forplot.liver <- cbind(liver_lipids, t(this.dataset[,which(colnames(this.dataset)%in%liver_lipids$identifier)])) 

forplot <- rbind(forplot.liver, forplot.plasma)
rm(dataset.key, this.key, this.dataset, forplot.plasma, forplot.liver, plasma_lipids, liver_lipids)

l <- length(forplot[1,])
forplot$max <- apply(forplot[,24:l],1,max) #ADD A COL FOR THE MAXIMUM VALUE
forplot$mean <- apply(forplot[,24:l],1,mean) #ADD A COL FOR THE AVERAGE VALUE
rownames(forplot) <- forplot$identifier

#SELECT TISSUE FOR PLOTTING
tis <- "liver"

#ACTUAL PLOT
#plot <- #only needed for save_plot
ggplotly( 
  ggplot() +
  
  geom_point(data=subset(forplot,tissue==tis&id.status!="identified"), aes(x=rt,y=mz,color="Unidentified",#size=mean,
                                                                           label=identifier,alpha=0.75)) + #making sure the hand-IDd ones are shown as UNK
  geom_point(data=subset(forplot,tissue==tis&(id.status=="identified"|category=="Fatty.Acyl")), aes(x=rt,y=mz,color=as.factor(category),#size=mean,
                                                                                                    label=identifier),alpha=0.75) + 
  
  scale_color_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
  ggtitle("Liver Lipidomics") + #make sure to change
  scale_alpha(guide = 'none') + #remove unuseful legend
  scale_size_continuous("Average Intensity", range = c(1.75, 5), breaks = c(10^4,10^6,10^8)) +
  scale_x_continuous(name="RT [min]") +
  scale_y_continuous(name="m/z") +
  theme(legend.position = "none")
)

# could be optimized works so-so
#save_plot("plot.pdf", plot, base_aspect_ratio = 1.3 # make room for figure legend

