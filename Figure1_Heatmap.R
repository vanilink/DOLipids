options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("plotly")
library("pheatmap")

setwd("E:/DO/R-Analyses/Summer_2019/") 
load("P_L_Lipids.RData")

#GOAL: make 384 mice vs 2,000 features heatmap for both tissues (for Fig 1/S1)

# SELECT TISSUE 
tis <- "plasma" #For Figure 1: "plasma" #For Figure S1: "liver"

#Select normalized dataset
if(tis=="plasma") {
  dataset.lipids <- dataset.plasma.lipids
  dataset <- dataset.plasma.lipids$norm
} else if (tis=="liver") {
  dataset.lipids <- dataset.liver.lipids
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

#remove gangliosides from IDd ones
data[which(data$class=="Ganglioside"),]$category <- "UNK" 

#annotation_row is info on lipids, rownames are feature IDs
ano <- data.frame("Category" = factor(data[,7]))
rownames(ano) <- data$identifier
  
Var1 = c(coon_red, coon_purp, coon_blue, coon_yel, coon_turq, "grey")
names(Var1) = unique(data$category)
ann_colors = list(Category = Var1)

start <- which(colnames(data)=="DO021")
end <- length(data[1,])
hm.data <- data[,start:end]

heat <- pheatmap(hm.data[,], 
         annotation_row = ano, 
         annotation_colors = ann_colors,
         scale = "row",
         breaks = c(-10,seq(-3,3,length.out=99),10),
         show_rownames = F, show_colnames = F,
         annotation_legend = F)

#export pdf with 4 x 4
# 
# library("fields")
# dist.data <- rdist(hm.data)
# rownames(dist.data) <- rownames(hm.data)
# colnames(dist.data) <- rownames(hm.data)
# 
# for (i in 1:6) {
#   class.lipids <- subset(data, category==names(table(data$category)[i]))$identifier
#   class <- subset(dist.data, rownames(dist.data) %in% class.lipids)
#   class <- subset(class, select=(colnames(class) %in% class.lipids))
#   class[upper.tri(class, diag=T)] <- NA
#   print(paste0("The mean distance of class ", names(table(data$category)[i]), " is ", as.integer(median(class, na.rm=T))))
# }
# 
# sterol.lipids <- subset(data, category=="Sterol.Lipid")$identifier
# sterol <- subset(dist.data, rownames(dist.data) %in% sterol.lipids)
# sterol <- subset(sterol, select=(colnames(sterol) %in% sterol.lipids))
# sterol[upper.tri(sterol, diag=T)] <- NA
# mean(sterol, na.rm=T)

#pdf(file=paste0(tis,"_HM_mice.pdf"), pointsize = 8, width = 3, height = 4, useDingbats=FALSE)
#heat
#dev.off()

new_data <- data.frame("order"=as.numeric(1:dim(hm.data)[1]), "category" = data[heat$tree_row$order,]$category)

pdf(file=paste0(tis,"_class_density.pdf"), pointsize = 8, 
    width = 3.3, height = 1.6, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), 
    cex.main = 1, font.main = 1, cex.lab = 1)

#ggplotly(
  ggplot(subset(new_data,category!="UNK"), aes(order, stat(count), fill=category)) + 
    geom_density(adjust=.5, alpha=0.8)+
    #facet_wrap(~category) + 
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
    theme(legend.position = "none") +
    scale_x_reverse() #+ coord_flip()
  
  #)
dev.off()
