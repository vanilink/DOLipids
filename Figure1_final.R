# GOAL: Create panels for figures 1 and S1

options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")

setwd("E:/DO/R-Analyses/GithubDOLipids/")

# SELECT TISSUE #For Figure 1: "plasma" #For Figure S1: "liver"
tis <- "plasma" 

# LOAD IN LIPID AND LIPID QTL DATA
lipids <- read.csv("lipids.csv")[,-1]
lipid_qtls <- read.csv("lipid_qtls.csv")[,-1]

lipids <- lipids[which(lipids$tissue==tis),]
rownames(lipids) <- lipids$identifier
lipid_qtls <- lipid_qtls[which(lipid_qtls$tissue==tis),]

# LOAD IN DATASET
if(tis=="plasma") {
  dataset <- readRDS("dataset.plasma.lipids.RDS")
  } else if (tis=="liver") {
    dataset <- readRDS("dataset.liver.lipids.RDS")
  }

#removes int.St. and other hiccups
data <- as.data.frame(t(dataset$norm[,which(colnames(dataset$norm)%in%lipids$identifier)]))
data$identifier <- rownames(data)

#remove gangliosides from IDd ones as they were later hand-identified
data[which(data$class=="Ganglioside"),]$category <- "UNK" 

#combine extended info from lipids table and quant. data from dataset
data <- merge(lipids,data,by = "identifier")
rownames(data) <- data$identifier


###############################################################################
# b. m/z vs. RT plot
theme_set(theme_classic())
pdf(file=paste0(tis,"_mzrt.pdf"), pointsize = 8, width = 4.5, height = 3.75, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), cex.main = 1, font.main = 1, cex.lab = 1)

ggplot() +
  
  geom_point(data=subset(data,category=="UNK"), #making sure unidentifed features are drawn first
             aes(x=rt,y=mz,color="Unidentified",alpha=0.75)) + 
  geom_point(data=subset(data,category!="UNK"),
             aes(x=rt,y=mz,color=as.factor(category)),alpha=0.75) + 
  
  scale_color_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
  scale_alpha(guide = 'none') + #remove unuseful legend
  scale_x_continuous(name="RT [min]") +
  scale_y_continuous(name="m/z") +
  theme(legend.position = "none")

dev.off()

###############################################################################
# c. DYNAMIC RANGE

ID <- subset(data, identifier=="PC_41.4_11.136_852.65076_plus")
#Plasma ID: "PC_41.4_11.136_852.65076_plus"
#Plasma UNK: "GM3.NGNA_d18.1_22.0_9.799_1253.81067_plus"
#Liver ID: "TG_16.0_18.1_20.0_1_20.103_906.84802_plus"
#Liver UNK: "UNK_5.124_465.30432_minus"

mrg <- data.frame("Sex" = ifelse(dataset$covar[,1]==0,coon_grey,"grey"), 
                  "Lip" = as.numeric(ID[,24:407]))
mrg <- mrg[order(mrg$Lip),] #rank-order
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
# d. Mice vs. Lipids Heatmap

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

#pdf(file=paste0(tis,"_HM_mice.pdf"), pointsize = 8, width = 3, height = 4, useDingbats=FALSE)
#heat
#dev.off()

###############################################################################
# d. Lipid Class Density Plot

new_data <- data.frame("order"=as.numeric(1:dim(hm.data)[1]), 
                       "category" = data[heat$tree_row$order,]$category)

pdf(file=paste0(tis,"_class_density.pdf"), pointsize = 8, 
    width = 3.3, height = 1.6, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), 
    cex.main = 1, font.main = 1, cex.lab = 1)

ggplot(subset(new_data,category!="UNK"), aes(order, stat(count), fill=category)) + 
  geom_density(adjust=.5, alpha=0.8)+
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
  theme(legend.position = "none") +
  scale_x_reverse()

dev.off()

###############################################################################
# e. Lipids vs. Genome Heatmap
# ATTENTION: This is a very large heatmap - make sure your computer is equipped!

# LOAD IN LOD DATA FOR HEATMAP

#Select normalized dataset
if(tis=="plasma") {
  comb <- read.csv("Pcombined.csv", nrows = 69006, 
                   colClasses = c("integer", "character", rep("numeric",1415)))
} else if (tis=="liver") {
  comb <- read.csv("combined.csv", nrows = 69006, 
                   colClasses = c("integer", "character", rep("numeric", 1148)))
}

rownames(comb) <- comb[,2] #69005 marker names as rownames
#colnames are already lipid feature IDs
comb <- t(comb[,-1:-2]) #1 is just relic numbers from csv, 2 is now stored as rownames
#transpose data to resemble heatmap: lipids in rows, markers in cols

# ensure identical feature IDs - needed?
cheat <- read.csv("ID_Cheatsheet_Liver.csv")
tobereplaced <- rownames(comb)[which(rownames(comb)%in%cheat$NEW_ID)]
replacement <- cheat$OLD_ID[which(tobereplaced%in%cheat$NEW_ID)]
rownames(comb)[which(rownames(comb)%in%tobereplaced)] <- replacement

#subset comb so there is no NAs (resulting from internal standard features)
comb <- subset(comb, (rownames(comb) %in% lipids$identifier))

#annotation_row is info on lipids, rownames are lipid feature IDs
ids <- data.frame("ID" = rownames(comb)) 
ids$order <- 1:nrow(comb)
#get lipid info by connecting lipid feature IDs (unique) to lipids table
id.data <- merge(x=ids, y=lipids, all.x=F, all.y=F, by.x="ID", by.y="identifier")
id.data <- id.data[order(id.data$order),]
ano <- data.frame("category" = factor(id.data$category))
rownames(id.data) <- id.data$identifier
rownames(ano) <- factor(id.data$ID)
rm(ids)

#annotation_col is info on genetic location (chr), rownames are marker IDs
#get access to this data from colnames of markers
ano_col <- data.frame("Chr" = factor(sapply(strsplit(colnames(comb),"_"), `[`, 1)))
rownames(ano_col) <- colnames(comb) #have ano_col have same marker names

Var1 = c(coon_red, coon_purp, coon_blue, coon_yel, coon_turq, "grey")
# Liver: c("grey", coon_yel, coon_turq, coon_purp, coon_red, coon_blue) 
names(Var1) = unique(ano$category)
Var2 = rep(c("grey", coon_grey),10)
names(Var2) = unique(ano_col[,])
ann_colors = list(category = Var1, Chr = Var2)

gaps = NULL
for (chr in as.character(c(2:19,"X"))){
  gaps <- append(gaps,which(ano_col$Chr==chr)[1]) #add gap at each first occurence of chr
}

gc() #garbage collector to make sure there is enough space for large heatmap

png(paste0(tis,"_HM_Chr.png"), width = 5.2, height = 5.24, unit="in", res=300) #can only be saved as png

heat <- pheatmap(comb[,], 
                 color = colorRampPalette(brewer.pal(n = 7, name = "Oranges"))(100),
                 cluster_cols = F, #we need chromosomes in order
                 gaps_col = gaps,
                 annotation_row = ano, 
                 annotation_col = ano_col, 
                 annotation_colors = ann_colors, #for chr and lipid classes
                 breaks = c(0,seq(min(comb),10,length.out=99),max(comb)), #make LOD contrast meaningful: 1 - 10
                 show_rownames = F, 
                 show_colnames = F,
                 annotation_legend = F)

dev.off()

saveRDS(heat, file = paste0(tis,"_genome_HM.RDS"))

###############################################################################
# e. Lipid Class Density Plot
new_data <- data.frame("order"=as.numeric(1:dim(comb)[1]), #or as.numeric(1:1405), #1085 for Liver !
                       "category" = id.data[heat$tree_row$order,]$category)

pdf(file=paste0(tis,"_Gen_class_density.pdf"), pointsize = 8, 
    width = 3.3, height = 1.6, useDingbats=FALSE)
par(mar=c(3,2.75,0.5,0.5), mgp = c(1.5, 0.5, 0), 
    cex.main = 1, font.main = 1, cex.lab = 1)

ggplot(subset(new_data,category!="UNK"), aes(order, stat(count), fill=category)) + 
  geom_density(adjust=.5, alpha=0.8)+
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_purp, coon_turq, coon_yel, coon_blue, coon_grey)) +
  theme(legend.position = "none") +
  scale_x_reverse()

dev.off()

