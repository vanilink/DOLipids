options(stringsAsFactors = F)
library("ggplot2")
library("qtl2")
library("tidyr")
library("reshape2")
library(dplyr)
library(RColorBrewer)
library(pheatmap)

setwd("E:/DO/R-Analyses/GithubDOLipids/")

# define colors
coon_blue <- "#2CA7DF"
coon_grey <- "#566977"
coon_purp <- "#955CA5"
coon_red <- "#EC6B63"
coon_turq <- "#63C29C"
coon_yel <- "#FFCB04"

# for liver eQTL and pQTL
load("E:/DO/R-Analyses/Svenson_DO850_for_eQTL_viewer_v4.Rdata")

# LOAD IN LIPID AND LIPID QTL DATA
lipids <- read.csv("TableS8.csv")
rownames(lipids) <- lipids$identifier

lipid_qtls <- read.csv("TableS9.csv")
lipid_qtls <- merge(lipid_qtls, lipids)

#remove gangliosides from IDd ones as they were later hand-identified
lipid_qtls[which(lipid_qtls$class=="Ganglioside"),]$category <- "UNK" 
lipids[which(lipids$class=="Ganglioside"),]$category <- "UNK" 

#load("Figure4.RData")

###############################################################################
# Load in B6 quantitative lipidomics data (4 females, 4 males)
B6_plasma <- read.csv("E:/FS/FS_Plasma/3-LipiDex-Identification/Results/FS_Plasma_B6_Sex.csv")
B6_plasma <- B6_plasma[-1,] # remove info on sex (1-4 is F, 5-8 is M)
rownames(B6_plasma) <- B6_plasma[,1] # give IDs as rownames
B6_plasma <- apply(B6_plasma[,2:9], c(1,2),as.numeric)  #remove ID column and save quant values as numeric
B6_plasma <- scale(B6_plasma, center=FALSE, scale=colSums(B6_plasma)) #sum-normalize each sample

ttest <- function(df, grp1, grp2){
  x = df[grp1]
  y = df[grp2]
  x = as.numeric(x)
  y = as.numeric(y)
  results = t.test(x,y)
  results$p.value
}

#calculate raw p-value M vs F for each lipid
rawpvalue = apply(B6_plasma, 1, ttest, grp1=c(1:4), grp2=c(5:8)) #1:4 is F, 5:8 is M

B6_plasma_log2 <- log2(B6_plasma) #log2 transform quan data
B6_F <- apply(B6_plasma_log2[,1:4],1,mean) #save average female value for each lipid
B6_M <- apply(B6_plasma_log2[,5:8],1,mean) #save average male value for each lipid
foldchange <- B6_F-B6_M # calculate FC as difference between log2 averages F-M

# for volcano plot: save FC and p-values in one dataframe
volcano = as.data.frame(cbind(foldchange,rawpvalue))
volcano$B6.ID <- rownames(volcano) #save IDs in separate column

#convert B6 identifications into info about analyte (m/z, polarity)
B6mz = sapply(strsplit(volcano$B6.ID,"*MZ: "), `[`, 2)
B6.polarity<-gsub(")","",sapply(strsplit(B6mz,"\\("), `[`, 2))

B6mz <- as.numeric(gsub(" ?\\(.\\)", "", B6mz))
volcano$B6mz <- B6mz
volcano$B6.polarity <- ifelse(B6.polarity=="+","plus","minus")
volcano$B6.UNK <- grepl("Unk",volcano$B6.ID)

p.cutoff = 0.05 #set p value cutoff
FC.cutoff = 1 # set fold change cutoff

volcano$sig <- ifelse(abs(volcano$foldchange)>FC.cutoff&
                        volcano$rawpvalue<p.cutoff,
                      yes = TRUE, no = FALSE)

info <- table(volcano$B6.UNK, volcano$sig)
cat(c("In the B6 dataset,", sum(info[3:4]), "out of", sum(info[1:4]), 
      "lipid features were different by sex, \nand", info[3], "of these were identified."))

# Now link B6 info to DO info to assess overlap and QTLs

volcano$inDO <- FALSE #instantiate inDO and DOmz to be not matching
volcano$DOmz <- 0

#calculate lipid masses plus and minus 10 ppm for all DO lipids
lipids$m10ppm <- lipids$mz - lipids$mz*10/1e6
lipids$p10ppm <- lipids$mz + lipids$mz*10/1e6

#reorder lipids so that in next step identified lipids will be picked last (preferentially)
lipids <- lipids[order(lipids$id.status,decreasing = T),]

# ! for loops rather slow:
for (B6mz in volcano$B6mz) { #go through all B6 features
  for (DOmz in lipids$mz) #for each DO feature check whether within 10 ppm
    
    if (volcano[which(volcano$B6mz==B6mz),]$B6.polarity == lipids[which(lipids$mz==DOmz),]$polarity #needs to be in same polarity
        && B6mz > lipids[which(lipids$mz==DOmz),]$m10ppm # B6mz needs to be above DOmz-10 pmm
        && B6mz < lipids[which(lipids$mz==DOmz),]$p10ppm) { #B6mz needs to be beow DOmz+10 ppm
      volcano[which(volcano$B6mz==B6mz),]$inDO <- TRUE #then set inDO to ture
      volcano[which(volcano$B6mz==B6mz),]$DOmz <- DOmz #and set DOmz to mz from DO lipids
      cat(c("\nB6:",volcano[which(volcano$B6mz==B6mz),]$B6.ID)) #print out to check whether it is doing it correctly and to see progress
      cat(c("\nDO:",lipids[which(lipids$mz==DOmz),]$identifier))
    }
} #still issues: should give preference to plasma DO ID, should give preference to closer RT
#warnings() In volcano[which(volcano$B6mz == B6mz), ]$B6.polarity ==  ... :
#longer object length is not a multiple of shorter object length

volcano_new <- inner_join(x = volcano, y=lipids, by=c("DOmz"="mz"))

#only call UNK if both in B6 and in DO unidentified:
volcano_new$UNK <- ifelse(volcano_new$B6.UNK==F,F,ifelse(volcano_new$id.status!="un-identified",F,T))

#add QTL info for chr6:91 +-2 Mbp and use ident to color accordingly
chr6.91Mbp <- subset(lipid_qtls[which(lipid_qtls$mz%in%volcano_new$DOmz),],qtl.chr==6&qtl.pos>89&qtl.pos<93)
volcano_new$ident  <- with(volcano_new, ifelse(identifier %in% chr6.91Mbp$identifier, "?", ifelse(UNK,"unidentified","identified")))

write.csv(volcano_new, "SupplementalTableS6.csv")

# Different approach, start from QTLs, then add info about whether or not in B6 dataset
###############################################################################
# S4a. Heatmap of allele effects at Chr6:91 Mbp locus

my.chr6.91Mbp <- subset(lipid_qtls, qtl.chr=="6" & qtl.pos>(91-2) & qtl.pos<(91+2))

ano <- data.frame("Class" = factor(my.chr6.91Mbp$category), 
                  "Tissue" = factor(my.chr6.91Mbp$tissue),
                  "inB6" = factor(my.chr6.91Mbp$identifier %in% volcano_new$identifier)
                    )
rownames(ano) <- my.chr6.91Mbp$identifier

Var1 = c(coon_yel, coon_turq, coon_purp, "grey")
names(Var1) = unique(my.chr6.91Mbp$category)
Var2 = c("black","white")
names(Var2) = unique(my.chr6.91Mbp$tissue)
Var3 = c(coon_red, "white")
names(Var3) = c(T,F)

ann_colors = list(Class = Var1, Tissue = Var2, inB6 = Var3)

cluster.my.chr6.91Mbp <- my.chr6.91Mbp[,28:35] #reduce to allele effects
rownames(cluster.my.chr6.91Mbp) <- my.chr6.91Mbp$identifier
colnames(cluster.my.chr6.91Mbp) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

heat <- pheatmap(cluster.my.chr6.91Mbp, 
                 annotation_row = ano, 
                 annotation_colors = ann_colors,
                 cluster_cols = F, 
                 show_rownames = F, 
                 breaks = seq(-2,2,length.out=100),
                 cutree_rows = 3, 
                 scale = "row",
                 main = "Allele Effect at Chr 6:91 +-2 Mbp",
                 annotation_legend = F,
                 filename = "Chr6_HM.pdf",
                 width = 2.8,
                 height = 3.5
                 )
dev.off()

table(cutree(heat$tree_row,h = 4)) 
sub.chr6.91Mbp <- cutree(heat$tree_row,h = 4)
sub.chr6.91Mbp <- names(sub.chr6.91Mbp[which(sub.chr6.91Mbp==3)])

chr691Mbp <- subset(lipid_qtls,identifier %in% sub.chr6.91Mbp & qtl.chr=="6")
chr691Mbp$inB6 <- chr691Mbp$identifier %in% volcano_new$identifier

###############################################################################
# S4b. FS Plasma boxplot of feature at Chr6:91 Mbp
fs <- read.csv("20170623_FS_Plasma_Log2.csv")

rownames(fs)[1:5] <- fs[1:5,6]

rownames(fs)[6:nrow(fs)] <- paste0(fs[6:nrow(fs),4],"_", fs[6:nrow(fs),1],"_",fs[6:nrow(fs),2],fs[6:nrow(fs),3])

colnames(fs)[7:ncol(fs)] <- fs[2, 7:ncol(fs)]

fs.q <- fs[-2:-5,-1:-6]
fs.sq <- fs[c(-2,-4,-5),-1:-6] #including sex

row <- which(grepl("_1130.", rownames(fs.q))) #1814

loi <- fs.sq[c(1,2,row+1),1:62]
rownames(loi)[3] <- "analyte"
loi <- data.frame(t(loi))

loi$analyte <- as.numeric(loi$analyte)
loi$Strain <- factor(loi$Strain,
                     levels = names(CCcolors))
loi$SexStrain <- factor(paste(loi$Strain, loi$Sex))

pdf(file=paste0("1130_boxplot.pdf"), pointsize = 6, 
    width = 4, height = 2.5, useDingbats=FALSE)

ggplot(loi, aes(x=SexStrain, y=analyte, fill = Strain)) + 
  stat_boxplot(geom ="errorbar") +
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  scale_fill_manual(values=CCcolors) +
  theme_classic() +
  ylab("Log2 Peak Area") +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

dev.off()

###############################################################################
# 4b. Volcano of B6 M/F
# add highlights based on which ones are at Chr 6:91 Mbp and whether or not id'd
volcano_new$col  <- with(volcano_new, ifelse(identifier %in% chr691Mbp$identifier, "?", ifelse(UNK,"unidentified","identified")))

pdf(file=paste0("Chr6_Volcano.pdf"), pointsize = 6, 
    width = 2, height = 3, useDingbats=FALSE)

ggplot(volcano_new, aes(x=foldchange,y=-1*log10(rawpvalue), 
                                 lable=B6.ID, lable2 = identifier,
                                 color=as.factor(col),
                                 alpha = 0.8,
                                 size=sig)) + 
           geom_point() +
           scale_color_manual("In B6", values=c(coon_red, coon_grey, "grey")) +
           scale_x_continuous(limits = c(-5,5), name="Fold Change F/M") +
           scale_y_continuous(limits = c(0,5), name="-log10(p-value)") +
           scale_size_manual("Significant", values = c(0.5,1.5)) +
           theme_classic() + 
           theme(legend.position="none")

dev.off()

###############################################################################
# 4c. Manhattan plot at Chr 6:91 Mbp
chr691Mbp$col  <- with(chr691Mbp, ifelse(identifier %in% volcano_new$identifier, "?", ifelse(class=="UNK","unidentified","identified")))

pdf(file=paste0("Chr6_Manhattan.pdf"), pointsize = 6, 
    width = 3, height = 3, useDingbats=FALSE)

ggplot(data=chr691Mbp, aes(x=qtl.pos,y=qtl.lod,
                                     color=col,
                                     label=identifier,alpha=0.85,size=1.5)) +
           geom_point() + 
           scale_x_continuous(limits = c(89,93), name="Chr 6 Pos [Mbp]") +
           scale_color_manual("In B6", values=c(coon_red, coon_grey, "grey")) +
           scale_y_continuous(name="LOD") + 
           theme_classic() +
           theme(legend.position="none")

dev.off()

###############################################################################
# 4d. MZ-RT at Chr 6:91 Mbp
pdf(file=paste0("Chr6_MZRT1.pdf"), pointsize = 6, 
    width = 2, height = 1.5, useDingbats=FALSE)

  ggplot(data=chr691Mbp) +
    geom_point(aes(x=rt,y=mz,color=as.factor(col), 
                   label=identifier),
               alpha=0.85) + 
    scale_color_manual("Identification",values = c(coon_red, coon_grey, "grey")) +
    scale_x_continuous(limits=c(0,20), name="RT [min]") +
    scale_y_continuous(limits=c(200,1600), name="m/z") +
    theme_classic() +
    theme(legend.position="none")
  
dev.off()

pdf(file=paste0("Chr6_MZRT2.pdf"), pointsize = 6, 
    width = 2, height = 1.5, useDingbats=FALSE)

ggplot(data=chr691Mbp) +
  geom_point(aes(x=rt,y=mz,color=as.factor(col), #size=sig, 
                 label=identifier),
             alpha=0.85) + 
  scale_color_manual("Identification",values = c(coon_red, coon_grey, "grey")) +
  scale_x_continuous(limits=c(12,17), name="RT [min]") +
  scale_y_continuous(limits=c(1000,1200), name="m/z") +
  theme_classic() +
  theme(legend.position="none")

dev.off()


###############################################################################
# S4c. Abhd1 allele effect cluster heatmap

my.abhd1 <- subset(lipid_qtls, qtl.chr=="5" & qtl.pos>(31-2) & qtl.pos<(31+2))

ano <- data.frame("Class" = factor(my.abhd1$category), 
                  "Plasma" = factor(my.abhd1$tissue),
                  "FA" = factor(my.abhd1$identifier %in% my.abhd1[grep(pattern = "14.0",my.abhd1$chain),]$identifier),
                  "LysoPC" = factor(my.abhd1$class=="LysoPC")
)

rownames(ano) <- my.abhd1$identifier

Var1 = c(coon_turq, coon_purp, "grey")
names(Var1) = unique(my.abhd1$category)
Var2 = c("black","white")
names(Var2) = unique(my.abhd1$tissue)
Var3 = c(coon_turq, "white")
names(Var3) = c(T,F)
Var4 = c(coon_turq, "white")
names(Var4) = c(T,F)

ann_colors = list(Class = Var1, Plasma = Var2, FA = Var3, LysoPC = Var4)

cluster.my.abhd1 <- my.abhd1[,28:35] #reduce to allele effects
rownames(cluster.my.abhd1) <- my.abhd1$identifier
colnames(cluster.my.abhd1) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

heat <- pheatmap(cluster.my.abhd1, 
                 annotation_row = ano, 
                 annotation_colors = ann_colors,
                 cluster_cols = F, 
                 show_rownames = F, 
                 clustering_method = "ward.D2",
                 breaks = seq(-2,2,length.out=100),
                 cutree_rows = 3, 
                 scale = "row",
                 main = "Allele Effect at Chr 5:31 +-2 Mbp",
                 annotation_legend = F,
                 filename = "Abhd1_HM.pdf",
                 width = 4,
                 height = 4
)
dev.off()

table(cutree(heat$tree_row,h = 8)) 
sub.abhd1 <- cutree(heat$tree_row,h = 8)
sub.abhd1 <- names(sub.abhd1[which(sub.abhd1==1)])

abhd1 <- subset(lipid_qtls,identifier %in% sub.abhd1 & qtl.chr=="5")

###############################################################################
# 4i. Abhd1 Manhattan plot
abhd1$col  <- with(abhd1, ifelse(identifier %in% abhd1[grep(pattern = "14.0",chain),]$identifier,"14.0",ifelse(class=="UNK","unidentified","identified")))

pdf(file=paste0("Abhd1_Manhattan.pdf"), pointsize = 6, 
    width = 2.8, height = 3.5, useDingbats=FALSE)

ggplot(data=abhd1, aes(x=qtl.pos,y=qtl.lod,
                       color=col,
                       label=identifier,alpha=0.85,size=1.5)) +
  geom_point() + 
  scale_x_continuous(limits = c(29,33), name="Chr 5 Pos [Mbp]") +
  scale_color_manual("FA 14:0", values=c(coon_turq, coon_grey, "grey")) +
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()

###############################################################################
# 4j. Manhattan plot of subset of LysoPCs across all chromosomes

dim(subset(lipids, class=="LysoPC")) #45
dim(subset(lipid_qtls, class=="LysoPC")) #118

chrLabels <- c("1"="1","2"="2","3"="3","4"="4","5"="5","6"="6","7"="7","8"="8","9"="9","10"="10","11"="11","12"="12","13"="13","14"="14","15"="15","16"="16","17"="17","18"="18","19"="19","20"="X")

lipid_qtls$chr = factor(lipid_qtls$qtl.chr, levels=chrLabels)

pdf(file=paste0("LysoPC_Manhattan.pdf"), pointsize = 5, 
    width = 3, height = 1.5, useDingbats=FALSE)

ggplot() + geom_point(data = subset(lipid_qtls, class=="LysoPC" & !is.na(qtl.chr)), 
                      aes(x = qtl.pos, y = qtl.lod), color = coon_turq, alpha = 0.85) +
  facet_wrap(facets = ~chr, strip.position = "bottom", nrow = 1, scales = "free_x") +
  theme_classic() +
  scale_x_continuous(name="Chromosome") +
  scale_y_continuous(name="LOD") +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(size = 5)
  )

dev.off()

###############################################################################
# 4k. Manhattan plot of subset of14:0 containing lipids
dim(lipids[grep(pattern = "14.0",lipids$chain),]) #30
dim(lipid_qtls[grep(pattern = "14.0",lipid_qtls$chain),]) #65

pdf(file=paste0("FA140_Manhattan.pdf"), pointsize = 5, 
    width = 3, height = 1.5, useDingbats=FALSE)

ggplot() + geom_point(data = lipid_qtls[grep(pattern = "14.0",lipid_qtls$chain),], 
                      aes(x = qtl.pos, y = qtl.lod), color = coon_turq, alpha = 0.85) +
  facet_wrap(facets = ~chr, strip.position = "bottom", nrow = 1, scales = "free_x") +
  theme_classic() +
  scale_x_continuous(name="Chromosome") +
  scale_y_continuous(name="LOD") +
  theme(legend.position="none",
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.spacing = unit(0, "lines"),
        strip.text.x = element_text(size = 5)
  )

dev.off()

###############################################################################
# S4d. FS Plasma LysoPC 14:0 boxplot

fs <- read.csv("20170623_FS_Plasma_Log2.csv")

rownames(fs)[1:5] <- fs[1:5,6]

rownames(fs)[6:nrow(fs)] <- paste0(fs[6:nrow(fs),4],"_", #roughly combine feature IDs
                                   fs[6:nrow(fs),1],"_",
                                   fs[6:nrow(fs),2],
                                   fs[6:nrow(fs),3])

colnames(fs)[7:ncol(fs)] <- fs[2, 7:ncol(fs)]

fs.q <- fs[-2:-5,-1:-6] #quant data

row <- which(grepl("LysoPC 14:0", rownames(fs.q))) #70
loi <- fs.q[c(1,row),1:62] # Lipid of Interest
rownames(loi)[2] <- "analyte"
loi <- data.frame(t(loi))

loi$analyte <- as.numeric(loi$analyte)
loi$Strain <- factor(loi$Strain, levels = names(CCcolors))

pdf(file=paste0("LysoPC140_boxplot.pdf"), pointsize = 6, 
    width = 4, height = 2.5, useDingbats=FALSE)

ggplot(loi, aes(x=Strain, y=analyte, fill = Strain)) + 
  stat_boxplot(geom ="errorbar") +
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  scale_fill_manual(values=CCcolors) +
  theme_classic() +
  ylab("Log2 Peak Area") +
  theme(legend.position = "none", 
        axis.title.x = element_blank())

dev.off()

###############################################################################
# S4e. Abhd1 Liver eQTL allele effect plot
eqtl <- dataset.mrna$lod.peaks$additive
eqtl_effect <- unlist(eqtl[which(eqtl$annot.id=="ENSMUSG00000006638"&eqtl$cis),11:18])

pdf(file=paste0("Abhd1_eQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),
                 y=eqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-2.2,2.2)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()

###############################################################################
# S4f. Heatmap Abhd2
my.abhd2 <- subset(lipid_qtls, qtl.chr == "7" & qtl.pos>79.3-2 & qtl.pos<79.3+2)

ano <- data.frame("Class" = factor(my.abhd2$category), 
                  "Tissue" = factor(my.abhd2$tissue),
                  "SumComp" = factor(my.abhd2$sum.comp)
)

rownames(ano) <- my.abhd2$identifier

Var1 = c(coon_turq, "grey")
names(Var1) = unique(my.abhd2$category)
Var2 = c("black","white")
names(Var2) = unique(my.abhd2$tissue)
Var3 = c("white", coon_turq)
names(Var3) = c(T,F)

ann_colors = list(Class = Var1, Tissue = Var2, SumComp = Var3)

cluster.my.abhd2 <- my.abhd2[,28:35] #reduce to allele effects
rownames(cluster.my.abhd2) <- my.abhd2$identifier
colnames(cluster.my.abhd2) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

heat <- pheatmap(cluster.my.abhd2, 
                 annotation_row = ano, 
                 annotation_colors = ann_colors,
                 cluster_cols = F, 
                 show_rownames = F, 
                 clustering_method = "ward.D2",
                 breaks = seq(-2,2,length.out=100),
                 cutree_rows = 3, 
                 scale = "row",
                 main = "Allele Effect at Chr 5:31.0 +-2 Mbp",
                 annotation_legend = F,
                 filename = "Abhd2_HM.pdf",
                 width = 4,
                 height = 4
)
dev.off()

table(cutree(heat$tree_row,h = 2.5)) 
sub.abhd2 <- cutree(heat$tree_row,h = 2.5)
sub.abhd2 <- names(sub.abhd2[which(sub.abhd2==2)])

abhd2 <- subset(lipid_qtls,identifier %in% sub.abhd2 & qtl.chr=="7")

###############################################################################
# S4g. Abhd2 Manhattan plot

pdf(file=paste0("Abhd2_Manhattan.pdf"), pointsize = 6, 
    width = 2.8, height = 3.5, useDingbats=FALSE)

ggplot() +
  geom_point(data=abhd2, 
             aes(x=qtl.pos,y=qtl.lod,color=category, alpha=0.85)) + 
  scale_x_continuous(name="Chr 7 Pos [Mbp]", limits = c(79.3-2, 79.3+2)) +
  scale_color_manual("Lipid Class",values = c(coon_turq, "grey")) + 
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position= "none") 

dev.off()

#####################################################################
# S4h. Abhd2 Liver eQTL allele effect plot
eqtl <- dataset.mrna$lod.peaks$additive
eqtl_effect <- unlist(eqtl[which(eqtl$annot.id=="ENSMUSG00000039202"&eqtl$cis),11:18])

pdf(file=paste0("Abhd2_eQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),
                 y=eqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-2.2,2.2)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()

###############################################################################
# S4i. Abhd2 Liver pQTL allele effect plot
pqtl <- dataset.protein$lod.peaks$additive
pqtl_effect <- unlist(pqtl[which(pqtl$annot.id=="ENSMUSP00000038361"&pqtl$qtl_chr==pqtl$gene_chr),6:13])

pdf(file=paste0("Abhd2_pQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),
                 y=pqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-2.2,2.2)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()