# GOAL: Create panels for figures 3 and S3

options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("qtl2") #has CC colors

setwd("E:/DO/R-Analyses/GithubDOLipids/")

# LOAD IN LIPID AND LIPID QTL DATA
lipids <- read.csv("lipids.csv")[,-1]
rownames(lipids) <- lipids$identifier

lipid_qtls <- read.csv("lipid_qtls.csv")[,-1]

#remove gangliosides from IDd ones as they were later hand-identified
lipid_qtls[which(lipid_qtls$class=="Ganglioside"),]$category <- "UNK" 
lipids[which(lipids$class=="Ganglioside"),]$category <- "UNK" 

# subset data at locus
my.B4galnt1 <- subset(lipid_qtls, qtl.chr=="10" & qtl.pos>(127.2-2) & qtl.pos<(127.+2))

###############################################################################
# b. Heatmap to cluster allele effects

cluster.my.B4galnt1 <- my.B4galnt1[,28:35] #reduce to allele effects
rownames(cluster.my.B4galnt1) <- my.B4galnt1$identifier
colnames(cluster.my.B4galnt1) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

ano <- data.frame("Class" = factor(my.B4galnt1$category), "Tissue" = factor(my.B4galnt1$tissue))
rownames(ano) <- my.B4galnt1$identifier

Var1 = c("grey", coon_purp)
names(Var1) = unique(my.B4galnt1$category)
Var2 = c("black","white")
names(Var2) = unique(my.B4galnt1$tissue)
ann_colors = list(Class = Var1, Tissue = Var2)

heat <- pheatmap(cluster.my.B4galnt1, 
                 annotation_row = ano, 
                 annotation_colors = ann_colors,
                 cluster_cols = F, 
                 show_rownames = F, 
                 clustering_method = "ward.D2",
                 breaks = seq(-2,2,length.out=100),
                 cutree_rows = 5, 
                 scale = "row",
                 main = "Allele Effect at Chr 10:127.2 +-2 Mbp",
                 annotation_legend = F,
                 filename = "B4galnt1_HM.pdf",
                 width = 2.8,
                 height = 3.5)

dev.off()

table(cutree(heat$tree_row,h = 3))
sub.B4galnt1 <- cutree(heat$tree_row,h = 3)
sub.B4galnt1 <- names(sub.B4galnt1[which(sub.B4galnt1<3)])

B4galnt1 <- subset(lipid_qtls,identifier %in% sub.B4galnt1 & qtl.chr=="10")

###############################################################################
# a. Manhattan plot at B4galnt1 locus

pdf(file=paste0("B4galnt1_Manhattan.pdf"), pointsize = 6, 
    width = 2.8, height = 3.5, useDingbats=FALSE)

ggplot() +
  geom_point(data=B4galnt1, 
             aes(x=qtl.pos,y=qtl.lod,fill=coon_grey, shape = (mz>1400 | mz<800),
                 alpha=0.85)) + 
  scale_x_continuous(name="Chr 10 Pos [Mbp]", limits = c(125.2, 129.2)) +
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position= "none") 

dev.off()

###############################################################################
# c. M/z plot of features at B4galnt1 locus

pdf(file=paste0("B4galnt1_MZRT.pdf"), pointsize = 6, 
    width = 2.5, height = 3.5, useDingbats=FALSE)

ggplot() +
  geom_point(data=B4galnt1, aes(x=rt,y=mz,fill=coon_grey, shape = (mz>1400 | mz<800),
                                size=qtl.lod),alpha=0.75) + 
  scale_alpha(guide = 'none') + #remove unuseful legend
  scale_shape(guide = 'none') + #remove unuseful legend
  scale_fill_discrete(guide = 'none') + #remove unuseful legend
  scale_size_continuous("LOD", range = c(1, 4), breaks = c(min(B4galnt1$qtl.lod),25,80)) +
  scale_x_continuous(name="RT [min]", limits = c(0,19)) +
  scale_y_continuous(name="m/z", limits = c(200,1600)) +
  theme_classic() +
  theme(legend.position = c(0.15,0.6))

dev.off()

###############################################################################
# f. Manhattan plot of one feature across mouse genome
Ganglioside <- sub.B4galnt1[1]

LOD_csv <- read.csv(paste0("E:/DO/R-Analyses/2018-10-11_DO_QTL/Plasma_Lipids/Plasma_Lipids_",Ganglioside,"_LOD_plot.csv"))
#SEPARATE X to CHR and BP
LOD_csv$CHR <- sapply(strsplit(LOD_csv$X,"_"), `[`, 1)
LOD_csv[LOD_csv$CHR == "X",]$CHR <- 20 #making sure X chromosome is numeric
LOD_csv$CHR <- as.numeric(LOD_csv$CHR)
LOD_csv$BP <- as.numeric(sapply(strsplit(LOD_csv$X,"_"), `[`, 2)) #in bp (not Mbp!)
Ganglioside <- colnames(LOD_csv)[2]
colnames(LOD_csv)[2] <- "LOD"

#FROM https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
#1 CALCULATE CUMULATIVE SNP POSITION
don <- LOD_csv %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(LOD_csv, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight 
  mutate(is_highlight=ifelse(LOD>6,"yes","no"))

#2 MAKE X AXIS WITH CHR DISPLAY
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#3 MAKE GGPLOT
plot <- ggplot(don, aes(x=BPcum, y=LOD)) + 
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("grey", coon_grey), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = c(1:19,"X"), breaks= axisdf$center ) +
  scale_y_continuous(#labels = scaleFUN,
                     expand = expand_scale(add=c(0,1)) ) + # remove space between plot area and x axis but add some above
  
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), color=coon_yel, size=1) +
  
  ggtitle(Ganglioside) + 
  
  xlab("Mouse Chromosome") +
  
  # Custom the theme:
  theme_classic() +
  theme(text = element_text(size=10),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

ggsave(paste0(Ganglioside,"_LOD_plot.pdf"), plot, height = 2.5, width = 2.8, units = "in")

###############################################################################
# S3a. B4galnt1 liver eqtl allele effect

load("Svenson_DO850_for_eQTL_viewer_v4.Rdata")
eqtl <- dataset.mrna$lod.peaks$additive
eqtl_effect <- unlist(eqtl[which(eqtl$annot.id=="ENSMUSG00000006731"&eqtl$cis),11:18])

pdf(file=paste0("B4galnt1_eQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),
                 y=eqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-3,3)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()

###############################################################################
# S3b. B4galnt1 liver pqtl allele effect
pqtl <- dataset.protein$lod.peaks$additive
pqtl_effect <- unlist(pqtl[which(pqtl$annot.id=="ENSMUSP00000006914"&pqtl$qtl_chr==pqtl$gene_chr),6:13])

pdf(file=paste0("B4galnt1_pQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),
                 y=pqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-3,3)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()
###############################################################################
# S3c. FS boxplot of max. LOD feature

fs <- read.csv("20170623_FS_Plasma_Log2.csv")

rownames(fs)[1:5] <- fs[1:5,6]

rownames(fs)[6:nrow(fs)] <- paste0(fs[6:nrow(fs),4],"_", #roughly combine feature IDs
                                   fs[6:nrow(fs),1],"_",
                                   fs[6:nrow(fs),2],
                                   fs[6:nrow(fs),3])

colnames(fs)[7:ncol(fs)] <- fs[2, 7:ncol(fs)]

fs.q <- fs[-2:-5,-1:-6] #quant data

row <- which(grepl("1253.", rownames(fs.q)))
loi <- fs.q[c(1,row),1:62] # Lipid of Interest
rownames(loi)[2] <- "analyte"
loi <- data.frame(t(loi))

loi$analyte <- as.numeric(loi$analyte)
loi$Strain <- factor(loi$Strain, levels = names(CCcolors))

pdf(file=paste0("Ganglioside_boxplot.pdf"), pointsize = 6, 
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
