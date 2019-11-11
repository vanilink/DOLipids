options(stringsAsFactors = F)
library("ggplot2")
library("qtl2")
library("plotly")
library("tidyr")
library("reshape2")
library("cowplot")
library(dplyr)
library(RColorBrewer)
library(pheatmap)

setwd("E:/DO/R-Analyses/Summer_2019/") 
load("P_L_Lipids.RData")

my.abhd1 <- subset(lipid_qtls, qtl.chr=="5" & qtl.pos>(31-2) & qtl.pos<(31+2))

ano <- data.frame("Class" = factor(my.abhd1$category), 
                  "Tissue" = factor(my.abhd1$tissue),
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

ann_colors = list(Class = Var1, Tissue = Var2, FA = Var3, LysoPC = Var4)

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
                 main = "Allele Effect at Chr 5:31.0 +-2 Mbp",
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

# PLOT i: Manhattan 
abhd1$col  <- with(abhd1, ifelse(identifier %in% abhd1[grep(pattern = "14.0",chain),]$identifier,"14.0",ifelse(class=="UNK","unidentified","identified")))

pdf(file=paste0("Abhd1_Manhattan.pdf"), pointsize = 6, 
    width = 3, height = 3, useDingbats=FALSE)

#ggplotly(
ggplot(data=abhd1, aes(x=qtl.pos,y=qtl.lod,
                           color=col,
                           label=identifier,alpha=0.85,size=1.5)) +
  geom_point() + 
  scale_x_continuous(limits = c(29,33), name="Chr 5 Pos [Mbp]") +
  scale_color_manual("FA 14:0", values=c(coon_turq, coon_grey, "grey")) +
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position="none")
#)
dev.off()

###############################################################################
#eQTL, pQTL allele effect plot
load("E:/DO/R-Analyses/Attie_QTL_Viewer.RData")
eqtl <- dataset.islet.rnaseq$lod.peaks
eqtl_effect <- unlist(eqtl[which(eqtl$annot.id=="ENSMUSG00000006638"),11:18])

pqtl <- dataset.islet.proteins$lod.peaks
pqtl_effect <- unlist(pqtl[which(pqtl$annot.id=="ENSMUSG00000006638"),11:18])

pdf(file=paste0("Abhd1_eQTL.pdf"), pointsize = 6, 
    width = 4, height = 2, useDingbats=FALSE)

ggplot() +
  geom_point(aes(x=factor(names(CCcolors),levels=names(CCcolors)),y=eqtl_effect,color=names(CCcolors)),size=6) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-2.2,2.2)) + 
  scale_x_discrete(name=NULL) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

dev.off()

###############################################################################
#manhattan plot of subset of data across all chromosomes
#LysoPCs
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

#14:0 containing lipids
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

####
lysopc <- subset(lipid_qtls, class=="LysoPC" & !is.na(qtl.chr))
topdata <- data.frame("PHE"=lysopc$identifier,"CHR"=lysopc$qtl.chr,
                      "MBP"=lysopc$qtl.pos,"LOD"=lysopc$qtl.lod)
  
FA140 <- subset(lipid_qtls[grep(pattern = "14.0",lipid_qtls$chain),],!is.na(qtl.chr))
topdata <- data.frame("PHE"=FA140$identifier,"CHR"=FA140$qtl.chr,
                         "MBP"=FA140$qtl.pos,"LOD"=FA140$qtl.lod)


#add dummy qtls so that each chr is present in dataset for plotting
#get lengths of mouse chr:
karyo <- read.delim("C:/circos-0.69-6/data/karyotype/karyotype.mouseX.txt", sep = " ", header = F)
karyo <- karyo[1:20,c(4,6)]

for (chr in c(1:20)) {
  newrow.1 <- data.frame("PHE" = "dummy",
                       "CHR" = chr,
                       "MBP" = 1,
                       "LOD" = 5.99)
  
  newrow.2 <- data.frame("PHE" = "dummy",
                        "CHR" = chr,
                        "MBP" = karyo[chr,2]/1e6,
                        "LOD" = 5.99)
  topdata <- rbind(topdata, newrow.1, newrow.2)
}



#SEPARATE X to CHR and BP
topdata[topdata$CHR == "X",]$CHR <- 20 #making sure X chromosome is numeric
topdata$CHR <- as.numeric(topdata$CHR)
topdata$BP <- as.numeric(topdata$MBP*1e6) #in bp (not Mbp!)

#FROM https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
#1 CALCULATE CUMULATIVE SNP POSITION
don <- topdata %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(topdata, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  #OPTIONAL Add highlight and annotation information
  mutate(is_highlight=ifelse(LOD>6,"yes","no"))

#2 MAKE X AXIS WITH CHR DISPLAY
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

#needed for y-axes to be same width
scaleFUN <- function(x) sprintf("%.0f", x)

#3 MAKE GGPLOT
plot <- ggplot(don, aes(x=BPcum, y=LOD)) + 
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("grey", coon_grey), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = c(1:19,"X"), breaks= axisdf$center ) +
  scale_y_continuous(labels = scaleFUN,
                     expand = expand_scale(add=c(0,1)) ) + # remove space between plot area and x axis but add some above
  
  #OPTIONAL Add highlighted points
  geom_point(data=subset(don, CHR == 5), color=coon_turq, size=1) +
  
  ggtitle("FA 14:0") + 
  
  xlab("Mouse Chromosome") +
  
  # Custom the theme:
  #theme_classic() +
  theme(text = element_text(size=10),
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_blank())

pdf(file=paste0("FA140_Manhattan.pdf"), pointsize = 5, 
    width = 3, height = 1.5, useDingbats=FALSE)
plot
dev.off()
#ggplotly(plot)

