# GOAL: Create panels for figures 2 and S2

options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("pheatmap")
library("qtl2") #has CC colors

setwd("E:/DO/R-Analyses/GithubDOLipids/")

# define colors
coon_blue <- "#2CA7DF"
coon_grey <- "#566977"
coon_purp <- "#955CA5"
coon_red <- "#EC6B63"
coon_turq <- "#63C29C"
coon_yel <- "#FFCB04"

# LOAD IN LIPID AND LIPID QTL DATA
lipids <- read.csv("TableS8.csv")
rownames(lipids) <- lipids$identifier

lipid_qtls <- read.csv("TableS9.csv")
lipid_qtls <- merge(lipid_qtls, lipids)

# subset data at locus
my.Apoa2 <- subset(lipid_qtls, qtl.chr=="1" & qtl.pos>(171.2-2) & qtl.pos<(171.2+2))
###############################################################################
# S2a. Heatmap to cluster allele effects

start <- which(colnames(my.Apoa2)=="A")
end <- which(colnames(my.Apoa2)=="H")

cluster.my.Apoa2 <- my.Apoa2[,start:end] #reduce to allele effects
rownames(cluster.my.Apoa2) <- my.Apoa2$identifier
colnames(cluster.my.Apoa2) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

ano <- data.frame("Class" = factor(my.Apoa2$category), "Tissue" = factor(my.Apoa2$tissue))
rownames(ano) <- my.Apoa2$identifier

Var1 = c(coon_red, coon_blue, coon_yel, coon_turq, "grey", coon_purp)
names(Var1) = unique(my.Apoa2$category)
Var2 = c("black","white")
names(Var2) = unique(my.Apoa2$tissue)
ann_colors = list(Class = Var1, Tissue = Var2)

heat <- pheatmap(cluster.my.Apoa2, 
         annotation_row = ano, 
         annotation_colors = ann_colors,
         cluster_cols = F, 
         show_rownames = F, 
         breaks = seq(-1,1,length.out=100),
         cutree_rows = 2, #determined below after running pheatmap once 
         main = "Allele Effect at Chr 1:171.2 +-2 Mbp",
         annotation_legend = F,
         filename = "Apoa2_HM.pdf",
         width = 4,
         height = 4)

table(cutree(heat$tree_row,h = 1.5))
sub.Apoa2 <- cutree(heat$tree_row,h = 1.5)
sub.Apoa2 <- names(sub.Apoa2[which(sub.Apoa2==1)])

# for further analyses work with reduced Apoa2 with matching allele effects
Apoa2 <- subset(lipid_qtls,identifier %in% sub.Apoa2 & qtl.chr=="1")

###############################################################################
# 2a. Manhattan plot of Apoa2 locus

pdf(file=paste0("Apoa2_Manhattan.pdf"), pointsize = 6, 
    width = 4, height = 3.25, useDingbats=FALSE)

ggplot() +
  geom_point(data=subset(Apoa2,id.status!="identified"), 
             aes(x=qtl.pos,y=qtl.lod,color="Unidentified",
                 alpha=0.85)) + #making sure the hand-IDd ones are shown as UNK
  geom_point(data=subset(Apoa2,id.status=="identified"), 
             aes(x=qtl.pos,y=qtl.lod,color=as.factor(category),
                 alpha=0.85)) + 
  scale_color_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, "grey")) + 
  scale_x_continuous(name="Chr 1 Pos [Mbp]") +
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position= "none") 

dev.off()

###############################################################################
# 2c. Bar plot of lipid classes at Apoa2 locus

pdf(file=paste0("Apoa2_Bar.pdf"), pointsize = 6, 
    width = 4, height = 3, useDingbats=FALSE)

ggplot(data=Apoa2, aes(x=subclass, fill=category)) + geom_bar(stat="count") + 
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, "grey")) +
  coord_flip() +
  scale_y_continuous(name="Count at Apoa2 Locus", expand = c(0.005, 0), 
                     breaks = seq(0, 130, 10), position="right") +
  scale_x_discrete(name="") + 
  theme_classic() +
  theme(legend.position= "none")   

dev.off()

###############################################################################
# 2d. m/z vs. RT plot of features mapping to Apoa2 locus

pdf(file=paste0("Apoa2_MZRT.pdf"), pointsize = 6, 
    width = 5, height = 3.5, useDingbats=FALSE)

ggplot() +
    
    geom_point(data=subset(Apoa2,id.status!="identified"), 
               aes(x=rt,y=mz,color="Unidentified",size=qtl.lod,alpha=0.75)) + #making sure the hand-IDd ones are shown as UNK
    geom_point(data=subset(Apoa2,(id.status=="identified"|category=="Fatty.Acyl")), 
               aes(x=rt,y=mz,color=as.factor(category),size=qtl.lod),alpha=0.75) + 
    
    scale_color_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, coon_grey)) +
    scale_alpha(guide = 'none') + #remove unuseful legend
    scale_size_continuous("LOD", range = c(1, 3), breaks = c(min(Apoa2$qtl.lod),10,14)) +
    scale_x_continuous(name="RT [min]") +
    scale_y_continuous(name="m/z") +
    theme_classic()
#)
dev.off()

###############################################################################
# 2e. Manhattan plots of CE features over all chromosomes
CEs <- subset(lipids, class == "CE" & tissue=="plasma")

#needed for y-axes to be same width
scaleFUN <- function(x) sprintf("%.0f", x)

for (i in 1:6) {
  CE <- CEs$identifier[i]
  LOD_csv <- read.csv(paste0("E:/DO/R-Analyses/2018-10-11_DO_QTL/Plasma_Lipids/Plasma_Lipids_",CE,"_LOD_plot.csv"))
  #SEPARATE X to CHR and BP
  LOD_csv$CHR <- sapply(strsplit(LOD_csv$X,"_"), `[`, 1)
  LOD_csv[LOD_csv$CHR == "X",]$CHR <- 20 #making sure X chromosome is numeric
  LOD_csv$CHR <- as.numeric(LOD_csv$CHR)
  LOD_csv$BP <- as.numeric(sapply(strsplit(LOD_csv$X,"_"), `[`, 2)) #in bp (not Mbp!)
  CE <- colnames(LOD_csv)[2]
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
    scale_y_continuous(labels = scaleFUN,
                       expand = expand_scale(add=c(0,1)) ) + 
    # remove space between plot area and x axis but add some above
    
    # Add highlighted points
    geom_point(data=subset(don, is_highlight=="yes"), color=coon_blue, size=1) +
    
    ggtitle(CE) + 
    
    xlab("Mouse Chromosome") +
    
    # Custom the theme:
    theme_classic() +
    theme(text = element_text(size=10),
          legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank())
  
  if (i < 6) {
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.y=element_blank())
  }
  
  
  if (i == 1) {p1 <- plot}
  else if (i == 2) {p2 <- plot}
  else if (i == 3) {p3 <- plot}
  else if (i == 4) {p4 <- plot}
  else if (i == 5) {p5 <- plot}
  else if (i == 6) {p6 <- plot}
  
  ggsave(paste0(CE,"_LOD_plot.pdf"), plot, height = 1, width = 3.5, units = "in")
  
}

###############################################################################
# S2b. FS boxplot of max. LOD feature

fs <- read.csv("20170623_FS_Plasma_Log2.csv")

rownames(fs)[1:5] <- fs[1:5,6]

rownames(fs)[6:nrow(fs)] <- paste0(fs[6:nrow(fs),4],"_", #roughly combine feature IDs
                                   fs[6:nrow(fs),1],"_",
                                   fs[6:nrow(fs),2],
                                   fs[6:nrow(fs),3])

colnames(fs)[7:ncol(fs)] <- fs[2, 7:ncol(fs)]

fs.q <- fs[-2:-5,-1:-6] #quant data

row <- which(grepl("CE 18:2", rownames(fs.q)))
loi <- fs.q[c(1,row+1),1:62] # Lipid of Interest
rownames(loi)[2] <- "analyte"
loi <- data.frame(t(loi))

loi$analyte <- as.numeric(loi$analyte)
loi$Strain <- factor(loi$Strain, levels = names(CCcolors))

pdf(file=paste0("CE182_boxplot.pdf"), pointsize = 6, 
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
# S2c. Apoa2 liver eqtl, pqtl allele effect

load("Svenson_DO850_for_eQTL_viewer_v4.Rdata")
eqtl <- dataset.mrna$lod.peaks$additive
eqtl_effect <- unlist(eqtl[which(eqtl$annot.id=="ENSMUSG00000005681"&eqtl$cis),11:18])

pdf(file=paste0("Apoa2_eQTL.pdf"), pointsize = 6, 
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

pqtl <- dataset.protein$lod.peaks$additive
pqtl_effect <- unlist(pqtl[which(pqtl$annot.id=="ENSMUSP00000106953"&pqtl$qtl_chr==pqtl$gene_chr),6:13])

pdf(file=paste0("Apoa2_pQTL.pdf"), pointsize = 6, 
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

###############################################################################
# 2f. RKMD of CE features (Referenced Kendrick Mass Defect)

data <- cbind(Apoa2, "KM" = Apoa2$mz*(14/14.01565)) #ADD A COL FOR THE KENDRICK MASS
data <- cbind(data, "KMD" = data$KM - floor(data$KM)) #ADD A COL FOR THE KENDRICK MASS DEFECT
data <- cbind(data, "KNM" = round(data$KM,0)) #ADD A COL FOR THE KENDRICK NOMINAL MASS

this.class = "CE"
Ref_CE <- (data[which(data$class=="CE"),]$KMD)[1]+0.0134
data <- cbind(data, "RKMD_CE" = (data$KMD-Ref_CE)/0.0134) #Reference to [CE+NH4]+

#select features within 3 SD of CE RT
rt_max = mean(data[which(data$class==this.class),]$rt) +
  3*sd(data[which(data$class==this.class),]$rt)
rt_min = mean(data[which(data$class==this.class),]$rt) -
  3*sd(data[which(data$class==this.class),]$rt)
#select features within 5 SD of CE m/z
mz_max = mean(data[which(data$class==this.class),]$mz) +
  5*sd(data[which(data$class==this.class),]$mz)
mz_min = mean(data[which(data$class==this.class),]$mz) -
  5*sd(data[which(data$class==this.class),]$mz)
RKMD_cutoff = 0.075

data.sub <- data[which((data$class==this.class|data$class=="UNK") &
                         data$rt > rt_min &
                         data$mz > mz_min &
                         data$mz < mz_max &
                         abs(data$RKMD_CE-round(data$RKMD_CE,0))<RKMD_cutoff
),]

ggplot(data.sub[order(data.sub$id.status,decreasing = T),]) + #draw unk's first
  geom_point(aes(x=KNM,y=RKMD_CE,color=id.status)) +
  scale_color_manual("Identified", values=c(coon_blue,"darkblue")) +
  theme(panel.grid.major.y = element_line(color = "grey80")) +
  scale_y_continuous(breaks = seq(-6, 0, 1),limits = c(-6.5,0)) +
  theme(legend.position="none")

