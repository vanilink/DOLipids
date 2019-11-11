options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("plotly")
library("pheatmap")
library("qtl2") #has CC colors

setwd("E:/DO/R-Analyses/Summer_2019/") 
load("P_L_Lipids.RData")

#GOAL: make Manhattan/LOD plot for Apoa2 locus (for Fig 2)

my.Apoa2 <- subset(lipid_qtls, qtl.chr=="1" & qtl.pos>(171.2-2) & qtl.pos<(171.2+2))

cluster.my.Apoa2 <- my.Apoa2[,28:35] #reduce to allele effects
rownames(cluster.my.Apoa2) <- my.Apoa2$identifier
colnames(cluster.my.Apoa2) <- c("AJ","B6","129","NOD","NZO","CAST","PWK","WSB")

ano <- data.frame("Class" = factor(my.Apoa2$category), "Tissue" = factor(my.Apoa2$tissue))
rownames(ano) <- my.Apoa2$identifier

Var1 = c(coon_red, coon_blue, coon_yel, coon_turq, "grey", coon_purp)
names(Var1) = unique(my.Apoa2$category)
Var2 = c("lightgrey",coon_grey)
names(Var2) = unique(my.Apoa2$tissue)
ann_colors = list(Class = Var1, Tissue = Var2)


heat <- pheatmap(cluster.my.Apoa2, 
         annotation_row = ano, 
         annotation_colors = ann_colors,
         cluster_cols = F, 
         show_rownames = F, 
         breaks = seq(-1,1,length.out=100),
         cutree_rows = 2, #scale = "row"
         main = "Allele Effect at Chr 1:171.2 +-2 Mbp",
         annotation_legend = F,
         filename = "Apoa2_HM.pdf",
         width = 4,
         height = 4)

table(cutree(heat$tree_row,h = 1.5))
sub.Apoa2 <- cutree(heat$tree_row,h = 1.5)
sub.Apoa2 <- names(sub.Apoa2[which(sub.Apoa2==1)])

#Apoa2_clusters <- read.delim("Apoa2_lipid_clusters.tab")
Apoa2 <- subset(lipid_qtls,identifier %in% sub.Apoa2 & qtl.chr=="1")

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

pdf(file=paste0("Apoa2_Bar.pdf"), pointsize = 6, 
    width = 4, height = 3, useDingbats=FALSE)

ggplot(data=Apoa2, aes(x=subclass, fill=category)) + geom_bar(stat="count") + 
  scale_fill_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, "grey")) +
  coord_flip() +
  scale_y_continuous(name="Count at Apoa2 Locus", expand = c(0.005, 0), breaks = seq(0, 130, 10), position="right") +
  scale_x_discrete(name="") + 
  theme_classic() +
  theme(legend.position= "none")   

dev.off()

pdf(file=paste0("Apoa2_MZRT.pdf"), pointsize = 6, 
    width = 5, height = 3.5, useDingbats=FALSE)
#ggplotly( 
  ggplot() +
    
    geom_point(data=subset(Apoa2,id.status!="identified"), aes(x=rt,y=mz,color="Unidentified",size=qtl.lod,label=identifier,alpha=0.75)) + #making sure the hand-IDd ones are shown as UNK
    geom_point(data=subset(Apoa2,(id.status=="identified"|category=="Fatty.Acyl")), aes(x=rt,y=mz,color=as.factor(category),size=qtl.lod,label=identifier),alpha=0.75) + 
    
    scale_color_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, coon_grey)) +
    scale_alpha(guide = 'none') + #remove unuseful legend
    scale_size_continuous("LOD", range = c(1, 3), breaks = c(min(Apoa2$qtl.lod),10,14)) +
    scale_x_continuous(name="RT [min]") +
    scale_y_continuous(name="m/z") +
    theme_classic()
#)
dev.off()


###############################################################################
#GOAL: make Manhattan/LOD plot for B4galnt1 locus (for Fig 3)

B4galnt1_clusters <- read.delim("B4galnt1_lipid_clusters.tab")
B4galnt1 <- subset(lipid_qtls,identifier %in% B4galnt1_clusters$X==T & qtl.chr == "10")

ggplot() +
  geom_point(data=subset(B4galnt1,id.status!="identified"), aes(x=qtl.pos,y=qtl.lod,color="Unidentified",label=identifier,alpha=0.85,size=1.5)) + #making sure the hand-IDd ones are shown as UNK
  geom_point(data=subset(B4galnt1,id.status=="identified"), aes(x=qtl.pos,y=qtl.lod,color=as.factor(category),label=identifier,alpha=0.85,size=1.5)) + 
  scale_color_manual("Lipid Class",values = c(coon_red, coon_turq, coon_yel, coon_blue, "grey")) + 
  ggtitle("B4galnt1 Locus") +
  scale_x_continuous(name="Chr 1 Pos [Mbp]") +
  scale_y_continuous(name="LOD") + 
  theme(legend.position="none")

#abhd1

Abhd1 <- subset(lipid_qtls, qtl.chr == "5" & qtl.pos>25 & qtl.pos<35)
gg <- ggplot() +
  geom_point(data=subset(Abhd1,id.status!="identified"), aes(label2=tissue,x=qtl.pos,y=qtl.lod,color="Unidentified",label=identifier,alpha=0.85,size=1.5)) + #making sure the hand-IDd ones are shown as UNK
  geom_point(data=subset(Abhd1,id.status=="identified"), aes(label2=tissue,x=qtl.pos,y=qtl.lod,color=as.factor(category),label=identifier,alpha=0.85,size=1.5)) + 
  scale_color_manual("Lipid Class",values = c(coon_purp, coon_turq,"grey")) + 
  ggtitle("Abhd1 Locus") +
  scale_x_continuous(name="Chr 5 Pos [Mbp]") +
  scale_y_continuous(name="LOD") 
ggplotly(gg)


#GOAL: show allele effect of max. LOD feature
#Apoa2
effect <- t(subset(Apoa2, qtl.lod==max(Apoa2$qtl.lod))[28:35])
  
ggplot() +
  geom_point(aes(x=names(CCcolors),y=effect,color=names(CCcolors)),size=10) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-1.1,1.1)) +
  scale_x_discrete(name=NULL) +
  ggtitle(subset(Apoa2, qtl.lod==max(Apoa2$qtl.lod))$identifier) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme_classic() +
  theme(legend.position="none")

#B4galnt1
effect <- t(subset(lipid_qtls, qtl.lod==max(lipid_qtls$qtl.lod,na.rm=T))[28:35])

ggplot() +
  geom_point(aes(x=names(CCcolors),y=effect,color=names(CCcolors)),size=10) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-1.1,2.2)) +
  scale_x_discrete(name=NULL) +
  ggtitle(subset(lipid_qtls, qtl.lod==max(lipid_qtls$qtl.lod,na.rm=T))$identifier) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme(legend.position="none")

#abhd1
effect <- t(subset(lipid_qtls, identifier=="LysoPC_14.0_1.031_468.30853_plus" & qtl.chr=="5",na.rm=T)[28:35])

ggplot() +
  geom_point(aes(x=names(CCcolors),y=effect,color=names(CCcolors)),size=10) +
  scale_color_manual(values=CCcolors) +
  scale_y_continuous(name="FS Allele Effect",limits=c(-1.1,1.1)) +
  scale_x_discrete(name=NULL) +
  ggtitle(subset(lipid_qtls, identifier=="LysoPC_14.0_1.031_468.30853_plus" & qtl.chr=="5",na.rm=T)$identifier) +
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme(legend.position="none")