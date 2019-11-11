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
#load("Figure4.RData")

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

# PLOT INITIAL VOLCANO PLOT
ggplotly(ggplot(
  volcano, aes(x=foldchange,y=-1*log10(rawpvalue), #use neg. log of pvalue for plotting
               label=B6.ID, color=B6.UNK, size=sig, alpha=0.5)) + 
  geom_point()+
  scale_color_manual("Un-Identified", values=c(coon_grey, "grey")) +
  scale_size_manual("Significant", values = c(1,3)) +
  scale_x_continuous(limits = c(-5,5), name="Fold Change F/M") +
  scale_y_continuous(limits = c(0,6), name="-log10(p-value)")
  )

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

#volcano_new <- volcano[which(volcano$inDO==TRUE),]
volcano_new <- inner_join(x = volcano, y=lipids, by=c("DOmz"="mz"))

#LoI_new <- volcano_new[which(volcano_new$sig==TRUE),]
#QoI <- lipid_qtls[which(lipid_qtls$identifier%in%LoI_new$identifier),c(1,25:27)]

#only call UNK if both in B6 and in DO unidentified:
volcano_new$UNK <- ifelse(volcano_new$B6.UNK==F,F,ifelse(volcano_new$id.status!="un-identified",F,T))

#add QTL info for chr6:91 +-2 Mbp and use ident to color accordingly
chr6.91Mbp <- subset(lipid_qtls[which(lipid_qtls$mz%in%volcano_new$DOmz),],qtl.chr==6&qtl.pos>89&qtl.pos<93)
volcano_new$ident  <- with(volcano_new, ifelse(identifier %in% chr6.91Mbp$identifier, "?", ifelse(UNK,"unidentified","identified")))

write.csv(volcano_new, "SupplementalTableS6.csv")

# PLOT B: Volcano 
# add highlights based on which ones are at Chr 6:90 Mbp

ggplotly(ggplot(volcano_new, aes(x=foldchange,y=-1*log10(rawpvalue), lable=B6.ID, lable2 = identifier,
                                 color=as.factor(ident), alpha = 0.8,
                                 size=sig)) + 
           geom_point() +
           scale_color_manual("Identification", values=c(coon_red,coon_grey, "grey")) +
           scale_x_continuous(limits = c(-5,5), name="Fold Change F/M") +
           scale_y_continuous(limits = c(0,6), name="-log10(p-value)") +
           scale_size_manual("Significant", values = c(1,3)) +
           theme_classic() + 
           theme(legend.position="none")
)

# PLOT C: Manhattan 
ggplotly(ggplot(data=chr6.91Mbp, aes(x=qtl.pos,y=qtl.lod,
                                     color=coon_red,
                                     label=identifier,alpha=0.85,size=1.5)) +
  geom_point() + 
  scale_x_continuous(limits = c(89,93), name="Chr 6 Pos [Mbp]") +
  scale_y_continuous(name="LOD") + 
  theme_classic() +
  theme(legend.position="none")
    )

# PLOT D: MZ-RT #need to add highlights
ggplotly( 
  ggplot(data=chr6.91Mbp) +
    geom_point(aes(x=rt,y=mz,color=as.factor(ident),label=identifier),alpha=0.85) + 
    scale_color_manual("Identification",values = c(coon_grey, "grey", coon_red)) +
    scale_x_continuous(name="RT [min]") +
    scale_y_continuous(name="m/z") +
    xlim(0, 21) +
    ylim(200,1600) +
    theme_classic() +
    theme(legend.position="none")
)

# Different approach, start from QTLs, then add info about whether or not in B6 dataset

my.chr6.91Mbp <- subset(lipid_qtls, qtl.chr=="6" & qtl.pos>(91-2) & qtl.pos<(91+2))

ano <- data.frame("Class" = factor(my.chr6.91Mbp$category), 
                  "Tissue" = factor(my.chr6.91Mbp$tissue),
                  "inB6" = factor(my.chr6.91Mbp$identifier %in% volcano_new$identifier)
                    )
rownames(ano) <- my.chr6.91Mbp$identifier

Var1 = c(coon_yel, coon_turq, coon_purp, "grey")
names(Var1) = unique(my.chr6.91Mbp$category)
Var2 = c("white","black")
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
                 show_rownames = T, 
                 clustering_method = "ward.D2",
                 breaks = seq(-2,2,length.out=100),
                 cutree_rows = 3, 
                 scale = "row",
                 main = "Allele Effect at Chr 6:91 +-2 Mbp",
                 annotation_legend = F,
                 filename = "Chr6_HM.pdf",
                 width = 4,
                 height = 4
                 )
dev.off()

table(cutree(heat$tree_row,h = 4)) 
sub.chr6.91Mbp <- cutree(heat$tree_row,h = 4)
sub.chr6.91Mbp <- names(sub.chr6.91Mbp[which(sub.chr6.91Mbp==3)])

chr691Mbp <- subset(lipid_qtls,identifier %in% sub.chr6.91Mbp & qtl.chr=="6")
chr691Mbp$inB6 <- chr691Mbp$identifier %in% volcano_new$identifier

# PLOT B: Volcano 
# add highlights based on which ones are at Chr 6:91 Mbp and whether or not id'd
volcano_new$col  <- with(volcano_new, ifelse(identifier %in% chr691Mbp$identifier, "?", ifelse(UNK,"unidentified","identified")))
volcano$col  <- with(volcano, ifelse(B6.ID %in% subset(volcano_new,col=="?"&sig==T)$B6.ID, "?", ifelse(B6.UNK,"unidentified","identified")))

pdf(file=paste0("Chr6_Volcano.pdf"), pointsize = 6, 
    width = 2, height = 3, useDingbats=FALSE)

#ggplotly(
  ggplot(volcano, aes(x=foldchange,y=-1*log10(rawpvalue), 
                                 lable=B6.ID,
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
#)
dev.off()

# PLOT C: Manhattan 
chr691Mbp$col  <- with(chr691Mbp, ifelse(identifier %in% volcano_new$identifier, "?", ifelse(class=="UNK","unidentified","identified")))

pdf(file=paste0("Chr6_Manhattan.pdf"), pointsize = 6, 
    width = 3, height = 3, useDingbats=FALSE)

#ggplotly(
  ggplot(data=chr691Mbp, aes(x=qtl.pos,y=qtl.lod,
                                     color=col,
                                     label=identifier,alpha=0.85,size=1.5)) +
           geom_point() + 
           scale_x_continuous(limits = c(89,93), name="Chr 6 Pos [Mbp]") +
           scale_color_manual("In B6", values=c(coon_red, coon_grey, "grey")) +
           scale_y_continuous(name="LOD") + 
           theme_classic() +
           theme(legend.position="none")
#)
dev.off()

# PLOT D: MZ-RT 
pdf(file=paste0("Chr6_MZRT2.pdf"), pointsize = 6, 
    width = 2, height = 1.5, useDingbats=FALSE)

#ggplotly( 
  ggplot(data=chr691Mbp) +
    geom_point(aes(x=rt,y=mz,color=as.factor(col), #size=sig, 
                   label=identifier),
               alpha=0.85) + 
    scale_color_manual("Identification",values = c(coon_red, coon_grey, "grey")) +
    #scale_size_manual("Significant", values = c(0.5,1.5)) +
    scale_x_continuous(limits=c(12,17), name="RT [min]") +
    scale_y_continuous(limits=c(1000,1200), name="m/z") +
    theme_classic() +
    theme(legend.position="none")
#)
dev.off()

pdf(file=paste0("Chr6_MZRT1.pdf"), pointsize = 6, 
    width = 2, height = 1.5, useDingbats=FALSE)

#ggplotly( 
  ggplot(data=chr691Mbp) +
    geom_point(aes(x=rt,y=mz,color=as.factor(col), 
                   label=identifier),
               alpha=0.85) + 
    scale_color_manual("Identification",values = c(coon_red, coon_grey, "grey")) +
    scale_x_continuous(limits=c(0,20), name="RT [min]") +
    scale_y_continuous(limits=c(200,1600), name="m/z") +
    theme_classic() +
    theme(legend.position="none")
#)
dev.off()

length(which((subset(volcano_new, sig==T & inDO==T)$DOmz)%in%(subset(lipid_qtls,qtl.lod>6)$mz)))
#127 sex-significant B6 features are found in DO and show QTL
length(which((subset(volcano_new, sig==T & inDO==T & B6.UNK == T)$DOmz)%in%(subset(lipid_qtls,qtl.lod>6)$mz)))
#79 of these were unidentified in the B6 dataset
