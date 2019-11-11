options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("cowplot")

setwd("E:/DO/R-Analyses/Feb_2019/") 
load("E:/DO/R-Analyses/Nov_2018/DO_Lipidomics.RData") #has Coon colors

#GOAL: make Manhattan plot (LOD) over whole genome for a selected lipid

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
    
    #OPTIONAL Add highlight and annotation information
    #mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
    mutate(is_highlight=ifelse(LOD>6,"yes","no"))
  #mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 
  
  #2 MAKE X AXIS WITH CHR DISPLAY
  axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  
  #3 MAKE GGPLOT
  plot <- ggplot(don, aes(x=BPcum, y=LOD)) + #should recalc to p-value: y=-log10(P)
    
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
    scale_color_manual(values = rep(c("grey", coon_grey), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = c(1:19,"X"), breaks= axisdf$center ) +
    scale_y_continuous(labels = scaleFUN,
      expand = expand_scale(add=c(0,1)) ) + # remove space between plot area and x axis but add some above
    
    #OPTIONAL Add highlighted points
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
  
  #ggsave(paste0(CE,"_LOD_plot.pdf"), plot, height = 1, width = 3.5, units = "in")
  
}
library("gridExtra")
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 6, ncol = 1, 
             widths = 3.5, heights = rep(1,times=6))
library("egg")
ggarrange(p1, p2, p3, p4, p5, p6, widths = 3.5, heights = rep(1,times=6))
