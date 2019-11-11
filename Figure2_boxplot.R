# founder strain boxplots

fs <- read.csv("E:/FS/FS_Plasma/3-LipiDex-Identification/Results/20170623_FS_Plasma_Log2.csv")

rownames(fs)[1:5] <- fs[1:5,6]

rownames(fs)[6:nrow(fs)] <- paste0(fs[6:nrow(fs),4],"_", fs[6:nrow(fs),1],"_",fs[6:nrow(fs),2],fs[6:nrow(fs),3])

colnames(fs)[7:ncol(fs)] <- fs[2, 7:ncol(fs)]

fs.q <- fs[-2:-5,-1:-6]

loi <- subset(fs.q, rownames(fs.q) == "CE 18:2_17.508_666.61725+" | rownames(fs.q) == "Strain")
rownames(loi)[2] <- "CE_18.2"
loi <- data.frame(t(loi))
loi$CE_18.2 <- as.numeric(loi$CE_18.2)
loi$Strain <- factor(loi$Strain,
                       levels = names(CCcolors))

#boxplot(data = loi, CE_18.2 ~ Strain, col = CCcolors)

ggplot(loi, aes(x=Strain, y=CE_18.2, fill = Strain)) + 
  stat_boxplot(geom ="errorbar") +
  geom_boxplot(outlier.shape = NA) + 
  geom_point() +
  scale_fill_manual(values=CCcolors) +
  ggtitle("Plasma CE 18:2") +
  theme_classic() +
  ylab("Log2 Peak Area") +
  theme(legend.position = "none", 
        axis.title.x = element_blank())
