options(stringsAsFactors = F)
library("ggplot2")
library("dplyr")
library("pheatmap")
library("RColorBrewer")
library("tidyr")

setwd("E:/DO/R-Analyses/GithubDOLipids/")

# define colors
coon_blue <- "#2CA7DF"
coon_grey <- "#566977"
coon_purp <- "#955CA5"
coon_red <- "#EC6B63"
coon_turq <- "#63C29C"
coon_yel <- "#FFCB04"

norm <- read.csv("E:/DO/NatureMetabolism/SD_Fig5.csv")

rownames(norm) <- norm$ID

norm <- norm[,grep("A1|A3|GFP", norm[1,])] #subset to Abhd1, Abhd3, and GFP

norm <- data.frame(t(norm[1:2,]), t(apply(norm[-1:-2,], c(1,2), as.numeric)))

#########

norm.long <- gather(norm, -Mutant, -Batch, key = "lipid", value="log2int")

GFP.long <- subset(norm.long, Mutant == "GFP")
mutant.long <- subset(norm.long, Mutant != "GFP")

GFP.avg <- GFP.long %>%
  group_by(lipid, Mutant, Batch) %>%
  summarise(
    n=n(),
    GFP.mean=mean(log2int), 
    GFP.sd=sd(log2int) 
    ) %>%
  mutate( GFP.se=GFP.sd/sqrt(n))  %>%
  mutate( GFP.ic=GFP.se * qt((1-0.05)/2 + .5, n-1))

norm.comb <- merge(mutant.long, GFP.avg, by= c("lipid", "Batch"))
norm.comb <- norm.comb[, !(colnames(norm.comb) %in% c("Mutant.y", "n") )]
colnames(norm.comb)[which(colnames(norm.comb)=="Mutant.x")] <- "Mutant"

Mutant.FC <- norm.comb %>%
  group_by(lipid, Mutant, Batch, GFP.mean, GFP.sd, GFP.se, GFP.ic) %>%
  mutate(
    FC = log2int - GFP.mean)

Mutant.avg <- Mutant.FC %>%
  group_by(lipid, Mutant) %>%
  summarise(
    n=n(),
    mean.FC=mean(FC),
    sd=sd(FC) 
  ) %>%
  mutate( se=sd/sqrt(n)) %>%
  mutate( ic=se * qt((1-0.05)/2 + 0.5, n-1))

#boxplot per lipidclass
Mutant.avg$class <- sapply(strsplit(Mutant.avg$lipid, "[.]"),'[[',1)

classes <- data.frame(table(Mutant.avg$class))[which(data.frame(table(Mutant.avg$class))$Freq>10),]$Var1 #more than 10 members of lipid class

p1 <- ggplot(subset(Mutant.avg,class %in% classes), 
       aes(x=class, y=abs(mean.FC), ymax=max(abs(mean.FC)), ymin=min(abs(mean.FC)), fill=Mutant)) + 
  geom_boxplot(position = 'dodge2', alpha=0.9) + 
  scale_fill_manual(values=c(coon_turq, coon_grey), labels = c("Abhd1", "Abhd3")) +
  ylab("Average Absolute Fold Change Mutant/GFP") +
  xlab("Lipid Classes with n > 10") +
  geom_hline(yintercept = 0.4, linetype = "dashed") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        legend.position = "none")

#ggsave("classes_FC.pdf", width = 8, height = 2.5, units = "in")

lipidclass <- "LysoPC"
forplot <- Mutant.avg[which(substring(Mutant.avg$lipid, 1, nchar(lipidclass)) == lipidclass),]
forplot$shortname <- sapply(strsplit(forplot$lipid, "_"),'[[',1) #for LysoPC

p2 <- ggplot(forplot, aes(x=lipid, y=mean.FC, ymin=mean.FC-ic, ymax=mean.FC+ic, fill=Mutant)) + 
  geom_bar(position = 'dodge', stat='identity', alpha=0.9) +
  scale_fill_manual(values=c(coon_turq, coon_grey), labels = c("Abhd1", "Abhd3")) +
  geom_errorbar(position = 'dodge', colour=coon_grey, alpha=0.5, size=0.7) +
  ylab("Average Fold Change Mutant/GFP") +
  scale_x_discrete(labels=subset(forplot,Mutant=="A1")$shortname) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.x = element_blank(),
        legend.position = "none")

#ggsave("LysoPC_FC.pdf", width = 4, height = 4, units = "in")

lipidclass <- "PC"
forplot <- Mutant.avg[which(substring(Mutant.avg$lipid, 1, nchar(lipidclass)) == lipidclass),]
forplot$shortname <- substring(forplot$lipid,first = 1, last=12) #for PC

p3 <- ggplot(forplot[1:38,], aes(x=lipid, y=mean.FC, ymin=mean.FC-ic, ymax=mean.FC+ic, fill=Mutant)) + 
  geom_bar(position = 'dodge', stat='identity', alpha=0.9) +
  scale_fill_manual(values=c(coon_turq, coon_grey), labels = c("Abhd1", "Abhd3")) +
  geom_errorbar(position = 'dodge', colour=coon_grey, alpha=0.5, size=0.7) +
  ylab("Average Fold Change Mutant/GFP") +
  scale_x_discrete(labels=subset(forplot,Mutant=="A1")$shortname) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5),
        axis.title.x = element_blank(),
        legend.position = c(.1, .9))

#ggsave("PC_FC.pdf", width = 4, height = 4, units = "in")
###############################################################################
library('gridExtra')
grid.arrange(p1, p2, p3, nrow = 1)
