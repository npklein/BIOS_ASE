library(ggplot2)
library(data.table)
library(reshape2)
library(ggforce)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(gridExtra)

data <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/concordanceGTEx/counts.matrix.inclPatrickASE.txt'))
snpInfo <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt')

#Merge SNP info into table
data_merged <- merge(data, snpInfo, by=c('VARIANT'))

#Select tissue whole blood
data1<-data_merged[data_merged$TISSUE == "WHLBLD",]

#Apply binominal test
data1$binom <- apply(data1, 1, function(x) binom.test(as.numeric(x[["SUMMAJOR"]]), as.numeric(x[["TOTAL"]]), 
                                                      p=0.5, alternative = "two.sided", conf.level = 0.95)$p.value)

data1$binomCorrected <- p.adjust(data1$binom, method = "bonferroni")
data1$ZSCORE_ASE <- -1*qnorm(data1$binomCorrected/2)
data2<-data1[data1$binomCorrected < 1e-5,]

#Check positive or negative ratios and annotate with them
data2$OURDIR <- ifelse(data2$RATIO<0.5, "NEG", "POS")
data2$GTEXDIR <- ifelse(data2$GTEXRATIO<0.5, "NEG", "POS")

#Calculate concordance
comp <- table(data2$OURDIR==data2$GTEXDIR)
concordance <- (comp["TRUE"]/(comp["TRUE"]+comp["FALSE"]))*100
concordance

#Calculate correlation
res <- cor.test(data2$RATIO, data2$GTEXRATIO, 
                method = "pearson")
correlation<-res$estimate


#ggplot(data=data2, aes(x=GTEXRATIO, y=RATIO, colour=ZSCORE_ASE))+
#  geom_point()+
#  ggtitle(paste0("GENE ratio vs GTEx ratio\nConcordance: ", concordance, "    Correlation: ", correlation))+
#  geom_vline(xintercept=0.5)+
#  geom_hline(yintercept=0.5)


######
data2$OURDIR <- ifelse(data2$PATRATIO<0.5, "NEG", "POS")
data2$GTEXDIR <- ifelse(data2$GTEXRATIO<0.5, "NEG", "POS")

#Calculate concordance
comp <- table(data2$OURDIR==data2$GTEXDIR)
concordance <- (comp["TRUE"]/(comp["TRUE"]+comp["FALSE"]))*100
concordance

#Calculate correlation
res <- cor.test(data2$PATRATIO, data2$GTEXRATIO, 
                method = "pearson")
correlation<-res$estimate

#Plot and save
ggplot(data=data2, aes(x=GTEXRATIO, y=PATRATIO,))+
  theme_bw()+
  geom_point(alpha=0.5,shape=16)+
  geom_vline(xintercept=0.5, lty=2,)+
  geom_hline(yintercept=0.5, lty=2,)+
  ylab('Observed variant ratio')+
  xlab('GTEx variant ratio')+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  annotate("text", x = 0.25, y = 1, label = paste0("Concordance: ",(signif(concordance,3)/100)),
           size=8,parse=TRUE)+
  annotate("text", x = 0.25, y = 0.95, label = paste0("Correlation: ",signif(correlation,3) ),
           size=8,parse=TRUE)+
  theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
        axis.title.x = element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=16))
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_d/variant_ratio_observed_vs_gtex.pdf', width=600, height=600)



