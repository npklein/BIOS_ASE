library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)

data <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/concordanceGTEx/counts.matrix.AlleleAdded.txt'))
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
data1$FDR<-p.adjust(data1$binom, method = "fdr")

#Select subset based on p-val
#data2<-data1[data1$binomCorrected < 1e-5,]
data2<-data1[data1$FDR < 0.05,]

#Check positive or negative ratios and annotate with them
data2$OURDIR <- ifelse(data2$RATIO<0.5, "NEG", "POS")
data2$GTEXDIR <- ifelse(data2$GTEXRATIO<0.5, "NEG", "POS")

#Binom test on GTEx ASE counts
data2$binomGTEx<- apply(data2, 1, function(x) binom.test(as.numeric(x[["GTEXSUMMAJOR"]]), as.numeric(x[["GTEXTOTAL"]]), p=0.5, alternative = "two.sided", conf.level = 0.95)$p.value)
data2$binomCorrectedGTEx <- p.adjust(data2$binomGTEx, method = "bonferroni")
data2$ZSCORE_ASE_GTEx<- -1*qnorm(data2$binomCorrectedGTEx/2)
data2$FDRGTEx<-p.adjust(data2$binomGTEx, method = "fdr")

#Filter GTEx on FDR
data5<-data2[data2$FDRGTEx < 0.05,]

#Calculate concordance
comp <- table(data5$OURDIR==data5$GTEXDIR)
concordance <- (comp["TRUE"]/(comp["TRUE"]+comp["FALSE"]))*100
concordance

#Calculate correlation
res <- cor.test(data5$RATIO, data5$GTEXRATIO, 
                method = "spearman")
correlation<-res$estimate


#Make plot
png(paste0("/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/concordance.Observed.vs.GTEx.ASE.png"), width = 800, height = 800)

myplot<-
ggplot(data=data5, aes(x=GTEXRATIO, y=RATIO, colour=ZSCORE_ASE))+
  geom_point(alpha=0.5,shape=16)+
  #ggtitle(paste0("GENE ratio vs GTEx ratio\nConcordance: ", concordance, "    Correlation: ", correlation))+
  ggtitle(paste0("Observed variant ratio vs GTEx ratio"))+
  ylab('Observed variant ratio')+
  xlab('GTEx variant ratio')+
  geom_vline(xintercept=0.5)+
  geom_hline(yintercept=0.5)+
  labs(color = paste0("Observed variant \nratio Z-score"))+
  scale_x_continuous(limits = c(0, 1))+
  scale_y_continuous(limits = c(0, 1))+
  #scale_color_gradientn(colours = rainbow(5))+
  annotate("text", x = 0.25, y = 1, label = paste0("Concordance: ",(signif(concordance,3)/100)),
           size=4,parse=TRUE)+
  annotate("text", x = 0.25, y = 0.95, label = paste0("Correlation: ",signif(correlation,3) ),
           size=4,parse=TRUE)+
  theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
        axis.title.x = element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=16))


#Print plot to output device
print(myplot)
dev.off()


