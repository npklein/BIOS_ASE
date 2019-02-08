library(ggplot2)
library(Hmisc)
library(reshape2)
library(dplyr)
library(ggpubr)
library(EnvStats)
library(scales)
library(ggrepel)



dataTab1<-read.table("/Users/freerkvandijk/Downloads/AllSamples.phenotypes.nGenesANDnOutliers.txt", sep="\t", header=TRUE)

ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=factor(dataTab1$biobank_id))) + theme_bw() + geom_point(alpha = 0.6)
ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=dataTab1$PF_READS)) + theme_bw() + geom_point(alpha = 0.6) + scale_color_gradientn(colours = rainbow(5))
ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=dataTab1$fastqc_clean.R1_clean_GC_mean)) + theme_bw() + geom_point(alpha = 0.6) + scale_color_gradientn(colours = rainbow(5))
ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=dataTab1$PF_HQ_ALIGNED_Q20_BASES)) + theme_bw() + geom_point(alpha = 0.6) + scale_color_gradientn(colours = rainbow(5))
ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=dataTab1$PERCENT_DUPLICATION)) + theme_bw() + geom_point(alpha = 0.6) + scale_color_gradientn(colours = rainbow(5))
ggplot(dataTab1,aes(x=dataTab1$NGENES, y=dataTab1$NOUTLIERS, colour=dataTab1$MEDIAN_5PRIME_TO_3PRIME_BIAS)) + theme_bw() + geom_point(alpha = 0.6) + scale_color_gradientn(colours = rainbow(5))

#Only keep sample having less than 1000 outlier ASE genes
genes_to_keep <- dataTab1[dataTab1$NOUTLIERS < 1000, ]
ggplot(genes_to_keep,aes(x=NGENES, y=NOUTLIERS, colour=factor(biobank_id))) + theme_bw() + geom_point(alpha = 0.6)

<<<<<<< HEAD
=======
#Replace the LLDeepNotInBios biobank IDs with LL
genes_to_keep1$biobank_id[genes_to_keep1$biobank_id == "LLDeepNotInBIOS"] <- "LL"

>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734
#Remove CODAM samples
genes_to_keep1<-genes_to_keep[genes_to_keep$biobank_id != "CODAM" | is.na(genes_to_keep$biobank_id),]

#Plot again
<<<<<<< HEAD
ggplot(genes_to_keep1,aes(x=NGENES, y=NOUTLIERS, colour=factor(biobank_id))) + theme_bw() + geom_point(alpha = 0.6)

#Write sample IDs
write.table(unique(as.character(genes_to_keep1$SAMPLE)), "/Users/freerkvandijk/Downloads/samples_NOUTLIERS500.depthFiltered.bonferroni.txt",sep='\t',quote=F,
=======
plot<-ggplot(genes_to_keep1,aes(x=NGENES, y=NOUTLIERS, colour=factor(biobank_id))) + theme_bw() + geom_point(alpha = 0.4) +
  ggtitle(paste0("Number of significant ASE genes per sample")) +
  ylab('Number of significant (P < 0.05) ASE genes')+
  xlab('Number of observed ASE genes per sample')+
  labs(colour= "Biobank")+
  theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
        axis.title.x = element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=16))

#Print plot to pdf
pdf("/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/suppl_ASEgenesVsSigASEgenes.pdf")
print(plot)
dev.off()


#Write sample IDs
write.table(unique(as.character(genes_to_keep1$SAMPLE)), "/Users/freerkvandijk/Downloads/samples_NOUTLIERS1000.depthFiltered.bonferroni.txt",sep='\t',quote=F,
>>>>>>> 2a80f51c302edfb1917d2925ee56cc7cdaab0734
            row.names=F, col.names=F)

