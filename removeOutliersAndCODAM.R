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

#Remove CODAM samples
genes_to_keep1<-genes_to_keep[genes_to_keep$biobank_id != "CODAM" | is.na(genes_to_keep$biobank_id),]

#Plot again
ggplot(genes_to_keep1,aes(x=NGENES, y=NOUTLIERS, colour=factor(biobank_id))) + theme_bw() + geom_point(alpha = 0.6)

#Write sample IDs
write.table(unique(as.character(genes_to_keep1$SAMPLE)), "/Users/freerkvandijk/Downloads/samples_NOUTLIERS500.depthFiltered.bonferroni.txt",sep='\t',quote=F,
            row.names=F, col.names=F)

