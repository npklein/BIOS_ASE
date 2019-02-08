

library(ggplot2)
library(data.table)
library(reshape2)
library(ggforce)
library(ggpubr)
library(ggsignif)

###### READ IN DATA #####
major_allele_counts <- data.frame(fread('/Users/freerkvandijk/Downloads/counts.matrix.majorAllelle.chrALL.txt.filtered.txt'))
minor_allele_counts <- data.frame(fread('/Users/freerkvandijk/Downloads/counts.matrix.minorAllelle.chrALL.txt.filtered.txt'))
snp_info <- data.frame(fread('/Users/freerkvandijk/Downloads/counts.chr22.addedCADD.txt'))
######

###### MELT AND FILTER #####
major_allele_counts_melt <- melt(major_allele_counts)
minor_allele_counts_melt <- melt(minor_allele_counts)
major_allele_counts_melt <- major_allele_counts_melt[!is.na(major_allele_counts_melt$value),]
minor_allele_counts_melt <- minor_allele_counts_melt[!is.na(minor_allele_counts_melt$value),]


allele_counts <- merge(major_allele_counts_melt, minor_allele_counts_melt, by=c('VARIANT','variable'))
colnames(allele_counts) <- c('VARIANT','variable','major','minor')
#write.table(allele_counts,file = "/Users/freerkvandijk/Downloads/MajorMinorCountsPerVariantPerSample.txt",row.names=FALSE,col.names=TRUE,quote=FALSE,sep = "\t")
allele_counts$logfc <- log2(allele_counts$minor / allele_counts$major)
allele_counts <- merge(allele_counts, snp_info, by='VARIANT')

allele_counts$SNPEFFIMPACT <- factor(allele_counts$SNPEFFIMPACT, levels=c('MODIFIER','LOW','MODERATE','HIGH'))
allele_counts$total <- allele_counts$minor+allele_counts$major
######


###### sum per variant #####
minor_summed <- aggregate(allele_counts$minor, by=list(VARIANT=allele_counts$VARIANT), FUN=sum)
major_summed <- aggregate(allele_counts$major, by=list(VARIANT=allele_counts$VARIANT), FUN=sum)
summed_counts <- merge(minor_summed, major_summed, by='VARIANT')
colnames(summed_counts) <- c('VARIANT','minor','major')
summed_counts <- merge(summed_counts, snp_info, by='VARIANT')

#Remove MODIFIER category from dataframe
summed_counts<-summed_counts[summed_counts$SNPEFFIMPACT != "MODIFIER" | is.na(summed_counts$SNPEFFIMPACT),]

summed_counts$SNPEFFIMPACT <- factor(summed_counts$SNPEFFIMPACT, levels=c('LOW','MODERATE','HIGH'))

summed_counts$minorRatio <- summed_counts$minor/(summed_counts$minor+summed_counts$major)
t.test(summed_counts[summed_counts$SNPEFFIMPACT=="HIGH",]$minorRatio, summed_counts[summed_counts$SNPEFFIMPACT=='MODERATE',]$minorRatio)
t.test(summed_counts[summed_counts$SNPEFFIMPACT=="HIGH",]$minorRatio, summed_counts[summed_counts$SNPEFFIMPACT=='LOW',]$minorRatio)
t.test(summed_counts[summed_counts$SNPEFFIMPACT=="MODERATE",]$minorRatio, summed_counts[summed_counts$SNPEFFIMPACT=='LOW',]$minorRatio)

comparisons <- list( c("HIGH", "MODERATE"), c("HIGH", "LOW"), c("MODERATE", "LOW") )
give.n <- function(x){
  return(c(y = -0.03, label = length(x)))
}
ggplot(summed_counts, aes(SNPEFFIMPACT, minor/(minor+major), fill = SNPEFFIMPACT))+
  geom_violin()+
  geom_boxplot(width=0.05)+
  theme_bw(base_size = 18)+
  xlab('Functional impact')+
  ylab('Fraction minor alleles\n(Minor allele counts / total counts)')+
  guides(fill=F)+
  scale_fill_brewer(palette="Dark2")+ 
  stat_compare_means(comparisons = comparisons, method='wilcox.test') + 
  stat_summary(fun.data = give.n, geom = "text")
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig4//minor_allele_fraction.manuscript.20190129.png',width=8, height=8)
####

##### sampling same sample size
summed_counts_HIGH <- summed_counts[summed_counts$SNPEFFIMPACT=="HIGH",]
summed_counts_MODERATE <- summed_counts[summed_counts$SNPEFFIMPACT=="MODERATE",]
summed_counts_MODERATE_subset <- summed_counts_MODERATE[summed_counts_MODERATE$MAF %in% summed_counts_HIGH$MAF,]
summed_counts_LOW <- summed_counts[summed_counts$SNPEFFIMPACT=="LOW",]
summed_counts_LOW_subset <- summed_counts_LOW[summed_counts_LOW$MAF %in% summed_counts_HIGH$MAF,]
summed_counts_MODIFIER <- summed_counts[summed_counts$SNPEFFIMPACT=="MODIFIER",]
summed_counts_MODIFIER_subset <- summed_counts_MODIFIER[summed_counts_MODIFIER$MAF %in% summed_counts_HIGH$MAF,]

summed_counts_subsetted <- rbind(summed_counts_HIGH, summed_counts_MODERATE_subset)
summed_counts_subsetted <- rbind(summed_counts_subsetted, summed_counts_LOW_subset)
summed_counts_subsetted <- rbind(summed_counts_subsetted, summed_counts_MODIFIER_subset)

ggplot(summed_counts_subsetted, aes(SNPEFFIMPACT, minor/(minor+major), fill = SNPEFFIMPACT))+
  geom_violin()+
  geom_boxplot(width=0.05)+
  theme_bw(base_size = 18)+
  xlab('Impact')+
  ylab('Minor allele counts / total counts')+
  guides(fill=F)+
  scale_fill_brewer(palette="Dark2")+ 
  stat_compare_means(comparisons = comparisons, method='wilcox.test') + 
  stat_summary(fun.data = give.n, geom = "text")
ggsave('figures/minor_allele_fraction_MAF_select.png',width=8, height=8)
#####

##### SNPeff annotation
summed_counts_annotation_subset <- summed_counts[!grepl('&', summed_counts$SNPEFFANNOTATION),]
ggplot(summed_counts_annotation_subset, aes(SNPEFFANNOTATION, minor/(minor+major)))+
  geom_sina()+
  theme_bw(base_size = 18)+
  xlab('Impact')+
  ylab('Minor allele counts / total counts')+
  guides(fill=F)
ggsave('figures/minor_allele_fraction.png',width=8, height=8)
####
