library(ggplot2)
library(data.table)
library(reshape2)
library(ggforce)
library(ggpubr)
library(ggsignif)

###### READ IN DATA #####
major_allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.majorAllelle.chrALL.txt.filtered.txt'))
minor_allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.minorAllelle.chrALL.txt.filtered.txt'))
snp_info <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.txt'))
######

###### MELT AND FILTER #####
major_allele_counts_melt <- melt(major_allele_counts)
minor_allele_counts_melt <- melt(minor_allele_counts)
major_allele_counts_melt <- major_allele_counts_melt[!is.na(major_allele_counts_melt$value),]
minor_allele_counts_melt <- minor_allele_counts_melt[!is.na(minor_allele_counts_melt$value),]


allele_counts <- merge(major_allele_counts_melt, minor_allele_counts_melt, by=c('VARIANT','variable'))
colnames(allele_counts) <- c('VARIANT','variable','major','minor')
allele_counts$logfc <- log2(allele_counts$minor / allele_counts$major)
allele_counts <- merge(allele_counts, snp_info, by='VARIANT')

allele_counts$SNPEFFIMPACT <- factor(allele_counts$SNPEFFIMPACT, levels=c('MODIFIER','LOW','MODERATE','HIGH'))
allele_counts$total <- allele_counts$minor+allele_counts$major



###### sum per variant #####
minor_summed <- aggregate(allele_counts$minor, by=list(VARIANT=allele_counts$VARIANT), FUN=sum)
major_summed <- aggregate(allele_counts$major, by=list(VARIANT=allele_counts$VARIANT), FUN=sum)
summed_counts <- merge(minor_summed, major_summed, by='VARIANT')
colnames(summed_counts) <- c('VARIANT','minor','major')
summed_counts <- merge(summed_counts, snp_info, by='VARIANT')

summed_counts$SNPEFFIMPACT <- factor(summed_counts$SNPEFFIMPACT, levels=c('MODIFIER','LOW','MODERATE','HIGH'))

summed_counts$minorRatio <- summed_counts$minor/(summed_counts$minor+summed_counts$major)


summed_counts_filtered<-summed_counts[summed_counts$SNPEFFIMPACT == "LOW" | summed_counts$SNPEFFIMPACT == "MODERATE" | summed_counts$SNPEFFIMPACT == "HIGH",]

#summed_counts1<-summed_counts[summed_counts$AF <= 0.01,]
summed_counts1<-summed_counts_filtered

t.test(summed_counts1[summed_counts1$SNPEFFIMPACT=="HIGH",]$minorRatio, summed_counts1[summed_counts1$SNPEFFIMPACT=='MODERATE',]$minorRatio)$p.val
t.test(summed_counts1[summed_counts1$SNPEFFIMPACT=="HIGH",]$minorRatio, summed_counts1[summed_counts1$SNPEFFIMPACT=='LOW',]$minorRatio)$p.val
t.test(summed_counts1[summed_counts1$SNPEFFIMPACT=="MODERATE",]$minorRatio, summed_counts1[summed_counts1$SNPEFFIMPACT=='LOW',]$minorRatio)$p.val

comparisons <- list( c("HIGH", "MODERATE"), c("HIGH", "LOW"), c("MODERATE", "LOW") )
give.n <- function(x){
  return(c(y = -0.03, label = length(x)))
}


ggplot(data=summed_counts1, aes(summed_counts1$SNPEFFIMPACT, summed_counts1$minor/(summed_counts1$minor+summed_counts1$major), fill = summed_counts1$SNPEFFIMPACT))+
  geom_violin()+
  geom_boxplot(width=0.05)+
  theme_bw(base_size = 18)+
  xlab('Functional impact')+
  ylab('Fraction minor alleles\n(Minor allele counts / total counts)')+
  guides(fill=F)+
  #ggtitle("Rare variants (AF <= 0.01)")+
  scale_fill_brewer(palette="Dark2")+ 
  stat_compare_means(comparisons = comparisons, method='wilcox.test')+
  stat_summary(fun.data = give.n, geom = "text")
  
  
#ggsave('/Users/freerkvandijk/Downloads/minor_allele_fraction.manuscript.20190303.png',width=8, height=8)

