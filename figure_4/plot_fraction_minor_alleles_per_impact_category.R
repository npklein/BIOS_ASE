library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)
library(ggsignif)

###### READ IN DATA #####
major_allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.majorAllelle.chrALL.txt.filtered.txt'))
minor_allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.minorAllelle.chrALL.txt.filtered.txt'))
snp_info <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.txt'))
#####

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


summed_counts_filtered <- summed_counts[summed_counts$SNPEFFIMPACT == "LOW" | summed_counts$SNPEFFIMPACT == "MODERATE" | summed_counts$SNPEFFIMPACT == "HIGH",]


give.n <- function(x){
  return(c(y = -0.03, label = length(x)))
}

summed_counts_filtered$frac_minor_allele <- summed_counts_filtered$minor/(summed_counts_filtered$minor+summed_counts_filtered$major)

# Calculate the p-values
comparison_df <- data.frame(compare_means(frac_minor_allele ~SNPEFFIMPACT,  
                             data = summed_counts_filtered,
                              method = "wilcox.test"))

# make datarfame with same colnames as summed_counts_filtered, set locations of text
comparison_df$frac_minor_allele <- NA
comparison_df$SNPEFFIMPACT <- NA
comparison_df[comparison_df$group2 == 'HIGH' & comparison_df$group1=='MODERATE',]$frac_minor_allele <- 1.13
comparison_df[comparison_df$group2 == 'HIGH' & comparison_df$group1=='MODERATE',]$SNPEFFIMPACT <- 'MODERATE'
comparison_df[comparison_df$group2 == 'HIGH' & comparison_df$group1=='LOW',]$frac_minor_allele <- 1.23
comparison_df[comparison_df$group2 == 'HIGH' & comparison_df$group1=='LOW',]$SNPEFFIMPACT <- 'MODERATE'
comparison_df[comparison_df$group2 == 'MODERATE' & comparison_df$group1=='LOW',]$frac_minor_allele <- 1.33
comparison_df[comparison_df$group2 == 'MODERATE' & comparison_df$group1=='LOW',]$SNPEFFIMPACT <- 'LOW'
 
p <- ggplot(data=summed_counts_filtered, aes(SNPEFFIMPACT, 
                                             frac_minor_allele,
                                             fill = SNPEFFIMPACT))+
  geom_violin()+
  geom_boxplot(width=0.05)+
  theme_bw(base_size = 18)+
  xlab('Functional impact')+
  ylab('fracion minor alleles\n(Minor allele counts / total counts)')+
  guides(fill=F)+
  scale_fill_brewer(palette="Dark2")+ 
  stat_summary(fun.data = give.n, geom = "text")+
  geom_text(data = comparison_df[!(comparison_df$group2 == 'HIGH' & comparison_df$group1=='LOW'),], 
            aes(x = SNPEFFIMPACT, y =frac_minor_allele , label = p.adj), size = 4,
            hjust=-2.2)+
  geom_text(data = comparison_df[(comparison_df$group2 == 'HIGH' & comparison_df$group1=='LOW'),], 
            aes(x = SNPEFFIMPACT, y =frac_minor_allele , label = p.adj), size = 4)+
  geom_segment(aes(x='HIGH', y=1.1, xend='MODERATE', yend=1.1), size=0.3) + # line 1
  geom_segment(aes(x='HIGH', y=1.1, xend='HIGH', yend=1.075), size=0.3) + # line 1
  geom_segment(aes(x='MODERATE', y=1.1, xend='MODERATE', yend=1.075), size=0.3) + # line 1
  geom_segment(aes(x='HIGH', y=1.2, xend='LOW', yend=1.2), size=0.3) + # line 2
  geom_segment(aes(x='HIGH', y=1.2, xend='HIGH', yend=1.175), size=0.3) + # line 2
  geom_segment(aes(x='LOW', y=1.2, xend='LOW', yend=1.175), size=0.3) + # line 2
  geom_segment(aes(x='MODERATE', y=1.3, xend='LOW', yend=1.3), size=0.3) + # line 3
  geom_segment(aes(x='MODERATE', y=1.3, xend='MODERATE', yend=1.275), size=0.3) + # line 3
  geom_segment(aes(x='LOW', y=1.3, xend='LOW', yend=1.275), size=0.3) # line 3

  #ggtitle("Rare variants (AF <= 0.01)")+  
#  stat_compare_means(aes(label = ..p.adj..),
#                    comparisons = comparisons, method='wilcox.test')+

  
#ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig4/minor_allele_fracion.manuscript.20190303.png',
#        plot=p,  width=8, height=8)


ggsave('minor_allele_fracion.manuscript.20190303.png',
        plot=p,  width=8, height=8)
