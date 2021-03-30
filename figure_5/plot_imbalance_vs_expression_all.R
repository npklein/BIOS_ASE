
library(data.table)
library(ggplot2)

#allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/geneExpressionAndMajorMinorAlleleCounts.GOI.list.txt'))
#allele_counts <- fread('geneExpressionAndMajorMinorAlleleCounts.allGenes.txt.gz')
allele_counts <- data.frame(fread('geneExpressionAndMajorMinorAlleleCounts.GOI.list.txt'))

object.size(allele_counts)

allele_counts$SAMPLE_GENE <- paste0(allele_counts$SAMPLE,'_',allele_counts$GENENAME)
allele_freqs <- fread('chrALL.AFsFromData.txt')
#allele_freqs <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.txt')

colnames(allele_freqs) <- c('VARIANT', 'allele_freq')
allele_counts_AF <- merge(allele_counts,allele_freqs, by='VARIANT' )
allele_counts_AF$majorAllele <- allele_counts_AF$REF
allele_counts_AF$minorAllele <- allele_counts_AF$ALT

allele_counts_AF[allele_counts_AF$allele_freq > 0.5,]$majorAllele <- allele_counts_AF[allele_counts_AF$allele_freq > 0.5,]$ALT
allele_counts_AF[allele_counts_AF$allele_freq > 0.5,]$minorAllele <- allele_counts_AF[allele_counts_AF$allele_freq > 0.5,]$REF
# cause allele counts 


#snpInfo <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt')
snpInfo <- fread('counts.chr22.addedCADD.addedVKGL.txt')
allele_counts_snpInfo <- merge(allele_counts, snpInfo, by=c('VARIANT','GENENAME'))
allele_counts_snpInfo <- allele_counts_snpInfo[!is.na(allele_counts_snpInfo$MINOR ),]
allele_counts_snpInfo_with_count_high_impact_variants <- allele_counts_snpInfo[allele_counts_snpInfo$CGDINHERITANCE=="AD" & allele_counts_snpInfo$CADDPHRED > 20 & 
                                                                                 allele_counts_snpInfo$SNPEFFIMPACT=="HIGH",]
high_impact_genes <- unique(allele_counts_snpInfo_with_count_high_impact_variants$GENENAME)

unique_expression <- allele_counts_AF[!duplicated(allele_counts_AF$SAMPLE_GENE),]
ggplot(data=unique_expression[unique_expression$GENENAME%in%high_impact_genes, ], aes(x=GENENAME, y=log10(GENEEXPRESSION+1)))+
  geom_violin(fill="grey90")+
  geom_boxplot(data=unique_expression[unique_expression$GENENAME%in%high_impact_genes, ], width=0.1, outlier.shape=NA)+
  geom_point(data=allele_counts_snpInfo_with_count_high_impact_variants[allele_counts_snpInfo_with_count_high_impact_variants$GENENAME!="ALOX5",],
              aes(x=GENENAME, y=log10(na.omit(GENEEXPRESSION)+1), fill = MINORRATIO),
             size=2.5,pch=21)+
  geom_point(data=allele_counts_snpInfo_with_count_high_impact_variants[allele_counts_snpInfo_with_count_high_impact_variants$GENENAME=="ALOX5"
                                                                        & allele_counts_snpInfo_with_count_high_impact_variants$VARIANT=="10_45891347_C_T",],
             aes(x=GENENAME, y=log10(na.omit(GENEEXPRESSION)+1), fill = MINORRATIO),
             size=2.5,pch=21, position = position_nudge(x = 0.1))+
  geom_point(data=allele_counts_snpInfo_with_count_high_impact_variants[allele_counts_snpInfo_with_count_high_impact_variants$GENENAME=="ALOX5"
                                                                        & allele_counts_snpInfo_with_count_high_impact_variants$VARIANT=="10_45938895_G_T",],
             aes(x=GENENAME, y=log10(na.omit(GENEEXPRESSION)+1), fill = MINORRATIO),
             size=2.5,pch=21, position = position_nudge(x = -0.1))+
  theme_bw(base_size=13)+
  xlab('')+
  ylab('log10(TPM+1)')+
  scale_fill_gradientn(colours = terrain.colors(10))+
  labs(fill="minor allele / major allele")
ggsave('~/Downloads/ase_samples_overall_expression_only_high_impact.pdf', width=8, height=5)
ggsave('~/Downloads/ase_samples_overall_expression_only_high_impact.png', width=8, height=5)









ggplot(data=allele_counts_AF, aes(x=GENENAME, y=log10(GENEEXPRESSION+1)))+
  geom_violin(fill="grey90")+
  geom_boxplot(data=allele_counts_AF[allele_counts_AF$GENENAME%in%high_impact_genes, ], width=0.1, outlier.shape=NA)+
  theme_bw(base_size=13)+
  xlab('')+
  ylab('log10(TPM+1)')+
  scale_fill_gradientn(colours = terrain.colors(10))+
  labs(fill="minor allele / major allele")



