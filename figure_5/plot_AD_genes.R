
library(data.table)
library(ggplot2)

allele_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/geneExpressionAndMajorMinorAlleleCounts.GOI.list.txt'))
allele_counts$SAMPLE_GENE <- paste0(allele_counts$SAMPLE,'_',allele_counts$GENENAME)


snpInfo <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt')
allele_counts <- merge(allele_counts, snpInfo, by=c('VARIANT','GENENAME'))


allele_counts_high_impact <- allele_counts[allele_counts$IMPACT == "HIGH",]

allele_counts_high_impact_unique <- allele_counts_high_impact[!duplicated(allele_counts_high_impact$SAMPLE_GENE),]

allele_counts_high <- allele_counts[allele_counts$IMPACT == "HIGH",]
allele_counts_low <- allele_counts[allele_counts$IMPACT == "LOW",]
allele_counts_moderate <- allele_counts[allele_counts$IMPACT == "MODERATE",]
allele_counts_modifier <- allele_counts[allele_counts$IMPACT == "MODIFIER",]

allele_counts_unique_1 <- allele_counts_high[!duplicated(allele_counts_high$SAMPLE_GENE),]
allele_counts_unique_2 <- allele_counts_low[!duplicated(allele_counts_low$SAMPLE_GENE),]
allele_counts_unique_3 <- allele_counts_moderate[!duplicated(allele_counts_moderate$SAMPLE_GENE),]
allele_counts_unique_4 <- allele_counts_modifier[!duplicated(allele_counts_modifier$SAMPLE_GENE) ,]


allele_counts_unique_minorRatio <- rbind(allele_counts_unique_1, allele_counts_unique_2)
allele_counts_unique_minorRatio <- rbind(allele_counts_unique_minorRatio, allele_counts_unique_3)
allele_counts_unique_minorRatio <- rbind(allele_counts_unique_minorRatio, allele_counts_unique_4)

allele_counts_unique_minorRatio_overlappingGene <- allele_counts_unique_minorRatio[allele_counts_unique_minorRatio$GENENAME %in% allele_counts_high_impact_unique$GENENAME,]
ggplot(data=allele_counts_high_impact_unique, aes(x=GENENAME, y=log2(GENEEXPRESSION+1)))+
  geom_violin()+
  geom_boxplot(width=0.05, outlier.shape=NA)+
  geom_jitter(data=allele_counts_unique_minorRatio_overlappingGene[allele_counts_unique_minorRatio_overlappingGene$IMPACT != "MODIFIER",],
              aes(x=GENENAME, y=log2(na.omit(GENEEXPRESSION)+1), size = MINORRATIO, colour=IMPACT, alpha=0.8), 
              position=position_jitter(w=0.4,h=0.1))+
  theme_bw(base_size=18)+
#  scale_y_continuous(limit=c(0,680))+
  guides(fill=F, alpha=F, colour=F)+
  scale_colour_brewer(palette="Dark2")+
  theme(legend.position="top",
        legend.box = "vertical")+
  xlab('')+
  ylab('log2(TPM+1)')
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig5/panel_b/ase_samples_overall_expression.png', width=8, height=8)

# plot this to stitch the legend together in .ai
ggplot(data=allele_counts_high_impact_unique, aes(x=GENENAME, y=log2(GENEEXPRESSION+1), fill = GENENAME))+
  geom_jitter(data=allele_counts_unique_minorRatio_overlappingGene,
              aes(x=GENENAME, y=log2(na.omit(GENEEXPRESSION)+1), size = MINORRATIO, colour=IMPACT, alpha=0.8),
              position=position_jitter(w=0.4,h=0.1))+
  theme_bw(base_size=18)+
  guides(fill=F, alpha=F)+
  scale_colour_brewer(palette="Dark2")+
  guides(colour = guide_legend(override.aes = list(size=5)))
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig5/panel_b/ase_samples_overall_expression_legend.png', width=8, height=8)



ggplot(data=allele_counts_high_impact_unique, aes(x=GENENAME, y=log2(GENEEXPRESSION+1), fill = GENENAME))+
  geom_violin()+
  geom_boxplot(width=0.05, outlier.shape=NA)+
  geom_jitter(data=allele_counts_unique_minorRatio[allele_counts_unique_minorRatio$IMPACT=="HIGH",],
              aes(x=GENENAME, y=log2(na.omit(GENEEXPRESSION)+1), size = MINORRATIO), 
              position=position_jitter(w=0,h=0.1))+
  theme_bw(base_size=18)+
#  scale_y_continuous(limit=c(0,680))+
  guides(fill=F, alpha=F)+
  scale_colour_brewer(palette="Dark2")+
  xlab('')+
  ylab('log2(TPM+1)')+ 
  guides(colour=guide_legend(title="Alt counts / ref counts"))
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig5/panel_a//ase_samples_overall_expression_only_high_impact.png', width=8, height=8)



