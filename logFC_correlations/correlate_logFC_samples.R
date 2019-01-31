
# correlate the logFC values to all pheno and run statistic data

library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)



logFC <- data.frame(fread('genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.removedCODAM.4outliersRemoved.SD3.table.txt', fill=T))
logFC <- logFC[!grepl(";", logFC$ENSEMBLID),]
logFC_melt <- melt(logFC)
colnames(logFC_melt) <- c("GENE", "sample", "logFC")

logFC$ENSEMBLID <- NULL
average_per_sample <- data.frame(apply(logFC, 2, function(x) sum(abs(x), na.rm=T) / length(x[!is.na(x)])))
average_per_sample$sample <- colnames(logFC)
colnames(average_per_sample)[1] <-'abs_logfc_summed'
# because hap A of one gene is not hap A of other gene, we can not look at direction
# change to absolute
annotations <- as.data.frame(fread("/Users/NPK/UMCG/projects/BIOS/number_of_outliers_per_sample/sample_genes_outliers_phenotypes.removedCODAM.txt"))

annotations$Sex <- ifelse(annotations$Sex == "Male", 0,  ifelse(annotations$Sex == "male", 0,  ifelse(annotations$Sex == "female", 1,  ifelse(annotations$Sex == "Female", 1, NA))))
annotations$Smoking <- ifelse(annotations$Smoking=="non-smoker", 0, ifelse(annotations$Smoking=='former smoker', 1, ifelse(annotations$Smoking=='current smoker', 2, NA)))
annotations$Lipids_BloodSampling_Fasting <- ifelse(annotations$Lipids_BloodSampling_Fasting=="yes", 1, ifelse(annotations$Lipids_BloodSampling_Fasting=="no", 0, NA))
annotations$LipidMed <- ifelse(annotations$LipidMed=="no", 0, ifelse(annotations$LipidMed=='yes but no statins', 1, 2))
annotations$biobank_id <- ifelse(annotations$biobank_id=="LL", 0, ifelse(annotations$biobank_id=='LLDeepNotInBIOS', 1, 
                               ifelse(annotations$biobank_id=='LLDeepNotInBIOS', 2, 
                                  ifelse(annotations$biobank_id=='LLS', 3, 
                                      ifelse(annotations$biobank_id=='NTR', 4, 
                                             ifelse(annotations$biobank_id=='PAN', 5,
                                                    ifelse(annotations$biobank_id=='RS', 6,NA)))))))
annotations$RNA_Extraction_Date <-  as.numeric(as.POSIXct(annotations$RNA_Extraction_Date, format="%Y-%m-%d"))
annotations$DNA_Extraction_Date<-  as.numeric(as.POSIXct(annotations$DNA_Extraction_Date, format="%Y-%m-%d"))
annotations$Sampling_Date<-  as.numeric(as.POSIXct(annotations$Sampling_Date, format="%Y-%m-%d"))


annotations_numeric <- annotations[, sapply(annotations, is.numeric)] 
annotation_cor <- data.frame(cor(annotations_numeric, use='pairwise.complete.obs'))
annotation_cor[is.na(annotation_cor)] <- 0
col<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("correlation_pheno_data.pdf",width=30, height=30)
heatmap(as.matrix(annotation_cor), col=col, symm=TRUE)
dev.off()

to_filter <- c('num_splice_annotated','star.num_splice_atac','star.num_splice_gcag','star.num_splice_gtag','star.num_splice_noncanonical',
               'star.num_splice_total','star.num_input','star.num_unique_mapped','PF_READS','PF_HQ_ALIGNED_BASES',
               'PF_HQ_ALIGNED_Q20_BASES','READS_ALIGNED_IN_PAIRS','PF_READS_ALIGNED','PF_HQ_ALIGNED_BASES','bam.genome_mapped',
               'bam.exon_total','bam.exon_mapped','bam.genome_mapped','CODING_BASES','UTR_BASES','PF_ALIGNED_BASES',
               'READ_PAIRS_EXAMINED','MEDIAN_INSERT_SIZE','bam.genome_insert_std','bam.genome_insert_mean',
               'WIDTH_OF_10_PERCENT','WIDTH_OF_20_PERCENT','WIDTH_OF_30_PERCENT','WIDTH_OF_40_PERCENT','WIDTH_OF_50_PERCENT',
               'WIDTH_OF_60_PERCENT','WIDTH_OF_70_PERCENT','INTERGENIC_BASES','prime_bias.MEDIAN_3PRIME_BIAS',
               'WIDTH_OF_80_PERCENT','WIDTH_OF_90_PERCENT',
               'MEDIAN_5PRIME_BIAS', 'MEDIAN_5PRIME_TO_3PRIME_BIAS','fastqc_clean.R1_clean_GC_mean',
               'fastqc_raw.R1_raw_GC_mean','bam.genome_total','PF_ALIGNED_BASES.1','PF_HQ_ALIGNED_READS','prime_bias.MEDIAN_5PRIME_BIAS',
               'fastqc_clean.R1_clean_GC_std','fastqc_clean.R2_clean_GC_mean','fastqc_clean.R2_clean_GC_std',
               'STANDARD_DEVIATION','MEDIAN_ABSOLUTE_DEVIATION', 'LUC','Eos','PCT_ADAPTER','PF_MISMATCH_RATE','PF_HQ_ERROR_RATE',
               'BAD_CYCLES','UNPAIRED_READS_EXAMINED','bam.genome_duplicates','READ_PAIR_OPTICAL_DUPLICATES','PERCENT_DUPLICATION',
               'READ_PAIRS','star.num_mapped_many','star.num_splice_annotated','PF_BASES','LipidsMed_Age','Anthropometry_Age',
               'TotChol','Baso','Granulocyte_Perc','star.avg_input_length','MCH','DNA_Extraction_Date','PCT_MRNA_BASES')

for(column in to_filter){
  annotations_numeric[column] <- NULL
  annotations[column] <- NULL
}

annotation_cor <- data.frame(cor(annotations_numeric, use='pairwise.complete.obs'))
annotation_cor[is.na(annotation_cor)] <- 0
col<- colorRampPalette(c("blue", "white", "red"))(20)
pdf("correlation_pheno_data_filtered.pdf",width=30, height=30)
heatmap(as.matrix(annotation_cor), col=col, symm=TRUE,cexRow=0.5,cexCol=0.5)
dev.off()



average_per_sample_with_annotation <- merge(average_per_sample, annotations, by="sample", by.y="SAMPLENAME")


rownames(average_per_sample_with_annotation) <- average_per_sample_with_annotation$sample
average_per_sample_with_annotation_numeric <- average_per_sample_with_annotation[, sapply(average_per_sample_with_annotation, is.numeric)]


correlations <- as.data.frame(t(cor(average_per_sample_with_annotation_numeric$abs_logfc_summed, average_per_sample_with_annotation_numeric , use="pairwise.complete.obs",
                                    method="spearman")))
pval <- data.frame(sapply(average_per_sample_with_annotation_numeric, function(x) {
  if(sum(complete.cases(average_per_sample_with_annotation_numeric$abs_logfc_summed, x))<3){
    return(NA)
  }
  cor.test(average_per_sample_with_annotation_numeric$abs_logfc_summed, x, method="spearman")$p.value
}))
cor_and_pval <- cbind(correlations,pval)
colnames(cor_and_pval) <- c("correlation","pvalue")


cor_and_pval$pheno <- rownames(cor_and_pval)
cor_and_pval <- cor_and_pval[!is.na(cor_and_pval$correlation),]

cor_and_pval %>% 
  mutate(pheno = factor(pheno, levels = pheno[order(correlation)])) %>%  # Order by correlation strength
  ggplot(aes(x = pheno, y = correlation,fill=-log10(pvalue+0.000000001))) +
  geom_bar(stat = "identity") +
  ylab("Correlation with summed abs(logFC)") +
  theme_pubr()+ 
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab('')+ 
  guides(fill=guide_legend(title="-log10(pvalue)"))+
  scale_fill_gradient(low = "white", high = "black")

ggsave("correlation_abs_logFC_samples.pdf", width=20, height=6)
