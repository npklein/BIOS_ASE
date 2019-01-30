
# correlate the logFC values to all pheno and run statistic data

library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)



do_correlations <- function(){
  
  logFC <- as.data.frame(fread('genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.removedCODAM.4outliersRemoved.SD3.table.txt',fill=T))
  logFC <- logFC[!grepl(";", logFC$ENSEMBLID),]
  logFC$NSAMPLES <- NULL
  logFC$MEDIAN <- NULL
  logFC$MAD <- NULL
  logFC_melt <- melt(logFC)
  colnames(logFC_melt) <- c("GENE", "sample", "logFC")
  
  # because hap A of one gene is not hap A of other gene, we can not look at direction
  # change to absolute
  logFC_melt$logFC_abs <- abs(logFC_melt$logFC)
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
  logFC_correlations <- data.frame()
  logFC_pval <- data.frame()
  i <- 0
  l <- length(unique(logFC_melt$GENE))
  
  
  
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
    annotations[column] <- NULL
  }
  
  # per gene calculate the correlations
  for( gene in unique(logFC_melt$GENE)){
    i <- i+1
    if(grepl(';', gene)){next}
    logFC_subset <- logFC_melt[logFC_melt$GENE == gene,]
    logFC_subset <- logFC_subset[!is.na(logFC_subset$logFC_abs),]
    # because we filtered some samples out of the table, some genes are NA for all
    print(paste0(i,'/',l,'  ', nrow(logFC_subset) ))
    if(nrow(logFC_subset) < 2){next}

      logFC_subset_with_info <- merge(logFC_subset, annotations, by.x='sample', by.y="SAMPLENAME")
    # extract all numeric columns, as these are possible to correlate with
    logFC_subset_with_info_numeric <- logFC_subset_with_info[, sapply(logFC_subset_with_info, is.numeric)]
  
    # pairwise.complete.obs makes sure that when there is an NA in two pairwise compared columns, this row does not get included
    correlationTMP <- as.data.frame(cor(logFC_subset_with_info_numeric$logFC_abs, logFC_subset_with_info_numeric , use="pairwise.complete.obs"),
                                    method="spearman")
  
    pvalTMP <- as.data.frame(t(sapply(logFC_subset_with_info_numeric, function(x) {
                                    if(sum(complete.cases(logFC_subset_with_info_numeric$logFC, x))<3){
                                         return(NA)
                                    }
                                    cor.test(logFC_subset_with_info_numeric$logFC, x, method="spearman")$p.value
    })))
  
    
    rownames(correlationTMP) <- gene
    rownames(pvalTMP) <- gene
    logFC_correlations <- rbind(logFC_correlations, correlationTMP)
    logFC_pval <- rbind(logFC_pval, pvalTMP)
  }
  write.table(logFC_correlations,	'logFC_correlations.txt', sep='\t',quote=F, col.names=NA)
  write.table(logFC_pval,	'logFC_pvals.txt', sep='\t',quote=F, col.names=NA)
  
}

if(!file.exists('logFC_correlations.txt')){
  do_correlations()
}

logFC_correlations <- read.table('logFC_correlations.txt', sep='\t', header=T, row.names=1)
logFC_pval <- read.table('logFC_pvals.txt', sep='\t', header=T, row.names=1)



logFC_correlations_remove_NA <- logFC_correlations
logFC_correlations_remove_NA <- logFC_correlations_remove_NA[,colSums(is.na(logFC_correlations_remove_NA))<nrow(logFC_correlations_remove_NA)]

logFC_correlations_remove_NA[is.na(logFC_correlations_remove_NA)] <- 0
pdf("correlations_logFC.pdf", width=10, height=10)
pheatmap(logFC_correlations_remove_NA, 
         show_rownames=F,
         cluster_rows=T,
         cluster_cols=T)
dev.off()


not_significants <- logFC_pval >= 0.05
not_significants[is.na(not_significants)] <- TRUE
logFC_correlations_significant <- logFC_correlations
logFC_correlations_significant[not_significants] <- 0

logFC_correlations_significant_remove_NA <- logFC_correlations_significant
logFC_correlations_significant_remove_NA <- logFC_correlations_significant_remove_NA[,colSums(is.na(logFC_correlations_significant_remove_NA))<nrow(logFC_correlations_significant_remove_NA)]

logFC_correlations_significant_remove_NA[is.na(logFC_correlations_significant_remove_NA)] <- 0
logFC_correlations_significant_remove_NA_at_least_one <- logFC_correlations_significant_remove_NA[rowSums(logFC_correlations_significant_remove_NA) > 0,]
pdf("correlations_logFC_significants.pdf", width=10, height=10)
pheatmap(logFC_correlations_significant_remove_NA_at_least_one, 
         show_rownames=F,
         cluster_rows=T,
         cluster_cols=T)
dev.off()



t <- logFC_correlations_significant[order(-abs(logFC_correlations_significant$MEAN_INSERT_SIZE)),]
