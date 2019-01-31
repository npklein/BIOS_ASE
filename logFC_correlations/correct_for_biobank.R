
# correlate the logFC values to all pheno and run statistic data

library(ggplot2)
library(reshape2)
library(data.table)
library(RColorBrewer)
library(viridis)
library(pheatmap)



logFC <- data.frame(fread('genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.removedCODAM.4outliersRemoved.SD3.table.txt', fill=T))
logFC <- logFC[!grepl(";", logFC$ENSEMBLID),]
logFC$MAD <- NULL
logFC$MEDIAN <- NULL
logFC$NSAMPLES <- NULL
logFC_melt <- melt(logFC)
colnames(logFC_melt) <- c("GENE", "sample", "logFC")

annotations <- as.data.frame(fread("/Users/NPK/UMCG/projects/BIOS/number_of_outliers_per_sample/sample_genes_outliers_phenotypes.removedCODAM.txt"))
annotations$biobank_id_coded <- ifelse(annotations$biobank_id=="LL", 0, ifelse(annotations$biobank_id=='LLDeepNotInBIOS', 1, 
                                                                               ifelse(annotations$biobank_id=='LLDeepNotInBIOS', 2, 
                                                                                      ifelse(annotations$biobank_id=='LLS', 3, 
                                                                                             ifelse(annotations$biobank_id=='NTR', 4, 
                                                                                                    ifelse(annotations$biobank_id=='PAN', 5,
                                                                                                           ifelse(annotations$biobank_id=='RS', 6,NA)))))))

logFC_with_annotation <- merge(logFC_melt, annotations, by="sample", by.y="SAMPLE")
logFC_with_annotation$abs_logFC <- abs(logFC_with_annotation$logFC)


abs_logFC_with_biobank_lm <- lm(abs_logFC~biobank_id_coded, data=logFC_with_annotation)
