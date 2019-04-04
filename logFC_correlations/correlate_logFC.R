
# correlate the logFC values to all pheno and run statistic data

library(reshape2)
library(data.table)

out_dir <- "/groups/umcg-bios/tmp03/projects/outlierGeneASE/correlateLogFC/"

  
  logFC <- as.data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.values.txt'))
  logFC <- logFC[!grepl(";", logFC$ENSEMBLID),]
  logFC$NSAMPLES <- NULL
  logFC$MEDIAN <- NULL
  logFC$MAD <- NULL
  logFC_melt <- melt(logFC, id.vars='ENSEMBLID')
  colnames(logFC_melt) <- c("GENE", "sample", "logFC")
  
  # because hap A of one gene is not hap A of other gene, we can not look at direction
  # change to absolute
  logFC_melt$logFC_abs <- abs(logFC_melt$logFC)
  annotations <- as.data.frame(fread("/groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/AllSamples.phenotypes.nGenesANDnOutliers.txt"))

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
  
  
  annotations_lipids <- annotations[c("SAMPLENAME", "LDLchol", "HDLchol", "Triglycerides", "TotChol")]
  
  
  # per gene calculate the correlations
  for( gene in unique(logFC_melt$GENE)){
    i <- i+1
    if(grepl(';', gene)){next}
    logFC_subset <- logFC_melt[logFC_melt$GENE == gene,]
    logFC_subset <- logFC_subset[!is.na(logFC_subset$logFC_abs),]
    # because we filtered some samples out of the table, some genes are NA for all
    print(paste0(i,'/',l,'  ', nrow(logFC_subset) ))
    if(nrow(logFC_subset) < 2){next}

      logFC_subset_with_info <- merge(logFC_subset, annotations_lipids, by.x='sample', by.y="SAMPLENAME")
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
  write.table(logFC_correlations,	paste0(out_dir,'logFC_correlations.txt'), sep='\t',quote=F, col.names=NA)
  write.table(logFC_pval,	paste0(out_dir, 'logFC_pvals.txt'), sep='\t',quote=F, col.names=NA)
  

