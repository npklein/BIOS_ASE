

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.txt

input_and_output_dir <- '/groups/umcg-bios/tmp03/projects/outlierGeneASE/binomialTest/'
#### READ IN A AND B COUNTS ####
aCount <- read.table(paste0(input_and_output_dir,'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.20190124.aCount.txt'), header=T, check.names=FALSE)
bCount <- read.table(paste0(input_and_output_dir,'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.20190124.bCount.txt'),header=T, check.names=FALSE)
annotations <- as.data.frame(fread("/groups/umcg-bios/tmp03/projects/outlierGeneASE/binomialTest/sample_genes_outliers_phenotypes.removedCODAM.txt"))


abcount_all_genes <- data.frame()
number_of_genes <- length(unique(aCount$ENSEMBLID))
x <- 0
for(gene in unique(aCount$ENSEMBLID)){
  if (x %% 100 == 0){
      print(paste0(x,'/',number_of_genes))
  }
  x <- x + 1
  a_subset <- aCount[aCount$ENSEMBLID==gene,]
  b_subset <- bCount[bCount$ENSEMBLID==gene,]
  aCount_melt <- melt(a_subset, id.vars='ENSEMBLID')
  colnames(aCount_melt)[c(2,3)] = c('sample','aCount')
  
  
  bCount_melt <- melt(b_subset, id.vars='ENSEMBLID')
  colnames(bCount_melt)[c(2,3)] = c('sample','bCount')
  #####
  
  ##### MERGE AND FILTER, CALC RATIO ####
  ab_count <- merge(aCount_melt, bCount_melt, by=c('ENSEMBLID','sample'))
  ab_count <- ab_count[!is.na(ab_count$aCount),]
  ab_count <- ab_count[!is.na(ab_count$bCount),]
  ab_count <- ab_count[ab_count$aCount > 0 & ab_count$bCount >0,]
  ab_count$ratio <- ab_count$aCount / (ab_count$aCount+ab_count$bCount)
  if(nrow(ab_count) == 0){
    next
  }
  ####
  
  ##### READ IN ANNOTATION
  ab_count$biobank <- annotations[match(as.character(ab_count$sample), annotations$SAMPLE),]$biobank_id
  #####
  
  #ab_count %>%
  #  group_by(ENSEMBLID, biobank) %>%
  #  summarise_at(vars(-ratio), funs(mean(., na.rm=TRUE)))
  
  # calcualte the mean ratio per biobank to use as new p
#  av_ratio_per_biobank <- ddply(ab_count, .(ENSEMBLID, biobank), summarize,  ratio=mean(ratio))
  
#  ab_count$binom_pvals_other_p <- apply(ab_count, 1, function(x) {
#    av_gene  <- av_ratio_per_biobank[x['ENSEMBLID']==av_ratio_per_biobank$ENSEMBLID, ] 
#    biobank_p <- av_gene[match(x['biobank'], av_gene$biobank),]$ratio
                                         
#    binom.test(as.numeric(x[c('aCount','bCount')]), alternative="two.sided", 
#                 p=biobank_p, 
#                 conf.level=0.95)$p.value
#    })
  
  ab_count$binom_pvals <- apply(ab_count, 1, function(x) {
    binom.test(as.numeric(x[c('aCount','bCount')]), alternative="two.sided", 
               p=0.5, 
               conf.level=0.95)$p.value
  })
 # ab_count_all <- rbind(ab_count_all, ab_count)
  ab_count$p_bonferroni_gene <- p.adjust(ab_count$binom_pvals, method='bonferroni')
  ab_count$p_fdr_gene <- p.adjust(ab_count$binom_pvals, method='fdr')
  abcount_all_genes <- rbind(abcount_all_genes, ab_count)
  write.table(ab_count, paste0(input_and_output_dir,'genes/',gene,'.txt'))
}  

abcount_all_genes$p_bonferroni <- p.adjust(abcount_all_genes$binom_pvals, method='bonferroni')
abcount_all_genes$p_fdr <- p.adjust(abcount_all_genes$binom_pvals, method='fdr')
write.table(abcount_all_genes, paste0(input_and_output_dir, 'binom_ASE_all_genes.txt'), quote = F, sep='\t')

