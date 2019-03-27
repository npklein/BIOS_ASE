
library(ggplot2)

##### read in data ####
input_dir <- '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/ASE_counts_merged_with_LCL/'

# Can do it for multiple differnet input files with different filters, finding all those files
file.names <- dir(input_dir, pattern =".txt")

used_in_manuscript = 'ASE.cov20'

numbers_and_pvalues_results <- list()
##### read in data ####
for(f in file.names){
  name <- gsub('.mergedWithLCL.txt','', gsub('eQTL_freeze2_snpASE.','',f))
  print(paste(f, name))
  eqtl_and_ASE <- read.table(paste0(input_dir,f),header=T,sep='\t')
  eqtl_and_ASE$LCL_ratio <- eqtl_and_ASE$LCL_ratio-0.5
  eqtl_and_ASE$fdr <- p.adjust(eqtl_and_ASE$binom_pval, method='fdr')
  
  #####
  
  
  ##### calculate concordance and correlation using ours VS LCL ASE with FDR pvalues < 0.05 #####
  # only have to select on FDR for our data, because LCL data is already filtered to only be FDR < 0.05
  eqtl_and_ASE_fdr05 <- eqtl_and_ASE[eqtl_and_ASE$fdr < 0.05,]
  concordance_LCL_logFC <- sum(sign(eqtl_and_ASE_fdr05$logFC)==sign(eqtl_and_ASE_fdr05$LCL_ratio))/nrow(eqtl_and_ASE_fdr05)
  cor_zscore_logFC_LCL <- cor(eqtl_and_ASE_fdr05$logFC, eqtl_and_ASE_fdr05$LCL_ratio, method='spearman')

  if(name==used_in_manuscript){
      # TODO: write these values to a file
      numbers_and_pvalues_results$n_overlap_ASE_LCL <- nrow(eqtl_and_ASE_fdr05)
      numbers_and_pvalues_results$concordance_with_LCL <- concordance_LCL_logFC*100
      numbers_and_pvalues_results$correlation_with_LCL <- cor_zscore_logFC_LCL
  }

  ggplot(eqtl_and_ASE_fdr05, aes(logFC, LCL_ratio, size=-log(fdr+0.00000000001, base=10),
                           colour=-log10(LCL_pval+0.000000000001)))+
    geom_point(alpha=0.5,shape=16)+
    geom_hline(yintercept=0, lty=2, colour='red')+
    geom_vline(xintercept=0, lty=2, colour='red')+
    theme_classic(base_size = 20)+
    ylab('Allelic ratio LCL cell line')+
    xlab('logFC blood')+
    annotate("text", x = -3, y = 0.45, label = paste0("SNPs: ",nrow(eqtl_and_ASE_fdr05)),
             size=8)+
    annotate("text", x = -3, y = 0.4, label = paste0("Concordance: ",signif(concordance_LCL_logFC,3) ),
             size=8,parse=TRUE)+
    annotate("text", x = -3, y = 0.35, label = paste0("Correlation: ",signif(cor_zscore_logFC_LCL,3) ),
             size=8,parse=TRUE)+
    geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
    scale_colour_gradient(low = "black", high = "darkgrey")+
    theme(legend.justification = c(1, -0.01), legend.position = c(1, 0))+ 
    labs(size="-log10 ( p-value blood )", 
         colour="-log10 ( p-value LCL cell line )")+
    theme(text = element_text(size=24)) 
  
  ggsave(paste0('figures/ASE_eQTL_comparison_ours_vs_LCL.',name,'.fdr05.pdf'), width=12, height=12)
  ######  
  
}

write.table(numbers_and_pvalues_results, file='LCL_ASE_comparison.txt',quote=F,
        sep='\t', row.names=F)
