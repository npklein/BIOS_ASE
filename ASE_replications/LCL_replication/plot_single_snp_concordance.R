
library(ggplot2)

##### read in data ####
input_dir <- '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'
input_file <- paste0(input_dir,'/ASE.mergedWithLCL.txt')

##### read in data ####
LCL_and_ASE <- read.table(input_file,header=T,sep='\t')
LCL_and_ASE$LCL_ratio <- LCL_and_ASE$LCL_ratio-0.5
LCL_and_ASE$fdr <- p.adjust(LCL_and_ASE$binom_pval, method='fdr')

#####


##### calculate concordance and correlation using ours VS LCL ASE with FDR pvalues < 0.05 #####
# only have to select on FDR for our data, because LCL data is already filtered to only be FDR < 0.05
LCL_and_ASE_fdr05 <- LCL_and_ASE[LCL_and_ASE$fdr < 0.05,]
concordance_LCL_logFC <- sum(sign(LCL_and_ASE_fdr05$logFC)==sign(LCL_and_ASE_fdr05$LCL_ratio))/nrow(LCL_and_ASE_fdr05)
cor_logFC_LCL <- cor(LCL_and_ASE_fdr05$logFC, LCL_and_ASE_fdr05$LCL_ratio, method='spearman')

# TODO: write these values to a file
numbers_and_pvalues_results <- list()
numbers_and_pvalues_results$n_overlap_ASE_LCL <- nrow(LCL_and_ASE_fdr05)
numbers_and_pvalues_results$concordance_with_LCL <- concordance_LCL_logFC*100
numbers_and_pvalues_results$correlation_with_LCL <- cor_logFC_LCL

ggplot(LCL_and_ASE_fdr05, aes(logFC, LCL_ratio, size=-log(fdr+0.00000000001, base=10),
                         colour=-log10(LCL_pval+0.000000000001)))+
  geom_point(alpha=0.5,shape=16)+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_classic(base_size = 20)+
  ylab('Allelic ratio LCL cell line')+
  xlab('logFC blood')+
  annotate("text", x = -3, y = 0.45, label = paste0("SNPs: ",nrow(LCL_and_ASE_fdr05)),
           size=8)+
  annotate("text", x = -3, y = 0.4, label = paste0("Concordance: ",signif(concordance_LCL_logFC,3) ),
           size=8,parse=TRUE)+
  annotate("text", x = -3, y = 0.35, label = paste0("Correlation: ",signif(cor_logFC_LCL,3) ),
           size=8,parse=TRUE)+
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  scale_colour_gradient(low = "black", high = "darkgrey")+
  theme(legend.justification = c(1, -0.01), legend.position = c(1, 0))+ 
  labs(size="-log10 ( p-value blood )", 
       colour="-log10 ( p-value LCL cell line )")+
  theme(text = element_text(size=24)) 

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_eQTL_comparison_ours_vs_LCL.fdr05.png', width=12, height=12)
######  

###### at differnet p-values #####
cor_at_different_pvalue <- data.frame('pval'=c(),'correlation'=c(),'concodrance'=c())
for(pval in c(0.05, 1e-02,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,1e-10,1e-11,1e-12,1e-13,1e-14,1e-15, 1e-16, 1e-17,1e-18,1e-19,1e-20)){
  LCL_and_ASE_subset_pval <- LCL_and_ASE[LCL_and_ASE$fdr < pval & LCL_and_ASE$LCL_pval < pval,]
  cor_current_pval <- cor(LCL_and_ASE_subset_pval$logFC, LCL_and_ASE_subset_pval$LCL_ratio)
  concordance_current_pval <- sum(sign(LCL_and_ASE_subset_pval$logFC)==sign(LCL_and_ASE_subset_pval$LCL_ratio))/nrow(LCL_and_ASE_subset_pval)
  cor_at_different_pvalue <- rbind(cor_at_different_pvalue, data.frame('pval'=pval,
                                                                       'correlation'=cor_current_pval,
                                                                       'concordance'=concordance_current_pval))
}
ggplot(cor_at_different_pvalue, aes(-log10(pval), concordance))+
  geom_point()+
  geom_line()+
  theme_bw(base_size=18)+
  xlab('-log10 ( p-value )')+
  ylab('concordance')
ggsave(paste0('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_LCL_concordance_different_pvalues.png'), width=8, height=8)
ggplot(cor_at_different_pvalue, aes(-log10(pval), correlation))+
  geom_point()+
  geom_line()+
  theme_bw(base_size=18)+
  xlab('-log10 ( p-value )')+
  ylab('Correlation')
ggsave(paste0('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_LCL_correlation_different_pvalues.png'), width=8, height=8)
#####


write.table(numbers_and_pvalues_results, file='LCL_ASE_comparison.txt',quote=F,
        sep='\t', row.names=F)
