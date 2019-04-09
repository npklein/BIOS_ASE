
library(ggplot2)
library(data.table)
##### read in data ####
# Possibly use ASE_eQTL.topEqtlOnly.cov30.top10_eQTL.txt instead
ase_eqtl_file <- '/groups/umcg-bios/tmp03/projects/outlierGeneASE/compareASEcounts/compareWithEqtlGen/data/ASE_eQTL.binom.all_eQTL.txt'

eqtl_and_ASE <- fread(ase_eqtl_file, header=T,sep='\t')

eqtl_and_ASE_full <- eqtl_and_ASE[!is.na(eqtl_and_ASE$logFC),]
eqtl_and_ASE_full$binom_fdr <- p.adjust(eqtl_and_ASE_full$binom, method='fdr')

#####

###### fdr 0.05 ######
eqtl_and_ASE_fdr_05 <- eqtl_and_ASE_full[eqtl_and_ASE_full$binom_fdr < 0.05 & eqtl_and_ASE_full$fdr_eqtlGen < 0.05,]
concordance_fdr <- sum(sign(eqtl_and_ASE_fdr_05$zScore_alt_swapped)==sign(eqtl_and_ASE_fdr_05$logFC))/nrow(eqtl_and_ASE_fdr_05)
cor_zscore_logFC_fdr_05 <- cor(eqtl_and_ASE_fdr_05$zScore_alt_swapped, eqtl_and_ASE_fdr_05$logFC, method='spearman')
ggplot(eqtl_and_ASE_fdr_05, aes(zScore_alt_swapped, logFC, size=-log(binom_fdr+0.0000000001, base=10),
                                                                     colour=-log10(fdr_eqtlGen+0.0000000001)))+
  geom_point(alpha=0.5,shape=16)+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_classic(base_size = 30)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -60, y = 1, label = paste0("Concordance: ",signif(concordance_fdr,3) ),
           size=8,parse=TRUE)+
  annotate("text", x = -60, y = 0.9, label = paste0("Correlation: ",signif(cor_zscore_logFC_fdr_05,3) ),
           size=8,parse=TRUE)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  scale_colour_gradient(low = "black", high = "darkgrey")+
  labs(size="-log10 ( p-value ASE )", 
       colour="-log10 ( p-value eQTL )")+
  theme(legend.justification = c(1, -0.01), legend.position = c(1, 0))+ 
  labs(size="-log10 ( p-value ASE )", 
       colour="-log10 ( p-value eQTL )")
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_eQTLgen_top10eQTL_ASEfdr05_topSNPs.pdf', width=14, height=14)
######

###### rank 1 #####
eqtl_and_ASE_fdr_05_rank1 <- eqtl_and_ASE_fdr_05[eqtl_and_ASE_fdr_05$eqtl_rank==0,]
concordance_fdr_rank1 <- sum(sign(eqtl_and_ASE_fdr_05_rank1$zScore_alt_swapped)==sign(eqtl_and_ASE_fdr_05_rank1$logFC))/nrow(eqtl_and_ASE_fdr_05_rank1)
cor_zscore_logFC_fdr_05_rank1 <- cor(eqtl_and_ASE_fdr_05_rank1$zScore_alt_swapped, eqtl_and_ASE_fdr_05_rank1$logFC, method='spearman')
ggplot(eqtl_and_ASE_fdr_05_rank1, aes(zScore_alt_swapped, logFC, size=-log(binom_fdr+0.0000000001, base=10),
                                colour=-log10(fdr_eqtlGen+0.0000000001)))+
  geom_point(alpha=0.5,shape=16)+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_classic(base_size = 30)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -60, y = 1, label = paste0("Concordance: ",signif(concordance_fdr_rank1,3) ),
           size=8,parse=TRUE)+
  annotate("text", x = -60, y = 0.9, label = paste0("Correlation: ",signif(cor_zscore_logFC_fdr_05_rank1,3) ),
           size=8,parse=TRUE)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  scale_colour_gradient(low = "black", high = "darkgrey")+
  labs(size="-log10 ( p-value ASE )", 
       colour="-log10 ( p-value eQTL )")+
  theme(legend.justification = c(1, -0.01), legend.position = c(1, 0))+ 
  labs(size="-log10 ( p-value ASE )", 
       colour="-log10 ( p-value eQTL )")
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_eQTLgen_top10eQTL_ASEfdr05_topSNPs_rank1.pdf', width=14, height=14)
######

