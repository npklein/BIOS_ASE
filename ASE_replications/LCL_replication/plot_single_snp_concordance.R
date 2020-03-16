
library(ggplot2)

##### read in data ####
input_dir <- '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/merged_count_data/'
#input_dir <- './/'
input_file <- paste0(input_dir,'/ASE.mergedWithLCL.txt')

##### read in data ####
LCL_and_ASE <- read.table(input_file,header=T,sep='\t',row.names=NULL)

LCL_and_ASE_min30 <- LCL_and_ASE[LCL_and_ASE$BIOS_N_samples>=30 & LCL_and_ASE$LCL_N_samples >= 30,]

LCL_and_ASE_min30$LCL_ratio <- LCL_and_ASE_min30$LCL_ratio-0.5
LCL_and_ASE_min30$LCL_fdr <- p.adjust(LCL_and_ASE_min30$LCL_pval, method='BH')
LCL_and_ASE_min30$BIOS_binom <- apply(LCL_and_ASE_min30, 1, function(x) binom.test(as.numeric(x[["minorCount"]]), as.numeric(x[["minorCount"]])+as.numeric(x[["majorCount"]]), 
                                                                                        p=0.5, alternative = "two.sided", conf.level = 0.95)$p.value)
LCL_and_ASE_min30$BIOS_fdr <- p.adjust(LCL_and_ASE_min30$BIOS_binom, method='BH')

LCL_and_ASE_min30$BIOS_ratio <- LCL_and_ASE_min30$minorCount / (LCL_and_ASE_min30$majorCount+LCL_and_ASE_min30$minorCount)
LCL_and_ASE_min30$BIOS_ratio_directed <- -1*(0.5-LCL_and_ASE_min30$BIOS_ratio)
LCL_and_ASE_min30$LCL_ratio_swapped_directed <- -1*(0.5-LCL_and_ASE_min30$LCL_ratio_swapped)
#####

##### calculate concordance and correlation using ours VS LCL ASE with FDR pvalues < 0.05 #####
# only have to select on FDR for our data, because LCL data is already filtered to only be FDR < 0.05
LCL_and_ASE_fdr05 <- LCL_and_ASE_min30[LCL_and_ASE_min30$BIOS_fdr < 0.05 & LCL_and_ASE_min30$LCL_fdr < 0.05,]





concordance_LCL_logFC <- sum(sign(LCL_and_ASE_fdr05$BIOS_ratio_directed)==sign(LCL_and_ASE_fdr05$LCL_ratio_swapped_directed))/nrow(LCL_and_ASE_fdr05)
cor_logFC_LCL <- cor(LCL_and_ASE_fdr05$BIOS_ratio_directed, LCL_and_ASE_fdr05$LCL_ratio_swapped_directed, method='spearman')


ggplot(LCL_and_ASE_fdr05[sign(LCL_and_ASE_fdr05$BIOS_ratio_directed) != LCL_and_ASE_fdr05$LCL_ratio_swapped_directed,], 
       aes(BIOS_ratio_directed, LCL_ratio_swapped_directed))+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_point(data=LCL_and_ASE_fdr05[sign(LCL_and_ASE_fdr05$BIOS_ratio_directed)!=sign(LCL_and_ASE_fdr05$LCL_ratio_swapped_directed),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('Ratio LCL')+
  xlab('Ratio BIOS')+
  annotate("text", x = -0.15, y = 0.25, label = paste0("Concordance: ",signif(concordance_LCL_logFC,3),
                                                     "\nCorrelation: ",signif(cor_logFC_LCL,3)),
           size=4, hjust=1)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  labs(size="-log10 ( p-value ASE )")

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/ASE_comparison_ours_vs_LCL.fdr05.png', width=6, height=6)
######  
