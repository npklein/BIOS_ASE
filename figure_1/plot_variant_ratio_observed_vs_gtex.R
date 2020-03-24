library(ggplot2)
library(data.table)
library(reshape2)
library(ggforce)
library(ggpubr)
library(ggsignif)
library(ggrepel)
library(gridExtra)

BIOS_ASE_counts <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/concordanceGTEx/counts.matrix.AlleleAdded.txt'))

# at least 30 samples in bios
BIOS_ASE_counts <- BIOS_ASE_counts[BIOS_ASE_counts$SAMPLECOUNT >= 30,]
BIOS_ASE_counts <- BIOS_ASE_counts[BIOS_ASE_counts$GTEXSAMPLECOUNT >= 30,]

BIOS_ASE_counts_deduplicated <- BIOS_ASE_counts[!duplicated(BIOS_ASE_counts$VARIANT),]

BIOS_ASE_counts_deduplicated$BIOS_ASE_binom <- apply(BIOS_ASE_counts_deduplicated, 1, function(x) binom.test(as.numeric(x[["SUMMAJOR"]]), as.numeric(x[["TOTAL"]]), 
                                                      p=0.5, alternative = "two.sided", conf.level = 0.95)$p.value)


BIOS_ASE_counts_deduplicated$BIOS_ASE_binom_FDR <- p.adjust(BIOS_ASE_counts_deduplicated$BIOS_ASE_binom, method = "BH")

BIOS_ASE_counts_deduplicated[,c('TISSUE','GTEXSAMPLECOUNT', 'GTEXSUMMAJOR','GTEXSUMMINOR','GTEXTOTAL','GTEXRATIO')] <- NULL
BIOS_ASE_counts_with_fdr <- merge(BIOS_ASE_counts, BIOS_ASE_counts_deduplicated, 
                                  by=colnames(BIOS_ASE_counts_deduplicated)[!colnames(BIOS_ASE_counts_deduplicated) %in% c('BIOS_ASE_binom', 
                                                                                             'BIOS_ASE_binom_FDR')],
                                  all=T)

BIOS_ASE_counts_with_fdr$GTEX_binom <- apply(BIOS_ASE_counts_with_fdr, 1, function(x) binom.test(as.numeric(x[["GTEXSUMMAJOR"]]), as.numeric(x[["GTEXTOTAL"]]), 
                                                                                             p=0.5, alternative = "two.sided", conf.level = 0.95)$p.value)

BIOS_ASE_counts_with_fdr$GTEX_fdr <- p.adjust(BIOS_ASE_counts_with_fdr$GTEX_binom, method = "BH")


BIOS_ASE_counts_with_fdr$RATIO_directed <- -1*(0.5-BIOS_ASE_counts_with_fdr$RATIO)
BIOS_ASE_counts_with_fdr$GTEXRATIO_directed <- -1*(0.5-BIOS_ASE_counts_with_fdr$GTEXRATIO)

BIOS_ASE_counts_sign <- BIOS_ASE_counts_with_fdr[BIOS_ASE_counts_with_fdr$GTEX_fdr < 0.05 & BIOS_ASE_counts_with_fdr$BIOS_ASE_binom_FDR < 0.05,]
p1 <- ggplot(BIOS_ASE_counts_sign[sign(BIOS_ASE_counts_sign$RATIO_directed)==sign(BIOS_ASE_counts_sign$GTEXRATIO_directed),], 
             aes(GTEXRATIO_directed, RATIO_directed))+
  geom_point(data=BIOS_ASE_counts_sign[sign(BIOS_ASE_counts_sign$RATIO_directed)!=sign(BIOS_ASE_counts_sign$GTEXRATIO_directed),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('ASE ratio BIOS')+
  xlab('ASE ratio GTEx')+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  facet_wrap(~TISSUE,ncol=6)+
  stat_cor(
    aes(label = paste(..rr.label..,  sep = "~`,`~")), 
    method = "spearman")

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/variant_ratio_observed_vs_gtex_all_tissue.pdf', width=24, height=24, plot=p1)

BIOS_ASE_counts_blood <- BIOS_ASE_counts_sign[BIOS_ASE_counts_sign$TISSUE=="WHLBLD",]
concordance <- sum(sign(BIOS_ASE_counts_blood$RATIO_directed)==sign(BIOS_ASE_counts_blood$GTEXRATIO_directed))/nrow(BIOS_ASE_counts_blood)
cor_zscore <- cor(BIOS_ASE_counts_blood$RATIO_directed, BIOS_ASE_counts_blood$GTEXRATIO_directed, method='spearman')
p2 <- ggplot(BIOS_ASE_counts_blood[sign(BIOS_ASE_counts_blood$RATIO_directed)==sign(BIOS_ASE_counts_blood$GTEXRATIO_directed),], 
             aes(GTEXRATIO_directed, RATIO_directed))+
  geom_point(data=BIOS_ASE_counts_blood[sign(BIOS_ASE_counts_blood$RATIO_directed)!=sign(BIOS_ASE_counts_blood$GTEXRATIO_directed),],
             alpha=0.5,colour='#fdbb84', colour='#fdbb84')+
  geom_point(alpha=0.5,colour='#9ecae1',colour='#9ecae1')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 18)+
  ylab('ASE ratio BIOS')+
  xlab('ASE ratio GTEx')+
  annotate("text", x = -0.1, y = 0.3, label = paste0("Concordance = ",signif(concordance,3),
                                                     "\nCorrelation = ",signif(cor_zscore,3)),
           size=5, hjust=1)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_d/variant_ratio_observed_vs_gtex_blood.pdf', width=6, height=6, plot=p2)

