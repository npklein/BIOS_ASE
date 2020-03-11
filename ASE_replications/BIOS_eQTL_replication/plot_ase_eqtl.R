
library(ggplot2)
library(data.table)
library(patchwork)
##### read in data ####
ase_and_eqtl <- read.table('/groups/umcg-bios/tmp03/projects/2020-03-05-ASE-BIOS-eqtl-comparison/ase_and_eqtl.txt',sep='\t',header=T)

#ase_and_eqtl[ase_and_eqtl$zscore_qtl != ase_and_eqtl$zscoreSwapped_qtl,]
ase_and_eqtl[ase_and_eqtl$ALT!=ase_and_eqtl$minor.allele,c('VARIANT','REF','ALT','SUMMAJOR','SUMMINOR','RATIO',
                                                           'zscore_qtl','zscoreSwapped_qtl','qtl_genotype',
                                                           'qtl_assessed_allele','major.allele', 'minor.allele')]

ase_and_eqtl$bion.pval <- apply(ase_and_eqtl, 1,function(x){
  t <- binom.test(as.numeric(x[['SUMMINOR']]), as.numeric(x[['SUMMAJOR']])+as.numeric(x[['SUMMINOR']]))
  return(t$p.value)
})

ase_and_eqtl$binom.fdr <- p.adjust(ase_and_eqtl$bion.pval,method='BH')
ase_and_eqtl_binom_signific <-ase_and_eqtl[ase_and_eqtl$binom.fdr<0.05,]
ase_and_eqtl_binom_signific_samp30 <- ase_and_eqtl_binom_signific[ase_and_eqtl_binom_signific$SAMPLECOUNT>=30,]
ase_and_eqtl_binom_signific_samp30$RATIO_dir <- -1*(0.5-ase_and_eqtl_binom_signific_samp30$RATIO)

#####

###### fdr 0.05 ######
concordance <- sum(sign(ase_and_eqtl_binom_signific_samp30$zscoreSwapped_qtl)==sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir))/nrow(ase_and_eqtl_binom_signific_samp30)
cor_zscore <- cor(ase_and_eqtl_binom_signific_samp30$zscoreSwapped_qtl, ase_and_eqtl_binom_signific_samp30$RATIO_dir, method='spearman')
p1 <- ggplot(ase_and_eqtl_binom_signific_samp30[sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir)==sign(ase_and_eqtl_binom_signific_samp30$zscoreSwapped_qtl),], 
       aes(zscoreSwapped_qtl, RATIO_dir))+
  geom_point(data=ase_and_eqtl_binom_signific_samp30[sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir)!=sign(ase_and_eqtl_binom_signific_samp30$zscoreSwapped_qtl),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -30, y = 0.39, label = paste0("Concordance: ",signif(concordance,3),
                                                     "\nCorrelation: ",signif(cor_zscore,3)),
           size=4, hjust=1)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  labs(size="-log10 ( p-value ASE )")+
  ggtitle('Swapped for assessed allele')
######
concordance <- sum(sign(ase_and_eqtl_binom_signific_samp30$zscore_qtl)==sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir))/nrow(ase_and_eqtl_binom_signific_samp30)
cor_zscore <- cor(ase_and_eqtl_binom_signific_samp30$zscore_qtl, ase_and_eqtl_binom_signific_samp30$RATIO_dir, method='spearman')
p2 <- ggplot(ase_and_eqtl_binom_signific_samp30[sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir)==sign(ase_and_eqtl_binom_signific_samp30$zscore_qtl),], 
             aes(zscore_qtl, RATIO_dir))+
  geom_point(data=ase_and_eqtl_binom_signific_samp30[sign(ase_and_eqtl_binom_signific_samp30$RATIO_dir)!=sign(ase_and_eqtl_binom_signific_samp30$zscore_qtl),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -60, y = 0.39, label = paste0("Concordance: ",signif(concordance,3),
                                                     "\nCorrelation: ",signif(cor_zscore,3)),
           size=4)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  labs(size="-log10 ( p-value ASE )")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ggtitle('Not swapped')

######
ase_and_eqtl_same_genotype <- ase_and_eqtl_binom_signific_samp30[paste0(ase_and_eqtl_binom_signific_samp30$REF,'/',ase_and_eqtl_binom_signific_samp30$ALT)==ase_and_eqtl_binom_signific_samp30$qtl_genotype,]
concordance <- sum(sign(ase_and_eqtl_same_genotype$zscoreSwapped_qtl)==sign(ase_and_eqtl_same_genotype$RATIO_dir))/nrow(ase_and_eqtl_same_genotype)
cor_zscore <- cor(ase_and_eqtl_same_genotype$zscoreSwapped_qtl, ase_and_eqtl_same_genotype$RATIO_dir, method='spearman')
p3 <- ggplot(ase_and_eqtl_same_genotype[sign(ase_and_eqtl_same_genotype$RATIO_dir)==sign(ase_and_eqtl_same_genotype$zscoreSwapped_qtl),], 
             aes(zscoreSwapped_qtl, RATIO_dir))+
  geom_point(data=ase_and_eqtl_same_genotype[sign(ase_and_eqtl_same_genotype$RATIO_dir)!=sign(ase_and_eqtl_same_genotype$zscoreSwapped_qtl),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -60, y = 0.39, label = paste0("Concordance: ",signif(concordance,3),
                                                     "\nCorrelation: ",signif(cor_zscore,3)),
           size=4)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  labs(size="-log10 ( p-value ASE )")+
  ggtitle('Swapped, same genotype')



######
ase_and_eqtl_different_genotype <- ase_and_eqtl_binom_signific_samp30[paste0(ase_and_eqtl_binom_signific_samp30$REF,
                                                                             '/',
                                                                             ase_and_eqtl_binom_signific_samp30$ALT)!=ase_and_eqtl_binom_signific_samp30$qtl_genotype,]
concordance <- sum(sign(ase_and_eqtl_different_genotype$zscoreSwapped_qtl)==sign(ase_and_eqtl_different_genotype$RATIO_dir))/nrow(ase_and_eqtl_different_genotype)
cor_zscore <- cor(ase_and_eqtl_different_genotype$zscoreSwapped_qtl, ase_and_eqtl_different_genotype$RATIO_dir, method='spearman')
p4 <- ggplot(ase_and_eqtl_different_genotype[sign(ase_and_eqtl_different_genotype$RATIO_dir)==sign(ase_and_eqtl_different_genotype$zscoreSwapped_qtl),], 
             aes(zscoreSwapped_qtl, RATIO_dir))+
  geom_point(data=ase_and_eqtl_different_genotype[sign(ase_and_eqtl_different_genotype$RATIO_dir)!=sign(ase_and_eqtl_different_genotype$zscoreSwapped_qtl),],
             alpha=0.5,shape=21, fill='darkred')+
  geom_point(alpha=0.5,shape=21, fill='darkblue')+
  geom_hline(yintercept=0, lty=2, colour='red')+
  geom_vline(xintercept=0, lty=2, colour='red')+
  theme_bw(base_size = 15)+
  ylab('logFC ASE')+
  xlab('Z-Scores eQTLs')+
  annotate("text", x = -40, y = 0.39, label = paste0("Concordance: ",signif(concordance,3),
                                                     "\nCorrelation: ",signif(cor_zscore,3)),
           size=4)+
  theme(legend.position="top") + 
  geom_smooth(method='lm', formula = y~x, show.legend = FALSE)+
  labs(size="-log10 ( p-value ASE )")+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y= element_blank())+
  ggtitle('Swapped, different genotype')

# Use below plot to visually check if swapping of alleles goes correctly
(p1|p2)/(p3|p4)

# if true, save the swapped figure
p1 <- p1+ggtitle('')
ggsave("ASE-eqtl-comparison-BIOS.pdf", width=6,height=6)





