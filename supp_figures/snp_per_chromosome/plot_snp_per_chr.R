library(data.table)
library(ggplot2)

maf_per_chr <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.final.3810samples.txt')

colnames(maf_per_chr) <- c('snp','maf')
maf_per_chr$chr <- sapply(strsplit(as.character(maf_per_chr$snp), "_"), "[[", 1)
maf_per_chr[maf_per_chr$maf>0.5,]$maf <- 1-maf_per_chr[maf_per_chr$maf>0.5,]$maf



ggplot(maf_per_chr, aes(maf))+
  geom_histogram(bins=20)+
  theme_bw(base_size=18)+
  xlab('Minor Allele Frequency')+
  ylab('Number of SNPs')
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl//MAFdistribution.png',width=8, height=8)


snps_per_chr <- data.frame(table(maf_per_chr$chr))
colnames(snps_per_chr) <- c('chr','nSNPs')
snps_per_chr_levels <- snps_per_chr[order(as.numeric(as.character(snps_per_chr$chr))),]$chr
snps_per_chr$chr <- factor(snps_per_chr$chr, levels = snps_per_chr_levels)

ggplot(snps_per_chr, aes(chr, nSNPs))+
  geom_bar(stat='identity')+
  theme_bw(base_size=18)+
  xlab('Chromosome')+
  ylab('Number of SNPs')
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/SNPs_per_chromosome.png',width=9, height=8)
