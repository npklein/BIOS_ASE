

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

#### Read in file with the ref/alt counts ####
ase_counts <- fread('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.txt')


ase_counts$Pval <- apply(ase_counts, 1, function(x) {
    binom.test(as.numeric(x['Ref_Counts']),
               as.numeric(x['Ref_Counts'])+as.numeric(x['Alt_Counts']), 
               alternative="two.sided", 
               p=0.5, 
               conf.level=0.95)$p.value
})

ase_counts$Bonf_corr_Pval <- p.adjust(ase_counts$Pval, method='bonferroni')
ase_counts$FDR <- p.adjust(ase_counts$Pval, method='fdr')

write.table(ase_counts, '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.binom.txt',
            sep='\t',
            row.names=F)

