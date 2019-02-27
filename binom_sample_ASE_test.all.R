library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)

#### Read in file with the ref/alt counts ####
ase_counts <- fread('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.all.txt')


ase_counts$Pval <- apply(ase_counts, 1, function(x) {
  binom.test(as.numeric(x['Ref_Counts']),
             as.numeric(x['Ref_Counts'])+as.numeric(x['Alt_Counts']), 
             alternative="two.sided", 
             p=0.5, 
             conf.level=0.95)$p.value
})

ase_counts$Bonf_corr_Pval <- p.adjust(ase_counts$Pval, method='bonferroni')
ase_counts$FDR <- p.adjust(ase_counts$Pval, method='fdr')

ase_counts$Ref_Counts <- paste0("", ase_counts$Ref_Counts, "")
ase_counts$Alt_Counts <- paste0("", ase_counts$Alt_Counts, "")
ase_counts$Chromosome <- paste0("", ase_counts$Chromosome, "")
ase_counts$Position <- paste0("", ase_counts$Position, "")
ase_counts$ID <- paste0("", ase_counts$ID, "")
ase_counts$Pval <- paste0("", ase_counts$Pval, "")
ase_counts$Bonf_corr_Pval <- paste0("", ase_counts$Bonf_corr_Pval, "")
ase_counts$FDR <- paste0("", ase_counts$FDR, "")

write.table(ase_counts, '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.all.binom.txt',
            quote=TRUE,
            sep='\t',
            row.names=F)

