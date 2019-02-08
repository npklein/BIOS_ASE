

library(ggplot2)
library(reshape2)
library(dplyr)
library(data.table)
library(bbmle)

#### Read in file with the ref/alt counts ####
ase_counts_snp_aggregate <- fread('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_ase.txt')
ase_counts_sample <- fread('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_sampleAse.txt')
ase_counts_sample$Alt_Counts <- as.numeric(ase_counts_sample$Alt_Counts)
ase_counts_sample$Ref_Counts <- as.numeric(ase_counts_sample$Ref_Counts)

#### First sum the ref and alt per snp ####
allele_counts_per_snp <- ase_counts_sample %>% 
  group_by(snp_id) %>% 
  summarise(Ref_Counts = sum(Ref_Counts),
            Alt_Counts = sum(Alt_Counts))
allele_counts_per_snp$total <- allele_counts_per_snp$Ref_Counts + allele_counts_per_snp$Alt_Counts

allele_counts_per_snp <- head(allele_counts_per_snp)
####

#### Calculate the log ratio ####
logL_function <- function(p, snp) -sum(dbinom(as.numeric(snp[['Ref_Counts']]),
                                          as.numeric(snp[['total']]),
                                          p, log=T))


logL_sumary <- apply(allele_counts_per_snp, 1, function(x)
                            summary(mle2(logL_function, 
                                 start=list(p=0.5), 
                                 data=list(snp=x))))
# can't get it nicely into the column in a different way, so apply and unlist
allele_counts_per_snp$Likelihood_ratio <- unlist(lapply(logL_sumary, function(x) x@m2logL))
allele_counts_per_snp$Likelihood_ratio_test_Pval <- unlist(lapply(logL_sumary, function(x) x@coef[4]))

allele_counts_per_snp$Likelihood_ratio_test_FDR <- p.adjust(allele_counts_per_snp$Likelihood_ratio_test_Pval, 
                                                            method='fdr')
allele_counts_per_snp$Likelihood_ratio_test_bonferroni <- p.adjust(allele_counts_per_snp$Likelihood_ratio_test_Pval, 
                                                                    method='bonferroni')
####

#### Calculate binomial p-value ####
allele_counts_per_snp$binom_Pval <- apply(allele_counts_per_snp, 1, function(x) {
    binom.test(as.numeric(x['Ref_Counts']),
               as.numeric(x['Ref_Counts'])+as.numeric(x['Alt_Counts']), 
               alternative="two.sided", 
               p=0.5, 
               conf.level=0.95)$p.value
})

allele_counts_per_snp$binom_FDR <- p.adjust(allele_counts_per_snp$binom_Pval,
                                                            method='fdr')
allele_counts_per_snp$binom_bonferroni <- p.adjust(allele_counts_per_snp$binom_Pval,
                                                            method='bonferroni')
####
ase_counts_snp_aggregate <- data.frame(ase_counts_snp_aggregate)
allele_counts_per_snp <- data.frame(allele_counts_per_snp)
ase_counts_snp_aggregate_binom <- merge(ase_counts_snp_aggregate, allele_counts_per_snp, by.x='SNP_ID',by.y='snp_id')

ase_counts_snp_aggregate_binom <- ase_counts_snp_aggregate_binom[c('SNP_ID','Fraction_alternative_allele',
                                                                   'Alternative_allele','Reference_allele',
                                                                   'Samples','Chr','Pos','Likelihood_ratio',
                                                                   'Likelihood_ratio_test_Pval','Likelihood_ratio_test_FDR',
                                                                   'Likelihood_ratio_test_bonferroni','binom_Pval',
                                                                   'binom_FDR','binom_bonferroni')]
write.table(ase_counts_snp_aggregate_binom, '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ase_snpAse.binom.txt',
            sep='\t',
            row.names=F)

