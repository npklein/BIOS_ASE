
library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)
library(ggbeeswarm)

outliers_per_clinvar <- fread('clinvar_overlapped_with_SNPs.txt')
table(outliers_per_clinvar[which(outliers_per_clinvar$outlier+outliers_per_clinvar$not_outlier==0),]$clinstat)

outliers_per_clinvar_at_least_1 <- outliers_per_clinvar[which(outliers_per_clinvar$outlier+outliers_per_clinvar$not_outlier>0),]
outliers_per_clinvar_at_least_1$fraction_outlier <- outliers_per_clinvar_at_least_1$outlier / (outliers_per_clinvar_at_least_1$not_outlier + outliers_per_clinvar_at_least_1$outlier)

total <- data.frame(clinstat=c('Benign','VUS','Pathogenic'),
                    outlier=c(sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='Benign',]$outlier_bonf),
                              sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='VUS',]$outlier_bonf),
                              sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='Pathogenic',]$outlier_bonf)),
                    not_outlier=c(sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='Benign',]$not_outlier_bonf),
                                  sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='VUS',]$not_outlier_bonf),
                                  sum(outliers_per_clinvar_at_least_1[outliers_per_clinvar_at_least_1$clinstat=='Pathogenic',]$not_outlier_bonf)))
total$fraction_outlier <- (total$outlier/(total$outlier+total$not_outlier))


ggplot(total, aes(clinstat, fraction_outlier, fill=clinstat))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('Pathogenic',
                           'VUS','Benign'))+
  guides(fill=FALSE)
outfile='/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig2//proportion_outlier_per_clinvar.png'
print(paste('write to:',outfile))
ggsave(outfile,width=8, height=8)

test_matrix <- data.frame(outlier=c(total[total$clinstat=='Pathogenic',]$outlier, 
                                    total[total$clinstat=='Benign',]$outlier),
                          not_outlier=c(total[total$clinstat=='Pathogenic',]$not_outlier, 
                                       total[total$clinstat=='Benign',]$not_outlier))
rownames(test_matrix) <- c('pathogenic','benign')
p_benign_path <- fisher.test(test_matrix)$p.value

test_matrix <- data.frame(outlier=c(total[total$clinstat=='Pathogenic',]$outlier, 
                                    total[total$clinstat=='VUS',]$outlier),
                          not_outlier=c(total[total$clinstat=='Pathogenic',]$not_outlier, 
                                        total[total$clinstat=='VUS',]$not_outlier))
rownames(test_matrix) <- c('pathogenic','VUS')
p_vus_path <- fisher.test(test_matrix)$p.value

test_matrix <- data.frame(outlier=c(total[total$clinstat=='Benign',]$outlier, 
                                    total[total$clinstat=='VUS',]$outlier),
                          not_outlier=c(total[total$clinstat=='Benign',]$not_outlier, 
                                        total[total$clinstat=='VUS',]$not_outlier))
rownames(test_matrix) <- c('pathogenic','VUS')
p_vus_benign <- fisher.test(test_matrix)$p.value

p_benign_path
p_vus_path
p_vus_benign

#####


##### VKGL #####
outliers_per_vkgl <- fread('vkgl_overlapped_with_SNPs.txt')
table(outliers_per_vkgl[which(outliers_per_vkgl $outlier+outliers_per_vkgl$not_outlier>0),]$clinstat)

outliers_per_vkgl_at_least_1 <- outliers_per_vkgl[which(outliers_per_vkgl$outlier+outliers_per_vkgl$not_outlier>0),]
outliers_per_vkgl$fraction_outlier <- outliers_per_vkgl_at_least_1$outlier / (outliers_per_vkgl_at_least_1$not_outlier + outliers_per_vkgl_at_least_1$outlier)

total_vkgl <- data.frame(clinstat=c('Benign','VUS','Pathogenic'),
                         outlier=c(sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='Benign',]$outlier_bonf),
                                   sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='VUS',]$outlier_bonf),
                                   sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='Pathogenic',]$outlier_bonf)),
                         not_outlier=c(sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='Benign',]$not_outlier_bonf),
                                       sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='VUS',]$not_outlier_bonf),
                                       sum(outliers_per_vkgl_at_least_1[outliers_per_vkgl_at_least_1$clinstat=='Pathogenic',]$not_outlier_bonf)))
total_vkgl$fraction_outlier <- (total_vkgl$outlier/(total_vkgl$outlier+total_vkgl$not_outlier))


ggplot(total_vkgl, aes(clinstat, fraction_outlier, fill=clinstat))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('Pathogenic','Benign'))+
  guides(fill=FALSE)

outfile='/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig2//proportion_outlier_per_vkgl.png'
print(paste('write to:',outfile))
ggsave(outfile,width=8, height=8)




test_matrix <- data.frame(outlier=c(total_vkgl[total_vkgl$clinstat=='Benign',]$outlier, 
                                    total_vkgl[total_vkgl$clinstat=='Pathogenic',]$outlier),
                          not_outlier=c(total_vkgl[total_vkgl$clinstat=='Benign',]$not_outlier, 
                                        total_vkgl[total_vkgl$clinstat=='Pathogenic',]$not_outlier))
rownames(test_matrix) <- c('Benign','Pathogenic')
p_benign_path <- fisher.test(test_matrix)$p.value
