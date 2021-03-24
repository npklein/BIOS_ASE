#TODO: write numbers to a file

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

outliers_per_clinvar_at_least_1 %>% group_by(clinstat) %>% summarize(outlier=sum(outlier_bonf), not_outlier=sum(not_outlier_bonf), median=median(fraction_outlier)*100)
ggplot(outliers_per_clinvar_at_least_1, aes(clinstat, fraction_outlier))+
  geom_violin()

ggplot(total, aes(clinstat, fraction_outlier, fill=clinstat))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('Pathogenic',
                           'VUS','Benign'))+
  scale_y_continuous(limit=c(0,1.1))+
  guides(fill=FALSE)

outfile='/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig2//proportion_outlier_per_clinvar.png'
print(paste('write to:',outfile))
ggsave(outfile,width=8, height=8)

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
  scale_y_continuous(limit=c(0,1.1))+
  guides(fill=FALSE)

outfile='/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig2//proportion_outlier_per_vkgl.png'
print(paste('write to:',outfile))
ggsave(outfile,width=8, height=8)
