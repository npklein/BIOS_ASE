#TODO: write numbers to a file

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

carriers_per_disease <- fread('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE/clinvar/clinvar_overlapped_with_SNPs.txt')




give.n <- function(x){
  return(c(y = -0.1, label = length(x)))
}

carriers_per_disease[carriers_per_disease$disease=='other',]$disease <- 'General'
carriers_per_disease <- carriers_per_disease[!carriers_per_disease$disease=='General',]

outliers_per_disease <- carriers_per_disease %>%
  group_by(disease,clinstat) %>%
   summarise(outlier = sum(type=='outlier'),
             not_outlier = sum(type=='not_outlier'),
             n=length(unique(snp)))

outliers_per_disease$fraction_outlier <- outliers_per_disease$outlier/(outliers_per_disease$outlier+outliers_per_disease$not_outlier)

outliers_per_disease$clinstat <- factor(outliers_per_disease$clinstat, levels=c('Pathogenic','VUS','Benign'))
print(head(outliers_per_disease))
ggplot(outliers_per_disease, aes(clinstat, fraction_outlier, fill=clinstat))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 24)+
  facet_wrap(~disease, nrow=3)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of samples that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label=n), vjust=-1)+
  scale_y_continuous(limit=c(0,1.1))+
  scale_x_discrete(limit=c('LOW', 'MODERATE','HIGH'), labels='Low','Moderate','High')

outfile = '/groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/fig3//proportion_outlier_per_clinvar_per_disease.pdf'
print(paste('write to:',outfile))
ggsave(outfile,width=25, height=20)





