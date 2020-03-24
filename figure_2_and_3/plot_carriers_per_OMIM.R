#TODO: write numbers to a file

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

carriers_per_disease <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt')


carriers_per_disease <- carriers_per_disease[!is.na(carriers_per_disease$logFC),]


give.n <- function(x){
  return(c(y = -0.1, label = length(x)))
}

carriers_per_disease_snpInfo[carriers_per_disease_snpInfo$disease=='other',]$disease <- 'General'
carriers_per_disease_snpInfo <- carriers_per_disease_snpInfo[!carriers_per_disease_snpInfo$disease=='General',]
outliers_per_disease <- carriers_per_disease_snpInfo %>%
  group_by(disease,CLNVRSIG) %>%
   summarise(outlier = sum(type=='outlier'),
             not_outlier = sum(type=='not_outlier'),
             n=n())

outliers_per_disease$fraction_outlier <- outliers_per_disease$outlier/(outliers_per_disease$outlier+outliers_per_disease$not_outlier)

ggplot(outliers_per_disease, aes(CLNVRSIG, fraction_outlier, fill=CLNVRSIG))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  facet_wrap(~disease, nrow=3)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label=n), vjust=-1)+
  scale_x_discrete(limit=c('Pathogenic',
                           'Uncertain_significance','Benign'))+
  scale_y_continuous(limit=c(0,1.1))

outfile = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig3//proportion_outlier_per_clinvar_per_disease.png'
print(paste('write to:',outfile))
ggsave(outfile,width=25, height=20)





