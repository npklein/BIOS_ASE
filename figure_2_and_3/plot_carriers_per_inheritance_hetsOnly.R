

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

carriers_per_inheritance <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_inheritance/carriers_per_inheritance.hetsOnly.txt')
colnames(carriers_per_inheritance) <- c('inheritance','type','impact','gene','snp','sample','logFC','MAF')
carriers_per_inheritance <- carriers_per_inheritance[!is.na(carriers_per_inheritance$logFC),]

outliers_per_inheritance <- carriers_per_inheritance %>%
  group_by(inheritance,impact) %>%
  summarise(outlier = sum(type=='outlier'),
            not_outlier = sum(type=='not_outlier'),
            n=n())

outliers_per_inheritance$fraction_outlier <- outliers_per_inheritance$outlier/(outliers_per_inheritance$outlier+outliers_per_inheritance$not_outlier)

outliers_per_inheritance$impact <- factor(outliers_per_inheritance$impact, levels =c('HIGH','MODERATE','LOW','MODIFIER'))
outliers_per_inheritance <- outliers_per_inheritance[!outliers_per_inheritance$impact=='MODIFIER',]
ggplot(outliers_per_inheritance, aes(impact, fraction_outlier, fill=impact))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  facet_wrap(~inheritance,nrow=2)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('HIGH','MODERATE','LOW'))+
  geom_text(aes(label=n), vjust=-1)+
  scale_y_continuous(limit=c(0,0.52))
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/proportion_outlier_per_inheritance.png',width=25, height=20)




