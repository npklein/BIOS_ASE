

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

carriers_per_disease <- fread('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt')
carriers_per_disease <- carriers_per_disease[!is.na(carriers_per_disease$logFC),]
carriers_per_disease$sample <- NULL

snpInfo <- fread('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt')
snpInfo$snp <- paste0(sapply(strsplit(snpInfo$VARIANT, "_"), "[[", 1),'_',sapply(strsplit(snpInfo$VARIANT, "_"), "[[", 2))
snpInfo$CGDAGEGROUP <- NULL
snpInfo$CGDINHERITANCE <- NULL
snpInfo$CGDMANIFESTATION <- NULL
snpInfo$GENEID <- NULL
snpInfo$EXACAF <- NULL
snpInfo$GONLAF <- NULL
snpInfo$REF <- NULL
snpInfo$ALT <- NULL


carriers_per_disease_snpInfo <- merge(carriers_per_disease, snpInfo, by='snp')
carriers_per_disease_snpInfo$GONLAF <- NULL
rm(carriers_per_disease)
gc()

give.n <- function(x){
  return(c(y = -0.1, label = length(x)))
}

###### outliers per impact factor #####
carriers_per_disease_snpInfo[carriers_per_disease_snpInfo$disease=="other",]$disease <- "General"
outliers_per_disease <- carriers_per_disease_snpInfo %>%
  group_by(disease,impact) %>%
   summarise(outlier = sum(type=='outlier'),
             not_outlier = sum(type=='not_outlier'),
             n=n())

outliers_per_disease$fraction_outlier <- outliers_per_disease$outlier/(outliers_per_disease$outlier+outliers_per_disease$not_outlier)

outliers_per_disease$impact <- factor(outliers_per_disease$impact, levels =c('HIGH','MODERATE','LOW','MODIFIER'))
outliers_per_disease <- outliers_per_disease[!outliers_per_disease$impact=='MODIFIER',]
ggplot(outliers_per_disease, aes(impact, fraction_outlier, fill=impact))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 15)+
  facet_wrap(~disease)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('HIGH','MODERATE','LOW'))+
  geom_text(aes(label=n), vjust=-1)+
  scale_y_continuous(limit=c(0,0.55))
ggsave('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_disease.pdf',width=15, height=10)


##### outliers per categorie
outliers_per_category <- carriers_per_disease_snpInfo %>%
  group_by(SNPEFFANNOTATION,impact) %>%
  summarise(outlier = sum(type=='outlier'),
            not_outlier = sum(type=='not_outlier'),
            n=n())
outliers_per_category$fraction_outlier <- outliers_per_category$outlier/(outliers_per_category$outlier+outliers_per_category$not_outlier)

print(outliers_per_category)
ggplot(outliers_per_category, aes(impact, fraction_outlier, fill=impact))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  facet_wrap(~SNPEFFANNOTATION)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_x_discrete(limit=c('HIGH','MODERATE','LOW','MODIFIER'))+
  geom_text(aes(label=n), vjust=-1)+
  scale_y_continuous(limit=c(0,0.60))
ggsave('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_type.pdf',width=25, height=20)
#####

print('saved at /groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_disease.')
print('saved at /groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_type.pdf')
