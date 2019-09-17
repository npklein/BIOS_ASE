#TODO: write numbers to a file

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

carriers_per_disease <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt')

carriers_per_disease <- carriers_per_disease[!is.na(carriers_per_disease$logFC),]
snpInfo <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.chr22.addedCADD.addedVKGL.txt')
snpInfo$snp <- paste0(sapply(strsplit(snpInfo$VARIANT, "_"), "[[", 1),'_',sapply(strsplit(snpInfo$VARIANT, "_"), "[[", 2))
carriers_per_disease_snpInfo <- merge(carriers_per_disease, snpInfo, by='snp')
clinvarTable <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvarSNPs_2019-Feb-19.txt')

carriers_per_disease_clinvarInfo <- merge(carriers_per_disease_snpInfo, clinvarTable, by='snp', fill=T)

# filter: keep only those SNPs that are high confidence. 
carriers_per_disease_snpInfo <- carriers_per_disease_clinvarInfo[which(carriers_per_disease_clinvarInfo$Selected == "TRUE"),]

carriers_per_disease_snpInfo_noDuplicates <- carriers_per_disease_snpInfo[!duplicated(carriers_per_disease_snpInfo$snp),]
print(table(carriers_per_disease_snpInfo_noDuplicates$ClinSimple))


give.n <- function(x){
  return(c(y = -0.1, label = length(x)))
}


###### outliers per impact factor #####
carriers_per_disease_snpInfo[grepl('Benign',carriers_per_disease_snpInfo$CLNVRSIG,ignore.case=TRUE),]$CLNVRSIG <- 'Benign'
carriers_per_disease_snpInfo[carriers_per_disease_snpInfo$CLNVRSIG=='Pathogenic,_other',]$CLNVRSIG <- 'Pathogenic'
carriers_per_disease_snpInfo[carriers_per_disease_snpInfo$CLNVRSIG=='Pathogenic/Likely_pathogenic',]$CLNVRSIG <- 'Pathogenic'
carriers_per_disease_snpInfo[carriers_per_disease_snpInfo$disease=='other',]$disease <- 'General'
carriers_per_disease_snpInfo <- carriers_per_disease_snpInfo[!carriers_per_disease_snpInfo$disease=='General',]
outliers_per_disease <- carriers_per_disease_snpInfo %>%
  group_by(disease,CLNVRSIG) %>%
   summarise(outlier = sum(type=='outlier'),
             not_outlier = sum(type=='not_outlier'),
             n=n())

outliers_per_disease$fraction_outlier <- outliers_per_disease$outlier/(outliers_per_disease$outlier+outliers_per_disease$not_outlier)

outliers_per_disease$CLNVRSIG <- factor(outliers_per_disease$CLNVRSIG, levels=c('Pathogenic',
                                                                                'Uncertain_significance','Benign'))
print(outliers_per_disease)
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

##### make one overview plot ####
outliers_total <- carriers_per_disease_snpInfo %>%
  group_by(CLNVRSIG) %>%
  summarise(outlier = sum(type=='outlier'),
            not_outlier = sum(type=='not_outlier'),
            n=n())
outliers_total$fraction_outlier <- outliers_total$outlier/(outliers_total$outlier+outliers_total$not_outlier)
print(outliers_total)
outliers_total$CLNVRSIG <- factor(outliers_total$CLNVRSIG, levels=c('Pathogenic',
                                                                                'Uncertain_significance','Benign'))
print(outliers_total)
ggplot(outliers_total, aes(CLNVRSIG, fraction_outlier, fill=CLNVRSIG))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 18)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of SNPs that is an outlier')+
  xlab('')+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label=n), vjust=-1)+
  scale_x_discrete(limit=c('Pathogenic',
                           'Uncertain_significance','Benign'))+
  scale_y_continuous(limit=c(0,1.1))+
  guides(fill=FALSE)

outfile='/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig2//proportion_outlier_per_clinvar.png'
print(paste('write to:',outfile))
ggsave(outfile,width=12, height=8)

#####

## Leaving in below temporarily, will remove if not needed later
##### clinvar more trusted
#clinvar_outliers_per_disease <- carriers_per_disease_clinvarInfo %>%
#  group_by(disease,ClinSimple) %>%
#  summarise(outlier = sum(type=='outlier'),
#            not_outlier = sum(type=='not_outlier'),
#            n=n())

#clinvar_outliers_per_disease$fraction_outlier <- clinvar_outliers_per_disease$outlier/(clinvar_outliers_per_disease$outlier+clinvar_outliers_per_disease$not_outlier)


#ggplot(clinvar_outliers_per_disease[clinvar_outliers_per_disease$ClinSimple!='not provided',], aes(ClinSimple, fraction_outlier, fill=ClinSimple))+
#  geom_bar(stat='identity')+
#  theme_bw(base_size = 18)+
#  facet_wrap(~disease)+
#  scale_fill_brewer(palette="Dark2")+
#  ylab('Fraction of SNPs that is an outlier')+
#  xlab('')+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#  geom_text(aes(label=n), vjust=-1)+
#  scale_y_continuous(limit=c(0,0.80))
#ggsave('figures/proportion_outlier_per_clinvar_trusted.png',width=25, height=20)



#####


