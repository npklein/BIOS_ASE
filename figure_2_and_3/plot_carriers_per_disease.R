

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

#carriers_per_disease <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt')
carriers_per_disease <- fread('carriers_per_disease.txt.gz')

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
  scale_y_continuous(limit=c(0,1.3))
#ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/proportion_outlier_per_disease.png',width=15, height=10)
ggsave('proportion_outlier_per_disease.png',width=15, height=10)


x <- 0
for(disease in unique(outliers_per_disease$disease)){
  df <- outliers_per_disease[outliers_per_disease$disease==disease,]

  if('HIGH' %in% df$impact){
    test_matrix <- data.frame(outlier=c(df[df$impact=='HIGH',]$outlier, 
                                        df[df$impact=='LOW',]$outlier),
                              not_outlier=c(df[df$impact=='HIGH',]$not_outlier, 
                                            df[df$impact=='LOW',]$not_outlier))
    rownames(test_matrix) <- c('HIGH','LOW')
    p_high_low <- fisher.test(test_matrix)$p.value*36
  
    test_matrix <- data.frame(outlier=c(df[df$impact=='HIGH',]$outlier, 
                                        df[df$impact=='MODERATE',]$outlier),
                              not_outlier=c(df[df$impact=='HIGH',]$not_outlier, 
                                            df[df$impact=='MODERATE',]$not_outlier))
    rownames(test_matrix) <- c('HIGH','MODERATE')
    p_high_modarate <- fisher.test(test_matrix)$p.value*36
    x <- x+1
  }
  test_matrix <- data.frame(outlier=c(df[df$impact=='LOW',]$outlier, 
                                      df[df$impact=='MODERATE',]$outlier),
                            not_outlier=c(df[df$impact=='LOW',]$not_outlier, 
                                          df[df$impact=='MODERATE',]$not_outlier))
  rownames(test_matrix) <- c('LOW','MODERATE')
  p_high_modarate <- fisher.test(test_matrix)$p.value*36
  x <- x+1
  if('HIGH' %in% df$impact){
    cat(paste(disease,signif(p_high_modarate,2),signif(p_high_low,2), signif(p_high_modarate,2),'\n'))
  }else{
    cat(paste(disease, signif(p_high_modarate,2),'\n'))
  }
  
}


print(x)



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
  scale_y_continuous(limit=c(0,1))
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/proportion_outlier_per_type.png',width=25, height=20)
#####

print('saved at /groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_disease.')
print('saved at /groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/suppl/proportion_outlier_per_type.pdf')
