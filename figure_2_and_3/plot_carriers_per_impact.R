#TODO: write numbers to a file

library(data.table)
library(ggplot2)
library(plotly)
require(dplyr)

#carriers_per_disease <- fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/omim_enrichment/carriers_per_disease/carriers_per_disease.txt')
carriers_per_disease <- fread('carriers_per_disease.txt.gz')


carriers_per_disease <- carriers_per_disease[!is.na(carriers_per_disease$logFC),]


give.n <- function(x){
  return(c(y = -0.1, label = length(x)))
}

carriers_per_disease[carriers_per_disease$disease=='other',]$disease <- 'General'
carriers_per_disease <- carriers_per_disease[!carriers_per_disease$disease=='General',]

outliers_per_disease <- carriers_per_disease %>%
  group_by(disease,impact) %>%
   summarise(outlier = sum(type=='outlier'),
             not_outlier = sum(type=='not_outlier'),
             n=length(unique(snp)))

outliers_per_disease$fraction_outlier <- outliers_per_disease$outlier/(outliers_per_disease$outlier+outliers_per_disease$not_outlier)
outliers_per_disease_noModifier <- outliers_per_disease[outliers_per_disease$impact!="MODIFIER",]

outliers_per_disease_noModifier$impact <- factor(outliers_per_disease_noModifier$impact, levels=c('LOW','MODERATE','HIGH'))



x <- 0
for(disease in unique(outliers_per_disease_noModifier$disease)){
  df <- outliers_per_disease_noModifier[outliers_per_disease_noModifier$disease==disease,]
  
  if('HIGH' %in% df$impact){
    
    test_matrix <- data.frame(outlier=c(df[df$impact=='HIGH',]$outlier, 
                                        df[df$impact=='LOW',]$outlier),
                              not_outlier=c(df[df$impact=='HIGH',]$not_outlier, 
                                            df[df$impact=='LOW',]$not_outlier))
    rownames(test_matrix) <- c('HIGH','LOW')
    p_high_low <- fisher.test(test_matrix,alternative='two.sided')$p.value*54
    
    test_matrix <- data.frame(outlier=c(df[df$impact=='HIGH',]$outlier, 
                                        df[df$impact=='MODERATE',]$outlier),
                              not_outlier=c(df[df$impact=='HIGH',]$not_outlier, 
                                            df[df$impact=='MODERATE',]$not_outlier))
    rownames(test_matrix) <- c('HIGH','MODERATE')
    p_high_modarate <- fisher.test(test_matrix,alternative='two.sided')$p.value*54
    x <- x+1
  }
  test_matrix <- data.frame(outlier=c(df[df$impact=='LOW',]$outlier, 
                                      df[df$impact=='MODERATE',]$outlier),
                            not_outlier=c(df[df$impact=='LOW',]$not_outlier, 
                                          df[df$impact=='MODERATE',]$not_outlier))
  rownames(test_matrix) <- c('LOW','MODERATE')
  p_low_modarate <- fisher.test(test_matrix,alternative='two.sided')$p.value*54
  x <- x+1
  
  if('HIGH' %in% df$impact){
    cat(paste(disease,signif(p_high_modarate,2),signif(p_high_low,2), signif(p_low_modarate,2),'\n'))
  }else{
    cat(paste(disease, signif(p_high_modarate,2),'\n'))
  }
  
}




ggplot(outliers_per_disease_noModifier, aes(impact, fraction_outlier, fill=impact))+
  geom_bar(stat='identity')+
  theme_bw(base_size = 20)+
  facet_wrap(~disease, nrow=3)+
  scale_fill_brewer(palette="Dark2")+
  ylab('Fraction of samples that is an outlier')+
  xlab('')+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label=n), vjust=-1)+
  scale_y_continuous(limit=c(0,1.1))+
  theme(strip.text.x = element_text(size = 7))
  #scale_x_discrete(limit=c('LOW', 'MODERATE','HIGH'), labels='Low','Moderate','High')+
  

#outfile = '/groups/umcg-bios/tmp04/projects/copy_from_tmp03/BIOS_manuscript/fig3//proportion_outlier_per_impact_per_disease.pdf'
#outfile = '/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig3//
outfile = 'proportion_outlier_per_impact_per_disease.pdf'
print(paste('write to:',outfile))
ggsave(outfile,width=12, height=9)





