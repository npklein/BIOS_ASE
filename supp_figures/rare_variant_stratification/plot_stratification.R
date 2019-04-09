library(dplyr)
library(ggplot2)
library(reshape2)
library(ggsignif)

input_dir <- '/groups/umcg-bios/tmp03/projects/outlierGeneASE/variant_MAF_stratification/'
stratificationTable <- read.table(paste0(input_dir, 'variant_stratification.binom.countOncePerSNP.txt'),
                                  header = T)
#omim_genes <- read.table('omim.genes.txt')
#cgd <- read.table('tmp.important.columns.only.txt', sep='\t', header=T)
#cgd_genes <- cgd[cgd$GenePresentInCGD=="YES",]$ENSEMBLID

  
rownames(stratificationTable) <- stratificationTable$gene
stratificationTable$gene <- NULL
  
stratificationTable <- stratificationTable[!is.na(stratificationTable$X0.0.1_outlier),]
stratificationSum <- colSums(stratificationTable)
  
not_outlier <- stratificationSum[grepl('not_outlier', names(stratificationSum))]
names(not_outlier) <-  c('0-0.1','0.1-1','1-5', '5-10','10-25','25-50')
outlier <- stratificationSum[!grepl('not_outlier', names(stratificationSum))]
names(outlier) <-   c('0-0.1','0.1-1','1-5', '5-10','10-25','25-50')
stratification <- as.data.frame(t(data.frame('outliers'=outlier,'not_outliers' = not_outlier)))
  
  
#### poster plot, use correct colours, put legend in the plot, combine maf bins 5-50 #####
# merge the mafs and remove unused columns
stratification_merge_high_mafs <- stratification
stratification_merge_high_mafs$type <- NULL
stratification_merge_high_mafs$`5-50` <- stratification_merge_high_mafs$`5-10`+stratification_merge_high_mafs$`10-25`+stratification_merge_high_mafs$`25-50`
stratification_merge_high_mafs$`5-10` <- NULL
stratification_merge_high_mafs$`10-25` <- NULL
stratification_merge_high_mafs$`25-50` <- NULL
# melt the numbers per group
stratification_merge_high_mafs$type <- rownames(stratification_merge_high_mafs)
stratification_merge_high_mafs_melt <- melt(stratification_merge_high_mafs)
stratification_merge_high_mafs$type <- NULL
# calculate the percentages and melt
stratification_merge_high_mafs_percentage <- (stratification_merge_high_mafs/rowSums(stratification_merge_high_mafs))*100
stratification_merge_high_mafs_percentage$type <- rownames(stratification_merge_high_mafs_percentage)
stratification_merge_high_mafs_percentage_melt <- melt(stratification_merge_high_mafs_percentage)
stratification_merge_high_mafs_percentage_melt$variable <- as.factor(stratification_merge_high_mafs_percentage_melt$variable)
stratification_merge_high_mafs_percentage_melt$totalSnps <- stratification_merge_high_mafs_melt[stratification_merge_high_mafs_melt$type == stratification_merge_high_mafs_percentage_melt$type,]$value
  
# calculate the pvalue by comparing the proportions
chi_square_table_0 <- data.frame('observed_0-0.1' = stratification_merge_high_mafs$`0-0.1`,
                                 rest=rowSums(stratification_merge_high_mafs[c(2,3,4)]))
chi_square_table_01 <- data.frame('observed_0.1-1' = stratification_merge_high_mafs$`0.1-1`,
                                  rest=rowSums(stratification_merge_high_mafs[c(1,3,4)]))
chi_square_table_1 <- data.frame('observed_1-5' = stratification_merge_high_mafs$`1-5`,
                                 rest=rowSums(stratification_merge_high_mafs[c(1,2,4)]))
chi_square_table_5 <- data.frame('observed_5-10' = stratification_merge_high_mafs$`5-50`,
                                 rest=rowSums(stratification_merge_high_mafs[c(1,2,3)]))

pval0  <- fisher.test(as.matrix(chi_square_table_0))$p.value
pval01  <- fisher.test(as.matrix(chi_square_table_01))$p.value
pval1  <- fisher.test(as.matrix(chi_square_table_1))$p.value
pval5  <- fisher.test(as.matrix(chi_square_table_5))$p.value
label.df <- data.frame(variable = c("0-0.1","0.1-1", "1-5","5-50"),
                       value = c(76,20, 20,48), 
                       type='outliers','outliers','outliers', 'outliers','outliers','outliers',
                       pval=c(paste('p =',signif(pval0,digits=2)), paste('p =',signif(pval01,digits=2)), 
                              paste('p =',signif(pval1,digits=2)), paste('p =',signif(pval5,digits=2))))
ggplot(stratification_merge_high_mafs_percentage_melt, aes(x=variable, y=value,fill=type))+
    geom_bar(stat='identity', position='dodge', colour="black")+
    theme_bw(base_size = 26)+
    ylab('Relative percentage of SNPs')+
    xlab('MAF bins')+
    scale_x_discrete(labels = c('0-0.1', '0.1-1','1-5','5-50'))+ 
    scale_fill_manual(values=c("#999999", "#E69F00"), labels=c('No gene ASE','Gene ASE'))+
    geom_text(data = label.df, label = label.df$pval)+
    labs(fill='')+ 
    theme(legend.position = c(0.8, 0.85),
          panel.grid.major.x = element_blank(), axis.ticks.x=element_blank())+
    geom_text(aes(label=format(totalSnps, big.mark=",", scientific=FALSE,trim=T),
                  y=-2), 
                  color="black", size=4,
                  position = position_dodge(width = 1))
ggsave(paste0('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/rare_variant_enrichment_binom.supFigure9.pdf'),width=10, height=8)



