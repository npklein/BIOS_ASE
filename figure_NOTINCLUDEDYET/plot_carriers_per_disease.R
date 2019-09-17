library(ggplot2)
library(dplyr)
library(ggpubr)

carriers_HIGH <- read.table('carriers_per_disease/carriers_per_disease_HIGH_filtered_variants.txt', header=T, sep='\t')
carriers_LOW <- read.table('carriers_per_disease/carriers_per_disease_LOW_same_genes_as_HIGH.txt', header=T, sep='\t')
carriers_LOW_same_gene <- read.table('carriers_per_disease/carriers_per_disease_LOW_same_genes_as_HIGH_mafWithinGeneSelect.txt',header=T,sep='\t')
carriers_LOW_mafSelect <- read.table('carriers_per_disease/carriers_per_disease_LOW_same_genes_as_HIGH_mafSelect.txt',header=T,sep='\t')


melt_carriers <- function(carriers, disease_levels = c()){
  carriers[carriers$disease=="other",][2:4] <- carriers[carriers$disease=="other",][c(2:4)] + carriers[carriers$disease=="General",][c(2:4)]
  carriers <- carriers[!carriers$disease=="General",]
  colnames(carriers[carriers$disease=="other",]) <- 'Other'
  
  carriers$total <- carriers$outlier+carriers$not_outlier
  carriers$outlier_perc <- (carriers$outlier/carriers$total)*100
  carriers$not_outlier_perc <- (carriers$not_outlier/carriers$total)*100
  
  carriers_outliers <- carriers[c('disease','outlier','outlier_perc')]
  carriers_outliers$type <- 'outlier'
  colnames(carriers_outliers) <- c('disease', 'count','perc','type')
  
  carriers_not_outliers <- carriers[c('disease','not_outlier','not_outlier_perc')]
  carriers_not_outliers$type <- 'not_outlier'
  colnames(carriers_not_outliers) <- c('disease', 'count','perc','type')
  
  
  carriers_melt <- rbind(carriers_outliers, carriers_not_outliers)
  carriers_melt$pos <- NA
  carriers_melt[carriers_melt$type=='outlier',]$pos <- NA
  carriers_melt[carriers_melt$type=='not_outlier',]$pos <- 75
  carriers_melt$total <- carriers_melt[carriers_melt$type=='not_outlier',]$count + carriers_melt[carriers_melt$type=='outlier',]$count
  carriers_melt[carriers_melt$count==0,]$pos <- NA
  if(length(disease_levels)==0){
    carriers_melt$disease <- factor(carriers_melt$disease, levels=unique(carriers_melt$disease[order(carriers_melt[carriers_melt$type=="outlier",]$perc)]), 
                                  ordered=TRUE)
  }
  else{
    carriers_melt$disease <- factor(carriers_melt$disease, levels=disease_levels, 
                                    ordered=TRUE)
  }
  
  carriers_melt <- carriers_melt[carriers_melt$count!=0,]
}
carriers_HIGH_melt <- melt_carriers(carriers_HIGH)
carriers_LOW_melt <- melt_carriers(carriers_LOW, levels(carriers_HIGH_melt$disease))

# carriers_LOW_same_gene: same genes and MAF is also select within the same gene
# carriers_LOW_mafSelect: same gene but MAF can be selected from any of the genes
carriers_LOW_melt_same_gene <- melt_carriers(carriers_LOW_same_gene, levels(carriers_HIGH_melt$disease))
carriers_LOW_melt_mafSelect <- melt_carriers(carriers_LOW_mafSelect, levels(carriers_HIGH_melt$disease))

plot_carriers <- function(carriers_melt, name,title, xaxis=T){
  averagePerc <- round((sum(carriers_melt[carriers_melt$type=="outlier",]$count) / ( sum(carriers_melt[carriers_melt$type=="outlier",]$count) + sum(carriers_melt[carriers_melt$type=="not_outlier",]$count)))*100)
#  carriers_melt <- rbind(data.frame('disease'='Total', 
#                                                   count=sum(carriers_melt[carriers_melt$type=='outlier',]$count), 
#                                                   perc=sum(carriers_melt[carriers_melt$type=='outlier',]$perc)/nrow(carriers_melt[carriers_melt$type=='outlier',]),
#                                                   type='outlier', pos=NA, 
#                                                   total=sum(carriers_melt[carriers_melt$type=='outlier',]$total)),carriers_melt)
#  carriers_melt <- rbind(carriers_melt, data.frame('disease'='Total', 
#                                                   count=sum(carriers_melt[carriers_melt$type=='not_outlier',]$count), 
#                                                   perc=sum(carriers_melt[carriers_melt$type=='not_outlier',]$perc)/nrow(carriers_melt[carriers_melt$type=='not_outlier',]),
#                                                   type='not_outlier', pos=85, 
#                                                   total=sum(carriers_melt[carriers_melt$type=='not_outlier',]$total)), carriers_melt)
  p <- ggplot(carriers_melt, aes(x=disease, y=perc, fill=type))+
    geom_bar(stat='identity',position='stack', colour="black")+
    theme_pubr(base_size = 12) + 
    coord_flip()+ 
    scale_fill_manual(values=c("#999999", "#E69F00"), labels=c('No gene ASE','Gene ASE'))+
    ylab('Percentage of carriers')+
    geom_text(data=carriers_melt, aes(x = disease, y = pos,label = paste("Total: ",count)), size=4,
              colour="white", fontface = "bold")+ 
    theme(legend.title=element_blank())+
    geom_hline(yintercept=averagePerc, lty=2, colour='white',size=1)+
    scale_y_continuous(breaks=c(seq(0, 100, by=25), averagePerc))+
    xlab('')+
    ggtitle(title)
  if(xaxis == F){
    p <- p + theme(axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank())
    ggsave(paste0('figures/carriers_per_disease_',name,'.pdf'), width=3, height=7,
           plot=p) 
  }
  else{
    ggsave(paste0('figures/carriers_per_disease_',name,'.pdf'), width=5, height=7,
           plot=p)   
  }
  
}

plot_carriers(carriers_HIGH_melt,'HIGH','High impact SNPs')
plot_carriers(carriers_LOW_melt,'LOW_noMafSelect','Low impact SNPs\nNo MAF select', xaxis=F)
plot_carriers(carriers_LOW_melt_same_gene,'LOW_sameGene_mafSelect','Low impact SNPs\nMAF select all genes',xaxis=F)
plot_carriers(carriers_LOW_melt_mafSelect,'LOW_mafSelect','Low impact SNPs\nMAF select within genes', xaxis=F)

# test
fisher_table <- data.frame(outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="outlier",]$count),
                                      sum(carriers_LOW_melt[carriers_LOW_melt$type=="outlier",]$count)),
                           not_outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="not_outlier",]$count),
                                                 sum(carriers_LOW_melt[carriers_LOW_melt$type=="not_outlier",]$count)))
rownames(fisher_table) <- c('HIGH','LOW')

fisher.test(fisher_table)$p.value                        


fisher_table <- data.frame(outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="outlier",]$count),
                                      sum(carriers_LOW_melt_same_gene[carriers_LOW_melt_same_gene$type=="outlier",]$count)),
                           not_outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="not_outlier",]$count),
                                          sum(carriers_LOW_melt_same_gene[carriers_LOW_melt_same_gene$type=="not_outlier",]$count)))
rownames(fisher_table) <- c('HIGH','LOW')

fisher.test(fisher_table)$p.value                        





fisher_table <- data.frame(outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="outlier",]$count),
                                      sum(carriers_LOW_melt_mafSelect[carriers_LOW_melt_mafSelect$type=="outlier",]$count)),
                           not_outliers=c(sum(carriers_HIGH_melt[carriers_HIGH_melt$type=="not_outlier",]$count),
                                          sum(carriers_LOW_melt_mafSelect[carriers_LOW_melt_mafSelect$type=="not_outlier",]$count)))
rownames(fisher_table) <- c('HIGH','LOW')

fisher.test(fisher_table)$p.value        
