# based on compare_allele_freqs.R, but only for WGS/RNAseq genotype comparison
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggExtra)
library(viridis)
library(grid)
library(gtable)

# READ DATA
# gonl wgs
gonl_af <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/GoNL.AF.txt'))
colnames(gonl_af) <- c('snp','ref','alt','AF_GoNL')
# gonl RNAseq genotypes
gonl_RNA_af <- data.frame(fread('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/exac_comparison/GoNL.RNAseqGenotypes.AF.txt'))
colnames(gonl_RNA_af) <- c('snp','ref','alt','AF_GoNL_RNAseq')


# FUNCTION TO MERGE AND PLOT
merge_AF <- function(AF1, AF2, name1, name2){
  AF <- merge(AF1, AF2, by='snp')
  colnames(AF) <- c('snp', 'ref.x','alt.x', name1, 'ref.y','alt.y',name2)
  AF <- AF[AF[[name1]]> 0 & AF[[name2]]>0, ]
  #AF$MAF1 <- ifelse(AF[[name1]] >= 0.5, AF[[name1]]-0.5, AF[[name1]])
  #AF$MAF2 <- ifelse(AF[[name2]] >= 0.5, AF[[name2]]-0.5, AF[[name2]])
  AF$AF1 <- AF[[name1]]
  AF$AF2 <- AF[[name2]]
  AF$AF2swapped <- ifelse(AF$ref.x == AF$ref.y, AF$AF2, 1-AF$AF2)
  AF <- AF[(AF$ref.x == AF$ref.y & AF$alt.x==AF$alt.y) | (AF$ref.x == AF$alt.y & AF$alt.x == AF$ref.y),]
  return(AF)
}

############################ PLOT AFs ###################################
AF <- merge_AF(gonl_RNA_af, gonl_af, 'GoNL RNA','GoNL WGS')
correlation <- cor(AF$AF1, AF$AF2,method='spearman')
correlation <- cor(AF$AF1, AF$AF2,method='pearson')
# GoNL samples gotten from the script that makes the A
model = lm(AF1 ~ AF2, data = AF)
label = paste0('Samples == 271\n',
               'SNPs == ', nrow(AF),
               '\nR^2 == ', signif(summary(model)$adj.r.squared,3))


p <- ggplot(AF, aes(AF1*100, AF2*100))+
#  geom_point(alpha=0.1)+
  geom_hex()+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  annotate("text", x = 10, y = 95, label = label,
           size=6,parse=TRUE, hjust=1)+
  geom_segment(aes(x = -1, y = 10, xend = 10, yend = 10, colour="red"))+ 
  geom_segment(aes(x = -1, y = -1, xend = 10, yend = 0, colour="red"))+
  geom_segment(aes(x = -1, y = -1, xend = -1, yend = 10, colour = "red"))+
  geom_segment(aes(x = 10, y = -1, xend = 10, yend = 10, colour="red"))+
  guides(colour=F, size=F)+ 
  scale_x_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))+ 
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))+
  scale_fill_viridis(trans = "log10")+
  labs(fill = "log10(count)")+
  theme(legend.position="top")+
  theme(plot.margin = unit(c(.5,6,.5,.5),"lines"),
          legend.background = element_rect(colour = "black"))



ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_a//WGS_vs_RNAseq_genotypes.pdf', width=6, height=6, plot = p)


AF_rare <- AF[AF$AF1 < 0.1 & AF$AF2 < 0.1,]
model = lm(AF1 ~ AF2, data = AF_rare)

label = paste0('Samples == 271\n',
               'SNPs == ',nrow(AF_rare),
               '\nR^2 == ', signif(summary(model)$adj.r.squared,3))
p <- ggplot(AF_rare, aes(AF1*100, AF2*100))+
#  geom_point(alpha=0.1)+
  geom_hex()+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  annotate("text", x = 1.1, y = 10, label = label,
           size=6,parse=TRUE,hjust=0)+
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))+ 
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))+
  scale_fill_viridis(trans = "log10")+
  labs(fill = "log10(count)")+
  theme(legend.position="top")

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_b/WGS_vs_RNAseq_genotypes_rare.pdf', width=6, height=6, plot = p)


