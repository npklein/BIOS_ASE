# based on compare_allele_freqs.R, but only for WGS/RNAseq genotype comparison
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggExtra)
library(viridis)

# READ DATA
# gonl wgs
gonl_af <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/GoNL.AF.txt'))
#gonl_af <- data.frame(fread('GoNL.AF.txt'))
colnames(gonl_af) <- c('snp','ref','alt','AF_GoNL')
# gonl RNAseq genotypes
gonl_RNA_af <- data.frame(fread('/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged/exac_comparison/GoNL.RNAseqGenotypes.AF.txt'))
#gonl_RNA_af <- data.frame(fread('GoNL.RNAseqGenotypes.AF.txt'))
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

p <- ggplot(AF, aes(AF1*100, AF2*100))+
  geom_hex()+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  geom_segment(aes(x = -3, y = 10, xend = 10, yend = 10, colour="red"), size=1.5)+ 
  geom_segment(aes(x = -3, y = -3, xend = 10, yend = -3, colour="red"), size=1.5)+
  geom_segment(aes(x = -3, y = -3, xend = -3, yend = 10, colour = "red"), size=1.5)+
  geom_segment(aes(x = 10, y = -3, xend = 10, yend = 10, colour="red"), size=1.5)+
  guides(colour=F, size=F)+ 
  scale_x_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))+ 
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))+
  scale_fill_viridis(trans = "log10")+
  labs(fill = "# SNPs")+
  theme(legend.position="top",
        legend.key = element_rect(fill = "white", colour = "black"))+
  stat_cor(aes(label = ..r.label..), size=5,geom='label')

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_a//WGS_vs_RNAseq_genotypes.pdf', width=6, height=7, plot = p)


AF_rare <- AF[AF$AF1 < 0.1 & AF$AF2 < 0.1,]
model = lm(AF1 ~ AF2, data = AF_rare)

p <- ggplot(AF_rare, aes(AF1*100, AF2*100))+
  geom_hex()+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))+ 
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))+
  scale_fill_viridis(trans = "log10")+
  labs(fill = "# SNPs")+
  theme(legend.position="top",
        legend.key = element_rect(fill = "white", colour = "black"))+
  stat_cor(aes(label = ..r.label..), size=5,geom='label')

ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_b/WGS_vs_RNAseq_genotypes_rare.pdf', width=6, height=7, plot = p)


