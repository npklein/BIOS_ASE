# based on compare_allele_freqs.R, but only for WGS/RNAseq genotype comparison
library(ggplot2)
library(ggpubr)
library(data.table)
library(ggExtra)

# READ DATA
# gonl wgs
gonl_af <- data.frame(fread('GoNL.AF.txt'))
colnames(gonl_af) <- c('snp','ref','alt','AF_GoNL')
# gonl RNAseq genotypes
gonl_RNA_af <- data.frame(fread('GoNL.RNAseqGenotypes.AF.txt'))
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
# GoNL samples gotten from the script that makes the AFs
samples_label = paste0('Samples == 271')
n_label <- paste0("SNPs == ",nrow(AF))
model = lm(AF1 ~ AF2, data = AF)
r2_label <- paste0("\nR^2 == ", signif(summary(model)$adj.r.squared,3))


p <- ggplot(AF, aes(AF1*100, AF2*100))+
  geom_point(alpha=0.1)+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  annotate("text", x = 10, y = 95, label = samples_label,
           size=6,parse=TRUE)+
  annotate("text", x = 10, y = 90, label = n_label,
           size=6,parse=TRUE)+
  annotate("text", x = 10, y = 85, label = r2_label,
           size=6,parse=TRUE)+ 
  geom_segment(aes(x = 0, y = 10, xend = 10, yend = 10, colour="red"))+ 
  geom_segment(aes(x = 0, y = 0, xend = 10, yend = 0, colour="red"))+
  geom_segment(aes(x = 0, y = 0, xend = 0, yend = 10, colour = "red"))+
  geom_segment(aes(x = 10, y = 0, xend = 10, yend = 10, colour="red"))+
  guides(colour=F, size=F)+ 
  scale_x_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))+ 
  scale_y_continuous(breaks = c(0,25,50,75,100),labels = paste0(c("0%", "25%", "50%", "75%", "100%")))



png('figures/WGS_vs_RNAseq_genotypes.png', width=600, height=600)
ggMarginal(p, type = "histogram", bins=1000)
dev.off()



AF_rare <- AF[AF$AF1 < 0.1 & AF$AF2 < 0.1,]
n_label <- paste0("N == ",nrow(AF_rare))
model = lm(AF1 ~ AF2, data = AF_rare)
r2_label <- paste0("\nR^2 == ", signif(summary(model)$adj.r.squared,3))

p <- ggplot(AF_rare, aes(AF1*100, AF2*100))+
  geom_point(alpha=0.1)+
  theme_bw(base_size = 18)+
  xlab(paste0('Allele Frequency WGS genotypes'))+
  ylab(paste0('Allele Frequency RNAseq genotypes'))+
  geom_abline(lty=2, colour='red')+
  annotate("text", x = 1.1, y = 10, label = samples_label,
           size=6,parse=TRUE)+
  annotate("text", x = 1, y = 9.6, label = n_label,
           size=6,parse=TRUE)+
  annotate("text", x = 1, y = 9.2, label = r2_label,
           size=6,parse=TRUE)+ 
  scale_x_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))+ 
  scale_y_continuous(breaks = c(0,2.5,5,7.5,10),labels = paste0(c("0%", "2.5%", "5%", "7.5%", "10%")))


png('figures/WGS_vs_RNAseq_genotypes_rare.png', width=600, height=600)
ggMarginal(p, type = "histogram", bins=40)
dev.off()


