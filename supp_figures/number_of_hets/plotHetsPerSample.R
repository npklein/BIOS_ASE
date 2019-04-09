library(data.table)
library(ggplot2)
hets_per_sample <- read.table('hets_per_sample.txt', sep='\t',header=T,stringsAsFactors=F)

print(paste0('Number of genes where at least 1 sample had 1 het: ',length(unique(hets_per_sample$gene))))

hets_per_sample_genes <- hets_per_sample[!hets_per_sample$gene=='noGene',]
hets_per_sample_genes <- hets_per_sample_genes[!hets_per_sample_genes$chr=="chr",]
samples_per_gene <- as.data.frame(table(hets_per_sample_genes$gene))
colnames(samples_per_gene) <- c('gene', 'numberOfSamples')
samples_per_gene$chr <- hets_per_sample_genes[match(samples_per_gene$gene,hets_per_sample_genes$gene),]$chr
samples_per_gene$chr <- factor(samples_per_gene$chr, levels=c('1','2','3','4','5','6','7','8','9',
                                                                 '10','11','12','13','14','15','16',
                                                                 '17','18','19','20','21','22'))
##### number of samples per gene #####
ggplot(samples_per_gene, aes(numberOfSamples))+
  geom_histogram(binwidth = 80)+
  theme_bw(base_size=18)+
  ylab("Number of genes")+
  xlab("Number of samples with a heterozygous SNP per gene")
  
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/numberSamplesPerGene.png', width=8, height=6)

genes_per_sample <- as.data.frame(table(hets_per_sample_genes$sample))
colnames(genes_per_sample) <- c('sample', 'numberOfGenes')
ggplot(genes_per_sample, aes(numberOfGenes))+
  geom_histogram()+
  theme_bw(base_size=18)+
  xlab('Number of genes with a heterozygous SNP per sample')+
  ylab("Number of samples")
  
  
  ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/numberOfGenesPerSamples.png',  width=8, height=6)
