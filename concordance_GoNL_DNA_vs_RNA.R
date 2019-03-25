library(ggplot2)
library(data.table)
library(reshape2)
library(ggpubr)



###### READ IN DATA #####
gonl_samples <- data.frame(fread('/groups/umcg-bios/tmp03/projects/outlierGeneASE/GoNL_WGS_vs_RNA_concordance/concordance.GoNL.out.sample'))

## Open output file
png(paste0("/groups/umcg-bios/tmp03/projects/BIOS_manuscript/suppl/concordance.GoNL.WGS.vs.RNA.png"), width = 800, height = 800)

## Plot samples
myplot<-
ggplot(gonl_samples, aes(x="", y=identicalCall, fill = identicalCall))+
  geom_violin()+
  geom_boxplot(width=0.10)+
  theme_bw(base_size = 18)+
  ggtitle(paste0("Genotype concordance DNA vs RNA-seq per sample")) +
  xlab('')+
  ylab('Identical call')+
  guides(fill=F)+
  scale_fill_brewer(palette="Dark2")+
  stat_boxplot(geom ='errorbar', width = 0.10) +
  scale_y_continuous(limits = c(0.95, 1))+
  theme_pubr()+
  theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
        axis.title.x = element_text(size=16),
        axis.text=element_text(size=16),
        axis.title.y = element_text(size=16)
  )

#Print plot to output device
print(myplot)
dev.off()
