library(ggplot2)
library(data.table)
library(ggpubr)

#Retrieve CMD args
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#Testing purpose
#dataTab<-data.frame(fread('/Users/freerkvandijk/Downloads/ase_sampleAse.all.binom.txt'))
#dataTab<-data.frame(fread('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/test.ase_sampleAse.all.binom.txt'))

#Read input file
dataTab<-data.frame(fread(args[1]))

dataTab$Ref_Counts <- as.numeric(dataTab$Ref_Counts)
dataTab$Alt_Counts <- as.numeric(dataTab$Alt_Counts)

#Extract unique SNPs
snps<-unique(dataTab$snp_id)
#snps<-c("rs1801287")

#snps<-c("rs307934")
#snps<-c("10:100004576")

#Iterate over all SNPs
for (SNP in snps){
  dataTab1<-c()
  dataTab2<-c()
  dataTab1<-dataTab[dataTab$snp_id == SNP,] #Select SNP of interest
  
  # Merge alt and ref counts and select max value
  mergedCounts<-c(dataTab1$Ref_Counts,dataTab1$Alt_Counts) #Merge ref and alt counts and take max value
  max<-max(mergedCounts)
  
  #Subset dataTab1 to create dataframe for Ref_Counts >= 10 & Alt_Count >= 10
  dataTab2<-dataTab1[dataTab1$Ref_Counts >= 10 & dataTab1$Alt_Counts >= 10,]
  
  #Calculate slope
  if (dim(dataTab2)[1] != 0) {
    lmResult<-coef(lm(Alt_Counts ~ Ref_Counts -1, data=dataTab2))
    #intercept<-lmResult["(Intercept)"]
    slope<-lmResult["Ref_Counts"]
  }
  
  #Set significance levels
  dataTab1$Significance <- "Non-significant"
  dataTab1$Significance[dataTab1$FDR < 0.05] <- "FDR < 0.05"
  dataTab1$Significance[dataTab1$Bonf_corr_Pval < 0.05 & dataTab1$FDR < 0.05] <- "Bonferroni P-val < 0.05"
  dataTab1$Significance[dataTab1$Ref_Counts < 10 | dataTab1$Alt_Counts < 10] <- "Excluded"
  
  #Set fixed colors to significance factors
  sign_levels<-c("Excluded", "Non-significant", "FDR < 0.05", "Bonferroni P-val < 0.05")
  dataTab1$Significance <- factor(dataTab1$Significance, levels=sign_levels)

  
#print(
  
  #Set output file path for png files
  #png(paste0("/Users/freerkvandijk/Downloads/testPlots/",dataTab1$Chromosome,"_",dataTab1$Position,".png"), width = 800, height = 600)
  png(paste0("/groups/umcg-bios/tmp03/projects/BIOS_manuscript/ASEbrowserplots/chr",dataTab1$Chromosome,"/png/",dataTab1$Chromosome,"_",dataTab1$Position,".png"), width = 800, height = 600)
  if (dim(dataTab2)[1] <= 1) { #If 0 or 1 samples to calculate from, don't print slope
    myplot<-
      ggplot(dataTab1, aes(x=Ref_Counts, y=Alt_Counts, colour=Significance))+
      geom_point(alpha=0.7)+
      xlab(paste0('Reference allele count'))+
      ylab(paste0('Alternative allele count'))+
      guides(fill=F)+
      scale_x_continuous(limits = c(0, max))+
      scale_y_continuous(limits = c(0, max))+
      geom_abline(linetype=2)+
      scale_color_manual(values=setNames(c("gray", "orange", "steelblue", "blue"),sign_levels))+
      theme_pubr()+
      scale_fill_brewer(palette="Dark2")+
      theme(legend.title=element_blank(), legend.position="right")+
      guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
      coord_fixed(ratio = 1)+
      ggtitle(paste0("Variant: ", dataTab1$snp_id, "               Gene: ", dataTab1$Gene_symbol, " (", dataTab1$Ensembl_ID, ")"))+
      theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
            axis.title.x = element_text(size=20),
            axis.text=element_text(size=20),
            axis.title.y = element_text(size=20),
            plot.title = element_text(size = 20, face = "bold"),
            legend.title=element_text(size=20), 
            legend.text=element_text(size=20)
      )
    
  }else{ #Create normal plot with slope
  myplot<-
      ggplot(dataTab1, aes(x=Ref_Counts, y=Alt_Counts, colour=Significance))+
      geom_point(alpha=0.7)+
      xlab(paste0('Reference allele count'))+
      ylab(paste0('Alternative allele count'))+
      guides(fill=F)+
      scale_x_continuous(limits = c(0, max))+
      scale_y_continuous(limits = c(0, max))+
      geom_abline(intercept = 0, slope = (slope), colour="red", lty=2)+
      geom_abline(linetype=2)+
      scale_color_manual(values=setNames(c("gray", "orange", "steelblue", "blue"),sign_levels))+
      theme_pubr()+
      scale_fill_brewer(palette="Dark2")+
      theme(legend.title=element_blank(), legend.position="right")+
      guides(colour = guide_legend(override.aes = list(alpha = 1, size=3)))+
      coord_fixed(ratio = 1)+
      ggtitle(paste0("Variant: ", dataTab1$snp_id, "               Gene: ", dataTab1$Gene_symbol, " (", dataTab1$Ensembl_ID, ")"))+
      theme(panel.grid.major = element_line(colour = "grey", size = 0.1),panel.grid.minor = element_line(colour = "grey", size = 0.1),
            axis.title.x = element_text(size=20),
            axis.text=element_text(size=20),
            axis.title.y = element_text(size=20),
            plot.title = element_text(size = 20, face = "bold"),
            legend.title=element_text(size=20), 
            legend.text=element_text(size=20)
            )
  }
#)
  
  #Print plot to output device
  print(myplot)
  dev.off()

}
