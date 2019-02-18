library(stringr)
library(data.table)
library(ggplot2)
library(plyr)

files = list.files(pattern="*.perSNP.txt", path="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/switchAndErrorsPerGene/",
                   full.names=T)

switch_errors <- data.frame()
for(f in files){
  dfTMP <- data.frame(fread(f))
  dfTMP$CHR <- str_extract_all(f, 'chr\\d+')[[1]]
  switch_errors <- rbind(switch_errors, dfTMP)
}


switch_errors$switchError <- (switch_errors$switch/(switch_errors$switch+switch_errors$no_switch))*100
switch_errors <- switch_errors[!is.na(switch_errors$switchError),]


switch_errors[switch_errors$switchError > 50, ]$switchError <- 100-switch_errors[switch_errors$switchError > 50, ]$switchError

switch_errors[switch_errors$maf_test > 0.5,]$maf_test <- 1-switch_errors[switch_errors$maf_test > 0.5,]$maf_test

#switch_errors$maf_breaks <- cut(switch_errors$maf_test, breaks=c(0, 0.05, 0.1, 0.03, 0.04,0.05, 0.06, 0.07,0.08, 0.09, 0.1, 0.2, 0.3, 0.4, 0.51))
switch_errors$maf_breaks <- cut(switch_errors$maf_test, breaks=c(0,0.001,0.01, 0.1, 0.2, 0.3, 0.4, 0.51))



switch_errors_no_na <- switch_errors[!is.na(switch_errors$switchError),]
switch_per_maf_bin <- ddply(switch_errors_no_na,~maf_breaks,summarise,mean_error=mean(switchError),sd=sd(switchError))
switch_per_maf_bin <- switch_per_maf_bin[!is.na(switch_per_maf_bin$maf_breaks),]

switch_per_maf_bin$maf <- c(0.005, 0.05, 0.1,
                            0.2, 0.3, 0.4, 0.5)

###### make maf bins and plot boxplots ######
switch_errors_no_na <- switch_errors_no_na[!is.na(switch_errors_no_na$maf_breaks),]

give.n <- function(x){
  return(c(y = max(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

switch_errors_no_na_subset <- switch_errors_no_na[switch_errors_no_na$switchError < 10,]
ggplot(switch_errors_no_na, aes(as.factor(maf_breaks), switchError, fill = maf_breaks))+
  geom_violin()+
  geom_jitter(data=switch_errors_no_na, aes(as.factor(maf_breaks), switchError), 
              alpha=0.1, width=0.2)+
  theme_bw(base_size=18)+
  xlab('MAF bins')+
  ylab('Percentage of samples for which SNP has a swap error')+
  stat_summary(fun.data = give.n, geom = "text", fun.y = median,
               position = position_dodge(width = 0.75))+
  scale_x_discrete(labels=c('< 0.1%','0.1-1%','1-10%','10-20%','20-30%','30-40%','40-50%'))+
  guides(fill=F)+
  scale_colour_brewer(palette="Dark2")
ggsave('/groups/umcg-bios/tmp03/projects/BIOS_manuscript/fig1/panel_c//switch_error_per_maf_bin.png', width=8, height = 8)

######
