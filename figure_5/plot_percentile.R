
library("plyr")   
library(ggplot2)
library(patchwork)

df <- read.table('expression_ase_vs_non_ase.txt',
                          header=T)
df$ceil <- round_any(df$percentile,10,f=ceiling)

df_subset <- df[df$major_count > 5 & df$minor_count>5,]


#df_subset$binom <- 
  #
df_subset$binom <- apply(df_subset, 1, function(x) {
        binom.test(as.numeric(x[['major_count']]),
                   as.numeric(x[['minor_count']])+as.numeric(x[['major_count']]))$p.value})
df_subset$fdr <- p.adjust(df_subset$binom, method='fdr')
df_subset_sign <- df_subset[df_subset$fdr < 0.05,]

df_subset$ratio_merged <- df_subset$ratio
df_subset[df_subset$ratio > 0.5,]$ratio_merged <- 1-df_subset[df_subset$ratio >0.5,]$ratio


count_per_percentile <- data.frame(t(table(df_subset$ceil)))
count_per_percentile$Var1 <- NULL
colnames(count_per_percentile) <- c('ceil','count')
count_per_percentile$ratio_merged <- 0.2
count_per_percentile$ceil <- factor(count_per_percentile$ceil, 
                                    levels=c(10,20,30,40,50,
                                             60,70,80,90,100))


df_subset$ceil <- factor(df_subset$ceil, 
                    levels=c(10,20,30,40,50,
                             60,70,80,90,100))
p1 <- ggplot(df_subset, aes(ratio_merged))+
  geom_histogram(fill='grey')+
  theme_bw()+
  facet_wrap(~ceil,nrow=2)+
  geom_text(data=count_per_percentile,
           aes(label=paste0('# genes: ',count), y=350))

p2 <- ggplot(df_subset, aes(ratio_merged))+
  geom_density(fill='grey')+
  theme_bw()+
  facet_wrap(~ceil,nrow=2)

p1/p2

ggsave('outliers_per_expression_percentile.pdf',width=10, height=10)
