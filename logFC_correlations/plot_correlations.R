library('ggplot2')
library('reshape2')
library(data.table)
library(pheatmap)
library(RColorBrewer)
library(viridis)

correlations <- as.data.frame(fread('logFC_correlations.txt'))
correlations$logFC <- NULL
correlations$logFC_abs <- NULL
pvals <- as.data.frame(fread('logFC_pvals.txt'))
pvals$logFC <- NULL
pvals$logFC_abs <- NULL

correlations_melt <- melt(correlations, id.vars='V1')
pvals_melt <- melt(pvals, id.vars='V1')

combined <- merge(correlations_melt, pvals_melt, by=c('V1','variable'))
colnames(combined) <- c('gene','pheno','correlation','pvalue')

rownames(pvals) <- pvals$V1
rownames(correlations) <- correlations$V1
pvals$V1 <- NULL
correlations$V1 <- NULL
pvals[is.na(pvals)] <- 1
are_significant <- pvals < (0.05/nrow(combined))
correlations_significant <- correlations
correlations_significant[!are_significant] <- NA
correlations_significant <- as.data.frame(correlations_significant)
rownames(correlations_significant) <- rownames(correlations)


ind <- apply(correlations_significant, 1, function(x) all(is.na(x)))
correlations_naRemoved <- correlations_significant[!ind,]

correlations_naRemoved <- correlations_naRemoved[,colSums(is.na(correlations_naRemoved))<nrow(correlations_naRemoved)]
correlations_naRemoved_noSD <- correlations_naRemoved[apply(correlations_naRemoved, MARGIN = 1, FUN = function(x) sd(x, na.rm=T) != 0),]
correlations_naRemoved_noSD <- correlations_naRemoved_noSD[apply(correlations_naRemoved_noSD, MARGIN = 1, FUN = function(x) !is.na(sd(x, na.rm=T))),]
correlations_naRemoved_noSD <- correlations_naRemoved_noSD[apply(correlations_naRemoved_noSD, MARGIN = 2, FUN = function(x) !is.na(sd(x, na.rm=T)))]

most_correlations <- correlations_naRemoved_noSD[apply(correlations_naRemoved_noSD, 2, function(x) length(x) - sum(is.na(x)) > 50)]

pheatmap(most_correlations, cluster_rows=F, cluster_cols = F)

