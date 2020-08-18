

clinvar <- read.table('clinvar_overlapped_with_SNPs.txt', sep='\t',header=T)

clinvar <- clinvar[2:nrow(clinvar),]
clinvar$het <- as.numeric(as.character(clinvar$het))
clinvar$ref <- as.numeric(as.character(clinvar$ref))
clinvar$alt <- as.numeric(as.character(clinvar$alt))
table(clinvar$clinstat)


hist(clinvar[clinvar$clinstat=="Benign",]$het)


pat <- clinvar[clinvar$clinstat=="Pathogenic",]
