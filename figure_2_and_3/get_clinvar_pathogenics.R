library(readr)
library(RCurl)

date <- format(Sys.time(), "%Y-%b-%d")


# INPUT: clinvar file as downloaded from clinvar
clinvar_file <- "/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/variant_summary_2018-Nov-16.txt.gz"

# OUTPUT: Clinvar file with certain criteia (see ReviewStatus lower in script)
clinvar_out_file <- paste0("/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvarSNPs_",date,".txt")


# If the file does not exist, download it from clinvar. However, this file updates without chaning the name. To exactly
# reproduce figures from manuscript, you need to use variant_summary_2018-Nov-16.txt.gz. Can contact niekdeklein@gmail.com 
# for this file. If you want to do a similar plot but does not have to be exact reproduction, you can use newer version of ClinVar
if(!file.exists(clinvar_file)){
    print(paste(clinvar_file,"does not exist, downloading new clinvar file. Change clinvar_file in this script and rerun"))
    url <- "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
    # wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz <- like wget
    outfile <- paste0('/groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/variant_summary_',date,'.txt.gz')
    download.file(url, destfile=outfile)
    print(paste("Downloaded file to",outfile))
    print("Change input file in the script and run again")
    q()
}

clinvar <- read_delim(clinvar_file, delim = "\t", quote = "")


clinvar$Selected <- 
  (
    clinvar$ReviewStatus == "criteria provided, multiple submitters, no conflicts"
    | clinvar$ReviewStatus == "reviewed by expert panel"
    | clinvar$ReviewStatus == "practice guideline"
  ) 
  #& (
    #clinvar$ClinicalSignificance == "Pathogenic"
    #| clinvar$ClinicalSignificance == "Pathogenic;not provided"
  #) 


clinvarSNPs <-  clinvar[clinvar$Type == "single nucleotide variant" &
                        clinvar$Assembly == "GRCh37",]
clinvarSNPs <- data.frame(clinvarSNPs)

clinvarSNPs$ClinSimple <- clinvarSNPs$ClinicalSignificance
clinvarSNPs[grepl('pathogenic', clinvarSNPs$ClinicalSignificance,ignore.case=TRUE),]$ClinSimple <- "Pathogenic"
clinvarSNPs[grepl('benign', clinvarSNPs$ClinicalSignificance,ignore.case=TRUE),]$ClinSimple <- "Benign"
clinvarSNPs[grepl('uncertain significance', clinvarSNPs$ClinicalSignificance,ignore.case=TRUE),]$ClinSimple <- "Uncertain significance"
clinvarSNPs$snp <- paste0(clinvarSNPs$Chromosome,'_',clinvarSNPs$Start)

write.table(clinvarSNPs, clinvar_out_file, quote = F, sep = "\t", row.names = F)
print(paste("written to",clinvar_out_file))
