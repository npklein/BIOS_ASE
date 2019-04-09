# calculate the rare variant enrichment by counting for each gene how often variants are seen 
# in an allele frequency bin of 0–1%, 1–5%, 5–10%, 10–25%

def main():
    # input files are made with ../allele_count_tables/combine_genes_and_samples.p
    not_outliers = stratify_snps('/groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables/alleleCountPerGroupPerGene.not_outliers.binom.txt')
    outliers = stratify_snps('/groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables/alleleCountPerGroupPerGene.outliers.binom.txt')

    with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variant_MAF_stratification/variant_stratification.binom.countOncePerSNP.txt','w') as out:
        out.write('gene\t0-0.1_outlier\t0.1-1_outlier\t1-5_outlier\t5-10_outlier\t10-25_outlier\t25-50_outlier'+
                  '\t0-0.1_not_outlier\t0.1-1_not_outlier\t1-5_not_outlier\t5-10_not_outlier\t10-25_not_outlier\t25-50_not_outlier\n')
        for gene in not_outliers:
            out.write(gene+'\t'+str(outliers[gene]['0-0.1'])+'\t'+str(outliers[gene]['0.1-1'])+'\t'+str(outliers[gene]['1-5'])+'\t'+str(outliers[gene]['5-10'])+'\t'+str(outliers[gene]['10-25'])+'\t'+str(outliers[gene]['25-50']))
            out.write('\t'+str(not_outliers[gene]['0-0.1'])+'\t'+str(not_outliers[gene]['0.1-1'])+'\t'+str(not_outliers[gene]['1-5'])+'\t'+str(not_outliers[gene]['5-10'])+'\t'+str(not_outliers[gene]['10-25'])+'\t'+str(not_outliers[gene]['25-50'])+'\n')

def stratify_snps(input_file_name):
    seen_snps = set([])
    stratification = {}
    print(input_file_name)
    with open(input_file_name) as input_file:
        print(input_file.readline().split('\t'))
        
        for line in input_file:
            if ';' in line:
                continue

            line = line.strip().split('\t')
            snp = line[1]
            if snp in seen_snps:
                continue
            seen_snps.add(snp)

            gene = line[0]
            if gene not in stratification:
                stratification[gene] = {'0-0.1':0,'0.1-1':0,'1-5':0,'5-10':0,'10-25':0,'25-50':0}

            maf = float(line[2])
            # since we are interested in minor allele frequency, have to swap if > 0.5
            # this also means that we have to swap the hom ref/alt counts or they dont match with the maf
            if maf > 0.5:
                maf = 1 - maf
                tmp = line[5]
                line[5] = line[7]
                line[7] = tmp

            # do the stratification[gene]. For each bin, count how many people are either het or hom alt for the
            # variant
            # instead of summing alternative alleles (number_of_carriers = int(line[6])+(int(line[7])*2))
            # count either as carrier or not carrier
            if int(line[6]) > 0 or int(line[7]) > 0:
                number_of_carriers = 1
            else:
                number_of_carriers = 0
            if maf < 0.001:
                stratification[gene]['0-0.1'] += number_of_carriers
            elif maf < 0.01:
                stratification[gene]['0.1-1'] += number_of_carriers
            elif maf < 0.05:
                stratification[gene]['1-5'] += number_of_carriers
            elif maf < 0.1:
                stratification[gene]['5-10'] += number_of_carriers
            elif maf < 0.25:
                stratification[gene]['10-25'] += number_of_carriers
            elif maf >= 0.25:
                stratification[gene]['25-50'] += number_of_carriers
            elif maf > 0.5:
                raise RuntimeError("is not MAF")
            else:
                raise RuntimeError("should have covered everything")
    return(stratification)


main()
