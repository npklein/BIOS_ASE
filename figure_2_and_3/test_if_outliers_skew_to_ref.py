
is_outlier = {}
with open('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE//logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt') as input_file:
    samples = header.strip().split('\t')[1:]
    for line in input_file:
        line = line.strip().split('\t')
        gene = line[0]
        if gene not in is_outlier:
            is_outlier[gene] = {}
        for index, element in enumerat(line[1:]):
            is_outlier[gene][samples[index]] = elemenet

with open('/groups/umcg-bios/tmp04/projects/copy_from_tmp03/outlierGeneASE//logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.values.txt') as input_file:
    samples = input_file.readline().strip().split('\t')
