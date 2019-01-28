# subset the logfold change table on outlier samples that have > 1000 ASE genes
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS1000.depthFiltered.bonferroni.txt') as input_file:
    samples = set(input_file.read().split('\n'))

with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.txt') as input_file:
    with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt','w') as out:
        header = input_file.readline().strip().split('\t')
        indexes = []
        out.write(header[0])
        for index, element in enumerate(header):
            if element in samples:
                out.write('\t'+element)
                indexes.append(index)
        out.write('\n')
        for line in input_file:
            line = line.strip().split('\t')
            print(line[0])
            out.write(line[0])
            for index in indexes:
                out.write('\t'+line[index])
            out.write('\n')
