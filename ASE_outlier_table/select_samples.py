with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS500.depthFiltered.binom.txt') as input_file:
    samples = set(input_file.read().split('\n'))

input_dir = '/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/'
with open(input_dir+'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.txt') as input_file:
    with open(input_dir+'genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.txt','w') as out:
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
