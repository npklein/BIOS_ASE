
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferonni.txt') as input_file:
    print(input_file.readline().split('\t').index('AD2DBFACXX-6-10'))
    for line in input_file:
        line = line.strip().split('\t')
        if line[0] == 'ENSG00000139842':
            print(line[1769])
