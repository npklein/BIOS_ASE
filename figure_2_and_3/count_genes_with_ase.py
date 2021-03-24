


total_genes = 0
genes_ase = 0
with open('/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.txt') as input_file:
    input_file.readline()
    for line in input_file:
        line = line.strip().split('\t')
        yes = float(line.count('YES'))
        no = float(line.count('NO'))
        print(line.count('NA'), line.count('YES'), line.count('NO'))
