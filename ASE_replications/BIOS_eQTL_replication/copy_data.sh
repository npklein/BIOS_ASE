
prm="/groups/umcg-bios/prm03/projects/BIOS_EGCUT_for_eQTLGen/BIOS_EGCUT/eqtlpipeline_bios_egcut_backup010517/"
expression="$prm/gene_read_counts_BIOS_and_LLD_passQC.tsv.SampleSelection.ProbesWithZeroVarianceRemoved.TMM.CPM.Log2Transformed.ProbesCentered.SamplesZTransformed.CovariatesRemoved.txt.gz "
rsync -vP $prm/$expression /groups/umcg-bios/tmp03/projects/2020-03-05-ASE-BIOS-eqtl-comparison/data/
