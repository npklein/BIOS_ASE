# Take output from phASER that we have over several chunks and merge per chromosome and sample

echo "SAFETY EXIT, check if you really want to run this (removes/overwrites some data)"
echo "If so, remove the exit from the code"
exit

BASE="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/"

#For each chromosome do
for CHR in {1..22}
do

    echo "Processing chr$CHR .."

    #Create output dir
    OUTPUTDIR="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/readbackedPhasing/allelic_counts_mergedPerSample/chr$CHR/"
    mkdir -p $OUTPUTDIR

    #Change dir
    cd ./chr$CHR/

    echo "Extracting *.tar.gz files .."

    #For every directory do
    for DIR in $( ls . )
    do

        #Change dir
        cd $DIR

        #Unarchive
        tar -xvzf allelic_counts.tar.gz

        cd -

    done

    cd -

    echo "Done extracting *.tar.gz files .."

    #Now cat all files per sample
    SAMPLEFILE="/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/samples.txt"

    echo "Merging all allelic count files per sample"

    #Change dir
    cd $BASE/chr$CHR/

    #For every sample do
    while read line
    do

        SAMPLE=$line

        echo "Merging files for sample $SAMPLE .."

        #Echo header in output file
        echo -e -n "contig\tposition\tvariantID\trefAllele\taltAllele\trefCount\taltCount\ttotalCount\n" > $OUTPUTDIR/BIOS_LLDeep_noRNAeditSites_phASER.$SAMPLE.chr$CHR.allelic_counts.txt

        #For every region open allelic counts file and do something with it
        for REGION in $( ls . )
        do

            #If exists, otherwise skip this file
            FILE="$BASE/chr$CHR/$REGION/allelic_counts/BIOS_LLDeep_noRNAeditSites_phASER.$SAMPLE.chr$CHR.$REGION.allelic_counts.txt"
            if [ -f $FILE ]
            then
                sed '1,1d' $FILE >> $OUTPUTDIR/BIOS_LLDeep_noRNAeditSites_phASER.$SAMPLE.chr$CHR.allelic_counts.txt
            fi

        done

    done<$SAMPLEFILE

    cd -

    cd "$BASE"

    #Now remove the untarred files to clean up storage
    cd ./chr$CHR/

    echo "Removing untarred files .."

    #For every directory do
    for DIR in $( ls . )
    do

        #Change dir
        cd $DIR

        #Unarchive
        rm -r ./allelic_counts/

        cd -

    done

    cd -

    cd "$BASE"

    echo "Done removing untarred files .."

done

