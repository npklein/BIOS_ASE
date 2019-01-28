#!/usr/bin/perl -w
use strict;
use warnings;
use diagnostics;
use List::Util qw(min max);
use List::Util qw(sum);
use File::Glob ':glob';
use File::Basename;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);
use lib '/home/umcg-fvandijk/perl_modules/';
use lib '/home/umcg-fvandijk/perl_modules/lib/perl5/';
use Statistics::Basic qw(:all);


my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";

#Read old MDB information
open(MDB, "< /groups/umcg-bios/tmp03/projects/bbmriSampleInfo/sampleSheetDirectlyFromMdb26-01-2016.txt") || die "Can't open input MDB file!\n";
my @MDBfile = <MDB>;
close(MDB);
my $MDBheader = $MDBfile[0];
chomp($MDBheader);

my %mdbInfo;
for (my $l=1; $l<=$#MDBfile; $l++){
    my $line = $MDBfile[$l];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[0];
    my $fcID = $array[1];
    $mdbInfo{ $sample } = $fcID;
    $mdbInfo{ $fcID } = $sample;
}


#Read BIOS phenotype information file
open(BIOS, "< ./BIOS_RNASeq_Metadata.180102.txt") || die "Can't open input BIOS file!\n";
my @BIOSfile = <BIOS>;
close(BIOS);

my $BIOSheader = $BIOSfile[0];
chomp($BIOSheader);
my $BIOSheaderConcat = "SAMPLENAMEPHENODB\t$BIOSheader";
my @BIOSh = split("\t", $BIOSheaderConcat);


#Loop over all lines in BIOS file
my %sampleInfo;
my %BIOSIDtoFlowcell;
for (my $i=1; $i<=$#BIOSfile; $i++){
    my $line = $BIOSfile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $biosID = $array[1];
    my $sample = $array[4]; #In the BIOS_RNA_pheno file column index 3 needs to be used
    $sampleInfo{ $sample } = "$line";
    $BIOSIDtoFlowcell{ $biosID } = $sample;
}


#Read BAM Linking file
open(BAM, "< /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/individual_bam_link.txt") || die "Can't open input BAMlinking file!\n";
my @BAMfile = <BAM>;
close(BAM);

my $BAMheader = $BAMfile[0];
chomp($BAMheader);

#Loop over all lines in BAM file
my %bams;
my %bamInfo;
for (my $m=1; $m<=$#BAMfile; $m++){
    my $line = $BAMfile[$m];
    chomp($line);
    my @array = split(",", $line);
    my $fcID = $array[0];
    my $sample = $array[1];
    my $bam = $array[2];
    $bamInfo{ $fcID } = $sample;
    $bamInfo{ $sample } = $fcID;
    $bams{ $fcID } = $bam;
    $bams{ $sample } = $bam;
}



#Read LLDeep phenotype information file
open(LLD, "< /groups/umcg-lld/tmp03/phenotype_data/LLD_data.txt") || die "Can't open input BIOS file!\n";
my @LLDfile = <LLD>;
close(LLD);

my $LLDheader = $LLDfile[0];
chomp($LLDheader);

#Loop over all lines in LLDeep file
for (my $j=1; $j<=$#LLDfile; $j++){
    my $line = $LLDfile[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[0];
    $sample =~ s/"//gs;
    $sampleInfo{ $sample } = "$line";
}


#Read corrected LLDeep phenotype table
open(LLDEEP, "< ./LLD_all_pheno_for_BIOS_correctDates.txt") || die "Can't open input correct LLDeep phenotype file!\n";
my @LLDEEPfile = <LLDEEP>;
close(LLDEEP);
my $LLDEEPheader = $LLDEEPfile[0];
chomp $LLDEEPheader;
my @LLDEEPheaderArray = split("\t", $LLDEEPheader);


my %LLDEEPinfo;
for (my $m=1; $m<=$#LLDEEPfile; $m++){
    my $line = $LLDEEPfile[$m];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[1];
    $sample =~ s/"//gs;
    $LLDEEPinfo{ $sample } = "$line";
}


#Read file containing related other LLDeep sample
open(LLDEEPALL, "< ./all_sample_ids_table.txt") || die "Can't open input all_sample_ids_table.txt!\n";
my @LLDEEPALLfile = <LLDEEPALL>;
close(LLDEEPALL);

my %LLDeepAll;
for (my $q=1; $q<=$#LLDEEPALLfile; $q++){
    my $line = $LLDEEPALLfile[$q];
    chomp($line);
    my @array = split("\t", $line);
    my $rnaseqID = $array[1];
    my $personID = $array[2];
    $rnaseqID =~ s/"//gs;
    $LLDeepAll{ $rnaseqID } = $personID;
}


#Read statistics table
open(STATS, "< ./BIOS_RNA.stats.N4001.txt") || die "Can't open statsfile!\n";
my @stats = <STATS>;
close(STATS);
my $statsHeader = $stats[0];
chomp($statsHeader);


my %stats;
for (my $o=1; $o<= $#stats; $o++){
    my $line = $stats[$o];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[0];
    $stats{ $sample } = $line;
}



#Loop over input file and add annotations
open(INPUT, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/AllSamples.txt") || die "Can't open inputfile!\n";
my @input = <INPUT>;
close(INPUT);
my $inputHeader = $input[0];
chomp($inputHeader);


open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/phenotypeTables/AllSamples.phenotypes.txt") || die "Can't open outputfile!\n";
print OUTPUT "SAMPLE\tmdbLine\tsampleFlowcell\t$BIOSheaderConcat\t$statsHeader\n";


my $count=0;
my $countLLDnotBIOS=0;
for (my $k=1; $k<= $#input; $k++){
    my $line = $input[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $sample = $array[0];
    
    my $mdbLine;
    my $pheno;
    my $naPrinter=0;
    if (exists $sampleInfo{ $sample }) {
        $mdbLine = $mdbInfo{ $sample };
        $pheno = $sampleInfo{ $sample };
        #print OUTPUT "$line\t$pheno\t" . $stats{ $sample } . "\n";
    }elsif (exists $sampleInfo{ $bamInfo{ $sample } }){
        $mdbLine = $mdbInfo{ $sampleInfo{ $bamInfo{ $sample } } };
        $pheno = $sampleInfo{ $bamInfo{ $sample } };
        #print OUTPUT "$line\t$pheno\t" . $stats{ $sample } . "\n";
    }elsif (exists $mdbInfo{ $sample }){
        #print "$sample\t" . $mdbInfo{ $sample } . "\n";
        if (not exists $sampleInfo{ $mdbInfo{ $sample } }) {
            $mdbLine = $mdbInfo{ $sample };
            if ($mdbLine =~ m/^NTR-(.+[CD])-.+/gs){ #See if this is a NTR sample, if so, try to find twin ID and match phenotype data of twin.
                my $twinID = $1;
                #Grep this full ID from RelationsDB and extract BIOS ID from second column
                #Assumption: Grep returns the hit on "twin" in last column of line, corresponding twin sib ID is in second column
                my $grepCMD = "grep $mdbLine /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/BIOS_Relations_Metadata.180207.txt";
                my $exec = `$grepCMD`;
                if (defined $exec && $exec ne '') { #If exec is defined and has a value attached (grep return in this case)
                    #grep command can return more than 1 value
                    my @greps = split("\n", $exec);
                    my $countLoops=0;
                    my $countNonMatch=0;
                    foreach my $greppedLine (@greps){
                        $countLoops++;
                        chomp($greppedLine);
                        my @array = split("\t", $greppedLine); #split line
                        my $sampleID = $array[0];
                        my $biosID = $array[1];
                        my $biobank = $array[2];
                        my $relationType = $array[5];
                        my $relationID = $array[6]; #In this case the ID that was used to grep
                        if ($mdbLine eq $sampleID && $relationType eq "has monozygotic twin" && $biobank eq "NTR") { #Check if the twin was annotated correct using match with it's own ID
                            #true, so extract BIOSID and use it to convert to flowcellID and extract phenotype information
                            if (exists $BIOSIDtoFlowcell{ $biosID }) {
                                #code
                                my $flowcellID = $BIOSIDtoFlowcell{ $biosID };  
                                $pheno = $sampleInfo{ $flowcellID };
                                #print OUTPUT "$line\t$pheno\t" . $stats{ $sample } . "\n";
                                last;
                            }else{
                                $countNonMatch++;
                                
                            }
                        }elsif ($mdbLine eq $relationID && $relationType eq "has monozygotic twin" && $biobank eq "NTR") { #Check if relationID matches with grepped sample, if the sample is monozygotic twin and is part of NTR cohort
                            #true, so extract BIOSID and use it to convert to flowcellID and extract phenotype information
                            if (exists $BIOSIDtoFlowcell{ $biosID }) {
                                #code
                                my $flowcellID = $BIOSIDtoFlowcell{ $biosID };  
                                $pheno = $sampleInfo{ $flowcellID };
                                #print OUTPUT "$line\t$pheno\t" . $stats{ $sample } . "\n";
                                last;
                            }else{
                                $countNonMatch++;
                                
                            }
                        }else{
                            $naPrinter=1;
                            #print "NTR SAMPLE $mdbLine | grep returned value, but no monozygotic twin hit\n";
                        }
                    }
                    if ($countLoops-1 != $countNonMatch) {
                        #print "NON-EXISTING $mdbLine\t$sample\tFLOWCELL ID\n";
                    }
                    
                }else{
                    $naPrinter=1;
                    #print "NTR SAMPLE $mdbLine | returned empty grep command";
                }
            }else{
                #Try to grep mdbLine in relationsDB file to retrieve BIOS id
                my $grepCMD = "grep $mdbLine /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/BIOS_Relations_Metadata.180207.txt";
                my $exec = `$grepCMD`;
                
                if (defined $exec && $exec ne '') { #If exec is defined and has a value attached (grep return in this case)
                    my @greps = split("\n", $exec);
                    foreach my $greppedLine (@greps){
                        chomp($greppedLine);
                        my @array = split("\t", $greppedLine); #split line
                        my $sampleID = $array[0]; #In this case the ID that was used to grep
                        my $biosID = $array[1];
                        my $biobank = $array[2];
                        my $relationType = $array[5];
                        my $relationID = $array[6];
                        if ($mdbLine eq $sampleID) { #Check grepped sample corresponds to sample ID in this line, if so it's a correct match with itself
                            if (exists $BIOSIDtoFlowcell{ $biosID }) {
                                my $flowcellID = $BIOSIDtoFlowcell{ $biosID };  
                                $pheno = $sampleInfo{ $flowcellID };
                                last;
                            }else{
                                $naPrinter=1;
                                #print "NON-EXISTING $mdbLine\t$sample\tFLOWCELL ID\n";
                            }
                        }else{
                            $naPrinter=1;
                            #print "SAMPLE $mdbLine | grep returned value, but no database hit\n";
                        }
                    }
                }else{
                    $naPrinter=1;
                    #print "SAMPLE $mdbLine | grep returned no value, no results returned yet\n";
                }
            }
            #print "$sample" . "$pheno\n";
        }   
    }else{
        $count++;
        #print "SAMPLE: $sample\t$mdbLine\n";
        if ($bams{ $sample } =~ m/.+lldeepNotInBIOS.+/gs) {
            $countLLDnotBIOS++;
            $naPrinter=2; #LLDeep sample, so different phenotype line to create
        }else{
            #print "CANNOT FIND SAMPLE: $sample\n";
            $naPrinter=1;
        }
    }
    
    if ($sample eq "AD19YAACXX-3-25" || $sample eq "AD19YAACXX-7-19") { #Exclude these two samples from pheno annotation, they are linked to different flowcells as compared to our LLDeep links
        $naPrinter = 3;
    }
    
    #If $pheno is empty fill it with NA values
    if ($naPrinter == 1 || $naPrinter == 2) {
        my $sampleID = "NA";
        my $biosID = "NA";
        my $relationID = "NA";
        my $relationType = "NA";
        my $relatedSampleFound = "NA";
        if (defined $mdbLine && $mdbLine ne '') {
            #
        }else{
            $mdbLine = "NA";
        }
        if ($naPrinter == 1) {
            if ($mdbLine ne 'NA') {
                #Check whether this sample can be found in relations db with respect to parent/child relation. If yes then check if data for fam relation indeed exists
                #Try to grep mdbLine in relationsDB file to retrieve BIOS id
                my $grepCMD = "grep $mdbLine /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/BIOS_Relations_Metadata.180207.txt";
                my $exec = `$grepCMD`;
                
                if (defined $exec && $exec ne '') { #If exec is defined and has a value attached (grep return in this case)
                    my @greps = split("\n", $exec);
                    my %sampleToBios; #Hash to save sample id link to BIOS id per sample for this corresponding grep result (family)
                    foreach my $gLine (@greps){
                        chomp($gLine);
                        my @array = split("\t", $gLine); #split line
                        $sampleID = $array[0]; #In this case the ID that was used to grep
                        $biosID = $array[1];
                        $sampleToBios{ $sampleID } = $biosID;
                    }
                    foreach my $greppedLine (@greps){
                        chomp($greppedLine);
                        my @array = split("\t", $greppedLine); #split line
                        $sampleID = $array[0]; #In this case the ID that was used to grep
                        $biosID = $array[1];
                        my $biobank = $array[2];
                        $relationType = $array[5];
                        $relationID = $array[6];
                        if ($mdbLine eq $sampleID) { #this is the sample itself
                            #check if relationType is one of the following, if so grep the relationID to see whether this individual is present in the database
                            if ($relationType eq "has child" || $relationType eq "has parent" || $relationType eq "has sib" || $relationType eq "has dizygotic twin") {
                                #grep on BIOS ID in metadata DB
                                
                                my $biosID = $sampleToBios{ $relationID };
                                my $grepCMD = "grep $biosID /groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/BIOS_RNASeq_Metadata.180102.txt";
                                my $exec = `$grepCMD`;
                                if (defined $exec && $exec ne '') { #If exec is defined and has a value attached sample is found
                                    chomp($exec);
                                    my @array = split("\t", $exec);
                                    my $sampleID = $array[1];
                                    if ($biosID eq $sampleID) {
                                        $relatedSampleFound = "YES";
                                        last;
                                    }else{
                                        $relatedSampleFound = "NO";
                                    }
                                }else{
                                    $relatedSampleFound = "NO";
                                }
                            }
                        }
                    }
                }
            }
            
            print "Cannot find: $mdbLine\t$sample\t$biosID\t$relationType\t$relationID\t$relatedSampleFound\n";

        }elsif($naPrinter == 2){
            my $LLDeepSample = $bamInfo{ $sample };
            $LLDeepSample =~ s/LL_//gs; #Remove the LL_ prefix from the sample name
            my $LLDeepinfo = $LLDEEPinfo{ $LLDeepSample };
            #Iterate over LLDeepInfo header and line, extract name and value and put in hash
            my %samplePhenoVals;
            if ($mdbLine eq "NA") {
                #Split sample on _
                chomp($sample);
                my @array = split("\t", $sample);
                foreach my $fID (@array){
                    if (exists $LLDeepAll{ $fID }) {
                        $LLDeepSample = $LLDeepAll{ $sample };
                        $LLDeepinfo = $LLDEEPinfo{ $LLDeepSample };
                        if (defined $LLDeepinfo && $LLDeepinfo ne ''){ #If value found break out of loop and continue
                            last;
                        }
                    }
                }
                #print "##UNDEF\n";
            }
            #print "LLDeepSample: $LLDeepSample\t" . $LLDeepAll{ $sample } . "\n";
            if (defined $LLDeepinfo && $LLDeepinfo ne '') {
                #print "##DEFINED\n";
                chomp($LLDeepinfo);
                my @LLDeepinfoArray = split("\t", $LLDeepinfo);
                $pheno = $LLDeepSample;
                for (my $p=0; $p<=$#LLDEEPheaderArray; $p++){
                    $samplePhenoVals{ $LLDEEPheaderArray[$p] } = $LLDeepinfoArray[$p];
                }
                $samplePhenoVals{ "biobank_id" } = "LL";
                for (my $i=1; $i<=$#BIOSh; $i++){ #Loop over pheno names in BIOS metaDB line
                    my $valToPrint = "NA";
                    my $phenoToSearch = $BIOSh[$i];
                    if ($phenoToSearch eq "Sampling_Age") { #Map RNA_BloodSampling_Age to Sampling_Age
                        $phenoToSearch = "RNA_BloodSampling_Age";
                    }
                    
                    if (exists $samplePhenoVals{ $phenoToSearch }) {
                        $valToPrint = $samplePhenoVals{ $phenoToSearch };
                        if ($phenoToSearch eq "LipidMed") { #Convert values to text according to BBMRI specs, described on wiki here: http://www.bbmriwiki.nl/wiki/BIOS_Phenotype
                            #0=no/1=statins/2=yes, but no statins
                            $valToPrint =~ s/0/no/gs;
                            $valToPrint =~ s/1/statins/gs;
                            $valToPrint =~ s/2/yes, but no statins/gs;
                        }
                        if ($phenoToSearch eq "Smoking") {
                            #0=non-smoker/1=former smoker/2=current smoker
                            $valToPrint =~ s/0/non-smoker/gs;
                            $valToPrint =~ s/1/former smoker/gs;
                            $valToPrint =~ s/2/current smoker/gs;
                        }
                        if ($phenoToSearch eq "LDLcholMethod") {
                            #1=Friedewald estimation/2=measured
                            $valToPrint =~ s/1/Friedewald estimation/gs;
                            $valToPrint =~ s/2/measured/gs;
                        }
                        if ($phenoToSearch eq "Sex") {
                            #0=male/1=female
                            $valToPrint =~ s/0/Male/gs;
                            $valToPrint =~ s/1/Female/gs;
                        }
                        if ($phenoToSearch eq "Lipids_BloodSampling_Fasting") {
                            #0=no/1=yes
                            $valToPrint =~ s/0/no/gs;
                            $valToPrint =~ s/1/yes/gs;
                        }
                        
                    }
                    $pheno .= "\t$valToPrint";
                }
                $naPrinter=0; #Set printer counter to zero to prevent it from overwriting with NAs
            }else{
                $pheno = "NA";
                for (my $i=1; $i<=$#BIOSh; $i++){ #Loop over pheno names in BIOS metaDB line
                    $pheno .= "\tNA";
                }
                print "LLDeepNotInBios sample: $mdbLine\t$sample\t$biosID\t$relationType\t$relationID\t$relatedSampleFound\n";
                $naPrinter=2;
            }
            
            
            #Iterate over BIOS metadata header and extract corresponding values and print to output
       }
        
        if ($naPrinter != 0) {
            $pheno = "NA";
            for (my $p=1; $p<=$#BIOSh; $p++){
                if ($naPrinter == 2  && $p == 2) {
                    $pheno .= "\tLLDeepNotInBIOS";
                }else{
                    $pheno .= "\tNA";
                }
            }
        }
    }
    
    if ($naPrinter == 3){ #To correct for both discordant LLDeep annotated samples
        $pheno = "NA";
        for (my $p=1; $p<=$#BIOSh; $p++){
            if ($p == 2 ) {
                $pheno .= "\tLL";
            }else{
                $pheno .= "\tNA";
            }
        }
    }
    
    if (defined $mdbLine && $mdbLine ne '') {
        #
    }else{
        $mdbLine = "NA";
    }
    
    print OUTPUT "$line\t$mdbLine\t$sample\t$pheno\t" . $stats{ $sample } . "\n";
    
}
close(OUTPUT);

print "\n\nCount: $count\n";
print "Of which $countLLDnotBIOS not in BIOS\n\n";



#/groups/umcg-bios/tmp03/projects/bbmriSampleInfo/sampleSheetDirectlyFromMdb26-01-2016.txt
#/groups/umcg-bios/tmp03/projects/bbmriSampleInfo/sample_id_conversion_18-11-2015.txt
#/groups/umcg-bios/tmp03/projects/phenotypes/BIOS_RNA_pheno.txt
#/groups/umcg-bios/tmp03/projects/phenotypes/freeze2_complete_GTE_Groningen_Mdb_07092016.txt
#/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing/individual_bam_link.txt


