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


#Read coupling table
open(LLD,"< /groups/umcg-lld/tmp03/projects/quantification/lldeepIDs_to_flowcellIDs.txt") || die "can’t open LLD file\n";
my %lld;
while (my $lin = <LLD>){
    chomp($lin);
    my @array=split("\t", $lin);
    my $LLDid = $array[0];
    my $FCid = $array[1];
    $lld{ $FCid } = $LLDid;
}
close(LLD);



my $file = "/groups/umcg-bios/tmp03/projects/outlierGeneASE/geneExpressionTables/GRCh37/TPM/counts_genes_BiosFreeze2Only.txt.gz";
#my $file = "./testInput.txt";

print "Processing COUNTS file: $file ..\n";
if ($file =~ /.gz$/) {
    open(COUNTS, "gunzip -c $file |") || die "can’t open pipe to $file\n";
}else {
    open(COUNTS, $file) || die "can’t open $file\n";
}

#Process counts file
my @COUNTSfile = <COUNTS>;
close(COUNTS);
#my $COUNTShead = $COUNTSfile[0];
my $COUNTShead = `zcat $file | head -n1`;
chomp($COUNTShead);
my @COUNTSheader = split("\t", $COUNTShead);

my %counts;
my @allGenes;
#Iterate over lines, extract geneID and add value to corresponding sample in hash
for (my $k=1; $k<=$#COUNTSfile; $k++){
    my $line = $COUNTSfile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    push(@allGenes, $geneID);
    for (my $l=1; $l<=$#array; $l++){ #Iterate over all values in line
        my $val = $array[$l];
        my $sample = $COUNTSheader[$l];
        $counts{ $geneID }{ $sample } = $val;
    }
}
print "Done processing COUNTS file\n";


print "Processing 2nd COUNTS file ..\n";
open(COUNTS2, "gunzip -c /groups/umcg-lld/tmp03/projects/quantification/results/LLdeep/kallisto/tpm_GENES_sampleAdded.txt.gz |") || die "can’t open pipe to /groups/umcg-lld/tmp03/projects/quantification/results/LLdeep/kallisto/tpm_GENES_sampleAdded.txt.gz\n";
my @COUNTS2file = <COUNTS2>;
close(COUNTS2);
my $COUNTS2head = $COUNTS2file[0];
chomp($COUNTS2head);
my @COUNTS2header = split("\t", $COUNTS2head);

#Iterate over lines, extract geneID and add value to corresponding sample in hash
for (my $r=1; $r<=$#COUNTS2file; $r++){
    my $line = $COUNTS2file[$r];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    for (my $l=1; $l<=$#array; $l++){ #Iterate over all values in line
        my $val = $array[$l];
        my $sample = $COUNTS2header[$l];
        $counts{ $geneID }{ $sample } = $val;
    }
}
print "Done processing 2nd COUNTS file\n";



#Iterate over phenotype data and extract biobank
print "Processing PHENO file ..\n";
open(PHENO, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/AllSamplesExcludingCODAMand4outliers.phenotypes.txt") || die "Can't open input PHENO file!\n";
my @PHENOfile = <PHENO>;
close(PHENO);
my $PHENOhead = $PHENOfile[0];
chomp($PHENOhead);
my @PHENOheader = split("\t", $PHENOhead);

my %pheno;
my %phenoMatch;
for (my $p=1; $p<=$#PHENOfile; $p++){
    my $line = $PHENOfile[$p];
    chomp($line);
    my @array=split("\t", $line);
    my $sample = $array[0];
    my $biobank = $array[5];
    my $lldID = $array[115];
    $pheno{ $sample } = $biobank;
    $phenoMatch{ $sample } = $lldID;
}
print "Done processing PHENO file ..\n";


#Interate over logFoldChange file to get the same sample and gene order in this output
print "Processing SAMPLES file ..\n";
open(SAMPLES, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.Bonferroni.samplesFILTERED.values.txt") || die "Can't open input SAMPLES file!\n";
#open(SAMPLES, "< ./testInputSamples.txt") || die "Can't open input SAMPLES file!\n";
my @SAMPLESfile = <SAMPLES>;
close(SAMPLES);
my $SAMPLEShead = $SAMPLESfile[0];
chomp($SAMPLEShead);
my @SAMPLESheader = split("\t", $SAMPLEShead);

open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneExpressionTables/AllSamplesExcludingCODAMand4outliers.geneCounts.AllGenes.TPM.txt") || die "Can't open OUTPUT file!\n";
print OUTPUT "$SAMPLEShead\n";
foreach my $gene (@allGenes){
    my $geneID = $gene;
    print OUTPUT "$geneID";
    for (my $m=1; $m<=$#SAMPLESheader; $m++){
        my $sample = $SAMPLESheader[$m];
        my $val = $counts{ $geneID }{ $sample };
        my $biobank = $pheno{ $sample };
        if (defined $val && $val ne '') {
            print OUTPUT "\t$val";
        }else{ #Try sample matching via partial sample matching (on flowcell)
            my @samples = split("_", $sample);
            my @reverse = reverse(@samples);
            if ($phenoMatch{ $sample }) { #Sample AC1C40ACXX-1-18_BD2D5MACXX-6-18 somehow is not available in pheno data, prevend throwing error while running the code
                
            }

            my $val = $counts{ $geneID }{ $phenoMatch{ $sample } }; #Try to match with lld ID
            if (defined $val && $val ne '') {
                print OUTPUT "\t$val";
            }else{ #Try partial matching
                #Split sampleID by underscore
                my @array = split("_", $sample);
                my $val = "NA";
                my $FCid = $array[1];
                my $FCidSecondAttempt = $array[0];
                if (exists $lld{ $FCid }) {
                    $val = $counts{ $geneID }{ $lld{ $FCid } };
                }else{
                    $val = $counts{ $geneID }{ $lld{ $FCidSecondAttempt } };
                }

                #print "\n\n$sample\t$FCid\n\n";
                if (defined $val && $val ne '') {
                    print OUTPUT "\t$val";
                }else{
                    print "Sample not found: $sample\t$biobank\n";
                }
            }
            
        }
    }
    print OUTPUT "\n";
}
close(OUTPUT);
print "Done processing SAMPLES file\n";



