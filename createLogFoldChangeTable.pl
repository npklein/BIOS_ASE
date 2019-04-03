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


##Read GTF file
print "Processing GTF file..\n";
open(GTF, "< /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf") || die "Can't open GTF file\n";
my %geneName2geneID;
my %geneID2geneName;
my $geneCount=0;
my %geneChr;
while (my $line=<GTF>){ #Read GTF file line by line
    chomp($line);
    if ($line !~ m/^#.+/gs) { #If not header line process further
        my @array = split("\t", $line);
        my $chr = $array[0];
        $chr =~ s/X/23/gs;
        $chr =~ s/Y/24/gs;
        $chr =~ s/MT/25/gs;
        my $feature = $array[2];
        my $start = $array[3];
        my $stop = $array[4];
        my $info = $array[8];
        my $geneID = "NA";
        my $gene = "NA";
        if ($feature eq "gene") {
            if (looks_like_number($chr)) { #Check if chromosome is a number
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})"; gene_name "(.+)"; gene_sourc.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $gene = $2;
                    $geneName2geneID{ $gene } = $geneID;
                    $geneID2geneName{ $geneID } = $gene;
                    $geneChr{ $geneID } = $chr;
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


#Read CGD file
print "Processing CGD file ..\n";
open(CGD, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/CGD.20171220.txt") || die "Can't open CGD file\n";
my %geneInheritance;
my %geneGroup;
my %geneCategory;
while (my $line = <CGD>) {
    next if $. == 1;
    chomp($line);
    my @array = split("\t", $line);
    my $gene = $array[0];
    my $inheritance = $array[4];
    my $group = $array[5];
    my $category = $array[8];
    my $geneID;
    if (exists $geneName2geneID{ $gene }) {
        $geneID = $geneName2geneID{ $gene };
        $geneInheritance{ $geneID } = $inheritance;
        $geneGroup{ $geneID } = $group;
        $geneCategory{ $geneID } = $category;
    }
}
close(CGD);
print "Done processing CGD file.\n";


my @geneAE_files = glob("/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr*.txt");
print "\n# Files to process: " . ($#geneAE_files+1) . "\n\n";

my %geneSampleLogfold;
my %geneSampleCov;
my %genes;
my $count = 0;
foreach my $file (@geneAE_files) { #Iterate over all geneAE files
    my $fileName = basename($file);
    #Read file into array
    open(ALL, "< $file") || die "Can't open file: $file\n";
    my @all=<ALL>;
    close(ALL);
    #If file contains more than 1 line the alternative sampleID can be retrieved from the lines
    if ($#all > 0) {
        #More than one line, so process the file
        for (my $i=1; $i<=$#all; $i++){ #Loop over lines in allelic_counts file
            my $line = $all[$i];
            my @array = split("\t", $line);
            my $gene = $array[3];
            my $aCount = $array[4];
            my $bCount = $array[5];
            my $totCount = $array[6];
            my $logFold = $array[7];
            my $sample = $array[13];
            $sample =~ s/-lib1//gs;
            if ($aCount >= 1 && $bCount >= 1) { # aCount and bCount both >= 1
                $geneSampleLogfold{ $gene }{ $sample } = $logFold;
                $geneSampleCov{ $gene }{ $sample } = $totCount;
                $genes{ $gene }++;
            }
        }
    }
    #Print counter every 1000 lines processed
    if ($count % 1 == 0){
        print "Processed $count files ..\n"
    }
    $count++;
    undef(@all);
}


#Process BINOM filtered file
print "Processing BINOM filtered file ..\n";
# input from figure_2_and_3/select_samples_from_outlierTable.py
my $binomFile = "/groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.txt";
#my $binomFile = "./testInput.txt";
open(BINOM, "< $binomFile") || die "Can't open file: $binomFile\n";
my @BINOMfile = <BINOM>;
close(BINOM);
my $BINOMhead = $BINOMfile[0];
chomp($BINOMhead);
my @BINOMheader = split("\t", $BINOMhead);
print "Done processing BINOM filtered file\n";

print "Creating output files ..\n";
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.values.txt") || die "Can't open outputfile\n";
print OUTPUT "$BINOMhead\n";

my %outliers;
my %nonOutliers;
#Loop through lines of file
for (my $i=1; $i<=$#BINOMfile; $i++){
    my $line = $BINOMfile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    print OUTPUT "$geneID";
    for (my $j=1; $j<=$#array; $j++){ #Loop over elements in array
        my $sample = $BINOMheader[$j];
        my $currentValue = $array[$j];
        my $logFoldChange = "NA";
        if (exists $geneSampleLogfold{ $geneID }{ $sample }) { #Check if gene/sample combination exists
            $logFoldChange = $geneSampleLogfold{ $geneID }{ $sample };
            if ($logFoldChange eq "Inf" || $logFoldChange eq "-Inf"){ #If logfold value is infinite put it on NA
                $logFoldChange = "NA";
            }
        }
        #Checks if BINOM input table is correct
        if (looks_like_number($logFoldChange)) { #If the logFold has a value it should be a YES or NO in the binom table
            if ($currentValue eq "NA") {
                print "ERROR: Sample $sample has a logFoldChange but NA in binom table for gene: $geneID!\n";
                exit("ERROR\n");
            }
        }else{ #logFold is NA
            if ($currentValue eq "YES" || $currentValue eq "NO") { #If currentValue is Yes or No it was counted, so something wrong
                print "ERROR: Sample $sample has a logFoldChange of NA but a YES/NO in binom table for gene: $geneID!\n";
                exit("ERROR\n");
            }
        }
        print OUTPUT"\t$logFoldChange";
        if ($currentValue eq "YES") {
            $outliers{ $geneID }++;
        }elsif($currentValue eq "NO"){
            $nonOutliers{ $geneID }++;
        }
    }
    print OUTPUT "\n";
}
close(OUTPUT);


open(STATS, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltere.BINOM.samplesFILTERED.stats.txt") || die "Can't open stats file\n";
print STATS "ENSEMBLID\tGENE\tCHR\tINHERITANCE\tGROUP\tCATEGORY\tNNONOUTLIERS\tNOUTLIERS\tNCATEGORIES\n";

#Generate statistics table
for (my $k=1; $k<=$#BINOMfile; $k++){
    my $line = $BINOMfile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    my $geneName = $geneID2geneName{ $geneID };
    my $chr = $geneChr{ $geneID };
    my $inheritance = "NA";
    my $group = "NA";
    my $category = "NA";
    my $outliers = 0;
    my $nonOutliers = 0;
    if (exists $geneInheritance{ $geneID }) {
        $inheritance = $geneInheritance{ $geneID };
        $group = $geneGroup{ $geneID };
        $category = $geneCategory{ $geneID };
    }
    if (exists $outliers{ $geneID }) {
        $outliers = $outliers{ $geneID };
    }
    if (exists $nonOutliers{ $geneID }) {
        $nonOutliers = $nonOutliers{ $geneID };
    }
    $inheritance =~ s/ //g;
    my @categoryArray = split(";",$category); #If multiple categories split them and write away per line
    my $numCategories = $#categoryArray+1;
    foreach my $cat (@categoryArray) {
        $cat =~ s/ //g;
        print STATS "$geneID\t$geneName\t$chr\t$inheritance\t$group\t$cat\t$nonOutliers\t$outliers\t$numCategories\n";
    }
}
close(STATS);
print "Done creating output files\n";


