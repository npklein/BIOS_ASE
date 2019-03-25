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


my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";

my $nSamples = "100"; #Minimum number of samples needed for outlier analysis


#Read sample table containing samples to include in analysis and create array
print "Reading sample inclusion file..\n";
open(SAMPLES, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS1000.depthFiltered.bonferroni.txt") || die "Can't open file: /groups/umcg-bios/tmp03/projects/outlierGeneASE/samples_NOUTLIERS500.depthFiltered.bonferroni.txt\n";
my @samples;
while (my $line = <SAMPLES>) {
    chomp($line);
    my $sample = $line;
    push(@samples, $sample);
}
print "Done reading sample inclusion file..\n";


my @geneAE_files = glob("/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing/results/phasing/geneAE_metaGenes/merged/chr*/*.chr*.geneAE.txt");
print "\n# Files to process: " . ($#geneAE_files+1) . "\n\n";

my %geneSampleLogfold;
my %geneSampleCov;
my %genes;
my $count = 0;
foreach my $file (@geneAE_files) { #Iterate over all geneAE files
    my $fileName = basename($file);
    #Extract part of filename to use for matching/coupling samples further downstream
    my $fID;
    if ($fileName =~ m/(.+).chr.+.geneAE.txt/gs) { #Extract samplename from allelic_counts filename
        $fID = $1;
        $fID =~ s/-lib1//gs;
    }
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
            if ($totCount >= 30 && $aCount >= 5 && $bCount >= 5) { #total coverage >= 30 && aCount and bCount both >= 5
                $geneSampleLogfold{ $gene }{ $fID } = $logFold;
                $geneSampleCov{ $gene }{ $fID } = $totCount;
                $genes{ $gene }++;
            }
        }
    }
    #Print counter every 1000 lines processed
    if ($count % 1000 == 0){
        print "Processed $count files ..\n"
    }
    $count++;
    undef(@all);
}

print "\n\n\nGenerating output tables ..\n\n";
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.AllSamples.SD3.table.txt") || die "Can't open outputfile!\n";
open(OUTLIER, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.AllSamples.table.outliers.SD3.min100samples.txt") || die "Can't open outlier outputfile!\n";
open(SAMPLECOUNT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/logFoldChangeTables/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing.logFoldChange.depthFiltered.AllSamples.table.sampleCount.SD3.min100samples.txt") || die "Can't open sampleCount outputfile!\n";

my $sampleString = join("\t", @samples);
my $outputHeader = "ENSEMBLID\t$sampleString";
print OUTPUT "$outputHeader\tNSAMPLES\tMEDIAN\tMAD\n";
print OUTLIER "$outputHeader\tNSAMPLES\tNOUTLIERSAMPLES\n";
print SAMPLECOUNT "SAMPLE\tNGENES\tNOUTLIERS\n";

#Iterate over all genes and about values per sample
my %nGenes;
my %nOutliers;
foreach my $gene (sort keys %genes){
    print OUTPUT "$gene";
    print OUTLIER "$gene";
    my $outlierCount = 0;
    #Iterate over all samples
    my @calcMeanAndSD;
    my $mean;
    my $sd;
    my $mad;
    my $median;
    foreach my $sample (@samples){
        my $logFoldChange = "NA";
        if (exists $geneSampleLogfold{ $gene }{ $sample }) { #If logfold value is present continue
            $logFoldChange = $geneSampleLogfold{ $gene }{ $sample };
            if (looks_like_number($logFoldChange) && $logFoldChange ne "Inf" && $logFoldChange ne "-Inf"){ #If logfold is a number and not infinite, continue
                push(@calcMeanAndSD, $logFoldChange); #Push logFold in array
            }
        }
        if ($logFoldChange eq "Inf" || $logFoldChange eq "-Inf"){ #If logfold value is infinite put it on NA
            $logFoldChange = "NA";
        }
        print OUTPUT "\t$logFoldChange";
    }

    if (@calcMeanAndSD) { #Calculate mean and SD 
        $mean = mean(@calcMeanAndSD);
        $sd = stddev(@calcMeanAndSD);
        $mad = Statistics::Robust::Scale::MAD(\@calcMeanAndSD);
        $median = median(@calcMeanAndSD);
        print OUTPUT "\t" . ($#calcMeanAndSD+1) . "\t$median\t$mad";
    }
    print OUTPUT "\n";
    
    #Process mean and SD outcome to create new binary table depicting outliers
    foreach my $sample (@samples){
        my $outlier = "NA";
        my $logFoldChange = "NA";
        if (($#calcMeanAndSD+1) >= $nSamples) {

            if (exists $geneSampleLogfold{ $gene }{ $sample }) { #If logfold value is present continue
                $logFoldChange = $geneSampleLogfold{ $gene }{ $sample };
                if (looks_like_number($logFoldChange) && $logFoldChange ne "Inf" && $logFoldChange ne "-Inf"){ #If logfold is a number and not infinite, continue
                    $outlier = "NO";
                    $nGenes{ $sample }++;
                    #Check if logFold is an outlier
                    my $lowLimit = ( $median-(3*$mad) ); #Outlier when it's 3 stdDev outside mean
                    my $upLimit = ( $median+(3*$mad) );
                    if ($logFoldChange <= $lowLimit || $logFoldChange >= $upLimit) {
                        $outlier = "YES";
                        $outlierCount++;
                        $nOutliers{ $sample }++; 
                    }
                }
            }
        }
        print OUTLIER "\t$outlier";
    }
    print OUTLIER "\t" . ($#calcMeanAndSD+1) . "\t$outlierCount\n";
    undef(@calcMeanAndSD);
}
close(OUTPUT);
close(OUTLIER);


#Print output sample statistics
foreach my $sample (@samples){
    my $nGenes = "NA";
    my $nOutliers = "NA";
    if (exists $nGenes{ $sample }) {
        $nGenes = $nGenes{ $sample };
    }
    if (exists $nOutliers{ $sample }) {
        $nOutliers = $nOutliers{ $sample };
    }
    print SAMPLECOUNT "$sample\t$nGenes\t$nOutliers\n";
}

print "Done generating output tables\n\n";

close(SAMPLECOUNT);

