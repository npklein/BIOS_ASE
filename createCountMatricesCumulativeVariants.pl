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


##Read GTF file
print "Processing GTF file..\n";
open(GTF, "< /apps/data/ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf") || die "Can't open GTF file\n";
my %geneName2geneID;
my %geneID2geneName;
my $geneCount=0;
my %geneChr;
my %geneLength;
my %geneStart;
my %geneStop;
while (my $line=<GTF>){ #Read GTF file line by line
    chomp($line);
    if ($line !~ m/^#.+/gs) { #If not header line process further
        my @array = split("\t", $line);
        my $chr = $array[0];
        #$chr =~ s/X/23/gs;
        #$chr =~ s/Y/24/gs;
        #$chr =~ s/MT/25/gs;
        my $feature = $array[2];
        my $start = $array[3];
        my $stop = $array[4];
        my $info = $array[8];
        my $geneID = "NA";
        my $gene = "NA";
        my $geneLength = ($stop-$start);
        if ($feature eq "gene") {
            if (looks_like_number($chr)) { #Check if chromosome is a number
                if ($info =~ m/gene_id "(ENSG[0-9]{1,})"; gene_name "(.+)"; gene_sourc.+/gs) { #Extract Ensembl gene ID
                    $geneID = $1;
                    $gene = $2;
                    $gene =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    $geneID =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
                    #if ($gene eq "PPT1" || $geneID eq "ENSG00000164458") {
                    #    print "$geneID\t$gene\n";
                    #}
                    $geneName2geneID{ $gene } = $geneID;
                    $geneID2geneName{ $geneID } = $gene;
                    $geneChr{ $geneID } = $chr;
                    $geneLength{ $geneID } = $geneLength;
                    $geneStart{ $geneID } = $start;
                    $geneStop{ $geneID } = $stop;
                    $geneCount++;
                }
            }
        }
    }
}
close(GTF);
print "#Genes detected: $geneCount\n";
print "Done processing GTF file.\n\n";


#Read annotation file
print "Processing annotation file..\n";
my %counts;
open(COUNTS, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/AnnotationTable.chrALL.txt") || die "Can't open annotation file\n";
my @COUNTSfile = <COUNTS>;
close(COUNTS);
for (my $k=1; $k<=$#COUNTSfile; $k++){
    my $line = $COUNTSfile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    my $chr = $array[1];
    my $pos = $array[2];
    my $ref = $array[4];
    my $alt = $array[5];
    my $impact = $array[7];
    my $gene = $array[8];
    my $geneID = "NA";
    if (exists $geneName2geneID{ $gene }) {
        $geneID = $geneName2geneID{ $gene };
    }
    my $key = "$chr\_$pos\_$ref\_$alt";
    $counts{ $var } = "$var\t$chr\t$pos\t$ref\t$alt\t$gene\t$geneID";
}
print "Done processing annotation file.\n\n";


my $postfix = "nonASEsamples.chrALL.txt.filtered.txt";


##Read counts file
print "Processing count matrix files..\n";
my %variants;
open(MAJOR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.majorAllelle.$postfix") || die "Can't open major count file.\n";
#open(MAJOR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/test.majorAllele.txt") || die "Can't open major count file.\n";
my @majorFile = <MAJOR>;
close(MAJOR);
my $headerMajor = $majorFile[0];
chomp($headerMajor);
my @headerMAJOR = split("\t", $headerMajor);

open(MINOR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.minorAllelle.$postfix") || die "Can't open minor count file.\n";
#open(MINOR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/test.minorAllele.txt") || die "Can't open minor count file.\n";
my @minorFile = <MINOR>;
close(MINOR);
my $headerMinor = $minorFile[0];
chomp($headerMinor);
my @headerMINOR = split("\t", $headerMinor);

my %majorCounts;
for (my $i=1; $i<=$#majorFile; $i++){
    my $line = $majorFile[$i];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    for (my $j=1; $j<=$#array; $j++){
        my $val = $array[$j];
        my $sample = $headerMAJOR[$j];
        if ($val ne "NA") { #if not NA push in hash
            $majorCounts{ $var }{ $sample } = $val;
            $variants{ $var }++;
        }
    }
}

my %minorCounts;
for (my $k=1; $k<=$#minorFile; $k++){
    my $line = $minorFile[$k];
    chomp($line);
    my @array = split("\t", $line);
    my $var = $array[0];
    for (my $j=1; $j<=$#array; $j++){
        my $val = $array[$j];
        my $sample = $headerMINOR[$j];
        if ($val ne "NA") { #if not NA push in hash
            $minorCounts{ $var }{ $sample } = $val;
            $variants{ $var }++;
        }
    }
}
print "Done processing count matrix files\n\n";


##Read counts file
print "Generating output file..\n";
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.cumulativeVariants.ALLcounts.$postfix") || die "Can't open output file.\n";
print OUTPUT "VARIANT\tCHR\tPOS\tREF\tALT\tGENENAME\tGENEID\tSAMPLECOUNT\tSUMMAJOR\tSUMMINOR\tTOTAL\tRATIO\n";
my %cumMajor;
my %cumMinor;
foreach my $key (sort keys %variants){ #Foreach variant
    my $variant = $key;
    my $sampleCount = 0;
    foreach my $sample (keys %{ $majorCounts{ $variant }}) {
        $sampleCount++;
        #Retrieve current couns for variant
        my $curMajor;
        my $curMinor;
        if (exists $cumMajor{ $variant }) {
            $curMajor = $cumMajor{ $variant };
            $curMinor = $cumMinor{ $variant };
        }else{ #Initialize values at 0
            $curMajor = 0;
            $curMinor = 0;
        }
        #Counts from sample
        my $major = $majorCounts{ $variant }{ $sample };
        my $minor = $minorCounts{ $variant }{ $sample };
        #Add counts to existing value
        $cumMajor{ $variant } = ($curMajor + $major);
        $cumMinor{ $variant } = ($curMinor + $minor);
    }
    #Write output line    
    my $line = $counts{ $variant };
    my $sumMajor = $cumMajor{ $variant };
    my $sumMinor = $cumMinor{ $variant };
    my $total = ($sumMajor+$sumMinor);
    my $ratio = ($sumMinor/$total);
    print OUTPUT "$line\t$sampleCount\t$sumMajor\t$sumMinor\t$total\t$ratio\n";
}
close(OUTPUT);
print "Done generating output file\n\n";


