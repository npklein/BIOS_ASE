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
my %varAF;
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
    my $gonlaf = $array[13];
    my $exacaf = $array[14];
    my $af = $array[16];
    my $geneID = "NA";
    if (exists $geneName2geneID{ $gene }) {
        $geneID = $geneName2geneID{ $gene };
    }
    my $key = "$chr\_$pos\_$ref\_$alt";
    $counts{ $var } = "$var\t$chr\t$pos\t$ref\t$alt\t$gene\t$geneID";
    $varAF{ $var } = $af;
}
print "Done processing annotation file.\n\n";


#Read cumulative counts for both ASE and non-ASE variants
print "Processing count tables ..\n";
my @countFiles = glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/counts.matrix.cumulativeVariants.ALLcounts*.chrALL.txt.filtered.txt');
#Loop through these files to extract annotation information per variant
my %varGeneName;
my %varGeneID;
my %varLine;
my %varSampleCount;
my %varMajor;
my %varMinor;
foreach my $countFile (@countFiles){
    open(COUNT, "< $countFile") || die "Can't open annotation file: $countFile\n";
    my @cFile = <COUNT>;
    close(COUNT);
    for (my $i=1; $i<=$#cFile; $i++){
        my $line = $cFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $variant = $array[0];
        my $chr = $array[1];
        my $pos = $array[2];
        my $ref = $array[3];
        my $alt = $array[4];
        my $geneName = $array[5];
        my $geneID = $array[6];
        my $sampleCount = $array[7];
        my $major = $array[8];
        my $minor = $array[9];
        my $varline = "$chr\t$pos\t$ref\t$alt";
        $varGeneName{ $variant } = $geneName;
        $varGeneID{ $variant } = $geneID;
        $varLine{ $variant } = $varline;
        if ($varMajor{ $variant }) { #Already exists, so sum with known
            my $hashVarSampleCount = $varSampleCount{ $variant };
            my $hashVarMajor = $varMajor{ $variant };
            my $hashVarMinor = $varMinor{ $variant };
            $varSampleCount{ $variant } = ($sampleCount + $hashVarSampleCount);
            $varMajor{ $variant } = ($major + $hashVarMajor);
            $varMinor{ $variant } = ($minor + $hashVarMinor);
        }else{
            $varSampleCount{ $variant } = $sampleCount;
            $varMajor{ $variant } = $major;
            $varMinor{ $variant } = $minor;
        }
    }
}
print "Done processing count tables.\n\n";


#Process GTEx files
print "Processing GTEx tables ..\n";
my @gtexFiles = glob('/groups/umcg-bios/tmp03/projects/GTEx_ASE/65436/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/ExpressionFiles/decrypted/phe000024.v1.GTEx_ASE_SNPs.expression-matrixfmt-ase.c1/GTEX*.ase_table.tsv');
my %gtexVars;
my %gtexVarTis;
my %gtexVarTisMaj;
my %gtexVarTisMin;
#Loop through the files
foreach my $gtexFile (@gtexFiles){
    print "Processing file: $gtexFile ..\n";
    open(GTEX, "< $gtexFile") || die "Can't open annotation file: $gtexFile\n";
    my @gFile = <GTEX>;
    close(GTEX);
    for (my $i=1; $i<=$#gFile; $i++){
        my $line = $gFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $chr = $array[0];
        my $pos = $array[1];
        my $variant = $array[2];
        $variant =~ s/_b37//g; #Remove the _b37 from variant ID to make it match with other IDs
        my $ref = $array[3];
        my $alt = $array[4];
        my $tissue = $array[7];
        my $refC = $array[8];
        my $altC = $array[9];
        if (exists $varMajor{ $variant }) { #If exists we observed it, otherwise we didn't
            my $af = $varAF{ $variant };
            my $majorC = $refC;
            my $minorC = $altC;
            my $totalC = ($majorC + $minorC);
            if ($totalC >= 20 && $majorC > 0 && $minorC > 0) { #Only count if on both alleles at least 1 read observerd and total coverage is at least 20

                if ($af > 0.5) { #Swap the counts for major and minor allele
                    $majorC = $altC;
                    $minorC = $refC;
                }
                #Check if counts for this variant already exists for this tissue. If so, sum the counts
                if (exists $gtexVarTisMaj{ $variant }{ $tissue }) {
                    my $hashMajorC = $gtexVarTisMaj{ $variant }{ $tissue };
                    my $hashMinorC = $gtexVarTisMin{ $variant }{ $tissue };
                    $gtexVarTisMaj{ $variant }{ $tissue } = ($majorC + $hashMajorC);
                    $gtexVarTisMin{ $variant }{ $tissue } = ($minorC + $hashMinorC);
                }else{
                    $gtexVarTisMaj{ $variant }{ $tissue } = $majorC;
                    $gtexVarTisMin{ $variant }{ $tissue } = $minorC;
                }
                $gtexVarTis{ $variant }{ $tissue }++;
            
            }
        }
        $gtexVars{ $variant }++;
    }
}
print "Done processing GTEx tables.\n\n";


#Create output table
print "Generating output table ..\n";
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/concordanceGTEx/counts.matrix.AlleleAdded.txt") || die "Can't open output file.\n";
print OUTPUT "VARIANT\tCHR\tPOS\tREF\tALT\tAF\tMINORALLELEGTEX\tGENENAME\tGENEID\tSAMPLECOUNT\tSUMMAJOR\tSUMMINOR\tTOTAL\tRATIO\tTISSUE\tGTEXSAMPLECOUNT\tGTEXSUMMAJOR\tGTEXSUMMINOR\tGTEXTOTAL\tGTEXRATIO\n";
my $inGTEx = 0;
my $notInGTEx = 0;
foreach my $variant (sort keys %varMajor){
    if (exists $gtexVars{ $variant }) { #If exists in GTEx
        $inGTEx++;
        #Loop over all tissues for specific variant
        foreach my $tissue (keys %{ $gtexVarTisMaj{ $variant }}) {
            my $gtexMajorC = $gtexVarTisMaj{ $variant }{ $tissue };
            my $gtexMinorC = $gtexVarTisMin{ $variant }{ $tissue };
            my $gtexTotal = ($gtexMajorC + $gtexMinorC);
            my $gtexRatio = ($gtexMinorC / $gtexTotal);
            
            my $line = $varLine{ $variant };
            my $geneName = $varGeneName{ $variant };
            my $geneID = $varGeneID{ $variant };
            my $varSampleCount = $varSampleCount{ $variant };
            my $major = $varMajor{ $variant };
            my $minor = $varMinor{ $variant };
            my $total = ($major + $minor);
            my $ratio = ($minor/$total);
            
            my $gtexObserved = $gtexVarTis{ $variant }{ $tissue };
            my $af = $varAF{ $variant };
            my @ar = split("\t", $line);
            my $ref = $ar[2];
            my $alt = $ar[3];
            my $minorAll = $alt;
            if ($af > 0.5) {
                $minorAll = $ref;
            }
            print OUTPUT "$variant\t$line\t$af\t$minorAll\t$geneName\t$geneID\t$varSampleCount\t$major\t$minor\t$total\t$ratio\t$tissue\t$gtexObserved\t$gtexMajorC\t$gtexMinorC\t$gtexTotal\t$gtexRatio\n";
        }
    }else{
        $notInGTEx++;
    }
}
close(OUTPUT);
print "Done generating output table.\n\n";

print "\n\nNumber of samples also observed in GTEx: $inGTEx\n";
print "Number of samples not observed in GTEx: $notInGTEx\n\n";
