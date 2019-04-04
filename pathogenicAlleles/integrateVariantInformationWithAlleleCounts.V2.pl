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

my $tabixPath="/apps/software/HTSlib/1.3.2-foss-2015b/bin/";





print "\nReading annotation tables ..\n";

my @annotation_files = glob("/groups/umcg-bios/tmp03/projects/genotypes_BIOS_LLDeep_Diagnostics_merged_phasing_noRnaEditing_annotation/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.annotation.table");
my $annoHeader;
my %variantAnnotation;
my %variantAnnotationImpact;
my %variantAll;
foreach my $annotationFile (@annotation_files){
    open(ANNO, "< $annotationFile") || die "Can't open annotation file: $annotationFile!\n";
    my @annotations = <ANNO>;
    close(ANNO);
    $annoHeader = $annotations[0];
    chomp($annoHeader);
    my @annoHeaderArray = split("\t", $annoHeader);
    for (my $l=1; $l<=$#annotations; $l++){
        my $line = $annotations[$l];
        chomp($line);
        my @array = split("\t", $line);
        my $chr = $array[0];
        my $pos = $array[1];
        my $ref = $array[3];
        my $alt = $array[4];
        my $annotation = $array[8];
        my $impact = $array[9];
        my $geneID = $array[11];
        my $variant = "$chr\t$pos\t$ref\t$alt";
        $variantAnnotation{ $variant } = $annotation;
        $variantAnnotationImpact{ $variant } = $impact;
        $variantAll{ $variant } = $line;
    }
    undef(@annotations);
}
print "Done reading annotation table\n\n";


#/groups/umcg-bios/tmp03/users/umcg-ndeklein/ESHG_2018/infoTables//groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.not_outliers.TEST.txt
#gene    snp     AF      ref     alt     0       1       2
#ENSG00000125826 20_389098       0.0004998750312421895   G       A       3294    3       0

print "\nReading Non-outliers table ..\n";
#/groups/umcg-bios/tmp03/users/umcg-ndeklein/ESHG_2018/infoTables/
open(NON, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables//groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.not_outliers.binom.txt") || die "Can't open non-outliers file!\n";
my @nons = <NON>;
close(NON);
my $nonsHeader = $nons[0];
chomp($nonsHeader);
my @nonsHeaderArray = split("\t", $nonsHeader);
print "Done reading Non-outliers table ..\n";


print "\nReading Outliers table ..\n";
#/groups/umcg-bios/tmp03/users/umcg-ndeklein/ESHG_2018/infoTables/
open(OUTL, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/infoTables//groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.outliers.binom.txt") || die "Can't open outlier file!\n";
my @outliers = <OUTL>;
close(OUTL);
my $outlHeader = $outliers[0];
chomp($outlHeader);
my @outlHeaderArray = split("\t", $outlHeader);

open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/pathogenicAlleles/alleleCountPerGroupPerGene.binom.annotated.alleleFiltered.removedCODAMandOutliers.ALL.txt") || die "Can't open outlier file!\n";
print OUTPUT "uniq\tchr\tpos\tref\talt\tgene_id\taf\tNonOutliersTotal\tNonOutliersHomRef\tNonOutliersHet\tNonOutliersHomAlt\tNonOutliersTotalAlleles\tNonOutlierRefAlleles\tNonOutliersAltAlleles\tOutliersTotal\tOutliersHomRef\tOutliersHet\tOutliersHomAlt\tOutliersTotalAlleles\tOutliersRefAlleles\tOutliersAltAlleles\t$annoHeader\n";


for (my $j=1; $j<=$#outliers; $j++){
    my $line = $outliers[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $geneID = $array[0];
    my @snp = split("_", $array[1]);
    my $chr = $snp[0];
    my $pos = $snp[1];
    my $af = $array[2];
    my $ref = $array[3];
    my $alt = $array[4];
    my $homRef = $array[5];
    my $het = $array[6];
    my $homAlt = $array[7];
    my $total = ($homRef + $het + $homAlt);
    my $uniq = "$geneID\_$chr\_$pos\_$ref\_$alt";
    my $variant = "$geneID\t$chr\t$pos\t$ref\t$alt";
    chomp($variant);

    my $nonLine = $nons[$j];
    chomp($nonLine);
    my @nonArray = split("\t", $nonLine);
    my $nonSNP = $nonArray[1];
    
    my $nonHomRef = $nonArray[5];
    my $nonHet = $nonArray[6];
    my $nonHomAlt = $nonArray[7];
    my $nonTotal = ($nonHomRef + $nonHet + $nonHomAlt);
    
    my $annotation = $variantAll{ "$chr\t$pos\t$ref\t$alt" };
    
    my $totAll = ($total*2);
    my $totRef = (($homRef*2) + $het);
    my $totAlt = ($het + ($homAlt*2));
    
    my $totNonAll = ($nonTotal*2);
    my $totNonRef = (($nonHomRef*2) + $nonHet);
    my $totNonAlt = ($nonHet + ($nonHomAlt*2));
    if ($array[1] eq $nonSNP) { #Check if SNP on current line is the exact same SNP as on the line from other file
        print OUTPUT "$uniq\t$chr\t$pos\t$ref\t$alt\t$geneID\t$af\t$nonTotal\t$nonHomRef\t$nonHet\t$nonHomAlt\t$totNonAll\t$totNonRef\t$totNonAlt\t$total\t$homRef\t$het\t$homAlt\t$totAll\t$totRef\t$totAlt\t$annotation\n";
    }else{
        die("SNPs " . $array[1] . " and $nonSNP on line $j don't match!\n");
    }
}
close(OUTPUT);
print "\nDone reading Outliers table ..\n";



