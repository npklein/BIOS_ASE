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


#Read CGD file
print "Processing CGD file..\n";
open(CGD, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/geneAndVariantLists/CGD.20171220.txt") || die "Can't open CGD file\n";
my @CGDfile = <CGD>;
my %inheritance;
my %ageGroup;
my %manifestation;
my %CGDgenes;
for (my $j=1; $j<=$#CGDfile; $j++){
    my $line = $CGDfile[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $gene = $array[0];
    my $inheritance = $array[4];
    $inheritance =~ s/ //gs;
    my $ageGroup = $array[5];
    my $manifestation = $array[7];
    $inheritance{ $gene } = $inheritance;
    $ageGroup{ $gene } = $ageGroup;
    $manifestation{ $gene } = $manifestation;
    $CGDgenes{ $gene } = $gene;
}
close(CGD);
print "Done processing CGD file.\n\n";


#Read CLINVAR file
print "Processing CLINVAR file..\n";
my %clinvar;
open(CLNVR, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/clinvar/clinvar_data/clinvar_20180128.criteriaProvided_multiple_submitters_no_conflicts.txt") || die "Can't open CLINVAR file\n";
my @CLNVRfile = <CLNVR>;
for (my $m=1; $m<=$#CLNVRfile; $m++){
    my $line = $CLNVRfile[$m];
    chomp($line);
    my @array = split("\t", $line);
    my $chr = $array[1];
    my $pos = $array[2];
    my $ref = $array[3];
    my $alt = $array[4];
    my $clnsig = $array[15];
    my $variant = "$chr\_$pos\_$ref\_$alt";
    $clinvar{ $variant } = $clnsig;
}
close(CLNVR);
print "Done processing CLINVAR file.\n\n";

#Genes with high pLI scores (pLI <B3> 0.9) are extremely LoF intolerant, whereby genes with low pLI scores (pLI <B2> 0.1) are LoF tolerant.

print "Processing annotation tables ..\n";
my @annotationFiles = glob('/groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/BIOS_LLDeep_noRNAeditSites_phASER.snpEff.closest.VEP.chr*.addedExACandGONLAlleleFrequency.addedpLI.annotation.table');
#Loop through these files to extract annotation information per variant
my %annotations;
my %impacts;
my %exacs;
my %gonls;
my %plis;
my %geneplis;
my %geneImpacts;
my %geneImpactHigh;
my %geneImpactModerate;
my %geneImpactModifier;
my %geneImpactLow;
my %genes;
my %varToGenes;
my $totalGenes=0;
my $tolGenes=0;
my $intolGenes=0;
my $nonGenes=0;
my $naGenes=0;
foreach my $anFile (@annotationFiles){
    open(ANNO, "< $anFile") || die "Can't open annotation file: $anFile\n";
    my @annotationFile = <ANNO>;
    close(ANNO);
    for (my $i=1; $i<=$#annotationFile; $i++){
        my $line = $annotationFile[$i];
        chomp($line);
        my @array = split("\t", $line);
        my $chr = $array[0];
        my $pos = $array[1];
        my $rsID = $array[2];
        my $ref = $array[3];
        my $alt = $array[4];
        my $annotation = $array[8];
        my $impact = $array[9];
        my $gene = $array[10];
        my $vepConsequence = $array[32];
        my $vepImpact = $array[33];
        my $vepFeatureType = $array[36];
        my $vepBiotype = $array[38];
        $gene =~ s/(?>\x0D\x0A?|[\x0A-\x0C\x85\x{2028}\x{2029}])//;
        my $gonlAF = $array[140];
        my $exacAF = $array[141];
        my $pLI = $array[142];
        my $key = "$chr\_$pos\_$ref\_$alt";
        $annotations{ $key } = "$chr\t$pos\t$rsID\t$ref\t$alt\t$annotation\t$impact\t$gene\t$vepConsequence\t$vepImpact\t$vepFeatureType\t$vepBiotype\t$gonlAF\t$exacAF\t$pLI";
        $impacts{ $key } = $impact;
        $exacs{ $key } = $exacAF;
        $gonls{ $key } = $gonlAF;
        $plis{ $key } = $pLI;
        $geneplis{ $gene } = $pLI;
        if (not exists $genes{ $gene }) {
            $totalGenes++;
            if (looks_like_number($pLI)) {
                if ($pLI >= 0.9) {
                    $intolGenes++;
                }elsif ($pLI <= 0.1){
                    $tolGenes++;
                }else{
                    $nonGenes++;
                }
            }else{
                $naGenes++;
            }
        }
        $varToGenes{ $key } = $gene;
        $genes{ $gene } = $gene;

        $geneImpacts{ $gene }{ $impact }++; #Counts number of impacts per gene
        #print "$key\t$annotation\t$impact\n";
    }
}
print "#Total genes: $totalGenes\n";
print "#Tolerant genes: $tolGenes\n";
print "#Intolerant genes: $intolGenes\n";
print "#Neither tolerant or intolerant genes: $nonGenes\n";
print "#NA genes: $naGenes\n";

print "Done processing annotation tables.\n\n";


#Read AF genotype data file
print "Processing AF file..\n";
open(AF, "< /groups/umcg-bios/tmp03/projects/outlierGeneASE/annotatedWith.snpEff.closest.VEP/chrALL.AFsFromData.txt") || die "Can't open AF file\n";
my @AFfile = <AF>;
my %afs;
for (my $j=0; $j<=$#AFfile; $j++){
    my $line = $AFfile[$j];
    chomp($line);
    my @array = split("\t", $line);
    my $key = $array[0];
    my $af = $array[1];
    $afs{ $key } = $af;
}
close(AF);
print "Done processing AF file.\n\n";


#my $hapFile = "./chr$CHR.TEST.txt";
my %varRefCount;
my %varAltCount;
my %varStatusCount;
my %varSampleRefCount;
my %varSampleAltCount;
print "\nGenerating output file ..\n";
    
open(OUTPUT, "> /groups/umcg-bios/tmp03/projects/outlierGeneASE/variantPenetranceAndPLIAnalysis/AnnotationTable.chrALL.txt") || die "Can't open output file for chrALL!\n";
print OUTPUT "VARIANT\tCHR\tPOS\tRSID\tREF\tALT\tSNPEFFANNOTATION\tSNPEFFIMPACT\tSNPEFFGENE\tVEPCONSEQUENCE\tVEPIMPACT\tVEPFEATURETYPE\tVEPBIOTYPE\tGONLAF\tEXACAF\tPLI\tAF\t";
print OUTPUT "PLISTATUS\tGENEID\tGENENAME\tCGDINHERITANCE\tCGDAGEGROUP\tCGDMANIFESTATION\tCLNVRSIG\n";
foreach my $key (sort keys %annotations){
    my $variant = $key; #Variant
    my $af = $afs{ $variant }; #Allele frequency from data
    my $exacAF = $exacs{ $variant };
    my $gonlAF = $gonls{ $variant };
    my $clinvar = "NA";
    if (exists $clinvar{ $variant }) { #Extract Clinvar status
        $clinvar = $clinvar{ $variant };
    }
    my $annotation = $annotations{ $variant }; #Row of annotations
    my $pLI = "NA";
    my $pliStat = "NA";
    #Set pLI status
    if (exists $plis{ $variant }) {
        $pLI = $plis{ $variant }; #pLI score
        if (looks_like_number($pLI)) {
            if ($pLI >= 0.9) {
                $pliStat = "INTOLERANT";
            }elsif ($pLI <= 0.1){
                $pliStat = "TOLERANT";
            }else{
                $pliStat = "NEITHER";
            }
        }else{
            $pliStat = "NA";
        }
    }
    #Set CGD status
    my $geneID = "NA";
    if (exists $geneName2geneID{ $varToGenes{ $key } }) {
        $geneID = $geneName2geneID{ $varToGenes{ $key } };
    }
    my $inheritance = "NA";
    my $ageGroup = "NA";
    my $manifestation = "NA";
    my $geneName = "NA";
    #print "geneID; $geneID\n";
    if (exists $varToGenes{ $key }) {
        $geneName = $varToGenes{ $key };
        #print "geneName: $geneName\n";
        if (exists $CGDgenes{ $geneName }) {
            $inheritance = $inheritance{ $geneName };
            $ageGroup = $ageGroup{ $geneName };
            $ageGroup =~ s|N/A|NA|g;
            $manifestation = $manifestation{ $geneName };
        }
    }
    print OUTPUT "$variant\t$annotation\t$af\t$pliStat\t$geneID\t$geneName\t$inheritance\t$ageGroup\t$manifestation\t$clinvar\n";
}
close(OUTPUT);
print "Done generating output file for chrALL.\n\n";
