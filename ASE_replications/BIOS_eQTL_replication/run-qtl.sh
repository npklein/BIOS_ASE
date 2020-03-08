#!/bin/bash
#SBATCH --job-name=run-qtl
#SBATCH --output=run-qtl.log
#SBATCH --error=run-qtl.log
#SBATCH --ntasks=1
#SBATCH --time=23:10:00
#SBATCH --mem=100g
#SBATCH --cpus-per-task=10
#SBATCH  --qos=dev

set -e
set -u

ml Java/11.0.2

java -Xmx90g -jar /groups/umcg-biogen/tmp03/tools/eqtl-mapping-pipeline-1.4.8-SNAPSHOT/eqtl-mapping-pipeline.jar \
        --mode metaqtl \
        --settings /groups/umcg-bios/tmp03/projects/2020-03-05-ASE-BIOS-eqtl-comparison/2020-02-18-bios-qtlsettings.xml
