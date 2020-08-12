#!/bin/bash -l
#SBATCH -J SATSUMA2_lycPyr4.1_IL
#SBATCH -o SATSUMA2_lycPyr4.1_IL.output
#SBATCH -e SATSUMA2_lycPyr4.1_IL.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 5-00:00:00
#SBATCH -A snic2017-1-533
#SBATCH -p core
#SBATCH -n 20

module load bioinfo-tools satsuma2/2016-12-07

cp /home/vpeona/b2016267/lycPyr_TheGenome/4.PhaseGenomics_scaffold/38_Scaffolds_default/lycPyr4.1_onlyChr.fasta $SNIC_TMP
cp /home/vpeona/b2016267/lycPyr_Illumina/lycPyr.fa $SNIC_TMP
cd $SNIC_TMP

QUERY=lycPyr.fa
TARGET=lycPyr4.1_onlyChr.fasta
DIR=/home/vpeona/b2016267/lycPyr_TheGenome/5.Syntheny_correction/1.Satsuma/output_IL

mkdir $DIR

SatsumaSynteny2 -q $QUERY -t $TARGET -o $DIR -threads 12
