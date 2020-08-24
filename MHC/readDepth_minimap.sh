#!/bin/bash
#SBATCH -J MINIMAP2_lycPyr7.4
#SBATCH -o MINIMAP2_lycPyr7.4.output
#SBATCH -e MINIMAP2_lycPyr7.4.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 7-00:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools minimap2 samtools

READS=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_PacBio/subreads/filtered_subreads.fastq

cp /proj/sllstore2017073/private/BirdReferences/lycPyr7.4.fasta $SNIC_TMP/lycPyr7.4.fa
REF=lycPyr7.4.fa

ALIGN=lycPyr7.4_minimap2

DIR=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Minimap/

cd $SNIC_TMP

minimap2 -ax map-pb  $REF $READS > $ALIGN.sam

samtools view -b $ALIGN.sam -o $ALIGN.bam

cp $ALIGN.bam $DIR
