#!/bin/bash -l
#SBATCH -J ANNIE_BLAST_round2_lycPyrZW
#SBATCH -o ANNIE_BLAST_round2_lycPyrZW.output
#SBATCH -e ANNIE_BLAST_round2_lycPyrZW.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02-00:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 8

# load required modules
ml bioinfo-tools blast/2.9.0+

DB=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/protein_chicken_uniprot.fasta
BLAST=lycPyr_rnd2.all.maker.proteins.blast.out
GFF=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd2.maker.output/lycPyr_rnd2.all.maker.noseq.gff3
OUTPUT=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Annie/lycPyrZW_rnd2.all.maker.proteins.annie

cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Annie

/home/vpeona/Annie/annie.py -b $BLAST -db $DB -g $GFF -o $OUTPUT
