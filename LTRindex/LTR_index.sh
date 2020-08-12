#!/bin/bash -l
#SBATCH -J LTR_index_lycPyr2
#SBATCH -o LTR_index_lycPyr2.output
#SBATCH -e LTR_index_lycPyr2.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A snic2018-3-568
#SBATCH -p core
#SBATCH -n 4

FASTA=lycPyr2.fasta
NAME="${FASTA%.*}"
DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/LTR_retriever/lycPyr2

cd $DIR

FINDER=${NAME}_ltr_finder_results.scn

#/home/vpeona/LTR_Finder/source/ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.9 $FASTA > $FINDER

ml bioinfo-tools
ml GenomeTools/1.5.9

# TGCA motif finder
HARVEST=${NAME}_ltr_harvest_results.scn

gt suffixerator -db $FASTA -indexname $FASTA -tis -suf -lcp -des -ssp -sds -dna

gt ltrharvest -index $FASTA -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 > $HARVEST

# non TGCA motif finder
NONMOTIF=${NAME}_ltr_harvest_results_nonTGCA.scn

gt ltrharvest -index $FASTA -similar 90 -vic 10 -seed 20 -seqids yes -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 > $NONMOTIF

conda activate LTR_retriever_Vale
ml bioinfo-tools RepeatMasker

mkdir /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/LTR_retriever/$NAME

/home/vpeona/LTR_retriever/LTR_retriever -genome $FASTA -infinder $FINDER -inharvest $HARVEST –nonTGCA $NONMOTIF
