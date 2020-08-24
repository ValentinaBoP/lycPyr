#!/bin/bash
#SBATCH -J EXONERATE_lycPyr_IL
#SBATCH -o EXONERATE_lycPyr_IL.output
#SBATCH -e EXONERATE_lycPyr_IL.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 4

# the job has been run over all the assembly versions, here the example with the Illumina assembly #

ml bioinfo-tools exonerate

cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Exonerate

QUERY=lycPyrZW_rnd2.all.maker.proteins.fasta
TARGET=lycPyr_IL.fasta

exonerate --model protein2genome --bestn 1 --showvulgar no --showquerygff yes --ryo ">%ti %tab-%tae %ps\n%tas\n" --query $QUERY --querytype protein --target $TARGET --targettype dn
