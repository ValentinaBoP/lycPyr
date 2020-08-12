#!/bin/bash -l
#SBATCH -J WINDOW_GAP_all_lycPyr
#SBATCH -o WINDOW_GAP_all_lycPyr.output
#SBATCH -e WINDOW_GAP_all_lycPyr.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00:20:00
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 1

# load modules
ml bioinfo-tools BEDTools/2.27.1

# copy files to analyse into a temporary directory
cp /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Annotation/*.gaps $SNIC_TMP
cp /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/*.out $SNIC_TMP

# go to temporary directory
cd $SNIC_TMP

# set output directory
DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Flanks

# list all the RepeatMasker output file available
ls *.out > list_genomes

# begin for loop
for RMSK in $( cat list_genomes )
do

        # get basename of the RepeatMasker file
        filename=$(basename -- "$RMSK")
        filename="${filename%.*}"
        filename="${filename%.*}"
        # reconstruct gap file
        GAPS=$filename.fasta.gaps

        # convert the RepeatMasker output file into a BED file
        /proj/sllstore2017073/private/scripts/RM2BED.sh -i $RMSK -o $RMSK.bed

        # sort the BED files
        bedtools sort -i $GAPS > ${GAPS}_sorted
        bedtools sort -i $RMSK.bed > $RMSK.bed_sorted

        # find which gaps are flanked by repeats and put the result into output directory
        bedtools window -w 200 -a ${GAPS}_sorted -b $RMSK.bed_sorted -bed > $DIR/${filename}_window100.bed

done
# end for loop
