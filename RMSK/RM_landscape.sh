#!/bin/bash -l
#SBATCH -J RM_landscape_lycPyr_merged_library
#SBATCH -o RM_landscape_lycPyr_merged_library.output
#SBATCH -e RM_landscape_lycPyr_merged_library.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 2

# produce the repeat landscape of Lycocorax assemblies using only the RepBase library

# load the module
module load bioinfo-tools RepeatMasker/4.0.7

# copy files to temporary directory
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/P6156_102.fasta.align $SNIC_TMP

# go to temporary directory
cd $SNIC_TMP

# set variables
DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape/
mkdir $DIR

ls *.align > list_align

for ALIGN in $( cat list_align )
do
  # Convert the .align file into a .tbl file with Kimura 2-parameter distances (excluding CpG sites)
  perl /proj/sllstore2017073/private/scripts/RM_divergence/calcDivergenceFromAlign.pl -noCpG $ALIGN > $ALIGN.k2p.noCpG

  # Add sizes of the individual TE copies
  perl /proj/sllstore2017073/private/scripts/RM_divergence/RM.addsize.div.pl $ALIGN.k2p.noCpG > $ALIGN.k2p.noCpG.size

  # Remove the term "kimura=" from the last column
  sed -i 's/kimura=//g' $ALIGN.k2p.noCpG.size

  # copy the files to $DIR
  cp -r $ALIGN.* $DIR
done
