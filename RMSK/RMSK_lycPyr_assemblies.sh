#!/bin/bash -l
#SBATCH -J RMSK_lycPyr_assemblies
#SBATCH -o RMSK_lycPyr_assemblies.output
#SBATCH -e RMSK_lycPyr_assemblies.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 48:00:00
#SBATCH -A snic2018-3-568
#SBATCH -p core
#SBATCH -n 16

# load modules
module load bioinfo-tools RepeatMasker/4.0.7

# copy files to temporary directory
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/lycPyr2_rm2.1_merged.lib $SNIC_TMP
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/lycPyr*.fasta $SNIC_TMP

# go to temporary directory
cd $SNIC_TMP

# set variables
LIB=lycPyr2_rm2.1_merged.lib
DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK

# list of the genomes to mask
ls *.fasta > list_genomes
ls *.fa >> list_genomes

#awk 'NR <= 7 {print $0}' list_genomes > temp
#mv temp list_genomes

# begin for loop
for genome in $( cat list_genomes )
do
                # shorten scaffold names at the first occurrence of a pipe
                #perl /proj/sllstore2017073/private/scripts/shortenScaffoldnamesPipe.pl $genome temp
                perl /proj/sllstore2017073/private/scripts/shortenScaffoldnames.pl $genome temp
                perl /proj/sllstore2017073/private/scripts/shortenScaffoldnamesUnderscore.pl temp $genome

                # run RepeatMasker with the custom library
                RepeatMasker -pa 16 -a -xsmall -gccalc -dir ./ -lib $LIB $genome

                # check if the output folder already exist and create it
                if [ ! -d "$DIR/lycPyr2_rm2.1_final" ]; then
                        mkdir $DIR/lycPyr2_rm2.1_final
                fi

                # copy the output files to the output directory
                cp $genome.* $DIR/lycPyr2_rm2.1_final/

                # remove the RepeatMasker output
                rm $genome.*

                # run RepeatMasker with RepBase Aves library
                #RepeatMasker -pa 16 -a -xsmall -gccalc -dir ./ -species aves $genome

                # create the output directory if doesn't exist
                #if [ ! -d "$DIR/Aves" ]; then
                #       mkdir $DIR/Aves
                #fi

                # copy output files
                #cp $genome.* $DIR/Aves/
done
# end for loop
