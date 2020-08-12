#!/bin/bash -l
#SBATCH -J QUADRON_lycPyr3
#SBATCH -o QUADRON_lycPyr3.output
#SBATCH -e QUADRON_lycPyr3.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 4

ml R_packages/3.5.0

mkdir /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr3
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr3

Rscript --vanilla /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Code/quadronPerChromosome.R /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/lycPyr3.2.fasta
