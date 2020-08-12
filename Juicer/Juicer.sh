#! /bin/bash -l
#SBATCH -A snic2019-3-671
#SBATCH -o JUICER_final_round.output
#SBATCH -e JUICER_final_round.error
#SBATCH -J JUICER_final_round
#SBATCH -p node -n 8
#SBATCH -t 24:00:00
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL

cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/6.Juicer/Intermediate/final_round

FASTA=lycPyr7.4.fasta
NAME=lycPyr7.4

module load bioinfo-tools bwa/0.7.17
bwa index $FASTA

ml python/2.7.9
/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/6.Juicer/Code/generate_site_positions.py Sau3AI $NAME

awk 'BEGIN{OFS="\t"}{print $1, $NF}' ${NAME}_Sau3AI.txt > $NAME.chrom.sizes

#grep 'chr' $NAME.chrom.sizes > temp
head -n 43 $NAME.chrom.sizes > temp
mv temp $NAME.chrom.sizes

set -e
singularity exec /crex/proj/sllstore2017073/private/Remi/Juicer/juicer-cpu.simg juicer.sh -t 8 -s Sau3AI -y ${NAME}_Sau3AI.txt -D /opt/juicer/ -z $FASTA -p $NAME.chrom.sizes
