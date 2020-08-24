#!/bin/bash -l
#SBATCH -J DEPTH_MHC_lycPyr7.4
#SBATCH -o DEPTH_MHC_lycPyr7.4.output
#SBATCH -e DEPTH_MHC_lycPyr7.4.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 05:00:00
#SBATCH -A snic2020-15-57
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools samtools

DIR=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Minimap
BAM=lycPyr7.4_minimap2.bam

cp $DIR/$BAM $SNIC_TMP
cd $SNIC_TMP

samtools sort -o lycPyr7.4_minimap2.sorted.bam -O BAM -@20 $BAM
rm $BAM
BAM=lycPyr7.4_minimap2.sorted.bam
cp $BAM $DIR
samtools depth -m 100 -Q 60  $BAM | awk '{cov[$1]+=$3;n[$1]++}END{for(i in n){print i"\t"cov[i]/n[i]}}' > $BAM.cov
cp $BAM.cov $DIR
