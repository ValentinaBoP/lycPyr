#!/bin/bash
#SBATCH -J arcs_test4
#SBATCH -o arcs_test4.output
#SBATCH -e arcs_test4.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 08:00:00
#SBATCH -A b2017197
#SBATCH -p core
#SBATCH -n 16

module load bioinfo-tools samtools/1.4

### copy files to temporary directory
cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/lycPyr2.2.1.fasta.masked-renamed.fa $SNIC_TMP
cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/P6156_102.align.sorted.bam $SNIC_TMP
cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/alignment.fof $SNIC_TMP
cd $SNIC_TMP

/home/vpeona/arcs/Arcs/arcs -f lycPyr2.2.1.fasta.masked-renamed.fa -a alignment.fof -b lycPyr3_test4 -e 10000 -m 20-100000 -s 95

cp lycPyr3_test4* /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/4.Arcs_test4/
