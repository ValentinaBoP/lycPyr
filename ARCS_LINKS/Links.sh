#!/bin/bash
#SBATCH -J LINKS_test4_par5
#SBATCH -o LINKS_test4_par5.output
#SBATCH -e LINKS_test4_par5.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A b2017197
#SBATCH -p core
#SBATCH -n 8

#cd /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/4.Arcs_test4/

#../makeTSVfile.py lycPyr2.2.1.fasta_original.gv lycPyr2.2.1.fasta.tigpair_checkpoint.tsv /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/lycPyr2.2.1.fasta

cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/lycPyr2.2.1.fasta $SNIC_TMP
cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/empty.fof $SNIC_TMP
cp /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/4.Arcs_test4/*.tsv $SNIC_TMP
cd $SNIC_TMP

#/home/vpeona/links_v1.8.5/LINKS -f /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/lycPyr2.2.1.fasta.masked-renamed.fa -s /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/empty.fof -k 20 -b lycPyr3.LINKS

/home/vpeona/links_v1.8.5/LINKS -f lycPyr2.2.1.fasta -s empty.fof -b lycPyr2.2.1.fasta -a 0.2

cp lycPyr2.2.1.fasta.* /proj/b2017197/nobackup/private/Valentina/Pilon/4.Mapping10X/4.Arcs_test4/5.LINKS_par5/
