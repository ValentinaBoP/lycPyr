
#!/bin/bash -l
#SBATCH -J INTERPROSCAN_round2_lycPyrZW
#SBATCH -o INTERPROSCAN_round2_lycPyrZW.output
#SBATCH -e INTERPROSCAN_round2_lycPyrZW.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02-00:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 4

# load required modules
ml bioinfo-tools InterProScan/5.30-69.0

cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/InterProScan/

interproscan.sh -i lycPyrZW_rnd2.all.maker.proteins.fasta -t p -dp -pa --goterms --iprlookup
