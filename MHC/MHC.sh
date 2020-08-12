#!/bin/bash -l
#SBATCH -J MHC_lycPyr6.1
#SBATCH -o MHC_lycPyr6.1.output
#SBATCH -e MHC_lycPyr6.1.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 00:30:00
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 4

module load bioinfo-tools
# load 2.4.0 version to match Reto's version
module load blast/2.4.0+
module load BEDTools

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/MHC/Vale

# remove ~ (invalid character)
sed -i 's/~//g' ex2_PROT_query.fas
sed -i 's/~//g' ex3_PROT_query.fas

REF=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/3.ErrorCorrectionPB/lycPyr6.1_simpleNames.fasta
QUERY=ex2_PROT_query.fas
OUT=lycPyr6.1.tBLASTN.ex2.out
makeblastdb -out $REF -dbtype nucl -in $REF -parse_seqids

# it takes about 5 minutes
tblastn -query $QUERY -db $REF -out $OUT -evalue 10e-3 -outfmt '6 qseqid qstart qend sseqid sstart send sstrand pident length'

QUERY=ex3_PROT_query.fas
OUT=lycPyr6.1.tBLASTN.ex3.out

# it takes about 5 minutes
tblastn -query $QUERY -db $REF -out $OUT -evalue 10e-3 -outfmt '6 qseqid qstart qend sseqid sstart send sstrand pident length'

ml bioinfo-tools
ml BEDTools/2.27.1

OUT="${OUT%.*}"
OUT="${OUT%.*}"

for file in $( ls $OUT*.out )
do
        awk '!/^#/' $file | awk '{print $4, $5, $6, $7}' > $file.bed
        awk '{if($2 > $3) {print $1, $3, $2, "a"NR, "100", "-"} else if($2 < $3) {print $1, $2, $3, "a"NR, "100", "+"}}' OFS="\t" $file.bed > temp
        mv temp $file.bed
        bedtools sort -i $file.bed > $file.sorted.bed
        bedtools merge -s -c 6 -o distinct -i $file.sorted.bed > $file.sorted.merged.bed
        awk '{print $1, $2, $3, "a"NR, "100", $4}' OFS="\t" $file.sorted.merged.bed > temp
        mv temp $file.sorted.merged.bed
        bedtools getfasta -fi $REF -bed $file.sorted.merged.bed -fo $file.sorted.merged.fasta -s
done
