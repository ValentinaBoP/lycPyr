#!/bin/bash -l
#SBATCH -J LASTZ_lycPyr4.1_PB_v2
#SBATCH -o LASTZ_lycPyr4.1_PB_v2.output
#SBATCH -e LASTZ_lycPyr4.1_PB_v2.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 4-00:00:00
#SBATCH -A snic2017-1-533
#SBATCH -p core
#SBATCH -n 8

#### EXAMPLE OF LASTZ ALIGNMENT - IT WAS DONE FOR ALL THE ASSEMBLY VERSIONS ####

cp /home/vpeona/b2016267/lycPyr_PacBio/analysis/falcon_assembly_01/lycPyr2.fa $SNIC_TMP
cp /home/vpeona/b2016267/lycPyr_TheGenome/4.PhaseGenomics_scaffold/38_Scaffolds_default/lycPyr4.1.fasta $SNIC_TMP
cd $SNIC_TMP

TARGET=lycPyr2.fa
REF=lycPyr4.1.fasta

module load bioinfo-tools lastz/1.04.00

lastz_32 $REF[multiple] $TARGET M=254 K=4500 L=3000 Y=15000 C=2 T=2 --matchcount=10000 --ambiguous=iupac --format=general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2
