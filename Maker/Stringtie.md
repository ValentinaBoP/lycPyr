# Assembly with StringTie

Follow ScilifeLab protocol

# Assembling transcripts based on RNA-seq data

## TRIM the reads with Trimmomatic

```bash
#!/bin/bash -l
#SBATCH -J TRIMMOMATIC_lycPyr_RNAseq
#SBATCH -o TRIMMOMATIC_lycPyr_RNAseq.output
#SBATCH -e TRIMMOMATIC_lycPyr_RNAseq.errorwritten
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p node
#SBATCH -n 8

ml bioinfo-tools trimmomatic/0.36

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Assembly/Trimmomatic/

R1=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R1_001.fastq.gz
R2=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R2_001.fastq.gz

R1PAIRED=lycPyr_R1_paired.fastq.gz
R2PAIRED=lycPyr_R2_paired.fastq.gz
R1UNPAIRED=lycPyr_R1_unpaired.fastq.gz
R2UNPAIRED=lycPyr_R2_unpaired.fastq.gz

java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE -threads 10 $R1 $R2 $R1PAIRED $R1UNPAIRED $R2PAIRED $R2UNPAIRED LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

## READ MAPPING splice aware

```bash
#!/bin/bash -l
#SBATCH -J TOPHAT_lycPyr_RNAseq
#SBATCH -o TOPHAT_lycPyr_RNAseq.output
#SBATCH -e TOPHAT_lycPyr_RNAseq.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 2-00:00:00
#SBATCH -A snic2019-8-79
#SBATCH -p node
#SBATCH -n 8

ml bioinfo-tools tophat/2.1.1


INDEX=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/lycPyr7.4
R1PAIRED=lycPyr_R1_paired.fastq.gz
R2PAIRED=lycPyr_R2_paired.fastq.gz

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data

bowtie2-build lycPyr7.4.fasta lycPyr7.4

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Assembly/Tophat/

tophat --library-type=fr-firststrand $INDEX ../Trimmomatic/$R1PAIRED ../Trimmomatic/$R2PAIRED -p 8
```

```bash
#!/bin/bash -l
#SBATCH -J STRINGTIE_lycPyr_RNAseq
#SBATCH -o STRINGTIE_lycPyr_RNAseq.output
#SBATCH -e STRINGTIE_lycPyr_RNAseq.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 01:00:00
#SBATCH -A snic2020-15-44
#SBATCH -p node
#SBATCH -n 2

ml bioinfo-tools StringTie/1.3.3

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Assembly/StringTie/

stringtie ../Tophat/tophat_out/accepted_hits.bam -o outdir/transcripts.gtf
```

for webapollo
```bash
/domus/h1/vpeona/GAAS/annotation/Tools/Converter/gxf_to_gff3.pl --gff transcripts.gtf -o transcripts.gff3
```
