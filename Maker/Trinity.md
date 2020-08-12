# Lycocorax transcriptome assembly

Data in: `/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13`

Paired end reads: `SH-2274-THY-410-25-11-13_S9_L001_R1_001.fastq.gz` and `SH-2274-THY-410-25-11-13_S9_L001_R2_001.fastq.gz`

############ Step 1: Removing erroneous k-mers from Illumina paired-end reads #################

```bash
#!/bin/bash -l 
#SBATCH -J rcorrect 
#SBATCH -o rcorrect_lycPyr_RNAseq.output
#SBATCH -e rcorrect_lycPyr_RNAseq.error
#SBATCH --mail-user valentina.peona90@gmail.com   # Email to send notifications to
#SBATCH --mail-type=ALL                           # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH -t 01-00:00:00                            # Runtime in DD-HH:MM:SS
#SBATCH -A snic2019-8-79                          # Project name
#SBATCH -p core                                   # Partition to submit to
#SBATCH -n 10                                     # Cores Used for the job

ml bioinfo-tools
ml perl/5.18.4
ml perl_modules/5.18.4
ml jellyfish/2.2.6

dir1=/domus/h1/vpeona/rcorrector
R1=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R1_001.fastq.gz
R2=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R2_001.fastq.gz
OUTDIR=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/

# Removing erroneous k-mers from Illumina paired-end reads (sample SH-2274-THY-410-25-11-13)
perl $dir1/run_rcorrector.pl -t 10 -od $OUTDIR -1 $R1 -2 $R2
```
############### Second: Discard read pairs for which one of the reads is deemed unfixable ####################
```bash
#!/bin/bash -l
#SBATCH -J REMOVE_UNFIX_lycPyr_RNAseq
#SBATCH -o REMOVE_UNFIX_lycPyr_RNAseq.output                                # File to which STDOUT will be written
#SBATCH -e REMOVE_UNFIX_lycPyr_RNAseq.error                                # File to which STDERR will be written
#SBATCH --mail-user valentina.peona90@gmail.com                # Email to send notifications to
#SBATCH --mail-type=ALL                                 # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH -t 01-00:00:00                                  # Runtime in DD-HH:MM:SS
#SBATCH -A snic2019-8-79                               # Project name
#ABATCH -M snowy                                        # snowy partition
#SBATCH -p core                                         # Partition to submit to
#SBATCH -n 10                                           # Cores Used for the job

ml bioinfo-tools
ml python/2.7.9

dir1=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Code
dir2=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/
R1=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/SH-2274-THY-410-25-11-13_S9_L001_R1_001.cor.fq.gz
R2=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/SH-2274-THY-410-25-11-13_S9_L001_R2_001.cor.fq.gz

cd $dir2

# Discard read pairs for which one of the reads is deemed unfixable
python $dir1/FilterUncorrectabledPEfastq.py -1 $R1 -2 $R2 -o fixed 2>&1 > REMOVE_UNFIX_lycPyr_RNAseq.std.output
```
############################### Third 2: Run Trimgalore ################################
```bash
#!/bin/bash -l 
#SBATCH -J TRIMGALORE_lycPyr_RNAseq 
#SBATCH -o TRIMGALORE_lycPyr_RNAseq.output          
#SBATCH -e TRIMGALORE_lycPyr_RNAseq.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL                        
#SBATCH -t 03-00:00:00                      
#SBATCH -A snic2019-8-79                     
#SBATCH -M snowy                              
#SBATCH -p core                                
#SBATCH -n 10                               

# load required modules 
ml bioinfo-tools
ml TrimGalore/0.4.4

dir1=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/
dir2=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Trimgalore/
R1=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/fixed_SH-2274-THY-410-25-11-13_S9_L001_R1_001.cor.fq
R2=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Rcorrect/fixed_SH-2274-THY-410-25-11-13_S9_L001_R2_001.cor.fq

# Trim adapter and low quality bases from fastq files
trim_galore --paired --retain_unpaired --phred33 --output_dir $dir2 --length 20 -q 20 --stringency 1 -e 0.1 $R1 $R2 
```
###################### Stpe 4: Run Trinity ####################################
```bash
#!/bin/bash -l
#SBATCH -J TRINITY_lycPyr_RNAseq 
#SBATCH -o TRINITY_lycPyr_RNAseq.output
#SBATCH -e TRINITY_lycPyr_RNAseq.errorwritten
#SBATCH --mail-user valentina.peona90@gmail.com 
#SBATCH --mail-type=ALL                         
#SBATCH -t 5-00:00:00               
#SBATCH -A snic2019-8-79  
#SBATCH -p node               
#SBATCH -n 20 

ml bioinfo-tools
ml trinity/2.8.2
ml trimmomatic/0.36 samtools/1.10 jellyfish/2.2.6
ml bowtie2/2.3.5.1
ml Salmon/0.9.1

dir1=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Trimgalore/
dir2=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Trinity/
R1=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R1_001.fastq.gz
R2=/proj/sllstore2017073/private/01_raw_data/SH-2274/191218_A00181_0135_AHHKFMDRXX/Sample_SH-2274-THY-410-25-11-13/SH-2274-THY-410-25-11-13_S9_L001_R2_001.fastq.gz

# Run Trinity
Trinity --seqType fq --max_memory 100G --SS_lib_type RF --left $R1 --right $R2 --CPU 20 --trimmomatic --full_cleanup --output /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Trinity/
```
