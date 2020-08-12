# Lycocorax ab-initio MAKER annotation

## SNAP

SNAP is pretty quick and easy to train. Issuing the following commands will perform the training. It is best to put some thought into what kind of gene models you use from MAKER. In this case, we use models with an AED of 0.25 or better and a length of 50 or more amino acids, which helps get rid of junky models.

```bash
#!/bin/bash -l 
#SBATCH -J SNAP_lycPyrZW 
#SBATCH -o SNAP_lycPyrZW.output          
#SBATCH -e SNAP_lycPyrZW.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL                        
#SBATCH -t 00:15:00                      
#SBATCH -A snic2019-8-79                                 
#SBATCH -p core                                
#SBATCH -n 1                               

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate
mkdir Snap
mkdir Snap/Round1
cd Snap/Round1

# load required modules 
ml bioinfo-tools maker/3.01.2-beta

LOG=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output/lycPyrZW_rnd1_master_datastore_index.log

# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d $LOG
rename 's/genome/lycPyr_rnd1.zff.length50_aed0.25/g' *

# gather some stats and validate
fathom lycPyr_rnd1.zff.length50_aed0.25.ann lycPyr_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom lycPyr_rnd1.zff.length50_aed0.25.ann lycPyr_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1

# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom lycPyr_rnd1.zff.length50_aed0.25.ann lycPyr_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1

# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..

# assembly the HMM
hmm-assembler.pl lycPyr_rnd1.zff.length50_aed0.25 params > lycPyr_rnd1.zff.length50_aed0.25.hmm
```

## Augustus

```bash
#!/bin/bash -l 
#SBATCH -J AUGUSTUS_lycPyrZW 
#SBATCH -o AUGUSTUS_lycPyrZW.output          
#SBATCH -e AUGUSTUS_lycPyrZW.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL                        
#SBATCH -t 48:00:00                      
#SBATCH -A snic2019-8-79                                 
#SBATCH -p core                                
#SBATCH -n 4                              

ml bioinfo-tools BEDTools/2.27.1
ml BUSCO/3.0.2b

GFF=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output/lycPyr_rnd1.all.maker.noseq.gff3
FASTA=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/lycPyr7.4_ZW.fasta
OUT=lycPyr_rnd1.all.maker.transcripts1000.fasta

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus

awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' $GFF | awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | bedtools getfasta -fi $FASTA -bed stdin -fo $OUT

source $BUSCO_SETUP

run_BUSCO.py -i $OUT -o lycPyr_rnd1_augustus -l /sw/bioinfo/BUSCO/v2_lineage_sets/aves_odb9 -m genome -c 4 --long -sp chicken -z --augustus_parameters='--progress=true'
```

### Post Augustus

```bash
cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus/run_lycPyr_rnd1_augustus/augustus_output/retraining_parameters

R
```

```R
files = list.files(pattern = "BUSCO")

new = sub(pattern = "BUSCO_lycPyr_rnd1_augustus_3470765843", replacement = "Lycocorax_pyrrhopterus", x = files)

file.rename(from = files, to = new)

q()
```

```bash
cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus/run_lycPyr_rnd1_augustus/augustus_output/retraining_parameters

sed -i 's/BUSCO_lycPyr_rnd1_augustus_3470765843/Lycocorax_pyrrhopterus/g' Lycocorax_pyrrhopterus_parameters.cfg

sed -i 's/BUSCO_lycPyr_rnd1_augustus_3470765843/Lycocorax_pyrrhopterus/g' Lycocorax_pyrrhopterus_parameters.cfg.orig1

mkdir /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus/augustus_config/species/Lycocorax_pyrrhopterus

cp ./Lycocorax_pyrrhopterus*  /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus/augustus_config/species/Lycocorax_pyrrhopterus/
```

## Run Maker with ab-initio gene prediction

```bash
#!/bin/bash -l 
#SBATCH -J MAKER_round2_lycPyrZW 
#SBATCH -o MAKER_round2_lycPyrZW.output          
#SBATCH -e MAKER_round2_lycPyrZW.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL
#SBATCH -t 05-00:00:00
#SBATCH -A snic2019-8-79                          
#SBATCH -p core                                
#SBATCH -n 4                              

# load required modules 
ml bioinfo-tools maker/3.01.2-beta

OUTDIR=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker

#cd $OUTDIR/lycPyrZW_rnd1.maker.output

# subset maker gff3 file
# transcript alignments
#awk '{ if ($2 == "est2genome") print $0 }' lycPyr_rnd1.all.maker.noseq.gff3 > lycPyr_rnd1.all.maker.est2genome.gff3

# protein alignments
#awk '{ if ($2 == "protein2genome") print $0 }' lycPyr_rnd1.all.maker.noseq.gff3 > lycPyr_rnd1.all.maker.protein2genome.gff3

# repeat alignments
#awk '{ if ($2 ~ "repeat") print $0 }' lycPyr_rnd1.all.maker.noseq.gff3 > lycPyr_rnd1.all.maker.repeats.gff3

cd $OUTDIR

# path to the augustus_config folder where I put the Lycocorax gene models in the previous chunk of code
export AUGUSTUS_CONFIG_PATH=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Augustus/augustus_config

maker -cpus 4 -base lycPyrZW_rnd2 round2_maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```

Control file for the second round - ab-initio

```bash
#-----Genome (these are always required)
genome=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/lycPyr7.4_ZW.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output/lycPyr_rnd1.all.maker.est2genome.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output/lycPyr_rnd1.all.maker.protein2genome.gff3  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output/lycPyr_rnd1.all.maker.repeats.gff3 #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Snap/Round1/lycPyr_rnd1.zff.length50_aed0.25.hmm #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=Lycocorax_pyrrhopterus #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=12 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=300000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=20000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```
