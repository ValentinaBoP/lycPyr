# Lycocorax evidence based MAKER

## File preparation

## Repeat annotation
```bash
#!/bin/bash -l 
#SBATCH -J REPEATMASKER_lycPyrZW 
#SBATCH -o REPEATMASKER_lycPyrZW.output          
#SBATCH -e REPEATMASKER_lycPyrZW.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL                        
#SBATCH -t 02:00:00                      
#SBATCH -A snic2019-8-79                                 
#SBATCH -p core                                
#SBATCH -n 4                               

ml bioinfo-tools RepeatMasker/4.0.8

FASTA=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/lycPyr7.4_ZW.fasta
LIB=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMDL/lycPyr2_rm2.0/lycPyr2_rm2.1_merged.lib
OUTDIR=/proj/sllstore2017073/private/Valentina/2020ChrW/Data/

cd $SNIC_TMP

RepeatMasker -pa 4 -a -xsmall -gccalc -dir ./ -lib $LIB $FASTA

cp ./lycPyr7.4_ZW.fasta.out $OUTDIR

cd $OUTDIR

RMOUT=lycPyr7.4_ZW.fasta.out
RMGFF=lycPyr7.4_ZW.fasta.out.gff3
RMGFF_FILTERED=lycPyr7.4_ZW_filtered.out.gff3
RMGFF_MAKER=lycPyr7.4_ZW_filtered_reformat_maker.out.gff3

# RM out into GFF3 file
/sw/apps/bioinfo/RepeatMasker/4.0.8/rackham/util/rmOutToGFF3.pl $RMOUT > $RMGFF

# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" $RMGFF > $RMGFF_FILTERED

# reformat to work with MAKER
cat $RMGFF_FILTERED | perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' > $RMGFF_MAKER
```

## Run Maker - evidence based round

### Prepare control file for Maker

```bash
ml bioinfo-tools maker

cd lycPyr_annotation/Intermediate/Maker/

maker -CTL

# modify the control file according to the code below
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
est=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Trinity/Trinity.Trinity.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Assembly/StringTie/outdir/transcripts_ZW.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/protein_chicken_uniprot.fasta  #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org=simple #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein=/sw/bioinfo/maker/3.01.2-beta-mpi/rackham/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Data/lycPyr7.4_ZW_filtered_reformat_maker.out.gff3 #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=8 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=0 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```

Be careful of the parameters:
`genome, est, est_gff, est2genome, protein, model_org, repeat_protein, rm_gff, softmask, protein2genome, min_contig, cpus`

### Run Maker round1

```bash
#!/bin/bash -l 
#SBATCH -J MAKER_round1_lycPyrZW 
#SBATCH -o MAKER_round1_lycPyrZW.output          
#SBATCH -e MAKER_round1_lycPyrZW.error          
#SBATCH --mail-user valentina.peona90@gmail.com  
#SBATCH --mail-type=ALL
#SBATCH -t 05-00:00:00
#SBATCH -A snic2019-8-79                          
#SBATCH -p core                                
#SBATCH -n 12                               

# load required modules 
ml bioinfo-tools maker/3.01.2-beta

#export LD_PRELOAD=$MPI_ROOT/lib/libmpi.so

OUTDIR=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker

cd $OUTDIR

maker -base lycPyrZW_rnd1 round1_maker_opts.ctl maker_bopts.ctl maker_exe.ctl
```

```bash
(base) [vpeona@rackham3 Code]$ sbatch -M snowy MAKER_round1_lycPyrZW.sh
Submitted batch job 1107521 on cluster snowy
```

### Post round1
```bash

# load required modules 
ml bioinfo-tools maker/3.01.2-beta

cd /crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_annotation/Intermediate/Maker/lycPyrZW_rnd1.maker.output

gff3_merge -s -d lycPyrZW_rnd1_master_datastore_index.log > lycPyr_rnd1.all.maker.gff3



fasta_merge -d lycPyrZW_rnd1_master_datastore_index.log
```
