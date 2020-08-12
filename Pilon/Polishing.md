# Error correction with Pilon (Illumina + 10XGenomics reads)

Correction of small indels outside repeats

- Map the reads (Illumina and 10XGenomics) to the genome with bwa or bowtie2. BAM files needed as input.
- The BAM files must be sorted in coordinate order and indexed
- I need to divide the scaffolds in chunks, otherwise Pilon will collapse.
- I have 1703 scaffolds. Divide them in 4 chunks of 400 scaffolds each.
- The first 36 scaffolds are really large respect to the others. Probably I should divide the large scaffolds in 4 chunks and have an additional fifth chunk with all the remaining little scaffolds.
- Run RepeatMasker
- Merge the vcf files produced by Pilon and remove all the modifications listed within repeats
- Turn RepeatMasker out file into a bed file
- Use bedtools to remove the variants within repeats
- Select PASS variants with GATK SelectVariants
- Apply the variants with GATK FastaAlternateReferenceMaker


## Map the reads with bwa

```
BWA mem -t <threads> <reference_path> <SE_reads_path> > filename_out
BWA mem -t <threads> <reference_path> <PE_F_reads_path> <PE_R_reads_path> > filename_out

samtools view -@ 16 -Sb P6156_101.align.sam > P6156_101.align.bam
```

### Map Illumina reads

`/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/BWA_lycPyr6.1_IL.sh`

```bash
#!/bin/bash -l
#SBATCH -J BWA_lycPyr6.1_IL
#SBATCH -o BWA_lycPyr6.1_IL.output
#SBATCH -e BWA_lycPyr6.1_IL.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02-00:00:00
#SBATCH -A snic2018-8-266
#SBATCH -p node
#SBATCH -C fat

# load modules
ml bioinfo-tools
ml bwa/0.7.17
ml samtools/1.8

# set variables
REF=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/3.ErrorCorrectionPB/lycPyr6.1.fasta
PE_F=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_Illumina/LpyrF013_L001/LpyrF013_L001_R1.fastq.gz
PE_R=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_Illumina/LpyrF013_L001/LpyrF013_L001_R2.fastq.gz
OUT=lycPyr6.1_LPyrF013_L001_mem
DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/

# copy files to temporary directory
cp $REF $SNIC_TMP
cp $PE_F $SNIC_TMP
cp $PE_R $SNIC_TMP
cd $SNIC_TMP

# re-set variables
REF=`basename $REF`
PE_F=`basename $PE_F`
PE_R=`basename $PE_R`

# run bwa
bwa index $REF
bwa mem -t 20 $REF $PE_F $PE_R > $OUT.sam

# convert sam into bam and index the bam file
samtools view -@ 20 -Sb $OUT.sam > $OUT.bam
samtools index -b -@ 20 $OUT.bam $OUT.bai

# copy the output back to $DIR
cp -t $DIR $OUT.bam $OUT.bai
```

### Map 10XGenomics reads

`/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/BWA_lycPyr6.1_10X.sh`

```bash
#!/bin/bash -l
#SBATCH -J BWA_lycPyr6.1_10X
#SBATCH -o BWA_lycPyr6.1_10X.output
#SBATCH -e BWA_lycPyr6.1_10X.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02-00:00:00
#SBATCH -A snic2018-8-266
#SBATCH -p node
#SBATCH -C fat

# load modules
ml bioinfo-tools
ml bwa/0.7.17
ml samtools/1.8
longranger/2.1.4

# set variables
REF=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/3.ErrorCorrectionPB/lycPyr6.1.fasta
PE_F=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_102/P6156_102_S8_L008_R1_001.fastq.gz
PE_R=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_102/P6156_102_S8_L008_R1_002.fastq.gz
SAMPLE=P6156_102
OUT=lycPyr6.1_P6156_102_mem
DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/

# copy files to temporary directory
cp $REF $SNIC_TMP
cp $PE_F $SNIC_TMP
cp $PE_R $SNIC_TMP
cd $SNIC_TMP

# re-set variables
REF=`basename $REF`
PE_F=`basename $PE_F`
PE_R=`basename $PE_R`

# run bwa
bwa index $REF

### longranger basic
longranger basic --id P6156_102 --fastqs=./

### map reads
unpigz -c ./P6156_102/outs/barcoded.fastq.gz | perl -ne 'chomp;$ct++;$ct=1 if($ct>4);if($ct==1){if(/(\@\S+)\sBX\:Z\:(\S{16})/){$flag=1;$head=$1."_".$2;print "$head\n";}else{$flag=0;}}else{print "$_\n" if($flag);}' > $SAMPLE_interleaved.fastq

bwa mem -t 20 $REF -p $SAMPLE_interleaved.fastq > $OUT.sam

rm $SAMPLE_interleaved.fastq

# convert sam into bam and index the bam file
samtools view -@ 20 -Sb $OUT.sam > $OUT.bam
samtools index -b -@ 20 $OUT.bam $OUT.bai

# copy the output back to $DIR
cp -t $DIR $OUT.bam $OUT.bai
```

## Divide the scaffolds in chunks

- 2 big chunks containing all the little scaffolds: [830 + 835]
- 7 chunks containing the big chromosomes. The structure of the chunks is this:


```
Name	Number of scaffolds		Cumulative sum
chunk1	[2]						2
chunk2	[3]						5
chunk3	[4]						9
chunk4	[5]						14
chunk5	[6]						20
chunk6	[8]						28
chunk7	[10]					38
chunk8	[830]					868
chunk9	[835]					1703
```

Code also saved in `/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/createChunks.R`

```r
n = c(2,3,4,5,6,8,10,830,835)

createChunks = function(n, fasta){
	
	# n contains the number of scaffolds that end up in each chunk
	# cumSum is the cumulative sum of scaffolds already placed in the chunks
	counter = 0
	cumSum = 0
	for(i in 1:length(n)){
		
		cumSum[i] = c(n[i] + counter)
		counter = cumSum[i]
		
	}
	
	begin = c(1, cumSum[1:(length(cumSum)-1)] + 1)
	
	cmd = paste("grep '>' ", fasta, " | cut -c2-", sep = "")
	names = system(command = cmd, intern = TRUE)
	
	for(i in 1:length(begin)){

	chunk = paste(names[begin[i]:cumSum[i]], collapse = ',')
	write.table(x = chunk, file = paste("lycPyr6.1_target_chunk", i, ".txt", sep = ""), col.names = F, row.names = F, quote = F)
	
	}

}
```

## Write the jobs for Pilon

Code also saved in `/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/writePilonJobs.R`

```r
targets = paste("lycPyr6.1_target_chunk", 1:9, ".txt", sep = "")

jobs = paste("PILON_lycPyr6.1_chunk", 1:9, sep = "")

for(i in 1:length(jobs)){
	
	sink(paste(jobs[i], ".sh", sep = ""))
	
	cat("#!/bin/bash -l")
	cat("\n")
	cat(paste("#SBATCH -J ", jobs[i], sep = ''))
	cat("\n")
	cat(paste("#SBATCH -o ", jobs[i], ".output", sep = ''))
	cat("\n")
	cat(paste("#SBATCH -e ", jobs[i], ".error", sep = ''))
	cat("\n")
	cat("#SBATCH --mail-user valentina.peona90@gmail.com")
	cat("\n")
	cat("#SBATCH --mail-type=ALL")
	cat("\n")
	cat("#SBATCH -t 12:00:00")
	cat("\n")
	cat("#SBATCH -A snic2018-8-266")
	cat("\n")
	cat("#SBATCH -p node")
	cat("\n")
	cat("#SBATCH -C mem1TB")
	cat("\n")
	cat("ml bioinfo-tools Pilon/1.22")
	cat("\n")
	cat("REF=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/3.ErrorCorrectionPB/lycPyr6.1.fasta")
	cat("\n")
	cat("ILLUMINA=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/lycPyr6.1_LPyrF013_L001_mem_merged_sorted")
	cat("\n")
	cat("GEN10X=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/lycPyr6.1_P6156_102_mem")
	cat("\n")
	cat(paste("OUT=lycPyr6.1_target_chunk", i, sep = ""))
	cat("\n")
	cat(paste("TARGET=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/lycPyr6.1_target_chunk", i, ".txt", sep = ""))
	cat("\n")
	cat("DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1/")
	cat("\n")
	cat("cp -t $SNIC_TMP $REF $ILLUMINA.bam $ILLUMINA.bai $GEN10X.bam $GEN10X.bai $TARGET")
	cat("\n")
	cat("cd $SNIC_TMP")
	cat("\n")
	cat("REF=`basename $REF`")
	cat("\n")
	cat("ILLUMINA=`basename $ILLUMINA`")
	cat("\n")
	cat("GEN10X=`basename $GEN10X`")
	cat("\n")
	cat("TARGET=`basename $TARGET`")
	cat("\n")
	cat("java -Xmx1000G -jar $PILON_HOME/pilon.jar --genome $REF --bam $ILLUMINA.bam --bam $GEN10X.bam --output $OUT --outdir $DIR --vcf --vcfqe --diploid --threads 8 --fix bases --targets $TARGET")
	cat("\n")
	cat(paste("cp $OUT* $DIR"))
	cat("\n")
	sink()
}
```

### Select only indels with PASS flag

```bash
#!/bin/bash -l
#SBATCH -J VCFTOOLS_INDELS_lycPyr6.1_chunks
#SBATCH -o VCFTOOLS_INDELS_lycPyr6.1_chunks.output
#SBATCH -e VCFTOOLS_INDELS_lycPyr6.1_chunks.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 04:00:00
#SBATCH -A snic2018-3-568
#SBATCH -p core
#SBATCH -n 8

ml bioinfo-tools vcftools

DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1

cp $DIR/*.vcf $SNIC_TMP
cd $SNIC_TMP

ls *.vcf > list_vcf

for VCF in $( cat list_vcf )
do
	OUT="${VCF%.*}"_onlyIndels_PASS
	vcftools --vcf $VCF --remove-filtered-all --keep-only-indels --recode --out $OUT
done

cp $OUT* $DIR
```

### Divide the original fasta according to the chunks
- It's easier to apply the VCF files on them separately


```r
chunks = list.files("chunk")

createChunksFasta = function(chunks, fasta){
	
	if(!require(seqinr)){
	
		message("Install seqinr!")
	
	}
	
	cmd = paste("grep '>' ", fasta, " | cut -c2-", sep = "")
	names = system(command = cmd, intern = TRUE)
	fasta = read.fasta(fasta, seqonly = TRUE)
	
	i = 1
	for(chunk in chunks){
	
		names_chunk = read.table(file = chunk, header = F, stringsAsFactors = F)[,1]
		names_chunk = unlist(strsplit(x = names_chunk, split = ","))
		boo = names %in% names_chunk
		
		new_names = names[boo]
		new_fasta = fasta[boo]
		
		write.fasta(sequences = new_fasta, names = new_names, file.out = sub(pattern = "txt", replacement = "fasta", x = chunk))
		
		i = i + 1
	}

}
```

### Run RepeatMasker
- Run RepeatMasker on the chunks separately


```bash
#!/bin/bash -l
#SBATCH -J RMSK_lycPyr6.1_chunks
#SBATCH -o RMSK_lycPyr6.1_chunks.output
#SBATCH -e RMSK_lycPyr6.1_chunks.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 20:00:00
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 20

# load modules
module load bioinfo-tools RepeatMasker/4.0.7

# set variables
	LIB=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Data/3BOP_rm2nr_lycPyr2_rm2.1_hcrow_rm3.1_ficAlb15_rm3.1_aves.lib
	CHUNKS=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/lycPyr6.1_target_chunk
	DIR=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/RMSK/

# copy files to temporary directory
cp -t $SNIC_TMP $LIB $CHUNKS*.fasta

# go tot temporary directory
cd $SNIC_TMP

# re-set variables
LIB=`basename $LIB`
CHUNKS=`basename $CHUNKS`


for chunk in $( ls *chunk*.fasta )
do
		cat $chunk | perl -ne 'chomp;if(/^\>/){$ct++;print ">$ct\n";}else{print "$_\n";}' > ${chunk}_renamed.fasta
		grep '>' $chunk | cut -c2- > first
		grep '>' ${chunk}_renamed.fasta | cut -c2- > second
		paste second first> ${chunk}_correspondence_table.txt
        RepeatMasker -pa 16 -a -xsmall -gccalc -dir ./ -lib $LIB ${chunk}_renamed.fasta
		awk '{print $5, $6-1, $7-1, $11, $2, $9}' OFS='\t' ${chunk}_renamed.fasta.out | awk 'NR > 3' | awk '{ sub(/C$/, "-", $6) }1' OFS='\t' > ${chunk}_renamed.bed
		sed -f <(printf 's/^%s$/%s/g\n' $(<${chunk}_correspondence_table.txt)) <(cut -f1 ${chunk}_renamed.bed) > ${chunk}_col1
		cut -f2- ${chunk}_renamed.bed | paste ${chunk}_col1 - > ${chunk}.out.bed
		cp -t $DIR ${chunk}_correspondence_table.txt ${chunk}.out.bed
done
```

### Select variants outside repeats

```bash
#!/bin/bash -l
#SBATCH -J BEDTOOLS_INTERSECT_lycPyr6.1_chunks
#SBATCH -o BEDTOOLS_INTERSECT_lycPyr6.1_chunks.output
#SBATCH -e BEDTOOLS_INTERSECT_lycPyr6.1_chunks.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2018-3-568
#SBATCH -p core
#SBATCH -n 8

module load bioinfo-tools BEDTools/2.26.0

DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1/

cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1/*PASS.vcf $SNIC_TMP
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/RMSK/*.fasta.out.bed $SNIC_TMP
cd $SNIC_TMP

ls *.fasta.out.bed | sed 's/.fasta.out.bed//g' > list_chunks

for chunk in $( cat list_chunks)
do
	grep -v "Simple_repeat\|Low_complexity\|Satellite\|Unknown" $chunk.fasta.out.bed > temp
	mv temp $chunk.fasta.out.bed
	bedtools intersect -a ${chunk}_onlyIndels_PASS.vcf -b $chunk.fasta.out.bed -v -wa -header > ${chunk}_intersect_final.vcf
	cp ${chunk}_intersect_final.vcf $DIR
done
```

### Apply VCF with GATK

```bash
#!/bin/bash -l
#SBATCH -J GATK_ALTERNATE_lycPyr6.1_chunks
#SBATCH -o GATK_ALTERNATE_lycPyr6.1_chunks.output
#SBATCH -e GATK_ALTERNATE_lycPyr6.1_chunks.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 06:00:00
#SBATCH -A snic2018-3-568
#SBATCH -p core
#SBATCH -n 8

# load modules
ml bioinfo-tools GATK/3.7

DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1/

cp /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/Round1/*intersect_final.vcf $SNIC_TMP
cp /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/4.ErrorCorrectionIL/Intermediate/Pilon/*.fasta $SNIC_TMP
cd $SNIC_TMP

ls *.fasta | sed 's/.fasta//g' | > list_chunks

for chunk in $( cat list_chunks )
do
	java -jar GenomeAnalysisTK.jar -T FastaAlternateReferenceMaker -R $chunk.fasta -o ${chunk}_Round1.fasta -V ${chunk}_intersect_final.vcf
	cp ${chunk}_Round1.fasta $DIR
done
```

```
Pilon version 1.22 Wed Mar 15 16:38:30 2017 -0400

    Usage: pilon --genome genome.fasta [--frags frags.bam] [--jumps jumps.bam] [--unpaired unpaired.bam]
                 [...other options...]
           pilon --help for option details


         INPUTS:
           --genome genome.fasta
              The input genome we are trying to improve, which must be the reference used
              for the bam alignments.  At least one of --frags or --jumps must also be given.
           --frags frags.bam
              A bam file consisting of fragment paired-end alignments, aligned to the --genome
              argument using bwa or bowtie2.  This argument may be specifed more than once.
           --jumps jumps.bam
              A bam file consisting of jump (mate pair) paired-end alignments, aligned to the
              --genome argument using bwa or bowtie2.  This argument may be specifed more than once.
           --unpaired unpaired.bam
              A bam file consisting of unpaired alignments, aligned to the --genome argument
              using bwa or bowtie2.  This argument may be specifed more than once.
           --bam any.bam
              A bam file of unknown type; Pilon will scan it and attempt to classify it as one
              of the above bam types.
         OUTPUTS:
           --output prefix
              Prefix for output files
           --outdir directory
              Use this directory for all output files.
           --changes
              If specified, a file listing changes in the <output>.fasta will be generated.
           --vcf
              If specified, a vcf file will be generated
           --vcfqe
               If specified, the VCF will contain a QE (quality-weighted evidence) field rather
               than the default QP (quality-weighted percentage of evidence) field.
           --tracks
               This options will cause many track files (*.bed, *.wig) suitable for viewing in
               a genome browser to be written.
         CONTROL:
           --variant
              Sets up heuristics for variant calling, as opposed to assembly improvement;
              equivalent to "--vcf --fix all,breaks".
           --chunksize
              Input FASTA elements larger than this will be processed in smaller pieces not to
              exceed this size (default 10000000).
           --diploid
              Sample is from diploid organism; will eventually affect calling of heterozygous SNPs
           --fix fixlist
              A comma-separated list of categories of issues to try to fix:
                "snps": try to fix individual base errors;
                "indels": try to fix small indels;
                "gaps": try to fill gaps;
                "local": try to detect and fix local misassemblies;
                "all": all of the above (default);
                "bases": shorthand for "snps" and "indels" (for back compatibility);
                "none": none of the above; new fasta file will not be written.
              The following are experimental fix types:
                "amb": fix ambiguous bases in fasta output (to most likely alternative);
                "breaks": allow local reassembly to open new gaps (with "local");
                "circles": try to close circlar elements when used with long corrected reads;
                "novel": assemble novel sequence from unaligned non-jump reads.
           --dumpreads
              Dump reads for local re-assemblies.
           --duplicates
              Use reads marked as duplicates in the input BAMs (ignored by default).
           --iupac
              Output IUPAC ambiguous base codes in the output FASTA file when appropriate.
           --nonpf
              Use reads which failed sequencer quality filtering (ignored by default).
           --targets targetlist
              Only process the specified target(s).  Targets are comma-separated, and each target
              is a fasta element name optionally followed by a base range.
              Example: "scaffold00001,scaffold00002:10000-20000" would result in processing all of
              scaffold00001 and coordinates 10000-20000 of scaffold00002.
              If "targetlist" is the name of a file, each line will be treated as a target
              specification.
           --threads
              Degree of parallelism to use for certain processing (default 1). Experimental.
           --verbose
              More verbose output.
           --debug
              Debugging output (implies verbose).
           --version
              Print version string and exit.
         HEURISTICS:
           --defaultqual qual
              Assumes bases are of this quality if quals are no present in input BAMs (default 15).
           --flank nbases
              Controls how much of the well-aligned reads will be used; this many bases at each
              end of the good reads will be ignored (default 10).
           --gapmargin
              Closed gaps must be within this number of bases of true size to be closed (100000)
           --K
              Kmer size used by internal assembler (default 47).
           --mindepth depth
              Variants (snps and indels) will only be called if there is coverage of good pairs
              at this depth or more; if this value is >= 1, it is an absolute depth, if it is a
              fraction < 1, then minimum depth is computed by multiplying this value by the mean
              coverage for the region, with a minumum value of 5 (default 0.1: min depth to call
              is 10% of mean coverage or 5, whichever is greater).
           --mingap
              Minimum size for unclosed gaps (default 10)
           --minmq
              Minimum alignment mapping quality for a read to count in pileups (default 0)
           --minqual
              Minimum base quality to consider for pileups (default 0)
           --nostrays
              Skip making a pass through the input BAM files to identify stray pairs, that is,
              those pairs in which both reads are aligned but not marked valid because they have
              inconsistent orientation or separation. Identifying stray pairs can help fill gaps
              and assemble larger insertions, especially of repeat content.  However, doing so
              sometimes consumes considerable memory.
```
