# MHC (Reto's code)

## Run tblastn on exon2 and exon3 - SLURM job

```bash
module load bioinfo-tools
# load 2.4.0 version to match Reto's version
module load blast/2.4.0+

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/MHC

# remove ~ (invalid character)
sed -i 's/~//g' ex2_PROT_query.fas
sed -i 's/~//g' ex3_PROT_query.fas

REF=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/6.Juicer/Intermediate/lycPyr7.2/lycPyr7.2.fasta
QUERY=ex2_PROT_query.fas
OUT=lycPyr7.2.tBLASTN.ex2.out

# it takes about 5 minutes
tblastn -query $QUERY -db $REF -out $OUT -evalue 10e-3 -outfmt 7

QUERY=ex3_PROT_query.fas
OUT=lycPyr7.2.tBLASTN.ex3.out

# it takes about 5 minutes
tblastn -query $QUERY -db $REF -out $OUT -evalue 10e-3 -outfmt 7

for f in *.tBLASTN.*.out
do
	echo -e "query_id\tsubject_id\t%_identity\talignment_length\tmismatches\tgap_opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbit_score" > $f.summary
	cat $f | grep -v ^# | sort -k2,2 -k9,9n >> $f.summary
done
```

## Merge intervals - Reto's code

```R
### Scripts to merge overlapping intervals
### Does not work in bed tools due to inverse hits in which s.start > s.end

setwd("~/Work/Data&Analyses/Other/Birds_of_Paradise/blast_results")

merge.intervals <- function(datafile){

  d <- read.table(datafile, h=T, stringsAsFactors = F)

  # determine strandedness and turn s.start/s.end for negative ones
  for (i in 1:nrow(d)){
    if(d$s.start[i] < d$s.end[i]) d$strand[i] <- "+"
    if(d$s.start[i] > d$s.end[i]) {
          d$strand[i] <- "-"
          start <- d$s.end[i]
          end   <- d$s.start[i]
          d$s.start[i] <- start
          d$s.end[i] <- end
          }
  }

  ## now need to resort the file because of turned coordinates
  d <- d[with(d, order(d$subject_id, d$s.start)), ]

  ## merge the stuff

  # initialize data frame
  d.merged <- d[1,]

  #initialize counter
  counter <- 1

  # run through the data
  # This is constructed such that first, rows with non-overlapping intervals are written to the
  #   d.merged data frame, and these are then compared to the next rows in the original data, and
  #   start and end are updated if necessary

  for(i in 2:(nrow(d))){
    # get new values
    scaf  <- d$subject_id[i]
    start <- d$s.start[i]
    end   <- d$s.end[i]

    # if the scaffold in row i is different from the ine in row i+1, write the row out.
    if(scaf!=d.merged$subject_id[counter]) {
      d.merged[counter+1,] <- d[i,]
      counter <- counter + 1
    }


    # if the scaffold is the same, then look further
    if(scaf==d.merged$subject_id[counter]){

      # if the interval is not overlapping, write the row out, set the counter+1
      if(start > d.merged$s.end[counter]){
        d.merged[counter+1,] <- d[i,]
        counter <- counter+1
      }

      # if the interval is overlapping, look further.
      # Note that no new start should be needed, as the tables are sorted by scaf and start

      # if the new interval is contained in the old, only add the strand information
      # that is, if the new start is smaller than the old end, and so is the new end, keep all
      # This includes intervals with same start and end
      if(start <= d.merged$s.end[counter] & end <= d.merged$s.end[counter]){
        d.merged$strand[counter] <- paste(d.merged$strand[counter], d$strand[i], sep=",")
      }

      # if the new interval extends the old one, update the end, and check whether a negative strand is contained
      # that is, if the new start is smaller than the old end, but the new end larger than the old end, then update etc
      if(start <= d.merged$s.end[counter] & end > d.merged$s.end[counter] ){
        d.merged$s.end[counter] <- end
        d.merged$strand[counter] <- paste(d.merged$strand[counter], d$strand[i], sep=",")
        # note that in this case we need not updating the counter; do not want to. Might still extend the interval.
      }
    }
  }

  #once run through, check whether strandednesses are equal, write out warning if not
  #also update the alignment length
  for(i in 1:nrow(d.merged)){
    strand <- unique(unlist(strsplit(d.merged$strand[i],split=",")))
    d.merged$strand.states[i] <- length(strand)

    d.merged$alignment_length[i] <- d$s.end[i] - d$s.start[i]
  }

  return(d.merged)

}



#### run through the files

files <- list.files(path=".", pattern=c("summary"))

for(k in 1:length(files)){
  filename <- files[k]
  outname <- paste("merged/",filename, ".merged", sep="")
  outname.short <- paste("merged/reduced/",filename, ".merged.reduced", sep="")
  out <- merge.intervals(filename)
  out.short <- out[,c(2,4,9,10,14)]
  write.table(out, outname, sep="\t", row.names=F, col.names=T, quote=F)
  write.table(out.short, outname.short, sep="\t", row.names=F, col.names=T, quote=F)
}
```

## Make BED files

```R
#### to check on the results, use bedtools (on EVE). For this, first the coordinates need to be turned, and I discard all but the interval info

make.bed <- function(datafile){

  d <- read.table(datafile, h=T, stringsAsFactors = F)

  # determine strandedness and turn s.start/s.end for negative ones
  for (i in 1:nrow(d)){
    if(d$s.start[i] > d$s.end[i]) {
      d$strand[i] <- "-"
      start <- d$s.end[i]
      end   <- d$s.start[i]
      d$s.start[i] <- start
      d$s.end[i] <- end
    }
    if(d$s.start[i] < d$s.end[i]) d$strand[i] <- "+"
  }

  ## now need to resort the file because of turned coordinates
  d <- d[with(d, order(d$subject_id, d$s.start)), ]

  bed <- d[,c(2,9,10)]
  return(bed)
}

files <- list.files(path=".", pattern=c("summary"))

for(k in 1:length(files)){
  filename <- files[k]
  outname <- paste("bed_files/",filename, ".bed", sep="")
  out <- make.bed(filename)
  write.table(out, outname, sep="\t", row.names=F, col.names=F, quote=F)
}

```

## Get fasta sequences from alignments

```bash
module load samtools

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/MHC/merged/reduced

for file in *.reduced; do

	## from the blast hit files produce intervals to be extracted

	# set interval around hits
	interval=0

	# produce intervals
	for k in $( cat $file | sed 1d | tr '\t' '@' )
	do
		scaf=$(echo $k | cut -d'@' -f1)
		start=$(echo { $(echo $k | cut -d'@' -f3) - $interval } | bc)
		end=$(echo { $(echo $k | cut -d'@' -f4) + $interval } | bc)
		if [ "$start" -lt "0" ]
			then start=0
		fi
		ival=$(echo $scaf:$start-$end)
		echo $ival
	done > ival.coord

	# obtain the assembly name from this file name
	assembly=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/6.Juicer/Intermediate/lycPyr7.2/lycPyr7.2.fasta

	# extract the intervals
	outfile=$file.fasta
	for coord in $( cat ival.coord )
	do
		samtools faidx $assembly $coord
		
	done > $outfile
	
done
```

## Results so far

```bash
-rw-rw-r-- 1 vpeona sllstore2017073  4192 Apr 23 11:40 lycPyr7.2.tBLASTN.ex2.out.summary.merged.reduced
-rw-rw-r-- 1 vpeona sllstore2017073 33329 Apr 23 11:53 lycPyr7.2.tBLASTN.ex2.out.summary.merged.reduced.fasta
-rw-rw-r-- 1 vpeona sllstore2017073  1179 Apr 23 11:40 lycPyr7.2.tBLASTN.ex3.out.summary.merged.reduced
-rw-rw-r-- 1 vpeona sllstore2017073  8396 Apr 23 11:55 lycPyr7.2.tBLASTN.ex3.out.summary.merged.reduced.fasta
```

`lycPyr7.2.tBLASTN.ex3.out.summary.merged.reduced.fasta` 39 sequences
`lycPyr7.2.tBLASTN.ex2.out.summary.merged.reduced.fasta` 134 sequences

## BLASTN to GenBank - online

## Filter out

```bash
awk '$2 == 1 {print $0}' lycPyr7.2_ex3_Genbank.tsv | cut -f1 > lycPyr7.2_ex3_Genbank_filtered
awk '$2 == 1 {print $0}' lycPyr7.2_ex2_Genbank.tsv | cut -f1 > lycPyr7.2_ex2_Genbank_filtered

FASTA=lycPyr7.2.tBLASTN.ex3.out.summary.merged.reduced.fasta
LIST=lycPyr7.2_ex3_Genbank_filtered
OUT=lycPyr7.2_ex3_Genbank_filtered.fasta

perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $FASTA list $LIST > $OUT

FASTA=lycPyr7.2.tBLASTN.ex2.out.summary.merged.reduced.fasta
LIST=lycPyr7.2_ex2_Genbank_filtered
OUT=lycPyr7.2_ex2_Genbank_filtered.fasta

perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl $FASTA list $LIST > $OUT
```

## MAFFT

```bash
module load bioinfo-tools
module load MAFFT

SEQ=lycPyr7.2_ex2_Genbank_filtered.fasta
ALIGN=ex2_query.fasta
OUT=lycPyr7.2_ex2_Genbank_filtered_MAFFT.fasta

mafft --thread 3 --threadit 0 --treeout --inputorder --anysymbol --ep 0.0 --add new_sequences --globalpair input > output

SEQ=lycPyr7.2_ex3_Genbank_filtered.fasta
ALIGN=ex3_query.fasta
OUT=lycPyr7.2_ex3_Genbank_filtered_MAFFT.fasta

mafft --inputorder --keeplength --mapout --anysymbol --add $SEQ --globalpair $ALIGN > temp
mv temp $OUT
```
