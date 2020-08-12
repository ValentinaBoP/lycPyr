# Figure 4 - Gap content panel - new alignment strategy

## 1) Extract gaps with 200 bp flanking regions (only intra-scaffold gaps)

```bash
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Annotation

ml bioinfo-tools BEDTools R_packages/3.5.0
```

```bash
FLANK=200

for GAPS in $( ls lycPyr*.gaps )
do
	FASTA="${GAPS%.*}"
	bedtools flank -i $GAPS -g $FASTA.fai -b $FLANK > /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/$GAPS.bothflanks$FLANK.bed
done
```

```bash
cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/

R
```

```R
library(data.table)
library(dplyr)

files = list.files(pattern = "bothflanks200.bed$")

for(i in 1:length(files)){

	flanks = fread(file = files[i])
	boo = grepl(pattern = "start|end", x = flanks$V4)
	flanks = flanks[!boo,]
	flanks$ID = paste(flanks$V1, "__", flanks$V4, sep = "")
	
	new = flanks %>% group_by(ID) %>% summarise(start = min(V2), end = max(V3))

	new$names = new$ID

	new$ID = sapply(strsplit(x = new$names, split = "__"), "[[", 1)

	filename = paste(files[i], ".onlyintra.bed", sep = "")
	write.table(x = new, file = filename, quote = F, col.names = F, row.names = F, sep = "\t")

}
```

```bash
cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/

for GAPS in $( ls *.onlyintra.bed )
do
	FASTA=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/"${GAPS%.*.*.*.*.*}"
	OUT="${GAPS%.*}".fasta
	echo $FASTA
	bedtools getfasta -fi $FASTA -bed $GAPS -fo $OUT -name
done
```

## 2) Blast alignment of the gaps + flankings
SLURM job

`/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Code/BLAST_all_gaps_vs_final.sh`

```bash
#!/bin/bash -l
#SBATCH -J BLAST_all_gaps_vs_final
#SBATCH -o BLAST_all_gaps_vs_final.output
#SBATCH -e BLAST_all_gaps_vs_final.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 02:00:00
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 8

module load bioinfo-tools blast

cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/*.fasta $SNIC_TMP
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/lycPyr7.4.fasta $SNIC_TMP
cd $SNIC_TMP

REF=lycPyr7.4.fasta
DIR=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks

makeblastdb -in $REF -out $REF -dbtype nucl -parse_seqids

for FLANKS in $( ls *onlyintra.fasta )
do
        # run Blast on the gaps + flanks
        blastn -db $REF -query $FLANKS -outfmt '6 qseqid qstart qend sseqid sstart send length pident sstrand evalue' -out $FLANKS.blast -task megablast -num_threads 8

        # copy the output
        cp $FLANKS.blast $DIR/
done
```

## 3) Filter alignment (closest coordinates)

R script in the Code folder: `filterClosestBlastHits.R`

SLURM jobs: `FILTER_CLOSEST_Illumina, FILTER_CLOSEST_SN1, FILTER_CLOSEST_SN2, FILTER_CLOSEST_DV`

```bash
#!/bin/bash -l
#SBATCH -J FILTER_CLOSEST_Illumina
#SBATCH -o FILTER_CLOSEST_Illumina.output
#SBATCH -e FILTER_CLOSEST_Illumina.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00
#SBATCH -A snic2018-8-266
#SBATCH -p core
#SBATCH -n 2

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/lycPyr

ml R_packages/3.5.0

Rscript --vanilla /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Code/filterClosestBlastHits.R lycPyr.fasta.gaps.bothflanks200.bed.onlyintra.fasta.blast
```
```R
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# usage:
# Rscript --vanilla filterClosestBlastHits.R filename

# be already in the working directory
# /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks

library(data.table)
library(dplyr)

#files = list.files(pattern = "onlyintra.fasta.blast")
#files = files[c(1,2,3,4)]

#files = files[3]

files = args[1]

for(j in 1:length(files)){

	data = fread(files[j])
	data = data[data$V8 >= 95,]
	data _ data[data$10 < 1e-20,]
	filename_closest = sub(pattern = ".fasta.blast", replacement = ".closest", x = files[j])
	filename_bed = paste(filename_closest, ".bed", sep = "")
	
	gaps = unique(data$V1)

	for(i in 1:length(gaps)){
	#for(i in 1:10){
		temp = data[data$V1 == gaps[i],]
		minimum = min(temp$V2)
		gap5 = filter(temp, between(temp$V2, left = minimum, right = (minimum + 200)))
		temp5 = gap5[, c(4,5,6,1)]
		temp5 = temp5[order(temp5[,1], temp5[,2]),]
		maximum = max(temp$V3)
		gap3 = filter(temp, between(temp$V2, left = (maximum - 200), right = maximum))
		temp3 = gap3[, c(4,5,6,1)]
		temp3 = temp3[order(temp3[,1], temp3[,2]),]

		write.table(x = temp5, file = "temp", quote = F, sep = "\t", col.names = F, row.names = F)
		system("ml bioinfo-tools BEDTools && awk '{if($3 < $2) {print $1, $3, $2, $4} else {print $0}}' OFS='\t' temp | bedtools sort -i stdin > temp5.bed")

		write.table(x = temp3, file = "temp", quote = F, sep = "\t", col.names = F, row.names = F)
		system("ml bioinfo-tools BEDTools && awk '{if($3 < $2) {print $1, $3, $2, $4} else {print $0}}' OFS='\t' temp | bedtools sort -i stdin > temp3.bed")

		cmd = paste("ml bioinfo-tools BEDTools && bedtools closest -a temp5.bed -b temp3.bed -t first -d > ", filename_closest, sep = "")
		
		system(cmd)
		
		closest = fread(filename_closest)

		boo = closest$V9 == "-1"
		closest = closest[!boo,]
	
		if(nrow(closest) == 0){
		
			next
			
		} else {
		
			closestDef = closest %>% group_by(V1, V4) %>% filter(V9 == min(V9)) %>% filter(row_number()==1)
			new = closestDef %>% summarise(chr = V1, start = sort(c(V2, V3, V6, V7))[2], end = sort(c(V2, V3, V6, V7))[3], name = V4, size = V9)

			if(nrow(new) > 1){
	
				minimum = min(new$size)
				new = new[new$size == minimum,]
		
			}
	
			write.table(x = new, file = filename_bed, append = TRUE, sep = "\t", col.names = F, row.names = F, quote = F)
		
		}
	
	}

}


```

### 4) Intersect gap coordinates with RepeatMasker output file (BED)

```bash
# collect projected coordinates
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/lycPyr*.onlyintra.closest.bed $SNIC_TMP

# collect RM BED files
cp /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/lycPyr7.4.fasta.out.bed $SNIC_TMP

# go to SNIC_TMP
cd $SNIC_TMP

ml bioinfo-tools BEDTools/2.27.1

RMSK=lycPyr7.4.fasta.out.bed
bedtools sort -i $RMSK > temp.rmsk
mv temp.rmsk $RMSK

for GAPS in $( ls *.closest.bed )
do
	OUT="${GAPS%.*.*}".intersect.bed
	awk '{print $3, $4, $5, $6, $7}' OFS='\t' $GAPS > temp.gaps
	mv temp.gaps $GAPS
	bedtools sort -i $GAPS > temp.gaps
	mv temp.gaps $GAPS
	bedtools intersect -a $RMSK -b $GAPS -wo > $OUT
done
```

### 5) Process the intersect output
```bash
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Gaps/Alignment/Flanks/intersect

ml R_packages/3.5.0

R
```

```R
simplifyRepeats = function(data){
  
  data$Gap[grepl(pattern = "LTR", x = data$Gap)] = "LTR"
  data$Gap[grepl(pattern = "DNA", x = data$Gap)] = "DNA"
  data$Gap[grepl(pattern = "Low_complexity", x = data$Gap)] = "Low complexity"
  data$Gap[grepl(pattern = "SINE", x = data$Gap)] = "SINE"
  data$Gap[grepl(pattern = "Simple_repeat", x = data$Gap)] = "Simple repeat"
  data$Gap[grepl(pattern = "Satellite", x = data$Gap)] = "Satellite"
  data$Gap[grepl(pattern = "LINE", x = data$Gap)] = "LINE"
  data$Gap[grepl(pattern = "Helitron", x = data$Gap)] = "Helitron"
  data$Gap[grepl(pattern = "_tRNA", x = data$Gap)] = "tRNA"
  data$Gap[grepl(pattern = "Unspecified", x = data$Gap)] = "Unknown"
  
  return(data)
}
```

```R
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

files = list.files(pattern = "intersect.bed")

totalSummary = data.frame()

for(file in files){

	data = fread(file)

	boo = data$V11 == 0
	data = data[!boo,]

	data$Gap = sapply(strsplit(x = data$V4, split = "__"), "[", 2)

	data = simplifyRepeats(data = data)

	summary = data %>% group_by(V10, V1) %>% summarise(Family = paste(unique(Gap), collapse = ","))

	#table(data$Gap)

	boo = grepl(pattern = ",", x = summary$Family)
	summary$Gap = "gap"
	summary$Gap[boo] = "Complex"
	summary$Gap[!boo] = summary$Family[!boo]

	temp = as.data.frame(table(summary$Gap))
	assembly = sub(pattern = ".fasta.*", replacement = "", x = file)
	temp$Assembly = assembly
	totalSummary = rbind(totalSummary, temp)

}

totalSummary$Var1 = as.character(totalSummary$Var1)

totalSummary = totalSummary %>% group_by(Var1, Assembly) %>% summarise(Count = sum(Freq))

DISCARD = c("tRNA", "SINE")
totalSummary = totalSummary[!totalSummary$Var1 %in% DISCARD,]
totalSummary$Var1 = as.character(totalSummary$Var1)
```

### 5a) Make the plots: number of gaps containing repeats - Figure 4e

```R
levels = unique(totalSummary$Var1)[c(1,3,5,2,6,7,8,4,9)]
totalSummary$Var1 = factor(x = totalSummary$Var1, levels = levels)
```

```R
levels = unique(totalSummary$Assembly)[c(1,3,4,5,2)]
totalSummary$Assembly = factor(x = totalSummary$Assembly, levels = levels)
```

```R
plot = ggplot(data = totalSummary) + geom_bar(aes(x = Var1, y = Count, fill = Assembly), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")
```

```R
#ggsave(filename = "/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure4e_oct.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)

ggsave(filename = "/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure4e_nov.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)
```

### 5b) Make the plots: percentage of gaps containing repeats - Figure 4f

```R
#gapIL = 14573
# not counting gaps smaller than 10
gapIL = 7869
gapILPB = 3718
#gapPB = 3422 * 2
gapSN1 = 21550
gapSN2 = 20131
gapDV = 6736

totalSummary$NGaps[totalSummary$Assembly == "lycPyr"] = gapIL
totalSummary$NGaps[totalSummary$Assembly == "lycPyr_ILPB"] = gapILPB
totalSummary$NGaps[totalSummary$Assembly == "lycPyr_SN1"] = gapSN1
totalSummary$NGaps[totalSummary$Assembly == "lycPyr_SN2"] = gapSN2
totalSummary$NGaps[totalSummary$Assembly == "lycPyr2.1"] = gapDV

totalSummary$Freq = (totalSummary$Count/totalSummary$NGaps) * 100
```

```R
plot = ggplot(data = totalSummary) + geom_bar(aes(x = Var1, y = Freq, fill = Assembly), stat = "identity", position = position_dodge2(width = 0.5)) + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_fill_brewer(palette = "Set1")
```

```R
#ggsave(filename = "/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure4f_oct.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)

ggsave(filename = "/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure4f_nov.pdf", plot = plot, device = "pdf", scale = 1, width = 15, height = 8, units = "cm", dpi = 300, limitsize = FALSE)
```
