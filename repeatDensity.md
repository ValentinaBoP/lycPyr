# ChrW repeat content respect to Z and autosome - Figure 3 g

## 1) bedtools window: 10kb windows for W, Z and chr2

```bash
ml bioinfo-tools BEDTools

# need the index file
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data
grep 'chrZ\|chrW\|chr2\b' lycPyr7.4.fasta.fai > /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape/lycPyr7.4_chrzw2.fai

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape/

INDEX=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape/lycPyr7.4_chrzw2.fai

bedtools makewindows -g $INDEX -w 50000 > lycPyr7.4_chrzw2.windows
WINDOWS=lycPyr7.4_chrzw2.windows

awk '{print $5, $6, $7}' OFS='\t' lycPyr7.4.fasta.align.k2p.noCpG.size > lycPyr7.4_chrzw2.rmsk
RMSK=lycPyr7.4_chrzw2.rmsk

bedtools coverage -a $WINDOWS -b $RMSK > lycPyr7.4_chrzw2.coverage


bedtools nuc -fi /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Data/lycPyr7.4.fasta -bed lycPyr7.4_chrzw2.windows > lycPyr7.4_chrzw2.nuc

cut -f 10 lycPyr7.4_chrzw2.nuc | awk 'NR>1 {print $0}' > temp
paste -d '\t' lycPyr7.4_chrzw2.coverage temp > lycPyr7.4_chrzw2.coverage.corrected

awk '{print $0, $6-$8}' OFS='\t' lycPyr7.4_chrzw2.coverage.corrected > temp
awk '{print $0, $5/$9}' OFS='\t' temp > lycPyr7.4_chrzw2.coverage.corrected

# almost ready to be plotted
```

Need the coverage for each type of repeat

```bash
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape/

WINDOWS=lycPyr7.4_chrzw2.windows

grep 'LINE' lycPyr7.4.fasta.align.k2p.noCpG.size | awk '{print $5, $6, $7}' OFS='\t' > lycPyr7.4_chrzw2_LINE.rmsk
RMSK=lycPyr7.4_chrzw2_LINE.rmsk
bedtools coverage -a $WINDOWS -b $RMSK > lycPyr7.4_chrzw2_LINE.coverage

grep 'LTR' lycPyr7.4.fasta.align.k2p.noCpG.size | awk '{print $5, $6, $7}' OFS='\t' > lycPyr7.4_chrzw2_LTR.rmsk
RMSK=lycPyr7.4_chrzw2_LTR.rmsk
bedtools coverage -a $WINDOWS -b $RMSK > lycPyr7.4_chrzw2_LTR.coverage

grep 'Satellite' lycPyr7.4.fasta.align.k2p.noCpG.size | awk '{print $5, $6, $7}' OFS='\t' > lycPyr7.4_chrzw2_Satellite.rmsk
RMSK=lycPyr7.4_chrzw2_Satellite.rmsk
bedtools coverage -a $WINDOWS -b $RMSK > lycPyr7.4_chrzw2_Satellite.coverage

grep -v 'LINE\|LTR\|Satellite' lycPyr7.4.fasta.align.k2p.noCpG.size | awk '{print $5, $6, $7}' OFS='\t' > lycPyr7.4_chrzw2_Other.rmsk
RMSK=lycPyr7.4_chrzw2_Other.rmsk
bedtools coverage -a $WINDOWS -b $RMSK > lycPyr7.4_chrzw2_Other.coverage

```

```bash
ml R_packages/3.5.0
```

```R
library(ggplot2)

data = read.table("lycPyr7.4_chrzw2.coverage.corrected", stringsAsFactors = FALSE, sep = "\t")

data$colors = 'black'
data$colors[data$V1 == 'chrZ'] = 'dark grey'

# modify the positions of the windows fro Z and W
data$start = data$V2
data$end = data$V3

# length chromosomes
len2 = max(data$V3[data$V1 == 'chr2'])
lenz = max(data$V3[data$V1 == 'chrZ'])

# correct positions for Z
data$start[data$V1 == 'chrZ'] = data$start[data$V1 == 'chrZ'] + len2
data$end[data$V1 == 'chrZ'] = data$end[data$V1 == 'chrZ'] + len2

# correct positions for W
data$start[data$V1 == 'chrW'] = data$start[data$V1 == 'chrW'] + lenz + len2
data$end[data$V1 == 'chrW'] = data$end[data$V1 == 'chrW'] + lenz + len2

# plot
ggplot(data = data[data$V1 == "chr2",], aes(x = start, y = V10)) + geom_area(stat = "identity", fill = "black") + theme_bw() + geom_area(data = data[data$V1 == "chrZ",], aes(x = start, y = V10), stat = "identity", fill = "dark grey") + geom_area(data = data[data$V1 == "chrW",], aes(x = start, y = V10), stat = "identity", fill = "black")




```
