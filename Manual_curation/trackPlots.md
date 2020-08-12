---
title: "R Notebook"
output:
  html_document: default
  html_notebook: default
---

Notes:  
what about bopSat!?!?

```{r, message=FALSE, warning=FALSE, echo=FALSE}
library(Biostrings)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
source(file = "annotateGaps.R")
source(file = "readRM.R")
source(file = "getFlankingRangeGaps.R")
source(file = "RM2GRanges.R")
source(file = "getGenomicRanges.R")
```

# lycPyr4.1

##Read input files
```{r}
setwd("~/lycPyr_Analysis/CoordindateLiftOver")
# read RepeatMasker outout file
RM = readRM(filepath = "lycPyr4.1.fasta.RM.out")

# read fasta file
FASTA = readDNAStringSet(filepath = "lycPyr4.1.fasta", format = "fasta")

# length chromosomes
lenChr = read.table(file = "lycPyr4.1.chr.length.txt", stringsAsFactors = F, header = F, col.names = c("chromosome", "len"))

```

##List of global variables
```{r}
# names of the chromosomes
CHR = paste("PGA_scaffold", 0:37, sep = '')

# simple repeat categories
SPLREP = c("Simple_repeat", "NA", "A-rich", "T-rich", "G-rich", "C-rich")
```

##List of parameters
```{r}
# window size to get the overlaps between the repeats and the gaps
flank = 100
```

##Create the ranges for gaps and RM hits
```{r}
GAPS = annotateGaps(FASTA = FASTA)
range_gaps = getFlankingRangeGaps(GAPS = GAPS, flank = flank)
range_rm = RM2GRanges(RM = RM[!grepl(pattern = '\\|arrow', x = RM$Sequence),])
```

##Find the overlapping ranges between the gaps and the RM hits
```{r}
#mtch <- findOverlaps(range_rm, range_gaps)
overRM = subsetByOverlaps(range_rm, range_gaps)
overGAPS = subsetByOverlaps(range_gaps, range_rm)
```

##Statistics/plots of the overlapping regions [using only the chromosome model]
```{r}
overRM = overRM[seqnames(overRM) %in% CHR]
overGAPS = overGAPS[seqnames(overGAPS) %in% CHR]
```

```{r}
overRepeats = as.data.frame(overRM@elementMetadata)
overRepeats$Subfamily = sapply(strsplit(x = overRepeats$mcols.Repeat, split = '\\$'), "[[", 1)
overRepeats$Family = sapply(strsplit(x = overRepeats$mcols.Repeat, split = '\\$'), "[[", 2)
overRepeatsCountFamily = overRepeats %>% group_by(Family) %>% tally()
```

###plot overRepeatsCountFamily
```{r}
ggplot(data = overRepeatsCountFamily, aes(x = Family, y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))
```
without simple repeats
```{r}
boo = overRepeatsCountFamily$Family == "Simple_repeat" | overRepeatsCountFamily$Family == "Low_complexity"
newData = overRepeatsCountFamily[!boo,]
newData[3,2] = newData[3,2] + newData[5,2]
newData[4,2] = newData[4,2] + newData[6,2] + newData[7,2]
newData = newData[-c(5,6,7),]
ggplot(data = newData, aes(x = Family, y = n)) + geom_bar(stat = "identity", color = "red", fill = "red") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 20)) + xlab("Repeat") + ylab("Number of generated gaps")
ggsave(filename = "lycPyr4.1_repeat_gaps.png", device = "png", width = 30, height = 21, units = "cm", dpi = 300, limitsize = FALSE)
```


```{r}
overRepeatsCountSubfamily = overRepeats %>% group_by(Subfamily) %>% tally()
boo = grepl(pattern = '\\(', x = overRepeatsCountSubfamily$Subfamily)
overRepeatsCountSubfamily$Subfamily[boo] = 'Simple_repeat'
overRepeatsCountSubfamily = overRepeatsCountSubfamily %>% group_by(Subfamily) %>% summarize(n = sum(n))
```

```{r}
ggplot(data = overRepeatsCountSubfamily, aes(x = Subfamily, y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))
ggplot(data = overRepeatsCountSubfamily[!overRepeatsCountSubfamily$Subfamily %in% SPLREP,], aes(x = Subfamily, y = n)) + geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90))
```

###How many gaps overlap with TEs? [only chromosomes]
We have `r length(overGAPS[seqnames(overGAPS) %in% CHR])` out of `r length(range_gaps[seqnames(range_gaps) %in% CHR])` caused by repetitive elements

#Plot chromosomes, gaps and RM track
```{r}
chrPos = getGenomicRanges(FASTA, lenChr)
chrPos = as.data.frame(chrPos)
```

###Plot of chromosome Z
```{r}
chrZ = chrPos[chrPos$seqnames == "PGA_scaffold5",]
chrZ$y = 5
chrZ$category = "Genome"

RMZ = as.data.frame(range_rm)
RMZ = RMZ[RMZ$seqnames == "PGA_scaffold5",]
RMZ$y = 3
RMZ = RMZ[,c(1:5,7,6)]
names(RMZ)[7] = "category"
RMZ$category = sapply(strsplit(x = RMZ$category, split = '\\$'), "[[", 2)



MERGED = rbind(chrZ, RMZ)
subMERGED = MERGED[MERGED$end <= 113437,]

ys = length(unique(subMERGED$category)):1
cat = unique(subMERGED$category)
for(i in 1:length(cat)){
  
  subMERGED$y[subMERGED$category == cat[i]] = ys[i]
  
}
```

```{r}
ggplot(data = subMERGED, aes(fill = category)) + geom_rect(data = subMERGED, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5)) + scale_y_continuous(labels = unique(subMERGED$category)[8:1], breaks = 1:8)
```


#lycPyr1

#lycPyr2
```{r}
PB = fread(file = "LASTZ_lycPyr4.1_PB_v2.output", stringsAsFactors = F)
names(PB)[1] = "name1"
names(PB)

#subset only chromosomes
PB = PB[grepl(pattern = "PGA_scaffold", x = PB$name1),]

#sort by coordinates
PB = PB[order(PB[,1], PB[,2]),]
```
Arrange dataframe for plotting
```{r}
PB_Z = PB[PB$name1 == "PGA_scaffold5", c(6,2,3,10)]
PB_Z$width = PB_Z$end1 - PB_Z$start1 + 1
PB_Z$category = "PacBio"
names(PB_Z) = c("seqnames", "start", "end", "strand", "width", "category")
PB_Z$y = "NA"
PB_Z = PB_Z[,c(1,2,3,5,4,7,6)]

MERGED = rbind(chrZ, RMZ, PB_Z)

#set custome colors
MERGED$cols = "black"
MERGED[MERGED$strand == "+", 8] = "red"
MERGED[MERGED$strand == "-", 8] = "blue"

subMERGED = MERGED[MERGED$end <= 895179,]

ys = length(unique(subMERGED$category)):1
cat = unique(subMERGED$category)
for(i in 1:length(cat)){
  
  subMERGED$y[subMERGED$category == cat[i]] = ys[i]
  
}

subMERGED$y = as.integer(subMERGED$y)
subMERGED$cols = as.factor(subMERGED$cols)
```

Add PB track to the plot
```{r}
ggplot(data = subMERGED, aes(fill = category)) + geom_rect(data = subMERGED, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5)) + scale_y_continuous(labels = unique(subMERGED$category)[length(cat):1], breaks = 1:length(cat)) + theme(legend.position="none")

ggplot(data = subMERGED) + geom_rect(data = subMERGED, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = cols)) + scale_y_continuous(labels = unique(subMERGED$category)[length(cat):1], breaks = 1:length(cat)) + theme(legend.position="none") + scale_fill_manual(values = as.character(levels(subMERGED$cols)))
```

#Wrap function
```{r}
library(Biostrings)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(ggplot2)
source(file = "annotateGaps.R")
source(file = "readRM.R")
source(file = "getFlankingRangeGaps.R")
source(file = "RM2GRanges.R")
source(file = "getGenomicRanges.R")
source(file = "partitionTrack.R")
source(file = "subsetByPartition.R")
source(file = "listPlots.R")
source(file = "addLASTZTrack.R")
source(file = "createRangesPlot.R")
source(file = "addRMTrack.R")
source(file = "generatePlots.R")
source(file = "addPiTrack.R")
source(file = "addGapLines.R")
source(file = "addRMDensity.R")
source(file = "addRMTrack.R")
source(file = "estimateRMDensity.R")
source(file = "addSatsumaTrack.R")

# read RepeatMasker outout file
#RM = readRM(filepath = "lycPyr4.1.fasta.out")
RM = list.files(pattern = ".RM.out")

# read fasta file
FASTA = readDNAStringSet(filepath = "lycPyr4.1.fasta", format = "fasta")

# length chromosomes
lenChr = read.table(file = "lycPyr4.1.chr.length.txt", stringsAsFactors = F, header = F, col.names = c("chromosome", "len"))

chrPos = getGenomicRanges(FASTA, lenChr)
chrPos = as.data.frame(chrPos)

chrPos$y = 1
chrPos$category = "Genome"
chrPos$cols = "black"

list_LASTZ = list.files(pattern = ".output")

categories = sub(pattern = ".output", replacement = "", x = sapply(strsplit(x = list_LASTZ, split = '-'), '[[', 2))

MERGED = chrPos

for(i in 1:length(list_LASTZ)){
  MERGED = addLASTZTrack(file = list_LASTZ[i], category = categories[i], pattern = "PGA_scaffold", track = MERGED)
}

filename = "satsuma_summary.chained_PB.out"
MERGED = addSatsumaTrack(file = filename, category = "Satsuma_PB", pattern = "PGA_scaffold", track = MERGED)

MERGED = addRMTrack(data = MERGED, file = RM, pattern = "PGA_scaffold")

#PATTERN = paste("PGA_scaffold", 0:37, sep = '')
PATTERN = "PGA_scaffold0"

```
Generate plots with the LASTZ and RM tracks
```{r}
#PI_file = list.files(pattern = "windowed.pi")[2]
#PI = addPiTrack(filename = PI_file, partition = 10000000, simplify = TRUE, track = MERGED)

generatePlots(data = MERGED, pattern = PATTERN, categories = categories)
```

Generate plots with PI track as well
```{r}
PI_file = list.files(pattern = "windowed.pi")[2]
PI = addPiTrack(filename = PI_file, partition = 10000000, simplify = TRUE, track = MERGED)

# test the plotting
p = "PGA_scaffold5"
categories = as.factor(c("Genome", "RepeatMasker", categories))
data = MERGED
chr2Plot = data[data$seqnames == p,]
chr2Plot = partitionTrack(track = chr2Plot)
ranges = createRangesPlot(track = chr2Plot)
dfORDER = as.data.frame(unique(chr2Plot[,c(6,7)]))
o = match(x = categories, table = dfORDER$category)
v = (nrow(dfORDER)+1)-dfORDER[o,1]
categories = factor(x = categories, levels = categories[v])
chr2Plot$fill = paste(chr2Plot$category, chr2Plot$cols, sep = '-')
chr2Plot$fill = as.factor(chr2Plot$fill)
colori = sapply(strsplit(x = levels(chr2Plot$fill), split = '-'), "[[", 2)
i = 1
temp = chr2Plot[chr2Plot$start >= ranges[i] & chr2Plot$end <= ranges[i+1],]
temp$fill = droplevels(temp$fill)
colori = sapply(strsplit(x = levels(temp$fill), split = '-'), "[[", 2)

#subset PI
subPI = PI[PI$seqnames == p,]
subPI = subPI[subPI$start >= ranges[i] & subPI$end <= ranges[i+1],]

# even if in this dummy plot the labels are not correct, it's not important for the purpose of the test
plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill))) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme(legend.position="none", axis.text = element_text(size = 50)) + scale_fill_manual(values = colori) + geom_point(data = subPI, aes(x = start, y = x, size = 0.1))

ggplot(data = subPI, aes(x = start, y = pi)) + geom_line()

newPI = subPI
newPI$new = newPI$pi * 10000
newPI$x = newPI$x + newPI$new

ggplot(data = newPI, aes(x = start, y = new)) + geom_line()

ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill))) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme(legend.position="none", axis.text = element_text(size = 50)) + scale_fill_manual(values = colori) + geom_line(data = newPI, aes(x = start, y = x)) + theme_bw()

ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill))) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + scale_fill_manual(values = colori) + geom_point(data = newPI, aes(x = start, y = x)) + theme_bw() + theme(legend.position="none", axis.text = element_text(size = 50), panel.grid = element_blank())

```



#Figure 3
###barplot, types of repeats vs gaps vs technologies within (PB, 10X IL) and between scaffold gaps (DV PG 10X scaffolding)
