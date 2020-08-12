# GC content per window in the input assemblies

```bash
ml bioinfo-tools BEDTools samtools

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/GCContent

ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_PacBio/analysis/falcon_assembly_01/lycPyr2.fa lycPyr2.fasta

ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_Illumina/lycPyr.fa lycPyr.fasta

ln -s /home/vpeona/sllstore2017073/private/01_raw_data/A.Suh_16_01_supernova2/assemblies/P6156_102.fasta lycPyr_SN2.fasta

ln -s /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/assemblies/Supernova1/P6156_102_pseudohap.fasta lycPyr_SN1.fasta


for genome in $( ls *.fasta )
do
	#samtools faidx $genome > $genome.fai
	#bedtools makewindows -g $genome.fai -w 1000 -i srcwinnum > $genome.windows
	bedtools nuc -fi $genome -bed $genome.windows > $genome.windows.nuc
done

ml R_packages/3.5.0

R
```

```R
IL = read.table("lycPyr.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = IL$V6, prob = 0.95)
# 0.573

SN1 = read.table("lycPyr_SN1.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = SN1$V6, prob = 0.95)
# 0.565

SN2 = read.table("lycPyr_SN2.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = SN2$V6, prob = 0.95)
# 0.575

PB = read.table("lycPyr2.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = PB$V6, prob = 0.95)
# 0.588

ILPB = read.table("lycPyr_ILPB.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = ILPB$V6, prob = 0.95)

DV = read.table("lycPyr_DV.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = DV$V6, prob = 0.95)

L7 = read.table("lycPyr7.4.fasta.windows.nuc", stringsAsFactors = FALSE)
quantile(x = L7$V6, prob = 0.95)


library(ggplot2)


data = data.frame(genome = c(rep("lycPyrIL", nrow(IL)), rep("lycPyrSN1", nrow(SN1)), rep("lycPyrSN2", nrow(SN2)), rep("lycPyrPB", nrow(PB)), rep("lycPyrILPB", nrow(ILPB)), rep("lycPyr2", nrow(DV)), rep("lycPyr6", nrow(L7))), GC = c(IL$V6, SN1$V6, SN2$V6, PB$V6, ILPB$V6, DV$V6, L7$V6))

ggplot(data, aes(GC, colour=genome, group=genome)) + geom_density()
```
