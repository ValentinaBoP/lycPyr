# calculate coverage ratios between male and female individuals

# 1) get the BAM files for the male and female individuals
# unfortunately I don't have those BAM files anywhere so I need to generate them again through LongRanger align

# run Longranger on reference female and on the male

#!/bin/bash -l
#SBATCH -J LONGRANGER_ALIGN_lycPyr104
#SBATCH -o LONGRANGER_ALIGN_lycPyr104.output
#SBATCH -e LONGRANGER_ALIGN_lycPyr104.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 72:00:00                           
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools longranger/2.2.2

REF=/proj/sllstore2017073/private/Valentina/2020ChrW/Data/refdata-lycPyr7.4_tot
READS=/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_104/
SAMPLE=P6156_104

cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/

longranger align --id=lycPyr104 --sample=$SAMPLE --reference=$REF --fastqs=$READS 



#!/bin/bash -l
#SBATCH -J LONGRANGER_ALIGN_lycPyr102
#SBATCH -o LONGRANGER_ALIGN_lycPyr102.output
#SBATCH -e LONGRANGER_ALIGN_lycPyr102.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 72:00:00                           
#SBATCH -A snic2019-8-79
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools longranger/2.2.2

REF=/proj/sllstore2017073/private/Valentina/2020ChrW/Data/refdata-lycPyr7.4_tot
READS=/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_10XG/H5JJHALXX/outs/fastq_path/A_Suh_16_01/Sample_P6156_102
SAMPLE=P10652_201

cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/

longranger align --id=lycPyr102 --sample=$SAMPLE --reference=$REF --fastqs=$READS

# 2) calculate coverage ratio for the W chromosome following Luohao's code

#!/bin/bash -l
#SBATCH -J COVERAGE_SEX_RATIO_female
#SBATCH -o COVERAGE_SEX_RATIO_female.output
#SBATCH -e COVERAGE_SEX_RATIO_female.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 12:00:00                           
#SBATCH -A snic2020-15-57
#SBATCH -p core
#SBATCH -n 20

ml bioinfo-tools samtools

DIR=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage
BAM=possorted_bam.bam

# female coverage

cd $DIR/lycPyr102/outs
cp $BAM $SNIC_TMP
cd $SNIC_TMP

samtools sort $BAM -O BAM -o $BAM.sorted -T aln -@ 20
samtools depth -m 100 -Q 60  $BAM.sorted | awk '{cov[$1]+=$3;n[$1]++}END{for(i in n){print i"\t"cov[i]/n[i]}}' > $BAM.sorted.cov

cp $BAM.sorted $DIR
cp $BAM.sorted.cov $DIR

#!/bin/bash -l
#SBATCH -J COVERAGE_SEX_RATIO_male
#SBATCH -o COVERAGE_SEX_RATIO_male.output
#SBATCH -e COVERAGE_SEX_RATIO_male.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 12:00:00                           
#SBATCH -A snic2020-15-57
#SBATCH -p core
#SBATCH -n 4

ml bioinfo-tools samtools

DIR=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage
BAM=possorted_bam.bam

# female coverage

cd $DIR/lycPyr104/outs
cp $BAM $SNIC_TMP
cd $SNIC_TMP

samtools sort $BAM -O BAM -o $BAM.sorted -T aln -@ 20
samtools depth -m 100 -Q 60 -a $BAM.sorted | awk '{cov[$1]+=$3;n[$1]++}END{for(i in n){print i"\t"cov[i]/n[i]}}' > $BAM.sorted.cov

cp $BAM.sorted $DIR
cp $BAM.sorted.cov $DIR

#!/bin/bash -l
#SBATCH -J COVERAGE_SEX_CHR
#SBATCH -o COVERAGE_SEX_CHR.output
#SBATCH -e COVERAGE_SEX_CHR.error
#SBATCH --mail-user valentina.peona90@gmail.com
#SBATCH --mail-type=ALL                           
#SBATCH -t 2:00:00                           
#SBATCH -A snic2020-15-57
#SBATCH -p core
#SBATCH -n 1

ml bioinfo-tools samtools

cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage

BAMF=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/female
BAMM=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/male

CONTIGSW=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Validation/References/lycPyr6_chrW.fasta.contigs

CONTIGSZ=/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Validation/References/lycPyr6_chrZ.fasta.contigs

#samtools index -b $BAMM.bam $BAMM.bai

#samtools index -b $BAMF.bam $BAMF.bai

awk '{print $1":"$2"-"$3}' $CONTIGSW > chrW.regions
awk '{print $1":"$2"-"$3}' $CONTIGSZ > chrZ.regions


for line in $( cat chrW.regions )
do
samtools coverage -r $line -Q 60 $BAMF.bam | grep -v "#" >> female.coverage.chrW.out
samtools coverage -r $line -Q 60 $BAMM.bam | grep -v "#" >> male.coverage.chrW.out
done

for line in $( cat chrZ.regions )
do
samtools coverage -r $line -Q 60 $BAMF.bam | grep -v "#" >> female.coverage.chrZ.out
samtools coverage -r $line -Q 60 $BAMM.bam | grep -v "#" >> male.coverage.chrZ.out
done

cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Minimap

ln -s /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/chrW.regions .
ln -s /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage/chrZ.regions .

for line in $( cat chrW.regions ); do samtools coverage -r $line -Q 30 lycPyr7.4_minimap2.sorted.bam | grep -v "#" >> lycPyr7.4_chrW.PB.cov; done

for line in $( cat chrZ.regions ); do samtools coverage -r $line -Q 30 lycPyr7.4_minimap2.sorted.bam | grep -v "#" >> lycPyr7.4_chrZ.PB.cov; done


cd /proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage

ml R_packages/3.6.0
R

####R#####

# PLOT COVERAGE MALE AND FEMALE on chrW and chrZ


library(ggplot2)

# get length of the chromosomes for the plot
index = read.table("lycPyr6.fasta.fai", stringsAsFactors = F, sep = "\t")
index = index[index$V1 %in% c("chrW", "chrZ"),]

# coordinates of the contigs for the plot - maybe not necessary - info also in the coverage tables
contigsw = read.table("/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Validation/References/lycPyr6_chrW.fasta.contigs", stringsAsFactors = F, sep = "\t")
contigsz = read.table("/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Validation/References/lycPyr6_chrZ.fasta.contigs", stringsAsFactors = F, sep = "\t")

# the 3rd last column is the mean read depth
coveragewf = read.table("female.coverage.chrW.out", stringsAsFactors = F, sep = "\t")
coveragewm = read.table("male.coverage.chrW.out", stringsAsFactors = F, sep = "\t")
coveragezf = read.table("female.coverage.chrZ.out", stringsAsFactors = F, sep = "\t")
coveragezm = read.table("male.coverage.chrZ.out", stringsAsFactors = F, sep = "\t")

# coverage PB
coveragewpb = read.table("../Minimap/lycPyr7.4_chrW.PB.cov", stringsAsFactors = F, sep = "\t")
coveragezpb = read.table("../Minimap/lycPyr7.4_chrZ.PB.cov", stringsAsFactors = F, sep = "\t")

# for plotting reasons we need one middle point for each contig
coveragewf$xcoord = (coveragewf$V2 + coveragewf$V3) / 2
coveragewm$xcoord = (coveragewm$V2 + coveragewm$V3) / 2
coveragezf$xcoord = (coveragezf$V2 + coveragezf$V3) / 2
coveragezm$xcoord = (coveragezm$V2 + coveragezm$V3) / 2
coveragewpb$xcoord = (coveragewpb$V2 + coveragewpb$V3) / 2
coveragezpb$xcoord = (coveragezpb$V2 + coveragezpb$V3) / 2

# pre-plot
boo = coveragewf$V7 > 200
coveragewf$V7[boo] = -2
boo = coveragewm$V7 > 200
coveragewm$V7[boo] = -2

# pb contigs
pbw = read.table("./PacBioContigs/chrW_intersect.bed", stringsAsFactors = F, sep = "\t")
pbz = read.table("./PacBioContigs/chrZ_intersect.bed", stringsAsFactors = F, sep = "\t")

coveragewf$pb = "pink"
coveragewf$pb[coveragewf$V2 %in% pbw$V2] = "black"
coveragewm$pb = "blue"
coveragewm$pb[coveragewm$V2 %in% pbw$V2] = "white"

coveragezf$pb = "pink"
coveragezf$pb[coveragezf$V2 %in% pbz$V2] = "black"
coveragezm$pb = "blue"
coveragezm$pb[coveragezm$V2 %in% pbz$V2] = "white"


#plot
chrw_plot = ggplot() + geom_rect(data = index[index$V1 == "chrW",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragewf, aes(x = xcoord, y = V7), color = "pink") + geom_point(data = coveragewm, aes(x = xcoord, y = V7), color = "blue", alpha = 0.5) + theme_bw() + xlab("Chromosome W") + ylab("Coverage")

#chrw_plot = ggplot() + geom_rect(data = index[index$V1 == "chrW",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragewf, aes(x = xcoord, y = V7, color = pb, stroke = .1), fill = "pink", alpha = 0.5) + geom_point(data = coveragewm, aes(x = xcoord, y = V7, color = pb), fill = "blue", alpha = 0.5) + theme_bw() + xlab("Chromosome W") + ylab("Coverage") + theme(legend.position = "none")

ggsave(filename = "SupplementaryFigureCoverageWFM.pdf", plot = chrw_plot, device = "pdf", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)
ggsave(filename = "SupplementaryFigureCoverageWFM.png", plot = chrw_plot, device = "png", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

chrw_pb = ggplot() + geom_rect(data = index[index$V1 == "chrW",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragewpb, aes(x = xcoord, y = V7), color = "green") + theme_bw() + xlab("Chromosome W") + ylab("Coverage")

ggsave(filename = "SupplementaryFigureCoverageWPB.pdf", plot = chrw_pb, device = "pdf", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

ggsave(filename = "SupplementaryFigureCoverageWPB.png", plot = chrw_pb, device = "png", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

chrz_plot = ggplot() + geom_rect(data = index[index$V1 == "chrZ",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragezf, aes(x = xcoord, y = V7), color = "pink") + geom_point(data = coveragezm, aes(x = xcoord, y = V7), color = "blue", alpha = 0.5) + theme_bw() + xlab("Chromosome Z") + ylab("Coverage")

#chrz_plot = ggplot() + geom_rect(data = index[index$V1 == "chrZ",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragezf, aes(x = xcoord, y = V7, color = pb, stroke = .1), fill = "pink", alpha = 0.5) + geom_point(data = coveragezm, aes(x = xcoord, y = V7, color = pb), fill = "blue", alpha = 0.5) + theme_bw() + xlab("Chromosome Z") + ylab("Coverage") + theme(legend.position = "none")

ggsave(filename = "SupplementaryFigureCoverageZFM.pdf", plot = chrz_plot, device = "pdf", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

ggsave(filename = "SupplementaryFigureCoverageZFM.png", plot = chrz_plot, device = "png", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

chrz_pb = ggplot() + geom_rect(data = index[index$V1 == "chrZ",], aes(xmin = 0, xmax = V2, ymin = -3, ymax = -1)) + geom_point(data = coveragezpb, aes(x = xcoord, y = V7), color = "green") + theme_bw() + xlab("Chromosome Z") + ylab("Coverage")

ggsave(filename = "SupplementaryFigureCoverageZPB.pdf", plot = chrz_pb, device = "pdf", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

ggsave(filename = "SupplementaryFigureCoverageZPB.png", plot = chrz_pb, device = "png", scale = 0.5, width = 30, height = 30, units = "cm", dpi = 300, limitsize = FALSE)

q()
n

# CALCULATE COVERAGE RATIO BETWEEN MALE AND FEMALE

setwd("/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Coverage")

# the 3rd last column is the mean read depth
coveragewf = read.table("female.coverage.chrW.out", stringsAsFactors = F, sep = "\t")
coveragewm = read.table("male.coverage.chrW.out", stringsAsFactors = F, sep = "\t")
coveragezf = read.table("female.coverage.chrZ.out", stringsAsFactors = F, sep = "\t")
coveragezm = read.table("male.coverage.chrZ.out", stringsAsFactors = F, sep = "\t")

# coverage ratio calculated as male/female (to not have 0 as denominator)
ratiow = data.frame(Chromosome = coveragewm$V1, Start = coveragewm$V2, End = coveragewm$V3, CoverageM = coveragewm$V7, CoverageF = coveragewf$V7, Ratio = coveragewm$V7 / coveragewf$V7)
ratioz = data.frame(Chromosome = coveragezm$V1, Start = coveragezm$V2, End = coveragezm$V3, CoverageM = coveragezm$V7, CoverageF = coveragezf$V7,Ratio = coveragezm$V7 / coveragezf$V7)

write.table(x = ratiow, file = "Supplementary_table_ratiow.txt", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE)
write.table(x = ratioz, file = "Supplementary_table_ratioz.txt", sep = "\t", quote = F, col.names = TRUE, row.names = FALSE)
