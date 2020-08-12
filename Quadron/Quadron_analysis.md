# Filter results

```bash
ml bioinfo-tools BEDTools

	for file in $( ls *.Quadron )
	do
	 grep -v 'NOTE\|HEADER\|NA' $file | tr -s ' ' '\t' | awk '$5 > 19 {print $0}' > $file.filter
	 chr="${file%.*}"
	 chr="${chr%.*}"
	 cat $file.filter | sed "s/DATA:/$chr/g" | awk '{print $1, $2-1, $2+$4-1, $1"_"NR, $5, $3}' OFS='\t' > $file.filter.bed
	done

# FASTA=/home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/6.Juicer/Data/lycPyr7.1.fasta
# BED=chr10.fasta.Quadron.filter.bed
# bedtools getfasta -fi $FASTA -bed $BED -fo chr10.fasta.Quadron.filter.fasta

```

# Base pairs covered in G4 in the assemblies

```bash
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr7.1/Output

cat *.bed > lycPyr7.1_Quadron.bed

awk '{print $3-$2}' lycPyr7.1_Quadron.bed | awk '{sum += $1} END {print sum}'
# 9177510

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr7.4/Output

cat *.bed > lycPyr7.4_Quadron.bed

bedtools merge -i lycPyr7.4_Quadron.bed | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'

# 9177208


cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_SN1/Output

cat *.bed > lycPyr_SN1_Quadron.bed

bedtools merge -i lycPyr_SN1_Quadron.bed | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'
# 6390589

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_SN2/Output

cat *.bed > lycPyr_SN2_Quadron.bed

bedtools merge -i lycPyr_SN2_Quadron.bed | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'
# 7342290


cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_DV/Output

cat *.bed > lycPyr_DV_Quadron.bed

bedtools merge -i lycPyr_DV_Quadron.bed | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'
# 8432759

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr1/Output

cat *.filter.bed > lycPyr_Quadron.bed

bedtools sort -i lycPyr_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'
# 3001975

bedtools sort -i lycPyr_Quadron.bed | bedtools merge -i stdin > total.merge.bed

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_PB/Output

cat *.bed > lycPyr_PB_Quadron.bed

bedtools merge -i lycPyr_PB_Quadron.bed | awk '{print $3-$2}'| awk '{sum += $1} END {print sum}'
# 8428895

cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_DV_arrow/Output

cat *.bed > lycPyr_DV_arrow_Quadron.bed

bedtools merge -i lycPyr_DV_arrow_Quadron.bed | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 8955400



cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr_ILPB/
mkdir Output
mv *.fasta* ./Output

cd Output

cat *.bed > lycPyr_ILPB_Quadron.bed

bedtools sort -i lycPyr_ILPB_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 8755646

bedtools sort -i lycPyr_ILPB_Quadron.bed | bedtools merge -i stdin > total.merge.bed


cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr3/

cd Output

cat *.bed > lycPyr3_Quadron.bed

bedtools sort -i lycPyr3_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 8755646

bedtools sort -i lycPyr3_Quadron.bed | bedtools merge -i stdin > total.merge.bed



cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr4/

cd Output

cat *filter.bed > lycPyr4_Quadron.bed

bedtools sort -i lycPyr4_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 8755646

bedtools sort -i lycPyr4_Quadron.bed | bedtools merge -i stdin > total.merge.bed


cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr5/

cd Output

cat *filter.bed > lycPyr5_Quadron.bed

bedtools sort -i lycPyr5_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 8755646

bedtools sort -i lycPyr5_Quadron.bed | bedtools merge -i stdin > total.merge.bed


cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr1/

cd Output

cat *filter.bed > lycPyr1_Quadron.bed

bedtools sort -i lycPyr1_Quadron.bed | bedtools merge -i stdin | awk '{print $3-$2}' | awk '{sum += $1} END {print sum}'
# 7501295

bedtools sort -i lycPyr1_Quadron.bed | bedtools merge -i stdin > total.merge.bed
```
# Quadron density

```bash
cd /home/vpeona/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/Quadron/lycPyr7.4/Output/density

ml bioinfo-tools BEDTools
cp ../Output/chr*filter.bed ./
mv ../Output/index.fai ./

cat *bed | bedtools sort -i stdin > total.bed
bedtools merge -i total.bed > total.merge.bed
bedtools coverage -b total.merge.bed -a index.fai > density.txt
```

## Results

```
chr1    1       119122530       9143    389128  119122529       0.0032666
chr10   1       20644987        3686    166448  20644986        0.0080624
chr11   1       21157310        3410    147622  21157309        0.0069774
chr12   1       23449977        3984    174276  23449976        0.0074318
chr13   1       18905119        3842    164771  18905118        0.0087157
chr14   1       16917316        4006    175423  16917315        0.0103694
chr15   1       13473298        3918    182743  13473297        0.0135633
chr17   1       11120114        3777    168125  11120113        0.0151190
chr18   1       12081340        3864    172068  12081339        0.0142425
chr19   1       11658745        3510    155667  11658744        0.0133520
chr1A   1       74173823        7869    381909  74173822        0.0051488
chr2    1       152213986       11717   540206  152213985       0.0035490
chr20   1       18848860        6273    281630  18848859        0.0149415
chr21   1       3173624 2874    173996  3173623 0.0548257
chr23   1       7387094 5026    244077  7387093 0.0330410
chr24   1       9868923 4157    188989  9868922 0.0191499
chr26   1       7472877 5140    252664  7472876 0.0338108
chr27   1       5379378 5827    305239  5379377 0.0567424
chr28   1       1989137 2230    125790  1989136 0.0632385
chr3    1       111593273       9685    441448  111593272       0.0039559
chr33   1       2347063 4636    253893  2347062 0.1081748
chr4    1       74280243        6920    299736  74280242        0.0040352
chr4A   1       20447179        3432    146586  20447178        0.0071690
chr5    1       64599582        8305    355061  64599581        0.0054963
chr6    1       36399617        4262    174155  36399616        0.0047845
chr7    1       39599063        4332    177694  39599062        0.0044873
chr8    1       31934939        3983    158785  31934938        0.0049721
chr9    1       26246584        3667    153253  26246583        0.0058390
chrUN1  1       2411798 1705    107649  2411797 0.0446344
chrUN2  1       2310713 2040    102913  2310712 0.0445374
chrUN3  1       2248212 4313    269805  2248211 0.1200088
chrUN4  1       2167051 1101    53051   2167050 0.0244807
chrUN5  1       1430982 1522    76256   1430981 0.0532893
chrUN6  1       1160716 3852    241182  1160715 0.2077874
chrW    1       19679164        3168    119832  19679163        0.0060893
chrZ    1       74491454        5994    275955  74491453        0.0037045
```

## Plot in RStudio (local)

```R
library(ggplot2)

G4 = read.table(text = "
chr1    1       119122530       9143    389128  119122529       0.0032666       0.403891
chr10   1       20644987        3686    166448  20644986        0.0080624       0.441167
                chr11   1       21157310        3410    147622  21157309        0.0069774       0.435586
                chr12   1       23449977        3984    174276  23449976        0.0074318       0.448957
                chr13   1       18905119        3842    164771  18905118        0.0087157       0.458794
                chr14   1       16917316        4006    175423  16917315        0.0103694       0.460750
                chr15   1       13473298        3918    182743  13473297        0.0135633       0.472187
                chr17   1       11120114        3777    168125  11120113        0.0151190       0.494422
                chr18   1       12081340        3864    172068  12081339        0.0142425       0.483155
                chr19   1       11658745        3510    155667  11658744        0.0133520       0.480880
                chr1A   1       74173823        7869    381909  74173822        0.0051488       0.410084
                chr2    1       152213986       11717   540206  152213985       0.0035490       0.400756
                chr20   1       18848860        6273    281630  18848859        0.0149415       0.479158
                chr21   1       3173624 2874    173996  3173623 0.0548257       0.517804
                chr23   1       7387094 5026    244077  7387093 0.0330410       0.521430
                chr24   1       9868923 4157    188989  9868922 0.0191499       0.495445
                chr26   1       7472877 5140    252664  7472876 0.0338108       0.527575
                chr27   1       5379378 5827    305239  5379377 0.0567424       0.532182
                chr28   1       1989137 2230    125790  1989136 0.0632385       0.531315
                chr3    1       111593273       9685    441448  111593272       0.0039559       0.407193
                chr33   1       2347063 4636    253893  2347062 0.1081748       0.567865
                chr4    1       74280243        6920    299736  74280242        0.0040352       0.402267
                chr4A   1       20447179        3432    146586  20447178        0.0071690       0.445173
                chr5    1       64599582        8305    355061  64599581        0.0054963       0.418710
                chr6    1       36399617        4262    174155  36399616        0.0047845       0.427004
                chr7    1       39599063        4332    177694  39599062        0.0044873       0.420810
                chr8    1       31934939        3983    158785  31934938        0.0049721       0.430614
                chr9    1       26246584        3667    153253  26246583        0.0058390       0.438469
                chrUN1  1       2411798 1705    107649  2411797 0.0446344       0.437013
                chrUN2  1       2310713 2040    102913  2310712 0.0445374       0.457090
                chrUN3  1       2248212 4313    269805  2248211 0.1200088       0.533989
                chrUN4  1       2167051 1101    53051   2167050 0.0244807       0.434508
                chrUN5  1       1430982 1522    76256   1430981 0.0532893       0.463384
                chrUN6  1       1160716 3852    241182  1160715 0.2077874       0.568701
                chrW    1       19679164        3168    119832  19679163        0.0060893       0.438118
                chrZ    1       74491454        5994    275955  74491453        0.0037045       0.401739", stringsAsFactors = F, h = F)

G4 = G4[order(G4[,3], decreasing = T),]
lev = as.character(G4$V1)
lev = lev[c(1:28,30,34,29,31:33,35,36)]
G4$V1 = factor(x = G4$V1, levels = lev)
G4$type = "macro"
boo = G4[,3] < 20000000
G4$type[boo] = "micro"

plot = ggplot(data = G4) + geom_bar(aes(x = V1, y = V7, fill = type), stat = "identity", colour = "black", size = 0.4) + theme_bw() + xlab("Chromosomes") + ylab("Density") + scale_fill_manual(values = c("macro" = "gray80", "micro" = "grey40")) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(panel.grid.minor = element_blank())

plot = ggplot(data = G4) + geom_point(aes(x = V3, y = V7, shape = factor(type), size = 12, colour = V8), stat = "identity", colour = "black", size = 0.4) + theme_bw() + xlab("Chromosomes") + ylab("Density") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(panel.grid.minor = element_blank()) + scale_colour_gradientn(colours = terrain.colors(10))

ggplot(data = G4) + geom_point(aes(x = V3, y = V7, shape = factor(type), size = 12, colour = V8), stat = "identity", colour = "black", size = 0.4) + theme_bw() + xlab("Chromosomes") + ylab("Density") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(panel.grid.minor = element_blank()) + scale_colour_gradient(low = "dark blue", high = "red3", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")

plot = ggplot(data = G4) + geom_point(aes(x = V3, y = V7, shape = factor(type), size = 12, colour = V8)) + theme_bw() + xlab("Chromosomes") + ylab("Density") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(panel.grid.minor = element_blank()) + scale_colour_gradient(low = "grey60", high = "black", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")

ggsave(filename = "G4_density.pdf", plot = plot, device = "pdf", scale = 0.5, width = 21, height = 21, units = "cm", dpi = 300)
ggsave(filename = "G4_density.png", plot = plot, device = "png", scale = 0.5, width = 21, height = 21, units = "cm", dpi = 300)
```

## Plot of the abundance per assembly
(on local)
```R
data = read.table(text =
"IL	7501295
SN1	6390589
SN2	7342290
PB	8428895
ILPB	8755646
PBDV	8432808
FINAL	9177208", stringsAsFactors = F, h = F, sep = "\t")

lev = data$V1[c(1:length(data$V1))]
data$V1 = factor(x = data$V1, levels = lev)

data$col = "col"
plot = ggplot(data = data) + geom_bar(aes(x = V1, y = V2, fill = col), stat = "identity", colour = "black", size = 0.4) + theme_bw() + theme(panel.grid.minor = element_blank()) + ylab("Bases occupied by G4 motifs") + xlab("Assemblies") + scale_fill_manual(values = c("col" = "grey60")) + theme(legend.position = "none")
ggsave(filename = "G4_assemblies_oct_Figure3f.png", plot = plot, device = "png", scale = 0.5, width = 21, height = 21, units = "cm", dpi = 300)
ggsave(filename = "G4_assemblies_oct_Figure3f.pdf", plot = plot, device = "pdf", scale = 0.5, width = 21, height = 21, units = "cm", dpi = 300)
```
