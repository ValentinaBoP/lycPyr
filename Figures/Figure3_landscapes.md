# Figure 3b + Supplementary figure S8 - repeat landscapes

## FIGURE 3b lycPyrILPB

```R
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)

setwd("/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK")

aves = read.landscape(file = "./Aves/RMLandscape/lycPyr6.4.2.fasta.align.k2p.noCpG.size")
aves = clean.landscape(data = aves)
aves = fixClasses(data= aves)
aves = cleanClasses(data= aves)

merged = read.landscape(file = "./lycPyr2_rm2.1_final/RMLandscape/lycPyr6.4.2.fasta.align.k2p.noCpG.size")
merged = clean.landscape(data = merged)
merged = fixClasses(data= merged)
merged = cleanClasses(data= merged)

DISCARD = "Simple_repeat"

aves_bps = getBps(data = aves, threshold = 0)
aves_bps$Divergence = as.integer(aves_bps$Divergence)
boo = aves_bps$Divergence <= 60
aves_bps = aves_bps[boo,]
aves_bps = aves_bps[aves_bps$Class != DISCARD,]
levels = c("LINE", "SINE", "LTR", "DNA", "RC", "Satellite", "Other", "Unknown")
aves_bps$Class = factor(x = aves_bps$Class, levels = levels)
aves_bps$Colors = aves_bps$Class
levels(aves_bps$Colors) = getColors()


merged_bps = getBps(data = merged, threshold = 0)
merged_bps$Divergence = as.integer(merged_bps$Divergence)
boo = merged_bps$Divergence <= 60
merged_bps = merged_bps[boo,]
merged_bps = merged_bps[merged_bps$Class != DISCARD,]
levels = c("LINE", "SINE", "LTR", "DNA", "RC", "Satellite", "Other", "Unknown")
merged_bps$Class = factor(x = merged_bps$Class, levels = levels)
merged_bps$Colors = merged_bps$Class
levels(merged_bps$Colors) = getColors()

merged_bps$Colors = merged_bps$Class
levels(merged_bps$Colors) = getColors()
aves_bps$Colors = aves_bps$Class
levels(aves_bps$Colors) = getColors()
# create the plots
p1 = ggplot(data = aves_bps, aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("% Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "Repbase Aves repeat library") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p2 = ggplot(data = merged_bps, aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("% Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "Custom and Repbase Aves repeat libraries") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + scale_fill_manual(name = "Subclasses", values = as.character(levels(merged_bps$Colors)))

figure_double = grid.arrange(p1, p2, nrow = 1)
ggsave(figure_double, device = "png", filename = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure_D_Landscape_1.0.png", scale = .75, width = 40, height = 21, unit = "cm", limitsize = FALSE, dpi = 300)
ggsave(figure_double, device = "pdf", filename = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Figure_D_Landscape_1.0.pdf", scale = .75, width = 40, height = 21, unit = "cm", limitsize = FALSE, dpi = 300)
```

```R
fixClasses = function(data){
  
  old = c("bopSat1SU_consensus1", "bopSat2_1.1")
  new = c("bopSat1", "bopSat2")
  
  data$Subfamily[data$Subfamily == old[1]] = new[1]
  data$Subfamily[data$Subfamily == old[2]] = new[2]
  
  old = c("bopSat1SU_consensus1#Unspecified", "bopSat2_1.1#Unspecified", "CR1-10_CPB-La_astRot#Unspecified", "TguERV7i_LTR-La_lycPyr#Unspecified", "TguERV7i_LTR-Lb_lycPyr#Unspecified", "UCON4#Unspecified", "X6B_LINE-La_ptiPar#Unspecified")
  new = c("bopSat1#Satellite", "bopSat2#Satellite", "CR1-10_CPB-La_astRot#LINE/CR1", "TguERV7i_LTR-La_lycPyr#LTR/ERV1", "TguERV7i_LTR-Lb_lycPyr#LTR/ERV1", "UCON4#DNA", "X6B_LINE-La_ptiPar#LINE/CR1")
  
  for(i in 1:length(old)){
    
    boo = data$Repeat == old[i]
    data$Repeat[boo] = new[i]
    
  }
  
  elements = strsplit(x = data$Repeat, split = "#")
  classes = strsplit(x = unlist(sapply(elements, "[[", 2)), split = "/")
  data$Class = unlist(sapply(classes, "[[", 1))
  
  return(data)
  
}

cleanClasses = function(data){
  
  data$Class = sub(pattern = "SINE\\?", replacement = "SINE", x = data$Class)
  data$Class = sub(pattern = "RC|RC\\?", replacement = "Other", x = data$Class)
  data$Class = sub(pattern = "LTR\\?", replacement = "LTR", x = data$Class)
  data$Class = sub(pattern = "LINE\\?", replacement = "LINE", x = data$Class)
  data$Class = sub(pattern = "DNA\\?", replacement = "DNA", x = data$Class)
  data$Class = sub(pattern = "ARTEFACT|snRNA|scRNA|tRNA|rRNA", replacement = "Other", x = data$Class)
  data$Class = factor(x = data$Class, levels = unique(data$Class)[c(1,2,6,3,9,4,5,8,7)])
  
  
  return(data)
  
}

getColors = function(){
  
  require(RColorBrewer)
  
  # create color palette
  col = character()
  colfunc <- colorRampPalette(c("gold", "red3"))
  col = c(col, colfunc(2)[c(2,1)])
  colfunc <- colorRampPalette(c("#6BAED6", "#08306B"))
  col = c(col, colfunc(5)[c(1,3)])
  colfunc <- colorRampPalette(c("lightgreen", "darkgreen"))
  col = c(col, colfunc(8)[c(2,6)])
  colfunc <- colorRampPalette(c("#BCBDDC", "#3F007D"))
  col = c(col, colfunc(5)[c(4,1)])

  #col = c("#ED1C24", "gold", "#6BAED6", "#2B3990", "#8DC63F", "#298B29", "#BCBDDC", "#662D91")
  col = c("#ED1C24", "gold", "#6BAED6", "#2B3990", "#298B29", "#BCBDDC", "#662D91")
  return(col)
}
```




## SUPPLEMENTARY FIGURE LANDSCAPE

```R
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
```

```R
setwd("/crex/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/RMSK/lycPyr2_rm2.1_final/RMLandscape")

files = list.files(pattern = "noCpG.size")
aves_bps = data.frame()

for(i in 1:length(files)){

	# read input files
	aves = read.landscape(file = files[i])
	aves = clean.landscape(data = aves)
	aves$Species = getSpecies(files[i])
    aves = fixClasses(data = aves)
	aves = cleanClasses(data = aves)
  	# get bps occupied by repeats
  	aves_bps = rbind(aves_bps, getBps(data = aves, threshold = 0))

}

assemblies = unique(aves_bps$Species)

DISCARD = "Simple_repeat"
aves_bps$Divergence = as.integer(aves_bps$Divergence)
boo = aves_bps$Divergence <= 60
aves_bps = aves_bps[boo,]
aves_bps = aves_bps[aves_bps$Class != DISCARD,]
boo = is.na(aves_bps$Class)
aves_bps = aves_bps[!boo,]
levels = c("LINE", "SINE", "LTR", "DNA", "Satellite", "Other", "Unknown")
aves_bps$Class = factor(x = aves_bps$Class, levels = levels)
aves_bps$Colors = aves_bps$Class
levels(aves_bps$Colors) = getColors()


for(assembly in assemblies){

# create plot specific
  data_plot = aves_bps[aves_bps$Species == assembly,]  
  
  land_aves_bps = ggplot(data = data_plot, aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Percentage of divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = assembly) + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + scale_fill_manual(name = "Subclasses", values = as.character(levels(data_plot$Colors))) + theme(plot.title = element_text(face="bold"))+ coord_cartesian(ylim = c(0,6))
  
  filename = paste("RMLandscape_merged_6Mb_axis_", assembly, ".pdf", sep = "")
  # save png
  ggsave(filename = filename, plot = land_aves_bps, device = "pdf", width = 40, height = 27, units = "cm", scale = .5)

}

# combined plot


# create the plots
p1 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr_IL") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p2 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr2",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr_PB") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p3 = ggplot(data = aves_bps[aves_bps$Species == "P6156_102_pseudohap",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr_10X") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p4 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr2.1",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr2.1") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA)) + scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p5 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr2.2",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr2.2") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p6 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr2.2.1",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr2.2.1") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p7 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr3.2",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr3.2") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p8 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr4.2",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr4.2") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p9 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr5.1",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr5.1") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p10 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr6.1",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr6.1") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

p11 = ggplot(data = aves_bps[aves_bps$Species == "lycPyr6.4.2",], aes(x = as.integer(Divergence), y = Mb, fill = Class)) + geom_bar(stat = "identity") + xlab("Divergence") + ylab("Base pairs occupied (Mb)") + ggtitle(label = "lycPyr6.4.2") + theme_bw() + theme(panel.border = element_rect(colour = 'darkgrey', fill = NA))+ scale_fill_manual(name = "Subclasses", values = as.character(levels(aves_bps$Colors)))

plot_list = list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11)

figure = grid.arrange(grobs = plot_list, nrow = 5, col = 3)

ggsave(figure, device = "png", filename = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Supplementary_Figure_Landscape.png", scale = 0.25, width = 250, height = 300, unit = "cm", limitsize = FALSE, dpi = 300)
ggsave(figure, device = "pdf", filename = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figures/Supplementary_Figure_Landscape.pdf", scale = 0.25, width = 250, height = 300, unit = "cm", limitsize = FALSE, dpi = 300)
```


## Custom functions

```R
read.landscape = function(file){
  
  require(data.table)
  data = fread(input = file)
  data = data[,c(1:17)]
  
  return(data)
}

clean.landscape = function(data){
  
  # keep only the columns with: name, size, divergence, ID number
  data = data[,c(11,8,17,16)]
  names(data) = c("Repeat", "Size", "Divergence", "ID")
  elements = strsplit(x = data$Repeat, split = "#")
  data$Subfamily = unlist(sapply(elements, "[[", 1))
  classes = strsplit(x = unlist(sapply(elements, "[[", 2)), split = "/")
  data$Class = unlist(sapply(classes, "[[", 1))
  #data$Family = unlist(sapply(classes, "[[", 2))
  
  return(data)
}

getSpecies = function(filename){
  
  species = strsplit(x = filename, split = ".fasta|.fa")[[1]][1]
  
  return(species)
}

getBps = function(data, threshold = 0){
  
  data$Divergence = data$Divergence * 100
  data$RoundDiv = floor(data$Divergence)
  ### FACTOR IS THE IDENTIFIER YOU CAN WORK ON TO USE GENERAL OR SPECIFIC CLASSES OF REPEATS
  data$Factor = paste(data$Class, data$RoundDiv, sep = "$")
  data_bps = aggregate(Size ~ Factor, data, sum)
  data_bps$Class = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 1)
  data_bps$Divergence = sapply(strsplit(data_bps$Factor, "\\$"), "[[", 2)
  data_bps$Mb = data_bps$Size / 1000000
  #data_bps_sub = data_bps[data_bps$Mb > 0.01,]
  data_bps_sub = data_bps[data_bps$Mb > threshold,]
  data_bps_sub = as.data.frame(data_bps_sub)
  data_bps_sub$Species = data$Species[1]
  
  return(data_bps_sub)
}

fixClasses = function(data){
  
  old = c("bopSat1SU_consensus1", "bopSat2_1.1")
  new = c("bopSat1", "bopSat2")
  
  data$Subfamily[data$Subfamily == old[1]] = new[1]
  data$Subfamily[data$Subfamily == old[2]] = new[2]
  
  old = c("bopSat1SU_consensus1#Unspecified", "bopSat2_1.1#Unspecified", "CR1-10_CPB-La_astRot#Unspecified", "TguERV7i_LTR-La_lycPyr#Unspecified", "TguERV7i_LTR-Lb_lycPyr#Unspecified", "UCON4#Unspecified", "X6B_LINE-La_ptiPar#Unspecified")
  new = c("bopSat1#Satellite", "bopSat2#Satellite", "CR1-10_CPB-La_astRot#LINE/CR1", "TguERV7i_LTR-La_lycPyr#LTR/ERV1", "TguERV7i_LTR-Lb_lycPyr#LTR/ERV1", "UCON4#DNA", "X6B_LINE-La_ptiPar#LINE/CR1")
  
  for(i in 1:length(old)){
    
    boo = data$Repeat == old[i]
    data$Repeat[boo] = new[i]
    
  }
  
  elements = strsplit(x = data$Repeat, split = "#")
  classes = strsplit(x = unlist(sapply(elements, "[[", 2)), split = "/")
  data$Class = unlist(sapply(classes, "[[", 1))
  
  return(data)
  
}

cleanClasses = function(data){
  
  data$Class = sub(pattern = "SINE\\?", replacement = "SINE", x = data$Class)
  data$Class = sub(pattern = "RC|RC\\?", replacement = "Other", x = data$Class)
  data$Class = sub(pattern = "LTR\\?", replacement = "LTR", x = data$Class)
  data$Class = sub(pattern = "LINE\\?", replacement = "LINE", x = data$Class)
  data$Class = sub(pattern = "DNA\\?", replacement = "DNA", x = data$Class)
  data$Class = sub(pattern = "ARTEFACT|snRNA|scRNA|tRNA|rRNA", replacement = "Other", x = data$Class)
  data$Class = factor(x = data$Class, levels = unique(data$Class)[c(1,2,6,3,9,4,5,8,7)])
  
  
  return(data)
  
}

getColors = function(){
  
  require(RColorBrewer)
  
  # create color palette
  col = character()
  colfunc <- colorRampPalette(c("gold", "red3"))
  col = c(col, colfunc(2)[c(2,1)])
  colfunc <- colorRampPalette(c("#6BAED6", "#08306B"))
  col = c(col, colfunc(5)[c(1,3)])
  colfunc <- colorRampPalette(c("lightgreen", "darkgreen"))
  col = c(col, colfunc(8)[c(2,6)])
  colfunc <- colorRampPalette(c("#BCBDDC", "#3F007D"))
  col = c(col, colfunc(5)[c(4,1)])

  col = c("#ED1C24", "gold", "#6BAED6", "#2B3990", "#298B29", "#BCBDDC", "#662D91")
  return(col)
}
```
