# GC content - Figure 3

Remake the panel about GC content: make a scatterplots for each assembly showing the windows over 3 sigma from the mean


```bash
cd /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Intermediate/GCContent

ml R_packages/3.6.0

R
```

```R
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

meanCol = 0.50
meanGenome = 0.4
threeSigma = 0.588
# calculated before on the basis of PacBio with the function quantile(data$CorGC, probs = 0.95)

files = list.files(pattern = "windows.nuc")

data_plot = data.frame()

for(file in files){
	print(file)
	data = read.table(file, stringsAsFactors = F, header = T, comment.char = "", col.names = c("Scaffold", "Start", "End", "Window", "AT", "GC", "A", "C", "G", "T", "N", "OtherChr", "Length"), sep = "\t")

	data$CorLength = data$Length - data$N
	data$CorGC = (data$G + data$C) / data$CorLength
	
	assembly = sub(pattern = ".fasta.*", replacement = "", x = file)

	temp = data[data$GC >= threeSigma,]
	temp$Assembly = assembly
	temp = temp[,c(16,15,14)]
	
	data_plot = rbind(data_plot, temp)
}

data_plot = data_plot[data_plot$CorLength == 1000,]

data_plot$Assembly = factor(x = data_plot$Assembly, levels = unique(data_plot$Assembly)[c(5,3,4,6,2,1,7)])

#plot = ggplot(data= data_plot, aes(Assembly, CorGC)) + geom_point() + theme_bw()

plot = ggplot(data= data_plot, aes(Assembly, CorGC)) + geom_jitter(aes(alpha = .5)) + theme_bw()

ggsave(plot = plot, filename = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/5.Manuscript/Results/Figure3_GCpanel_Dec.pdf", device = "pdf", units = "cm", height = 21, width = 30, dpi = 300, limitsize = F)
```
