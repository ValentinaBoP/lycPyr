#Rscript --vanilla summaryGenePresence.R input.table

args=(commandArgs(TRUE))

table = args[1]

index = read.table(file = "/proj/uppstore2018073/private/Valentina/lycPyr/Intermediate/Exonerate/lycPyrZW_proteins.fasta.fai", sep = "\t", stringsAsFactors = F, header = FALSE)

table = read.table(file = table, sep = "\t", stringsAsFactors = F, header = FALSE)

# coverage of the protein from alignment
table$V7 = table$V3 - table$V2 + 1

# length of the query protein
for(i in unique(table$V1)){
  
  table$V8[table$V1 == i] = index$V2[index$V1 == i]
  
}

# proportion of covered protein (alignment / length)
table$V9 = table$V7 / table$V8

# similarity filter
table = table[table$V6 >= 0.9,]

# complete
complete = unique(table$V1[table$V9 >= 0.95])

# partial 
partial = unique(table$V1[table$V9 < 0.95])
partial = partial[!(partial %in% complete)]

# absent
absent = index$V1[!(index$V1 %in% table$V1)]

filename = paste(args[1], ".summary", sep = "")

sink(filename)
cat(paste("Complete", length(complete), sep = "\t"))
cat("\n")
cat(paste("Partial", length(partial), sep = "\t"))
cat("\n")
cat(paste("Absent", length(absent), sep = "\t"))
cat("\n")
sink()
