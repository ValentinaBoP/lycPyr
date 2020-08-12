#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(seqinr)
load("/home/vpeona/Quadron/Quadron.lib")

fasta = args[1]

seqs = read.fasta(fasta, seqonly = TRUE)
names = unlist(attributes(read.fasta(fasta)))

for(i in 1:length(names)){

        filename = paste(names[i], ".fasta", sep = "")
        write.fasta(sequences = seqs[i], names = names[i], file.out = filename)

        out = paste(filename, ".Quadron", sep = "")
        Quadron(FastaFile = filename, OutFile = out, nCPU = 16)

}
