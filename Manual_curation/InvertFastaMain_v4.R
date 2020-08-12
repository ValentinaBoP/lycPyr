# running the code in
# /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Code
# ml bioinfo-tools
# ml R_packages/3.5.0

#read all the input files

source("/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Code/InvertFastaFunctions.R")

general = read.table(file = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Data/Contig_Gap_Coordinates.txt", sep ="\t", stringsAsFactors = FALSE, header = FALSE)
chr_list = read.table(file = "/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Data/chr2modify.contigs", sep ="\t", stringsAsFactors = FALSE, header = FALSE)

# for testing
#chr_list[4,] = c("PGA_scaffold1", "Contig3,Contig4")
#chr_list[5,] = c("PGA_scaffold1", "Contig6,Contig7,Contig8")
#chr_list[6,] = c("PGA_scaffold1", "Contig11,Contig12,Contig13-Contig12")

chromosomes = unique(chr_list[,1])

setwd("/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Intermediate/")

for(chr in chromosomes){

        # get all the contigs to modify relative to the particular chromosome (chr)
        contigs = chr_list[chr_list[,1] == chr, 2]
        # get contigs to flip multiple times
        multiple_contigs = contigs[grepl(pattern = "-", x = contigs)]
        # discard them from the original contig list
        contigs = contigs[!grepl(pattern = "-", x = contigs)]
        # subset the contigs that must be flipped together
        double_contigs = contigs[grepl(pattern = ",", x = contigs)]
        # discard them from the original contig list
        contigs = contigs[!grepl(pattern = ",", x = contigs)]

        # get the coordinates of the contigs to be flipped
        coordinates = getContigs(chr, general, contigs)

        # deal with double contigs (contigs to be flipped together)
        if(length(double_contigs) > 0){

                coordinates = rbind(coordinates, getDoubleContigs(chr, general, double_contigs))

        }


        # deal with contigs that must be flipped twice
        if(length(multiple_contigs) > 0){

                coordinates = rbind(coordinates, getMultipleContigs(chr, general, multiple_contigs))

        }


        # write a coordinate file relative to the current chromosome (chr)
        write.table(x = coordinates, file = paste("/proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Intermediate/", chr, ".coordinates", sep = ""), sep = "\t", col.names = FALSE, row.names = FALSE, quote = F)

        # get the fasta sequence of the chromosome of interest
        cmd = paste("perl /proj/sllstore2017073/private/scripts/extractFromFasta.pl /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Data/lycPyr4.1_OnlyChr.fasta single '", chr, "' > ", chr, ".fasta", sep = "")
        system(cmd)

        # flip contigs using the Python script
        cmd = paste("python /proj/sllstore2017073/private/BirdsOfParadise/lycPyr_TheGenome/6.TheGenome/1.ManualCuration/Code/invert_fasta_v2.py", chr)
        system(cmd)

}
