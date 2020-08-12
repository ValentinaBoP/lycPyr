getGenomicRanges = function(FASTA, lenChr){
  # get the ranges for the actual sequences in the chromosomes excluding the gaps
  
  source(file = "annotateGaps.R")
  library(GenomicRanges)  
  
  GAPS = annotateGaps(FASTA = FASTA)
  
  GAPS = GAPS[GAPS$scaffold %in% lenChr$chromosome,]
  
  temp = lenChr$len
  names(temp) = lenChr$chromosome
  
  GAPS$scaffold = as.character(GAPS$scaffold)
  gr <- with(GAPS, GRanges(seqnames = GAPS$scaffold, IRanges(start, end), seqlengths = temp))
  
  gaps <- gaps(gr)
  
  gaps = gaps[strand(gaps) == "*"]
  
  return(gaps)
}
