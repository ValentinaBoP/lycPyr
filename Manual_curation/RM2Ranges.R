RM2GRanges = function(RM){
  
  library(GenomicRanges)
  
  RM$Sequence = as.character(RM$Sequence)
  RM$Repeat = as.character(RM$Repeat)
  RM = RM[,c(5,6,7,9,10,11)]
  RM = RM[order(RM[,1], RM[,2]),]
  range_rm = GRanges(seqnames = RM$Sequence, ranges = IRanges(RM$Begin, RM$End), strand = RM$Strand, mcols = data.frame(Repeat = paste(RM$Repeat, RM$Family, sep = '$'), stringsAsFactors = FALSE))
  
  return(range_rm)
  
}
