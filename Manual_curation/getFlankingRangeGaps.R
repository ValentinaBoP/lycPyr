getFlankingRangeGaps = function(GAPS, flank = 500){
  
  library(GenomicRanges)
  
  out_gaps = data.frame(start = integer(), end = integer(), scaffold = character(), width = integer())
  
  NEW = data.frame(start = GAPS$start - flank, end = GAPS$start, scaffold = GAPS$scaffold, width = GAPS$width)
  out_gaps = rbind(out_gaps, NEW)
  
  NEW = data.frame(start = GAPS$end, end = GAPS$end + flank, scaffold = GAPS$scaffold, width = GAPS$width)
  out_gaps = rbind(out_gaps, NEW)
  out_gaps = out_gaps[order(out_gaps[,3], out_gaps[,1]),]
  
  range_gaps = GRanges(out_gaps$scaffold, IRanges(out_gaps$start, out_gaps$end), mcols = data.frame(width = out_gaps$width))
  
  return(range_gaps)
}
