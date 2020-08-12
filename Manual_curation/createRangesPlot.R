createRangesPlot = function(track, partition = 10000000){
  
  MAX = max(track$end)
  
  ranges = as.integer(c(seq(from = 1, to = MAX, by = partition), MAX))
  
  return(ranges)
  
}
