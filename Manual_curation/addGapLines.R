addGapLines = function(data, pi){
  
  gapLines = data[data$category == "Genome",]
  
  gapLines = data.frame(x = c(gapLines$start, gapLines$end), xend = c(gapLines$start, gapLines$end), y = c(gapLines$y, gapLines$y))
  gapLines$yend = max(pi$newPi)
  
  return(gapLines)
  
}
