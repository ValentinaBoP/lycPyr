generatePlots = function(data, piTOT = data.frame(), pattern, categories, GC){
  
  library(ggplot2)
  #PATTERN = "PGA_scaffold5"
  categories = as.factor(c("Genome", "RepeatMasker", categories))
  
  for(p in pattern){
    
    chr2Plot = data[data$seqnames == p,]
    
    chr2Plot = partitionTrack(track = chr2Plot)
    
    ranges = createRangesPlot(track = chr2Plot)
    
    GC2plot = GC[GC$chrom == p,]
    
    if(length(piTOT) > 0){
      
      pi2plot = piTOT[piTOT$seqnames == p,]
      
    } else {
      
      pi2plot = data.frame()
      
    }
    
    listPlots(chr2Plot = chr2Plot, pi2plot = pi2plot, ranges = ranges, categories = categories, pattern = p, GC2plot = GC2plot)
    #listPlots(chr2Plot = chr2Plot, pi2plot = pi2plot, ranges = ranges, categories = categories, pattern = p)
    
  }
  
}
