estimateRMDensity = function(RMTrack, window){
  
  ranges = createRangesPlot(track = RMTrack, partition = window)
  
  RMDensity = data.frame(start = ranges[1:(length(ranges) - 1)], density = 0)
  
  for(i in 1:(length(ranges) - 1)){
    
    RMTemp = RMTrack[RMTrack$start >= ranges[i] & RMTrack$end <= ranges[i+1],]
    RMTemp = RMTemp %>% arrange(start) %>% group_by(g = cumsum(cummax(lag(end, default = first(end))) < start)) %>% summarise(start = first(start), end = max(end), width = first(width), strand = first(strand), query = first(query), y = first(y), category = first(category), cols = first(cols), orientation = first(orientation), fill = first(fill))
    RMTemp = RMTemp[,-1]
    RMTemp$start = as.integer(RMTemp$start)
    RMTemp$end = as.integer(RMTemp$end)
    RMTemp$width = (RMTemp$end - RMTemp$start) + 1
    RMDensity[i,2] = sum(RMTemp$width) / window
    
  }
  
  
  return(RMDensity)
  
}
