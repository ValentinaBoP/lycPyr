partitionTrack = function(track, partition = 10000000){
  
  MAX = max(track$end)
  
  ranges = createRangesPlot(track, partition)
  
  for(i in 2:length(ranges)){
    
    track = subsetByPartition(track = track, range = ranges[i])
    
  }
  
  return(track)

}
