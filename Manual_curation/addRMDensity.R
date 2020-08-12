addRMDensity = function(track, window = 10000){
  
  track = track[track$category == "RepeatMasker",]
  
  RMTrack = partitionTrack(track = track, partition = window)
  
  density = estimateRMDensity(RMTrack, window)
  
  return(density)
  
}
