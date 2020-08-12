subsetByPartition = function(track, range){
  
  boo = track$start < range & track$end > range
  boo = as.logical(boo)
  
  if(sum(boo) > 0){
    
    for(i in 1:sum(boo)){
      
      temp = track[which(boo)[i],]
      
      left = temp
      left[,3] = as.integer(range)
      
      right = temp
      right[,2] = as.integer(range + 1)
      
      track[which(boo)[i],] = left
      track = rbind(track, right)
      
    }
    
  }
  
  track = track[order(track[,1], track[,2]),]
  
  return(track)

}
