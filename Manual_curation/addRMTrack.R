addRMTrack = function(data, file, pattern){
  
  # add RM output of the type
  # names(RM) = c("Score", "Divergence", "Deletion", "Insertion", "Sequence", "Begin", "End", "SequenceLeft", "Strand", "Repeat", "Family", "RepeatBegin", "RepeatEnd", "RepeatLeft", "ID")
  
  library(data.table)
  
  newTrack = readRM(file)
    
  newTrack = newTrack[grepl(pattern = pattern, x = newTrack$Sequence), c(5,6,7,9,11)]
  
  newTrack$width = newTrack$End - newTrack$Begin + 1
  
  newTrack$category = "RepeatMasker"
  
  names(newTrack) = c("seqnames", "start", "end", "strand", "query", "width", "category")
  newTrack$y = "NA"
  newTrack = newTrack[,c(1,2,3,6,4,5,8,7)]
  newTrack$cols = "orange"
  
  newTrack$orientation = "purple"
  
  data$y = as.integer(data$y)
  LEN = max(data$y)
  
  newTrack$y = LEN
  data[data$y == LEN, 7] = LEN + 1
  
  MERGED = rbind(data, newTrack)
  
  MERGED = MERGED[order(MERGED[,1], MERGED[,2]),]
  row.names(MERGED) = 1:nrow(MERGED)
  return(MERGED)
  
}
