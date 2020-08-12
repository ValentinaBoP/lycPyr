addSatsumaTrack = function(file, category, pattern, track){
  
  # add Satsuma output of the type
  # "name1" "start1" "end1" "name2" "start2" "end2" "identity" "strand"
  
  library(data.table)
  
  newTrack = fread(file, stringsAsFactors = F, header = F, col.names = c("name1", "start1", "end1", "name2", "start2", "end2", "identity", "strand"))
  
  newTrack = newTrack[grepl(pattern = pattern, x = newTrack$name1), c(1,2,3,8,4)]
  
  newTrack$width = newTrack$end1 - newTrack$start1 + 1
  
  newTrack$category = category
  
  names(newTrack) = c("seqnames", "start", "end", "strand", "query", "width", "category")
  newTrack$y = 1
  newTrack = newTrack[,c(1,2,3,6,4,5,8,7)]
  newTrack$cols = "black"
  newTrack$seqnames = as.character(newTrack$seqnames)
  newTrack$seqnames = sapply(strsplit(x = newTrack$seqnames, split = "\\|"), "[[", 1)
  
  newTrack$start = newTrack$start + 1
  newTrack$end = newTrack$end + 1
  newTrack$orientation = "black"
  
  if(grepl(pattern = "edges", x = newTrack$query[1])){
    
    newTrack$query = sapply(strsplit(x = newTrack$query, split = "_edges"), "[[", 1)
    
  }
  
  prova = newTrack[newTrack$query != "Genome",] %>% group_by(seqnames, query) %>% summarize(n_unique = n_distinct(strand))
  prova_reverse = newTrack[newTrack$query != "Genome",] %>% group_by(query) %>% summarize(n_unique = n_distinct(seqnames))
  
  for(i in 1:nrow(prova_reverse)){
    
    element = prova_reverse$query[i]
    
    if(prova_reverse$n_unique[i] == 1){
      
      newTrack[newTrack$query == element, 10] = "green"  
      
    } else {
      
      newTrack[newTrack$query == element, 10] = "orange"
      
    }
    
  }
  
  for(i in 1:nrow(prova)){
    
    if(prova$n_unique[i] == 2){
      
      newTrack[newTrack$seqnames == prova$seqnames[i] & newTrack$query == prova$query[i] & newTrack$strand == "+", 10] = "red"
      newTrack[newTrack$seqnames == prova$seqnames[i] & newTrack$query == prova$query[i] & newTrack$strand == "-", 10] = "blue"
      
    }
    
  }
  
  track$y = as.integer(track$y)
  track$y = track$y + 1
  
  MERGED = rbind(track, newTrack)
  
  #set custom colors
  MERGED$cols = "black"
  #MERGED[MERGED$strand == "+", 8] = "red"
  #MERGED[MERGED$strand == "-", 8] = "blue"
  
  #MERGED$orientation = "black"
  

  
  MERGED = MERGED[order(MERGED[,1], MERGED[,2])]
  return(MERGED)
  
  
  
  
}
