addLASTZTrack_out = function(file, category, pattern, track){
  
  # add LASTZ output of the type
  # "name1" "start1" "end1" "length1" "strand1" "name2" "start2" "end2" "length2" "strand2"
  
  library(data.table)
  
  newTrack = fread(file, stringsAsFactors = F)
  names(newTrack)[1] = "name1"
  
  newTrack = newTrack[grepl(pattern = pattern, x = newTrack$name2), c(6,7,8,9,5,1)]
  
  #newTrack$width = newTrack$end1 - newTrack$start1 + 1
  
  newTrack$category = category
  
  names(newTrack) = c("seqnames", "start", "end", "width", "strand", "query", "category")
  newTrack$y = 1
  newTrack = newTrack[,c(1:6, 8, 7)]
  newTrack$cols = "black"
  newTrack$orientation = "black"
  
  track$y = as.integer(track$y)
  track$y = track$y + 1
  
  if(grepl(pattern = "edges", x = newTrack$query[1])){
    
    newTrack$query = sapply(strsplit(x = newTrack$query, split = "_edges"), "[[", 1)
    
  }
  
  if(grepl(pattern = "__", x = newTrack$seqnames[1])){
    
    newTrack$seqnames = sapply(strsplit(x = newTrack$seqnames, split = "__"), "[[", 1)
    
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
  
  MERGED = rbind(track, newTrack)
  
  #set custom colors
  MERGED$cols = "black"
  
  MERGED = MERGED[order(MERGED[,1], MERGED[,2])]
  return(MERGED)
  
}
