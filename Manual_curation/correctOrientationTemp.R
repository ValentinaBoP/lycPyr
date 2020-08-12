correctOrientationTemp = function(track, cats){
  
  cats = as.character(cats)
  cats = cats[!cats %in% c("Genome", "RepeatMasker", "Repeat density")]
  track$strand = as.character(track$strand)
  
  
  for(cat in cats){
    
    orientation = track[track$category == cat,] %>% group_by(query) %>% summarise(strandmin= sum(strand == "-"), strandplus = sum(strand == "+"))
    
    if(nrow(orientation) > 0){
      
      for(i in 1:nrow(orientation)){
        
        if((orientation$strandmin[i] > orientation$strandplus[i])){
          
          track$strand[track$query == orientation$query[i] & track$strand == "-"] = "C"
          track$strand[track$query == orientation$query[i] & track$strand == "+"] = "-"
          track$strand[track$query == orientation$query[i] & track$strand == "C"] = "+"
          
          ### ATTENTION TO THE COLORS!!!
          track$orientation[track$query == orientation$query[i] & track$strand == "-" & track$orientation == "red"] = "blue"
          track$orientation[track$query == orientation$query[i] & track$strand == "+" & track$orientation == "blue"] = "red"
          
        } else {
          
          next
          
        }
        
      }
      
    } else {
      
      next
      
    }
    
  }
  
  return(track)
  
}
