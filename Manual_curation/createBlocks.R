createBlocks = function(track, cats){
  
  cats = as.character(cats)
  cats = cats[!cats %in% c("Genome", "RepeatMasker", "Repeat density")]
  
  i = 1
  e = 1
  s = 1
  COND = TRUE
  blocks = data.frame(query = character(), xmin = integer(), xmax = integer(), category = character())
  
  
  for(cat in cats){
    
    track_cat = track[track$category == cat,]
    
    if(nrow(track_cat) > 0){
      
      track_cat$query = as.character(track_cat$query)
      row.names(track_cat) = 1:nrow(track_cat)
      
      i = 1
      e = 1
      s = 1
      
      while(COND){
        
        e = i
        s = i
        
        while(track_cat[i,6] == track_cat[e,6]){
          
          if((e + 1) > nrow(track_cat)){
            
            break
            
          } else {
            
            e = e + 1
            
          }
          
        }
        
        blocks = rbind(blocks, data.frame(query = track_cat[i,6], xmin = min(track_cat[c(s,(e - 1)), c(2,3)]), xmax = max(track_cat[c(s, (e - 1)), c(2,3)]), category = cat))
        
        if(e == nrow(track_cat)){
          
          N = nrow(blocks)
          blocks[N,] = c(query = track_cat[i,6], xmin = min(track_cat[c(s,(e)), c(2,3)]), xmax = max(track_cat[c(s, (e)), c(2,3)]), category = cat)
          break
          
        }
        
        i = e
        s = e
        
      }
      
      
    } else {
      
      next
      
    }
    
  }
  
  # add colors to blocks
  blocks$colors = rep_len(x = c("black", "gray60", "gray90"), length.out = nrow(blocks))
  
  # function to adjust the y coordinates!
  dfY = unique(track[,c(7,8)])
  o = match(x = blocks$category, table = dfY$category)
  blocks$ymin = dfY$y[o]
  blocks$ymax = blocks$ymin + 0.5
  
  blocks[,1] = as.character(blocks[,1])
  blocks[,4] = as.character(blocks[,4])
  blocks[,2] = as.integer(blocks[,2])
  blocks[,3] = as.integer(blocks[,3])
  
  return(blocks)
}
