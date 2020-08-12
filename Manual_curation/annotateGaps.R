annotateGaps = function(filepath = "0", FASTA){
  
  library(Biostrings)
  
  if(filepath == "0" & class(FASTA) == "DNAStringSet"){
    
    print("Analyzing DNAStringSet")
    
  } else if(filepath == "0" & class(FASTA) != "DNAStringSet") {
    
    print("DNAStringSet object is needed!")
    break()
    
  } else {
    
    FASTA = readDNAStringSet(filepath = filepath, format = "fasta")
    
  }
  
  GAPS = data.frame(start = integer(), end = integer(), width = integer())
  
  for(i in 1:length(FASTA)){
    
    y = maskMotif(FASTA[[i]],'N')
    z = as(gaps(y),"Views")
    NEW = as.data.frame(ranges(z))
    
    if(nrow(NEW) > 0){
      
      NEW$scaffold = names(FASTA[i])
      GAPS = rbind(GAPS, NEW)
      
    } else {
      
      next()
      
    }
    
  }
  
  return(GAPS)
  
}
