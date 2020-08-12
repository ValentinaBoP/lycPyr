readRM = function(filepath){
  
  RM = read.table(file = filepath, skip = 3, fill = TRUE, stringsAsFactors = FALSE)
  RM = RM[!is.na(RM[,2]),]
  names(RM) = c("Score", "Divergence", "Deletion", "Insertion", "Sequence", "Begin", "End", "SequenceLeft", "Strand", "Repeat", "Family", "RepeatBegin", "RepeatEnd", "RepeatLeft", "ID")
  RM$Strand = sub(pattern = "C", replacement = '-', x = RM$Strand)
  return(RM)
  
}
