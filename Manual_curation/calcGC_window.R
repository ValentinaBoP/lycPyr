calcGC_window = function(file = NULL, fasta = NULL, window, pattern = character()){
  
  # Valentina Peona 01/07/2018
  
  # usage by directly reading the fasta file
  # source("GC_content.R")
  # GCcontent_df = calcGC_window(file = "genome.fasta", window = 10000)
  
  # usage with the fasta file already imported as DNAStringSet object
  # require(Biostrings)
  # source("GC_content.R")
  # fasta = readDNAStringSet(filepath = "path/to/genome.fasta")
  # GCcontent_df = calcGC_window(file = character(), fasta = fasta, window = 10000)
  
  # example: calculate GC content on only two chromosomes
  # source("GC_content.R")
  # GCcontent_df = calcGC_window(file = character(), fasta = fasta, window = 10000, pattern = c("chr1", "chr2"))
  
  # required libraries
  require(Biostrings)
  require(BSgenome)
  require(seqinr)
  
  # check window size
  if(length(window) == 0 | window <= 0){
    
    window = 10000
    
  }
  
  # check if the correct data have been provided
  if(length(fasta) == 0 & length(file) == 0){
    
    print("you must provide either the filepath to the fasta file or a DNAStringSet object")
    stop()
    
  } else if (length(fasta) == 0 & length(file) > 0){
    
    print("reading fasta file")
    fasta = readDNAStringSet(filepath = file)
    print("fasta file imported")
    
  }
  
  # if pattern has been given, subset the sequences and names by it
  if(length(pattern) > 0){
    
    boo = which(fasta@ranges@NAMES %in% pattern)
    fasta = fasta[boo]
    
  }
  
  chrWidth = data.frame(chrom = fasta@ranges@NAMES, width = fasta@ranges@width)
  
  GCcontent = createWindows(data = chrWidth, window = window)
  
  print("calculating GC content")
  list_seqs = strsplit(as.character((getSeq(fasta, as(GCcontent, "GRanges")))), split = "")
  GCcontent$gc = sapply(list_seqs, GC)
  
  return(GCcontent)
  
}

createWindows = function(data, window){
  
  newData = data.frame(chrom = character(), start = integer(), end = integer())
  
  for(i in 1:nrow(data)){
    
    start = as.integer(seq(from = 1, to = data$width[i], by = window))
    end = as.integer(c(start[2:length(start)] - 1, data$width[i]))
    
    newData = rbind(newData, data.frame(chrom = data$chrom[i], start = start, end = end))
    
  }
  
  return(newData)
  
}

correctGCcontent = function(GC){
  
  rr <- range(GC$gc)
  svals <- (GC$gc-rr[1])/diff(rr)
  # [1] 0.2752527 0.0000000 0.9149839 0.3680242 1.0000000 0.2660587
  
  ## Play around with ends of the color range
  f <- colorRamp(c("white", "yellow", "red"))
  GC$colors <- rgb(f(svals)/255)
  
  return(GC)
  
}
