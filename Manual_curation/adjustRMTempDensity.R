adjustRMTempDensity = function(temp_RM_density, window = 10000){
  
  temp_RM_density$end = c(temp_RM_density$start[2:nrow(temp_RM_density)] - 1, temp_RM_density$start[nrow(temp_RM_density)] + window - 1)
  temp_RM_density$ymin = floor(temp_RM_density$y)
  
  rr <- range(temp_RM_density$density)
  svals <- (temp_RM_density$density-rr[1])/diff(rr)
  # [1] 0.2752527 0.0000000 0.9149839 0.3680242 1.0000000 0.2660587
  
  ## Play around with ends of the color range
  f <- colorRamp(c("white", "yellow", "red"))
  temp_RM_density$colors <- rgb(f(svals)/255)
  
  ## Check that it works
  #image(seq_along(svals), 1, as.matrix(seq_along(svals)), col=colors,
  #      axes=FALSE, xlab="", ylab="")
  
  return(temp_RM_density)
  
}
