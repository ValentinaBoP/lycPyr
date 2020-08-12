listPlots = function(chr2Plot, pi2plot, ranges, categories, pattern, GC2plot){
  
  #listPlots = list()
  chr2Plot$y = as.integer(chr2Plot$y)
  
  # order the categories according to the order in chr2Plot$y
  dfORDER = as.data.frame(unique(chr2Plot[,c(7,8)]))
  o = match(x = categories, table = dfORDER$category)
  v = (nrow(dfORDER)+1)-dfORDER[o,1]
  o = (nrow(dfORDER)+1)-dfORDER[,1]
  dfORDER = dfORDER[order(dfORDER$y),]
  
  #categories = factor(x = categories, levels = categories[v])
  categories = factor(x = categories, levels = dfORDER$category)
  categories = categories[!is.na(categories)]
  
  chr2Plot$fill = paste(chr2Plot$category, chr2Plot$orientation, sep = '-')
  chr2Plot$fill = as.factor(chr2Plot$fill)
  colori = sapply(strsplit(x = levels(chr2Plot$fill), split = '-'), "[[", 2)
  
  temp_levels = levels(categories)
  categories = as.character(categories)
  categories = c(categories[1], "Repeat density", categories[2:length(categories)])
  temp_levels[length(temp_levels)] = "Repeat density"
  temp_levels = c(temp_levels, "Genome")
  categories = factor(x = categories, levels = temp_levels)
  
  chr2Plot = correctOrientationTemp(track = chr2Plot, cats = categories)
  
  for(i in 1:(length(ranges) -1)){
    
    temp = chr2Plot[chr2Plot$start >= ranges[i] & chr2Plot$end <= ranges[i+1],]
    temp$fill = droplevels(temp$fill)
    temp$fill = as.character(temp$fill)
    #colori = sapply(strsplit(x = levels(temp$fill), split = '-'), "[[", 2)
    temp$orientation = sapply(strsplit(x = temp$fill, split = '-'), "[[", 2)
    temp$orientation[temp$orientation == "orange" & temp$strand == "-"] = "purple"
    
    temp_RM_density = addRMDensity(temp)
    temp_RM_density = temp_RM_density[temp_RM_density$start >= ranges[i] & temp_RM_density$start <= ranges[i+1],]
    temp$y[temp$category == "Genome"] = temp$y[temp$category == "Genome"] + 1.5
    #temp$y[temp$category == "Genome"] = temp$y[temp$category == "Genome"] + 0.5
    #temp_RM_density$y = temp$y[temp$category == "RepeatMasker"][1] + temp_RM_density$density + 1
    temp_RM_density$y = temp$y[temp$category == "RepeatMasker"][1] + 1.5
    temp_RM_density = adjustRMTempDensity(temp_RM_density)
    
    #temp = correctOrientationTemp(track = temp, cats = categories)
    
    temp_GC_content = GC2plot[GC2plot$start >= ranges[i] & GC2plot$start <= ranges[i+1],]
    temp_GC_content = correctGCcontent(temp_GC_content)
    temp$y[temp$category == "Genome"] = temp$y[temp$category == "Genome"] + 0.5
    temp_GC_content$y = temp_RM_density$y[1] + 0.5
    
    
    
    temp$fill = paste(temp$category, temp$orientation, sep = "-")
    temp$fill = factor(temp$fill)
    colori = sapply(strsplit(x = levels(temp$fill), split = '-'), "[[", 2)
    
    temp = temp[temp$category != "RepeatMasker",]
    categories = categories[-3]
    categories = droplevels(categories)
    temp_RM_density$ymin = temp$y[temp$category == "IL"][1] + 0.5
    temp_RM_density$y = temp_RM_density$ymin + 0.5
    temp_GC_content$y = temp_RM_density$y + 1
    temp$y[temp$category == "Genome"] = temp_GC_content$y[1] + 0.5
    
    blocks = createBlocks(temp, categories)
    annotation = createAnnotation(blocks)
    #plot = ggplot(data = chr2Plot[chr2Plot$start > ranges[i] & chr2Plot$end <= ranges[i+1],]) + geom_rect(data = chr2Plot[chr2Plot$start > ranges[i] & chr2Plot$end <= ranges[i+1],], mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = cols)) + scale_y_continuous(labels = unique(chr2Plot[chr2Plot$start > ranges[i] & chr2Plot$end <= ranges[i+1],]$category)[length(cat):1], breaks = 1:length(cat)) + theme(legend.position="none") + scale_fill_manual(values = as.character(unique(chr2Plot[chr2Plot$start > ranges[i] & chr2Plot$end <= ranges[i+1],]$cols)))
    #plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = cols)) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme(legend.position="none") + scale_fill_manual(values = temp$cols)
    #plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5)) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme(legend.position="none") + scale_color_manual(values = temp$cols)
    #plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill), alpha = 0.3)) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme_bw() + theme(legend.position="none", axis.text = element_text(size = 30), panel.grid = element_blank()) + scale_fill_manual(values = colori) + geom_line(data = temp_RM_density, aes(x = start, y = y))
    breaks = 1:length(levels(categories))
    #breaks[length(breaks)] = breaks[length(breaks)] + 1.5
    #plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill), alpha = 0.3)) + scale_y_continuous(labels = levels(categories), breaks = breaks) + theme_bw() + theme(legend.position="none", axis.text = element_text(size = 30), panel.grid = element_blank()) + scale_fill_manual(values = colori) + geom_line(data = temp_RM_density, aes(x = start, y = y))
    
    plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill), alpha = 0.3)) + scale_y_continuous(labels = levels(categories), breaks = breaks) + theme_bw() + theme(legend.position="none", axis.text = element_text(size = 30), panel.grid = element_blank()) + scale_fill_manual(values = colori)
    plot = plot + geom_rect(data = temp_RM_density, mapping = aes(xmin = start, xmax = end, ymin = ymin, ymax = y), fill = temp_RM_density$colors)
    plot = plot + geom_rect(data = temp_GC_content, mapping = aes(xmin = start, xmax = end, ymin = (y - 0.3), ymax = y), fill = temp_GC_content$colors)
    plot = plot + geom_rect(data = blocks, mapping = aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), color = rep_len(x = c("black", "gray60"), length.out = nrow(blocks)), fill = NA)
    
    # add annotation layer
    plot = plot + annotate(geom = "text", x = annotation$x_coord, y = annotation$y_coord, label = annotation$query, size = 1)
    
    # add gap blocks
    plot + geom_rect(data = gaps, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y + 0.5, fill = factor(cols))) + scale_fill_manual(values = gaps_colors)
    
    #if(length(pi2plot) > 0){
      
    #  temp_pi = pi2plot[pi2plot$start >= ranges[i] & pi2plot$end <= ranges[i+1],]
    #  
    #  gapLines = addGapLines(data = temp, pi = temp_pi)
    #  
    #  plot = plot + geom_point(data = temp_pi, aes(x = start, y = newPi)) + geom_segment(data = gapLines, aes(x = x, xend = xend, y = y, yend = yend))
    #  
    #}
    #plot = ggplot(data = temp) + geom_rect(data = temp, mapping = aes(xmin = start, xmax = end, ymin = y, ymax = y+0.5, fill = factor(fill))) + scale_y_continuous(labels = levels(categories)[length(levels(categories)):1], breaks = 1:length(levels(categories))) + theme(legend.position="none", axis.text = element_text(size = 50))
    ggsave(filename = paste(pattern, "_plot_", i, ".pdf", sep = ''), plot = plot, device = "pdf", width = 155, height = 20, units = "cm", limitsize = FALSE, dpi = 300)
    
  }
  
}

createAnnotation = function(blocks){
  
  boo = which(table(blocks$query) > 1)
  duplicatedNames = names(boo)
  
  annotation = blocks[blocks$query %in% duplicatedNames,]
  annotation$x_coord = (annotation$xmax + annotation$xmin)/2
  annotation$y_coord = (annotation$ymax + annotation$ymin)/2
  
  # simplify names 10X
  annotation$query = sapply(strsplit(x = annotation$query, split = ","), "[[", 1)
  
  annotation[, 2:7] = NULL
  
  return(annotation)
}

addGapTrack = function(data){
  
  temp_gaps = data[data$query == "Genome",]
  
  gaps = temp_gaps
  gaps$start = temp_gaps$end
  gaps$end = c(temp_gaps$start[2:length(temp_gaps$start)],temp_gaps$end[nrow(temp_gaps)])
  gaps$start = gaps$start + 1
  gaps$end = gaps$end - 1
  gaps = gaps[-(nrow(gaps)),]
  
  gaps$width = gaps$end - gaps$start + 1
  
  # PGA gaps
  gaps[gaps$width == 5000, 9] = "darkgreen"
  gaps[gaps$width == 5000, 10] = "darkgreen"
  # Dovetail gaps
  gaps[gaps$width == 100, 9] = "darkblue"
  gaps[gaps$width == 100, 10] = "darkblue"
  # Arcs gaps
  gaps[gaps$width == 10, 9] = "gold"
  gaps[gaps$width == 10, 10] = "gold"
  
  gaps = gaps[!(gaps$cols == "black"),]
  
  data = rbind(data, gaps)
  data = data[order(data[,1], data[,2]),]
  
  
  return(data)

}
