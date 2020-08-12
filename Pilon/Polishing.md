### BASH ###
cat 6.Arcs_test6/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds 5.Arcs_test5/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds 4.Arcs_test4/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds 3.Arcs_test3/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 3.Arcs_test3/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds 4.Arcs_test4/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 4.Arcs_test4/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds 5.Arcs_test5/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 5.Arcs_test5/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds 6.Arcs_test6/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 6.Arcs_test6/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds 7.Arcs_test7/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 7.Arcs_test7/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds 8.Arcs_test8/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds 8.Arcs_test8/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | grep '_' | sort -u > scaffolds_tests

grep '_' 3.Arcs_test3/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test31
grep '_' 3.Arcs_test3/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test32

grep '_' 4.Arcs_test4/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test41
grep '_' 4.Arcs_test4/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test42
grep '_' 4.Arcs_test4/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test45

grep '_' 5.Arcs_test5/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test51
grep '_' 5.Arcs_test5/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test52
grep '_' 5.Arcs_test5/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test55

grep '_' 6.Arcs_test6/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test61
grep '_' 6.Arcs_test6/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test62
grep '_' 6.Arcs_test6/5.LINKS_par5/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test65

grep '_' 7.Arcs_test7/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test71
grep '_' 7.Arcs_test7/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test72

grep '_' 8.Arcs_test8/1.LINKS_default/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test81
grep '_' 8.Arcs_test8/2.LINKS_par2/lycPyr2.2.1.fasta.scaffolds | sort -u > scaffolds_test82


### R ###
summary = read.table('scaffolds_tests', stringsAsFactors = F, header = F)
names(summary) = 'Scaffolds'
summary$ID = scaffold2ID(summary$Scaffolds)
boo = duplicated(summary$ID)

summary = summary[!boo,]
tests = list.files(pattern = 'scaffolds_test')

tests = tests[1:12]

for(i in 1:length(tests)){
  
  new = read.table(file = tests[i], stringsAsFactors = F, header = F)
  new$ID = scaffold2ID(new[,1])
  summary$new = '-'
  summary = createSummary(summary, new, col = (2 + i))
  names(summary)[(2 + i)] = sapply(strsplit(tests[i], split = '_'), '[[', 2)
  
}

summary = summary[order(summary[,2]),]

write.table(x = summary, file = 'summary_tests.stats', sep = '\t', quote = F, row.names = F, col.names = T)

###R###
disentagleArcsNames = function(DATA){
  
  NEW = data.frame(Scaffold = sapply(strsplit(x = DATA[,1], split = ','), '[[', 1))
  
  m = gregexpr(pattern = 'f[0-9]+z|f[0-9]+Z|r[0-9]+z|f[0-9]+Z', text = DATA$V1)
  
  contigs = regmatches(x = DATA$V1, m = m)
  
  for(i in 1:length(contigs)){
    
    contigs[[i]] = gsub(pattern = 'f|z|r', replacement = '', x = paste(contigs[[i]], collapse = ','), ignore.case = TRUE)
    
  }
  
  NEW$Contigs = unlist(contigs)
  
  NEW$Arcs = sub(pattern = 'scaffold', replacement = 'arcs', x = NEW$Scaffold)
  
}

# x data.frame con piu' scaffolds
# y data.frame da controllare
# col numero della colonna da riempire di x con dati presenza/assenza
#
createSummary = function(x, y, col){
  
  reg = sapply(strsplit(x$ID, '_'), '[[', 1)
  reg = paste('_', reg, '_|^', reg, '_', sep = '')
  o = integer()
  
  for(i in 1:nrow(x)){
    o = c(o, which(grepl(pattern = reg[i], x = y$ID)))
    if(length(which(grepl(pattern = reg[i], x = y$ID))) == 0){
      o = c(o, 0)
    } else {
      next
    }
  }
  
  for(i in 1:nrow(x)){
    if(o[i] != 0){
      x[i,col] = y$ID[o[i]]
    } else {
      next
    }
  }
  
  return(x)
  
}




