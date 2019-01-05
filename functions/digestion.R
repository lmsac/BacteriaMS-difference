trypsin = function(sequence, max.miss = 0) {
  miss.number = 0:max.miss

  k.r.index = gregexpr('(?<=[KR])(?!P)', sequence, perl = T, ignore.case = TRUE)[[1]]
  
  do.call(c, lapply(miss.number[which(miss.number < length(k.r.index) - 1)], function(miss) {
    indexs = if(miss > 0)
      cbind(c(1, k.r.index[1:(length(k.r.index) - miss)]), 
            c(k.r.index[-(1:miss)] - 1, nchar(sequence)))
    else
      cbind(c(1, k.r.index), 
            c(k.r.index - 1, nchar(sequence)))
    
    lapply(1:nrow(indexs), function(i) {
      start = indexs[i, 1]
      end = indexs[i, 2]
      peptide = substr(sequence, start, end)
      list(sequence = peptide,
           start = start,
           end = end,
           miss = miss)
    })
  }))
}

digest.proteins = function(sequence, digestion = trypsin, ...) {
  peptides = lapply(sequence, function(seq) {
    digestion(seq, ...)
  })
  
  peptides = do.call(c, lapply(1:length(peptides), function(i) {
    s = peptides[[i]]
    lapply(s[which(!duplicated(sapply(s, function(p)
      p$sequence)))], function(p) {
        list(sequence = p$sequence,
             from = list(list(
               protein = i, miss = p$miss
             )))
      })
  }))
  
  peptides.sequence = sapply(peptides, function(p) p$sequence)
  duplicated.indexs = which(duplicated(peptides.sequence))
  if(length(duplicated.indexs) > 0) {
    lapply(duplicated.indexs, function(i) {
      first.index = which(peptides.sequence == peptides[[i]]$sequence)[1]
      peptides[[first.index]]$from <<-
        append(peptides[[first.index]]$from, peptides[[i]]$from)
    })
    peptides = peptides[-duplicated.indexs]
  }
  peptides
}