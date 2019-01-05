fasta.file = 'uniprot-Klebsiella_pneumoniae-reviewed.fasta'

proteins = local({
  proteins = read.fasta(fasta.file)
  proteins.sequence = sapply(proteins, function(entry) entry$sequence)
  valid.indexs = grep('[^ARNDCEQGHILKMFPSTWYV]', proteins.sequence, invert = TRUE)
  proteins[valid.indexs]
})

proteins.sequence = sapply(proteins, function(entry) entry$sequence)


undigested.tolerance = 2000
undigested.peaklists = lapply(list.files(pattern = '.txt'), function(file) {
  read.table(file, header = T)
})
undigested.combined.peaklist = combine.peaklists(undigested.peaklists, tolerance = undigested.tolerance)

undigested.match = match.peaks.with.sequence(
  # undigested.combined.peaklist[undigested.combined.peaklist[, 5] >= 1, 1], 
  undigested.combined.peaklist[, 1],
  proteins.sequence, 
  tolerance = undigested.tolerance,
  monoisotopic = FALSE,
  iaa = FALSE,
  include.n.terminal.methionine.removal = TRUE
)

undigested.match.result = do.call(rbind, lapply(undigested.match, function(matches) {
  do.call(rbind, lapply(matches, function(match) {
    cbind(
      mz = match$mz,
      theoretical.mz = match$theoretical.mz,
      charge = match$charge,
      n.terminal.methionine.removal = if (is.null(match$n.terminal.methionine.removal))
        FALSE
      else
        match$n.terminal.methionine.removal,
      t(unlist(proteins[[match$index]]))
    )
  }))
}))

write.csv(undigested.match.result, 'match_undigested.csv', row.names = F)


undigested.match.proteins = proteins[unique(unlist(sapply(undigested.match, function(matches) {
  sapply(matches, function(match) match$index)
})))]
undigested.match.proteins.sequence = sapply(undigested.match.proteins, function(entry) entry$sequence)

peptides = digest.proteins(undigested.match.proteins.sequence, digestion = trypsin, max.miss = 1)
peptides.sequence = sapply(peptides, function(p) p$sequence)


digested.tolerance = 500
digested.peaklists = lapply(list.files(pattern = '.txt'), function(file) {
  read.table(file, header = T)
})
digested.combined.peaklist = combine.peaklists(digested.peaklists, tolerance = digested.tolerance)

digested.match = match.peaks.with.sequence(
  # digested.combined.peaklist[digested.combined.peaklist[, 5] >= 1, 1], 
  digested.combined.peaklist[, 1], 
  peptides.sequence,
  tolerance = digested.tolerance,
  charge = 1,
  monoisotopic = TRUE,
  iaa = TRUE,
  include.n.terminal.methionine.removal = FALSE
)

digested.match.result = do.call(rbind, lapply(digested.match, function(matches) {
  do.call(rbind, lapply(matches, function(match) {
    cbind(
      mz = match$mz,
      theoretical.mz = match$theoretical.mz,
      charge = match$charge,
      n.terminal.methionine.removal = if (is.null(match$n.terminal.methionine.removal))
        FALSE
      else
        match$n.terminal.methionine.removal,
      sequence = peptides[[match$index]]$sequence,
      protein = do.call(paste, lapply(undigested.match.proteins, function(entry) entry$name)[sapply(peptides[[match$index]]$from, function(p) p$protein)]),
      miss = peptides[[match$index]]$from[[1]]$miss
    )
  }))
}))

write.csv(digested.match.result, 'match_digested.csv', row.names = F)

protein.match.result = local({
  protein.match = unlist(strsplit(as.character(digested.match.result[, 'protein']), ' '))
  protein.match.count = table(protein.match)
  result = undigested.match.result[as.character(undigested.match.result[, 'name']) %in% protein.match, ]
  result = cbind(result[, 1:5], peptide.count = protein.match.count[result[, 'name']], result[, 6:7])
})

write.csv(protein.match.result, 'match_proteins.csv', row.names = F)
