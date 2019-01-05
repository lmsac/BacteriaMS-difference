peptide.mass = function(sequence, monoisotopic = FALSE) {
  aa.residues.mass.monoisotopic = c(
    A =	71.037114,
    R =	156.101111,
    N =	114.042927,
    D =	115.026943,
    C =	103.009185,
    E =	129.042593,
    Q =	128.058578,
    G =	57.021464,
    H =	137.058912,
    I =	113.084064,
    L =	113.084064,
    K =	128.094963,
    M =	131.040485,
    F =	147.068414,
    P =	97.052764,
    S =	87.032028,
    T =	101.047679,
    # U =	150.95363,
    W =	186.079313,
    Y =	163.06332,
    V =	99.068414
  )
  aa.residues.mass.average = c(
    A =	71.0779,
    R =	156.1857,
    N =	114.1026,
    D =	115.0874,
    C =	103.1429,
    E =	129.114,
    Q =	128.1292,
    G =	57.0513,
    H =	137.1393,
    I =	113.1576,
    L =	113.1576,
    K =	128.1723,
    M =	131.1961,
    F =	147.1739,
    P =	97.1152,
    S =	87.0773,
    T =	101.1039,
    U =	150.0379,
    W =	186.2099,
    Y =	163.1733,
    V =	99.1311
  )
  n.terminus.mass.monoisotopic = c(
    H = 1.007825
  )
  n.terminus.mass.average = c(
    H = 1.008
  )
  c.terminus.mass.monoisotopic = c(
    OH = 17.00274
  )
  c.terminus.mass.average = c(
    OH = 17.007
  )
  
  if(monoisotopic) 
    aa.residues.mass = aa.residues.mass.monoisotopic
  else
    aa.residues.mass = aa.residues.mass.average
  if(monoisotopic) 
    n.terminus.mass = n.terminus.mass.monoisotopic
  else
    n.terminus.mass = n.terminus.mass.average
  c.terminus.mass = if(monoisotopic) 
    c.terminus.mass = c.terminus.mass.monoisotopic
  else
    c.terminus.mass = c.terminus.mass.average
  
  result = sapply(strsplit(sequence, ''), function(aa.seq) {
    undefined.aa = aa.seq[!(aa.seq %in% names(aa.residues.mass))]
    if (length(undefined.aa) > 0) {
      stop(do.call(paste, as.list(c('undefined residue: ', undefined.aa))))
    }
    as.numeric(sum(aa.residues.mass[aa.seq]) + n.terminus.mass['H'] + c.terminus.mass['OH'])
  })
  names(result) = NULL
  result
}

match.peaks.with.sequence = function(mz, sequence, 
                                     tolerance = 2000, charge = 1:3, 
                                     monoisotopic = FALSE,
                                     iaa = FALSE,
                                     include.n.terminal.methionine.removal = FALSE) {
  mw = sapply(sequence, function(seq) {
    mass = peptide.mass(seq, monoisotopic = monoisotopic)
    # mass = Peptides::mw(seq, monoisotopic = monoisotopic)
    if(iaa) {
      carbamidomethyl = if(monoisotopic)
        57.021464
      else
        57.0513
      mass = mass + carbamidomethyl * sum(strsplit(seq, '')[[1]] == 'C')
    }
    mass
  })
  lapply(mz, function(z) {
    mz.range = c((1 - tolerance * 1e-6) * z, (1 + tolerance * 1e-6) * z)
    
    do.call(c, lapply(charge, function(ch) {
      result = {
        mass.range = mz.range * ch - ch
        matched.indexs = which(mw <= mass.range[2] &
                                 mw >= mass.range[1])
        lapply(matched.indexs, function(i) {
          list(
            index = i,
            theoretical.mz = (mw[i] + ch) / ch,
            mz = z,
            charge = ch
          )
        })
      }
      
      if(include.n.terminal.methionine.removal) {
        mw.methionine = 131.1986
        result = c(result, {
          mass.range = mz.range * ch - ch + mw.methionine
          matched.indexs = which(mw <= mass.range[2] &
                                   mw >= mass.range[1] &
                                   sapply(substring(sequence, 1, 1), function(char) char == 'M' || char == 'm'))
          lapply(matched.indexs, function(i) {
            list(index = i,
                 theoretical.mz = (mw[i] - mw.methionine + ch) / ch,
                 mz = z,
                 charge = ch,
                 n.terminal.methionine.removal = TRUE)
          })
        })
      }
      result
    }))
  })
}
