find.common.peaks.multi = function(peaklists, tolerance, method = "complete") {
  peaklists = lapply(peaklists, function(peaklist) {
    if(class(peaklist) == 'data.frame')
      as.matrix(peaklist)
    else
      peaklist
  })
  
  all.peaks = do.call(rbind, lapply(1:length(peaklists), function(i) {
    cbind(peaklists[[i]][, 1], peaklists[[i]][, 2], i)
  }))
  
  distance = apply(all.peaks, 1, function(x) {
    apply(all.peaks, 1, function(y) {
      if (x[3] == y[3])
        return(1)
      d = abs(x[1] - y[1]) / (x[1] + y[1]) * 2
      # d = abs(x[1] - y[1]) / max(x[1], y[1])
    })
  })
  
  hc = hclust(as.dist(distance), method)
  cluster.indexs = cutree(hc, h = tolerance * 1e-6)
  
  peak.clusters = lapply(1:max(cluster.indexs), function(i) {
    x = all.peaks[which(cluster.indexs == i), ]
    if(is.vector(x))
      x = matrix(x, nrow = 1)
    x
  })
  
  return(list(peak.clusters = peak.clusters, peaklists = peaklists, tolerance = tolerance))
}

align.peaks.multi = function(common.peaks.multi) {
  peak.clusters = common.peaks.multi$peak.clusters
  peaklists = common.peaks.multi$peaklists
  
  mzs = do.call(rbind, lapply(peak.clusters, function(cluster) {
    x = rep(NA, length = length(peaklists))
    if (is.vector(cluster))
      x[cluster[3]] = cluster[1]
    else
      apply(cluster, 1, function(peak) {
        x[peak[3]] <<- peak[1]
      })
    return(x)
  }))
  colnames(mzs) = unlist(lapply(1:length(peaklists), function(i) {
    paste0('mz', i)
  }))
  
  intensities = do.call(rbind, lapply(peak.clusters, function(cluster) {
    x = rep(0, length = length(peaklists))
    if (is.vector(cluster))
      x[cluster[3]] = cluster[2]
    else
      apply(cluster, 1, function(peak) {
        x[peak[3]] <<- peak[2]
      })
    return(x)
  }))
  colnames(intensities) = unlist(lapply(1:length(peaklists), function(i) {
    paste0('int', i)
  }))
  
  mz.mean = apply(cbind(mzs, intensities), 1, function(row) {
    mz = row[1:(length(row) / 2)]
    inten = row[(length(row) / 2 + 1):length(row)]
    sum(mz * inten, na.rm = TRUE) / sum(inten)
  })
  peak.order = order(mz.mean)
  mzs = mzs[peak.order, ]
  intensities = intensities[peak.order, ]
  
  return(cbind(mzs, intensities))
}

combine.peaklists = function(peaklists, tolerance) {
  aligned.peaks = align.peaks.multi(find.common.peaks.multi(peaklists, tolerance))
  
  mzs = aligned.peaks[, 1:length(peaklists)]
  intensities = aligned.peaks[, (length(peaklists) + 1):(length(peaklists) * 2)]
  
  mz.mean = apply(mzs, 1, function(row) {
    mean(row, na.rm = TRUE)
  })
  mz.sd = apply(mzs, 1, function(row) {
    s = sd(row, na.rm = TRUE)
    ifelse(is.na(s), 0, s)
  })
  
  intensities[is.na(mzs)] = NA
  
  i.mean = apply(intensities, 1, function(row) {
    mean(row, na.rm = TRUE)
  })
  i.sd = apply(intensities, 1, function(row) {
    s = sd(row, na.rm = TRUE)
    ifelse(is.na(s), 0, s)
  })
  
  peak.percentage = apply(mzs, 1, function(row) {
    sum(!is.na(row)) / length(row)
  })
  
  return(cbind(mz.mean = mz.mean, i.mean = i.mean, mz.sd = mz.sd, i.sd = i.sd, peak.percentage = peak.percentage))
}