# normalize.intensity.internal.standard = function(
#   data, tolerance = 2000, 
#   internal.standard.mz = 12383, 
#   internal.standard.intensity = 100
# ) {
#   internal.standard.index = which(data[, 1] >= internal.standard.mz * (1 - tolerance * 1e-6) &
#                                     data[, 1] <= internal.standard.mz * (1 + tolerance * 1e-6))
#   if(length(internal.standard.index) == 0)
#     stop('internal standard peak not found.')
#   internal.standard.index = internal.standard.index[which.min(abs(
#     data[internal.standard.index, 1] - internal.standard.mz))]
#   
#   normalized.intensity = data[, 2] / data[internal.standard.index, 2] * internal.standard.intensity
#   normalized.intensity[normalized.intensity < 0] = 0
#   new.data = data
#   new.data[, 2] = normalized.intensity
#   new.data
# }

# peaklists.1 = lapply(list.files(pattern = '.txt'), function(file) {
#   normalize.intensity.internal.standard(read.table(file, header = T))
# })
# peaklists.2 = lapply(list.files(pattern = '.txt'), function(file) {
#   normalize.intensity.internal.standard(read.table(file, header = T))
# })

normalize.area = function(peaklist) {
  cbind(peaklist[, 1], peaklist[, 'Area'] / sum(peaklist[, 'Area']) * 100)
}

peaklists.1 = lapply(list.files(pattern = '.txt'), function(file) {
  normalize.area(read.table(file, header = T))
})
peaklists.2 = lapply(list.files(pattern = '.txt'), function(file) {
  normalize.area(read.table(file, header = T))
})

aligned.peaks = align.peaks.multi(find.common.peaks.multi(c(peaklists.1, peaklists.2), tolerance = 2000))

mz = apply(aligned.peaks[, 1:(length(peaklists.1) + length(peaklists.2))], 1, function(x) {
  mean(x, na.rm = T)
})
intensity = aligned.peaks[, (length(peaklists.1) + length(peaklists.2) + 1):
                            (length(peaklists.1) * 2 + length(peaklists.2) * 2)]

# write.csv(cbind(mz = mz, intensity), 'peak_intensity_table.csv', row.names = F)

write.csv(cbind(mz = mz, intensity), 'peak_area_table.csv', row.names = F)
