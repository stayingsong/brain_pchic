# Implement bedtools in R.
# Inputs: The two bed files and all options.
# Outputs: The resulting intersection file.
bedTools.2in <- function(functionstring="bedtools intersect", bed1, bed2, opt.string="") {
  
  # Create temp files.
  a.file = tempfile()
  b.file = tempfile()
  out.file = tempfile()
  
  # Write BED data to temp files.
  write.table(bed1, file=a.file, quote=F, sep="\t", col.names=F, row.names=F)
  write.table(bed2, file=b.file, quote=F, sep="\t", col.names=F, row.names=F)
  
  # Create and run bedtools command.
  command = paste(functionstring, "-a", a.file, "-b", b.file, opt.string, ">", out.file, sep=" ")
  try(system(command))
  
  # Read in resulting table.
  result = read.table(out.file, header=F, sep="\t")
  unlink(a.file)
  unlink(b.file)
  unlink(out.file)
  
  # Return results.
  return(result)
  
}

# Sorts BED file according to lexographical order.
# Inputs:
# Outputs:
sortBed <- function(bed) {
  
  # Define chromosome order.
  chrs <- paste0("chr", c(1:22, "X", "Y"))
  
  #
  sorted.bed <- c()
  for (chr in chrs) {
    
    #
    subset.bed <- bed[bed[, 1] == chr, ]
    sorted.bed <- rbind(sorted.bed, subset.bed[order(subset.bed[, 2]), ])
    
  }
  
  #
  return(sorted.bed)
  
}

#
# Inputs:
# Outputs:
bedTools.2in.compact <- function(functionstring="bedtools intersect", bed1, bed2, opt.string="") {
  
  # Create temp files.
  a.file = tempfile()
  b.file = tempfile()
  out.file = tempfile()
  
  # Write BED data to temp files.
  write.table(bed1, file=a.file, quote=F, sep="\t", col.names=F, row.names=F)
  write.table(bed2, file=b.file, quote=F, sep="\t", col.names=F, row.names=F)
  
  # Create and run bedtools command.
  command = paste(functionstring, "-a", a.file, "-b", b.file, opt.string, ">", out.file, sep=" ")
  try(system(command))
  
  # Read in resulting table.
  result = read.table(out.file, header=F, sep="\t", col.names=paste0("V", 1:8))
  unlink(a.file)
  unlink(b.file)
  unlink(out.file)
  
  # Return results.
  return(result)
  
}

# Given two sets of interactions, checks how many interactions overlap between the two sets.
# Inputs:
# Outputs:
checkInteractionOverlap <- function(interaction.set.1, interaction.set.2) {
  
  #
  loci1 <- unique(interaction.set.1[, c(1:3, 5:7)])
  loci2 <- unique(interaction.set.2[, c(1:3, 5:7)])
  combined.loci <- rbind(loci1, loci2)
  unique.loci <- unique(combined.loci)
  
  #
  print(paste0("Interaction set 1 has ", length(loci1[, 1]), " unique interaction loci."))
  print(paste0("Interaction set 2 has ", length(loci2[, 1]), " unique interaction loci."))
  print(paste0("The two interaction sets have ", length(unique.loci[, 1]), " total unique loci."))
  print(paste0("The two interaction sets have ", (length(loci1[, 1]) + length(loci2[, 1]) - length(unique.loci[, 1])), " common loci."))
  print(paste0("Interaction set 1 has ", (-length(loci2[, 1]) + length(unique.loci[, 1])), " specific loci."))
  print(paste0("Interaction set 2 has ", (-length(loci1[, 1]) + length(unique.loci[, 1])), " specific loci."))
  
}

#
# Inputs:
# Outputs:
convertToWashUFormat <- function(interactions.con, cell.type, center=F) {
  
  #
  field1 <- cbind(interactions.con[[cell.type]]$bait_chr, interactions.con[[cell.type]]$bait_start, interactions.con[[cell.type]]$bait_end)
  field1_old <- data.frame(field1, stringsAsFactors=F)
  field1_old[, 1] <- as.character(field1_old[, 1])
  field1_old[, 2] <- as.numeric(field1_old[, 2])
  field1_old[, 3] <- as.numeric(field1_old[, 3])
  for (i in 1:length(field1_old[, 1])) {
    field1[i, 2] <- floor((field1_old[i, 2] + field1_old[i, 3])/2)
    field1[i, 3] <- floor((field1_old[i, 2] + field1_old[i, 3])/2) + 1
  }
  field1[, 2] <- as.character(field1[, 2])
  field1[, 3] <- as.character(field1[, 3])
  field1_collapsed <- apply(field1, 1, paste0, collapse=",")
  
  #
  field2 <- cbind(interactions.con[[cell.type]]$otherEnd_chr, interactions.con[[cell.type]]$otherEnd_start, interactions.con[[cell.type]]$otherEnd_end)
  field2_old <- data.frame(field2, stringsAsFactors=F)
  field2_old[, 1] <- as.character(field2_old[, 1])
  field2_old[, 2] <- as.numeric(field2_old[, 2])
  field2_old[, 3] <- as.numeric(field2_old[, 3])
  for (i in 1:length(field2_old[, 1])) {
    field2[i, 2] <- floor((field2_old[i, 2] + field2_old[i, 3])/2)
    field2[i, 3] <- floor((field2_old[i, 2] + field2_old[i, 3])/2) + 1
  }
  field2[, 2] <- as.character(field2[, 2])
  field2[, 3] <- as.character(field2[, 3])
  field2_collapsed <- apply(field2, 1, paste0, collapse=",")
  
  #
  field3 <- interactions.con[[cell.type]]$score
  
  #
  return(cbind(field1_collapsed, field2_collapsed, field3))
  
}

#
# Inputs:
# Outputs:
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), labels, title='') {
  
  scale = (length(lut)-1)/(max-min)
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title, cex.main=0.85)
  axis(2, ticks, labels=labels, las=2)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
  
}

# Filter upper and lower percentile of values to remove outliers.
# Inputs:
# Outputs:
filterDistribution <- function(vector, percent) {
  
  #
  high_cutoff <- quantile(vector, 1-percent)
  low_cutoff <- quantile(vector, percent)
  
  #
  return(vector[(vector > low_cutoff) & (vector < high_cutoff)])
  
}

