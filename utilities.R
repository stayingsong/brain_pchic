# Implement bedtools intersect command in R.
# Notes: Set 'num.cols' equal to the # of expected columns in the output file.
#
bedTools.2in <- function(functionstring="bedtools intersect", bed1, bed2, opt.string="", num.cols=NULL) {
  
  # Create temporary files for use with bedtools intersect.
  a.file <- tempfile()
  b.file <- tempfile()
  out.file <- tempfile()
  
  # Write BED input data to temporary files.
  write.table(bed1, file=a.file, sep="\t", row.names=F, col.names=F, quote=F)
  write.table(bed2, file=b.file, sep="\t", row.names=F, col.names=F, quote=F)
  
  # Run bedtools intersect command.
  command <- paste(functionstring, "-a", a.file, "-b", b.file, opt.string, ">", out.file, sep=" ")
  try(system(command))
  
  # Read in output.
  if (!is.null(num.cols)) {
    result <- read.table(out.file, header=F, sep="\t", stringsAsFactors=F, col.names=paste0("V", 1:num.cols))
  } else {
    result <- read.table(out.file, header=F, sep="\t", stringsAsFactors=F)
  }
  
  # Clean up temporary files.
  unlink(a.file)
  unlink(b.file)
  unlink(out.file)
  
  # Return intersection results.
  return(result)
  
}

# Sort BED intervals according to lexographical order.
# Notes:
#
sortBed <- function(bed) {
  
  # Define order of chromosomes.
  chrs <- paste0("chr", c(1:22, "X", "Y"))
  
  # Compile sorted BED intervals chromosome by chromosome.
  sorted.bed <- c()
  for (chr in chrs) {
    
    # Sort by start coordinate.
    subset.bed <- bed[bed[, 1] == chr, ]
    sorted.bed <- rbind(sorted.bed, subset.bed[order(subset.bed[, 2], subset.bed[, 3]), ])
    
  }
  
  # Return sorted BED intervals.
  return(sorted.bed)
  
}

# Sort interactions according to lexographical order of the lhs bin first, then the rhs bin.
# Notes: Assumes all interactions are cis.
#
sortInteractions <- function(bed) {
  
  # Define order of chromosomes.
  chrs <- paste0("chr", c(1:22, "X", "Y"))
  
  # Compile sorted BED intervals chromosome by chromosome.
  sorted.bed <- c()
  for (chr in chrs) {
    
    # Sort by start coordinate.
    subset.bed <- bed[bed[, 1] == chr, ]
    sorted.bed <- rbind(sorted.bed, subset.bed[order(subset.bed[, 2], subset.bed[, 3], subset.bed[, 6], subset.bed[, 7]), ])
    
  }
  
  # Return sorted BED intervals.
  return(sorted.bed)
  
}

# Converts interactions to a format compatible with the WashU Epigenome Browser.
# Notes: Set 'center' to TRUE to plot interaction arcs between the centers of bins.
#
convertToWashUFormat <- function(interactions, center=T) {
  
  # Process the lhs bin (first field).
  field1 <- cbind(interactions$bait_chr, interactions$bait_start, interactions$bait_end)
  field1.old <- data.frame(field1, stringsAsFactors=F)
  field1.old[, 1] <- as.character(field1.old[, 1])
  field1.old[, 2] <- as.numeric(field1.old[, 2])
  field1.old[, 3] <- as.numeric(field1.old[, 3])
  if (center == T) {
    for (i in 1:length(field1.old[, 1])) {
      field1[i, 2] <- floor((field1.old[i, 2] + field1.old[i, 3])/2)
      field1[i, 3] <- floor((field1.old[i, 2] + field1.old[i, 3])/2) + 1
    }
  }
  field1[, 2] <- as.character(field1[, 2])
  field1[, 3] <- as.character(field1[, 3])
  field1.collapsed <- apply(field1, 1, paste0, collapse=",")
  
  # Process the rhs bin (second field).
  field2 <- cbind(interactions$otherEnd_chr, interactions$otherEnd_start, interactions$otherEnd_end)
  field2.old <- data.frame(field2, stringsAsFactors=F)
  field2.old[, 1] <- as.character(field2.old[, 1])
  field2.old[, 2] <- as.numeric(field2.old[, 2])
  field2.old[, 3] <- as.numeric(field2.old[, 3])
  if (center == T) {
    for (i in 1:length(field2.old[, 1])) {
      field2[i, 2] <- floor((field2.old[i, 2] + field2.old[i, 3])/2)
      field2[i, 3] <- floor((field2.old[i, 2] + field2.old[i, 3])/2) + 1
    }
  }
  field2[, 2] <- as.character(field2[, 2])
  field2[, 3] <- as.character(field2[, 3])
  field2.collapsed <- apply(field2, 1, paste0, collapse=",")
  
  # Process the interaction score (third field).
  field3 <- interactions$score
  
  # Stitch together and return formatted interactions.
  return(cbind(field1.collapsed, field2.collapsed, field3))
  
}

# Converts interactions from the old WashU format to longrange format.
# Notes:
#
updateWashUFormat <- function(interactions) {
  
  # Convert each interaction in the old format to two entries in the new format.
  updated <- c()
  for (i in 1:length(interactions[, 1])) {

    lhs <- unlist(strsplit(interactions[i, 1], split=","))
    rhs <- unlist(strsplit(interactions[i, 2], split=","))
    updated <- rbind(updated, c(lhs, paste0(rhs[1], ":", rhs[2], "-", rhs[3], ",", interactions[i, 3])))
    updated <- rbind(updated, c(rhs, paste0(lhs[1], ":", lhs[2], "-", lhs[3], ",", interactions[i, 3])))
    
  }
  
  # Return the interactions in longrange format.
  return(updated)
  
}

# Update rsids that have been merged in dbSNP using the RsMergeArch table.
# Notes:
#
update.rsid <- function(rsid, rs.merge) {
  
  # Continue looking for rsid translations until no more can be found.
  repeat {
    
    # Translate the rsid using the RsMergeArch table.
    hit <- rs.merge[rs.merge[, 1] == rsid, 2]
    
    # If a hit is found, update the rsid and repeat. Otherwise, break.
    if (length(hit) == 1) {
      rsid <- hit
    } else if (length(hit) == 0) {
      break
    }
    
  }
  
  # Return the updated rsid.
  return(rsid)
  
}


# End ---------------------------------------------------------------------

