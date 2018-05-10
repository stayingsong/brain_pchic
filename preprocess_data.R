# Retain only cis interactions.
# Inputs:
# Outputs:
removeTransInteractions <- function(interactions) {
  
  # Return interactions where the lhs and rhs chromosomes match.
  return(interactions[interactions[, 1] == interactions[, 5], ])
  
}

# 
# Inputs:
# Outputs:
removeLowResInteractions <- function(interactions, min.bin.size) {
  
  # Return interactions where either the lhs or rhs bins are above a certain size.
  lhs.size <- abs(interactions[, 2] - interactions[, 3])
  rhs.size <- abs(interactions[, 6] - interactions[, 7])
  return(interactions[(lhs.size <= min.bin.size) & (rhs.size <= min.bin.size), ])
  
}

#
# Inputs: Individual data frames containing cell type specific expression data (needs mean_tpm field).
# Outputs: Takes the mean TPM for each gene for each cell type and summarizes it in a new data frame.
summarizeExpressionResults <- function(expression.results, gene.lengths=NULL, type) {
  
  #
  cell.types <- names(expression.results)
  expression.summary <- data.frame(matrix(NA, length(expression.results[[cell.types[1]]][, 1]), length(cell.types)), stringsAsFactors=F)
  colnames(expression.summary) <- cell.types
  rownames(expression.summary) <- rownames(expression.results[[cell.types[1]]])
  if (type == "TPM") {
    
    #
    for (i in 1:length(cell.types)) {
      
      expression.summary[, i] <- rowMeans(expression.results[[cell.types[i]]])
      
    }
    
  } else if (type == "RPKM") {
    
    #
    all.results <- rep(NA, length(expression.results[[cell.types[1]]]))
    all.lengths <- rep(NA, length(expression.results[[cell.types[1]]]))
    groups <- c()
    for (i in 1:length(cell.types)) {
      
      all.results <- cbind(all.results, expression.results[[cell.types[i]]])
      all.lengths <- cbind(all.lengths, gene.lengths[[cell.types[i]]])
      groups <- c(groups, rep(cell.types[i], length(expression.results[[cell.types[i]]])))
      
    }
    all.results <- all.results[, -1]
    all.lengths <- all.lengths[, -1]
    rownames(all.results) <- rownames(expression.results[[cell.types[1]]])
    
    #
    test <- (all.lengths == 0) & (all.results == 0)
    all.lengths[(all.lengths == 0) & (all.results == 0)] <- 1
    
    #
    y <- DGEList(counts=all.results, group=groups)
    y <- calcNormFactors(y, method="TMM")
    
    #
    all.rpkms <- rep(NA, length(expression.results[[cell.types[1]]]))
    for (i in 1:length(all.results)) {
      RPKM <- rpkm(y, gene.length=all.lengths[, i])
      all.rpkms <- cbind(all.rpkms, RPKM[, i])
    }
    all.rpkms <- all.rpkms[, -1]
    
    #
    for (i in 1:length(cell.types)) {
      
      expression.summary[, i] <- rowMeans(as.matrix(all.rpkms[, which(groups==cell.types[i])]))
      
    }
    
  }
  
  # Return expression summary results.
  return(expression.summary)
  
}

#
# Inputs:
# Outputs:
filterInteractionsBySignificance <- function(interactions, score.cutoff) {
  
  # Apply a score cutoff.
  filtered.interactions <- interactions[interactions$score >= score.cutoff, ]
  
  # Returns filtered interactions.
  return(filtered.interactions)
  
}

#
# Inputs:
# Outputs:
getNonSignificantInteractions <- function(interactions, low.score, high.score, downsampling.depth) {
  
  #
  interactions <- interactions[(interactions$score > low.score) & (interactions$score < high.score), ]
  
  #
  if (length(interactions[, 1]) >= downsampling.depth) {
    sampled.indices <- sample(length(interactions[, 1]), replace=F)[1:downsampling.depth]
  }
  else {
    sample.indices <- 1:length(interactions[, 1])
  }
  
  #
  filtered.interactions <- interactions[sampled.indices, ]
  
  # Returns filtered interactions.
  return(filtered.interactions)
  
}

# Downsamples interactions to specified depth.
# Inputs: Interactions and the # of reads to downsample to.
# Outputs: Downsampled interactions (with the highest scores).
downsampleInteractions <- function(interactions, downsampling.depth) {
  
  # Make sure there are enough interactions to downsample at the desired depth.
  if (length(interactions[, 1]) < downsampling.depth) {
    num.keep <- length(interactions[, 1])
  } else {
    num.keep <- downsampling.depth
  }
  
  # Sort by score, then take the top N interactions.
  temp <- interactions[order(-interactions$score), ]
  interactions.downsampled <- temp[1:num.keep, ]
  
  # Return downsampled interactions in same format as input.
  return(interactions.downsampled)
  
}

# Expands interactions to a minimum resolution.
# Inputs: Interactions and the minimum resolution to expand them to.  Can also consider neighboring fragments but not implemented for now.
# Outputs: Interactions with expanded lhs and rhs bins.
expandInteractions <- function(interactions, interaction.resolution, neighboring.fragments=FALSE) {
  
  # Extend interaction coordinates such that they cover minimum width defined by interaction resolution, unless they already exceed it.
  interactions.expanded <- interactions
  lhs.widths <- abs(interactions.expanded[, 2] - interactions.expanded[, 3])
  lhs.centers <- round((interactions.expanded[, 2] + interactions.expanded[, 3])/2)
  rhs.widths <- abs(interactions.expanded[, 6] - interactions.expanded[, 7])
  rhs.centers <- round((interactions.expanded[, 6] + interactions.expanded[, 7])/2)
  interactions.expanded[lhs.widths < interaction.resolution, 2] <- lhs.centers[lhs.widths < interaction.resolution] - interaction.resolution/2
  interactions.expanded[lhs.widths < interaction.resolution, 3] <- lhs.centers[lhs.widths < interaction.resolution] + interaction.resolution/2
  interactions.expanded[rhs.widths < interaction.resolution, 6] <- rhs.centers[rhs.widths < interaction.resolution] - interaction.resolution/2
  interactions.expanded[rhs.widths < interaction.resolution, 7] <- rhs.centers[rhs.widths < interaction.resolution] + interaction.resolution/2
  
  # Fix negative coordinates. 
  interactions.expanded[interactions.expanded[, 2] < 1, 2] <- 1
  interactions.expanded[interactions.expanded[, 6] < 1, 6] <- 1
  
  # Return expanded interactions.
  return(interactions.expanded)
  
}

# Expands peaks to a minimum width.
# Inputs: Peaks and the minimum width to expand them to.
# Outputs: Peaks with expanded bins (same format as input).
expandATACPeaks <- function(peaks, min.width) {
  
  # Extend peak coordinates such that they cover specified minimum width, unless they already exceed it.
  peaks.expanded <- peaks
  widths <- abs(peaks.expanded[, 2] - peaks.expanded[, 3])
  centers <- round((peaks.expanded[, 2] + peaks.expanded[, 3])/2)
  peaks.expanded[widths < min.width, 2] <- centers[widths < min.width] - min.width/2
  peaks.expanded[widths < min.width, 3] <- centers[widths < min.width] + min.width/2
  
  # Fix negative coordinates. 
  peaks.expanded[peaks.expanded[, 2] < 1, 2] <- 1
  peaks.expanded[peaks.expanded[, 6] < 1, 6] <- 1
  
  # Return expanded peaks.
  return(peaks.expanded)
  
}

