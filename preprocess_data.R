# Filter out trans interactions between different chromosomes.
# Notes:
#
removeTransInteractions <- function(interactions) {
  
  # Return only interactions where the lhs and rhs chromosomes match.
  return(interactions[interactions[, 1] == interactions[, 5], ])
  
}

# Filter out interactions where either end is not on a list of allowed chromosomes.
# Notes:
#
filterInteractionsByChromosome <- function(interactions, allowed.chrs) {
  
  # Return interactions where both ends are on the list of allowed chromosomes.
  return(interactions[(interactions[, 1] %in% allowed.chrs) & (interactions[, 5] %in% allowed.chrs), ])
  
}

# Reorder interactions based on the positions of their lhs and rhs bins.
# Notes: Some interaction calling pipelines such as CHiCAGO report interactions where the rhs bin is to the left of the lhs bin.
# 
reorderInteractionBins <- function(interactions) {
  
  # Find interactions where the rhs bin is to the left of the lhs bin.
  incorrect.start <- interactions$bait_start > interactions$otherEnd_start

  # Flip the lhs and rhs bins of interactions that need to be reordered.
  temp.start <- interactions[incorrect.start, 1:4]
  interactions[incorrect.start, 1:4] <- interactions[incorrect.start, 5:8]
  interactions[incorrect.start, 5:8] <- temp.start
  
  # Return reordered interactions.
  return(interactions)
  
}

# Remove interactions where the lhs and rhs bins are the same. This can occur in CHiCAGO when both ends are baited fragments.
# Notes: Assigns new interaction IDs (assumes interactions have 10 columns).
#
deduplicateInteractions <- function(interactions, cell.type) {

  # Find out which interactions have duplicated coordinates.
  is.duplicate <- duplicated(interactions[, 1:8]) | duplicated(interactions[, 1:8], fromLast=T)
  duplicates <- interactions[is.duplicate, 1:10]
  singletons <- interactions[!is.duplicate, 1:10]

  # Collapse duplicate interactions by the first 8 columns, summarizing the last 2 columns.
  reduced <- duplicates %>% group_by_(.dots=c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name")) %>% 
    summarize(N_reads=sum(N_reads), score=max(score))

  # Reformat collapsed interactions.
  reduced <- data.frame(reduced)
  reduced[, 1] <- as.character(reduced[, 1])
  reduced[, 2] <- as.integer(reduced[, 2])
  reduced[, 3] <- as.integer(reduced[, 3])
  reduced[, 4] <- as.character(reduced[, 4])
  reduced[, 5] <- as.character(reduced[, 5])
  reduced[, 6] <- as.integer(reduced[, 6])
  reduced[, 7] <- as.integer(reduced[, 7])
  reduced[, 8] <- as.character(reduced[, 8])
  reduced[, 9] <- as.integer(reduced[, 9])
  reduced[, 10] <- as.numeric(reduced[, 10])
  colnames(reduced) <- colnames(interactions)[1:10]
  
  # Recombine the collapsed interactions with the singleton interactions and assign interaction IDs.
  interactions <- rbind(singletons, reduced)
  interactions <- cbind(interactions, paste0(cell.type, "_", 1:length(interactions[, 1])))
  interactions[, 11] <- as.character(interactions[, 11])
  colnames(interactions)[11] <- "interaction_ID"

  # Return deduplicated interactions.
  return(interactions)
  
}

# Remove interactions where the lhs or rhs bin is larger than the minimum size (low resolution interactions).
# Notes:
#
removeLowResInteractions <- function(interactions, min.bin.size) {
  
  # Return interactions where the lhs or rhs bin is not above the minimum size.
  lhs.size <- abs(interactions[, 2] - interactions[, 3])
  rhs.size <- abs(interactions[, 6] - interactions[, 7])
  return(interactions[(lhs.size <= min.bin.size) & (rhs.size <= min.bin.size), ])
  
}

# Filter interactions by a specified score cutoff.
# Notes:
#
filterInteractionsBySignificance <- function(interactions, score.cutoff) {
  
  # Return interactions with scores equal to or above the specified score cutoff.
  filtered.interactions <- interactions[interactions$score >= score.cutoff, ]
  return(filtered.interactions)
  
}

# Expand lhs and rhs bins to the minimum resolution.
# Notes: 'neighboring.fragments' option is not implemented yet.
#
expandInteractions <- function(interactions, interaction.resolution, neighboring.fragments=F) {
  
  # Calculate interaction bin widths and centers.
  lhs.widths <- abs(interactions[, 2] - interactions[, 3])
  lhs.centers <- round((interactions[, 2] + interactions[, 3])/2)
  rhs.widths <- abs(interactions[, 6] - interactions[, 7])
  rhs.centers <- round((interactions[, 6] + interactions[, 7])/2)
  
  # For bins whose widths are smaller than the minimum resolution, expand them from the center to match the minimum resolution.
  interactions[lhs.widths < interaction.resolution, 2] <- lhs.centers[lhs.widths < interaction.resolution] - interaction.resolution/2
  interactions[lhs.widths < interaction.resolution, 3] <- lhs.centers[lhs.widths < interaction.resolution] + interaction.resolution/2
  interactions[rhs.widths < interaction.resolution, 6] <- rhs.centers[rhs.widths < interaction.resolution] - interaction.resolution/2
  interactions[rhs.widths < interaction.resolution, 7] <- rhs.centers[rhs.widths < interaction.resolution] + interaction.resolution/2
  
  # Fix negative coordinates. 
  interactions[interactions[, 2] < 1, 2] <- 1
  interactions[interactions[, 6] < 1, 6] <- 1
  
  # Return interaction with expanded bins.
  return(interactions)
  
}

# Expand features to a specified minimum width.
# Notes:
#
expandFeatures <- function(features, min.width) {
  
  # Calculate feature interval widths and centers.
  widths <- abs(features[, 2] - features[, 3])
  centers <- round((features[, 2] + features[, 3])/2)
  
  # Expand feature intervals such that they cover the specified minimum width, unless they already exceed it.
  features[widths < min.width, 2] <- centers[widths < min.width] - min.width/2
  features[widths < min.width, 3] <- centers[widths < min.width] + min.width/2
  
  # Fix negative coordinates. 
  features[features[, 2] < 1, 2] <- 1

  # Return expanded features.
  return(features)
  
}

# Expand features by a specified width.
# Notes:
#
expandFeaturesByWidth <- function(features, by.width) {
  
  # Expand feature intervals by the specified width.
  features[, 2] <- features[, 2] - by.width
  features[, 3] <- features[, 3] + by.width
  
  # Fix negative coordinates. 
  features[features[, 2] < 1, 2] <- 1
  
  # Return expanded features.
  return(features)
  
}

# Calculate mean TPM or TMM-normalized FPKM expression values for each gene and cell type.
# Notes: RNA-seq results for each cell type should have the same # and ordering of genes.
#
summarizeExpressionResults <- function(expression.results, gene.lengths=NULL, type) {
  
  # Create results matrix where rows are genes and columns are cell types.
  cell.types <- names(expression.results)
  num.genes <- length(expression.results[[cell.types[1]]][, 1])
  expression.summary <- data.frame(matrix(NA, num.genes, length(cell.types)), stringsAsFactors=F)
  colnames(expression.summary) <- cell.types
  rownames(expression.summary) <- rownames(expression.results[[cell.types[1]]])
  
  # Calculate the specified metric for gene expression.
  if (type == "TPM") {
    
    # Calculate the mean TPM across all replicates.
    for (i in 1:length(cell.types))
      expression.summary[, i] <- rowMeans(expression.results[[cell.types[i]]])
    
  } else if (type == "RPKM") {
    
    # Count number of replicates across all cell types.
    num.replicates <- 0
    for (cell.type in cell.types)
      num.replicates <- num.replicates + length(names(expression.results[[cell.type]]))
    
    # Compile matrices containing gene lengths and expected counts for each replicate across all cell types.
    all.results <- data.frame(matrix(NA, num.genes, num.replicates))
    all.lengths <- data.frame(matrix(NA, num.genes, num.replicates))
    rownames(all.results) <- rownames(expression.results[[cell.types[1]]])
    rownames(all.lengths) <- rownames(expression.results[[cell.types[1]]])
    groups <- c()
    col <- 1
    for (cell.type in cell.types) {
      
      all.results[, col:(col + length(names(expression.results[[cell.type]])) - 1)] <- expression.results[[cell.type]]
      all.lengths[, col:(col + length(names(expression.results[[cell.type]])) - 1)] <- gene.lengths[[cell.type]]
      colnames(all.results)[col:(col + length(names(expression.results[[cell.type]])) - 1)] <- names(expression.results[[cell.type]])
      colnames(all.lengths)[col:(col + length(names(expression.results[[cell.type]])) - 1)] <- names(expression.results[[cell.type]])
      groups <- c(groups, rep(cell.type, length(expression.results[[cell.type]])))
      col <- col + length(names(expression.results[[cell.type]]))
      
    }
    
    # Calculate TMM normalization factors for each replicate using edgeR.
    y <- DGEList(counts=all.results, group=groups)
    y <- calcNormFactors(y, method="TMM")
    
    # Calculate the TMM-normalized RPKM for each replicate (accounting for its reported gene lengths) using edgeR.
    all.rpkms <- data.frame(matrix(NA, num.genes, num.replicates))
    for (i in 1:length(all.results)) {
      
      RPKM <- rpkm(y, gene.length=all.lengths[, i])
      all.rpkms[, i] <- RPKM[, i]
      
    }

    # Calculate the mean TMM-normalized RPKM across all replicates.
    for (i in 1:length(cell.types))
      expression.summary[, i] <- rowMeans(as.matrix(all.rpkms[, which(groups==cell.types[i])]))
    
  }
  
  # Return summarized expression results.
  return(expression.summary)
  
}

# Randomly sample a specified number of interactions with scores in a specified range.
# Notes:
#
sampleInteractionsByScore <- function(interactions, low.score, high.score, sampling.depth) {
  
  # Retain only interactions with scores in the specified range.
  interactions <- interactions[(interactions$score >= low.score) & (interactions$score <= high.score), ]
  
  # If there are fewer interactions in the specified range of scores, just use all of them.
  if (length(interactions[, 1]) >= sampling.depth) {
    sampled.indices <- sample(length(interactions[, 1]), replace=F)[1:sampling.depth]
  }
  else {
    print(paste0("Not enough interactions in the specified range of scores, using only ", length(interactions[, 1]), 
                 " interactions (tried to sample ", sampling.depth, " interactions)"))
    sampled.indices <- 1:length(interactions[, 1])
  }
  
  # Return sampled interactions.
  return(interactions[sampled.indices, ])
  
}

# Downsample interactions to a specified depth.
# Notes:
#
downsampleInteractions <- function(interactions, downsampling.depth) {
  
  # Determine how man interactions to keep. If there are fewer interactions than the desired depth, keep all of them.
  if (length(interactions[, 1]) < downsampling.depth) {
    num.keep <- length(interactions[, 1])
  } else {
    num.keep <- downsampling.depth
  }
  
  # First sort interactions by score, then retain the top 'num.keep' interactions.
  interactions <- interactions[order(-interactions$score), ]
  interactions.downsampled <- interactions[1:num.keep, ]
  
  # Return downsampled interactions.
  return(interactions.downsampled)
  
}
