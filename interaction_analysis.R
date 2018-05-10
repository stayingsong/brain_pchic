# Analyzes number of interactions per gene (for genes with multiple promoters, assume they are all functioning together).
# Inputs: Expression results need to be arranged by Ensembl gene ID.
# Outputs:
analyzePromoterHubs <- function(interaction.data, interaction.annotations, promoters, expression.summary) {
  
  # Create a data frame for each cell type which has a row for each unique gene ID and columns for gene ID, number of interactions, expression value, and gene type.
  hub.results <- list()
  gene.ids <- unique(promoters$gene_id)
  cell.types <- names(interaction.data)
  for (cell.type in cell.types) {
    
    print(paste0("Analyzing promoter hubs for ", length(gene.ids), " genes for cell type: ", cell.type))
    
    # Get promoter to promoter and promoter to other interactions (lists of IDs) for current cell type.
    annotation <- interaction.annotations[[cell.type]]
    Pr_Pr <- as.character(interaction.data[[cell.type]][((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
                                                          ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)), 11])
    Pr_NP <- as.character(interaction.data[[cell.type]][(((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
                                                           (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)) |
                                                          ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) &
                                                             ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0))), 11])
    
    # For each promoter, determine how many interactions it overlaps with, and also retrieve the corresponding expression value.
    hub.results[[cell.type]] <- data.frame(matrix("", length(gene.ids), 9), stringsAsFactors=F)
    colnames(hub.results[[cell.type]]) <- c("gene_id", "gene_type", "expression_value", "num_interactions", "num_interactions_PA", "num_interactions_PP", "num_interactions_PO", "num_interactions_PE", "num_interactions_PR")
    hub.results[[cell.type]][, 1] <- as.character(hub.results[[cell.type]][, 1])
    hub.results[[cell.type]][, 2] <- as.character(hub.results[[cell.type]][, 2])
    hub.results[[cell.type]][, 3] <- as.numeric(hub.results[[cell.type]][, 3])
    hub.results[[cell.type]][, 4] <- as.integer(hub.results[[cell.type]][, 4])
    hub.results[[cell.type]][, 5] <- as.integer(hub.results[[cell.type]][, 5]) #PA
    hub.results[[cell.type]][, 6] <- as.integer(hub.results[[cell.type]][, 6]) #PP
    hub.results[[cell.type]][, 7] <- as.integer(hub.results[[cell.type]][, 7]) #PO
    hub.results[[cell.type]][, 8] <- as.integer(hub.results[[cell.type]][, 8]) #PE
    hub.results[[cell.type]][, 9] <- as.integer(hub.results[[cell.type]][, 9]) #PR
    hub.results[[cell.type]][, 1] <- gene.ids
    interaction.loci <- rbind(interaction.data[[cell.type]][, c(1:3, 11)], setNames(interaction.data[[cell.type]][, c(5:7, 11)], names(interaction.data[[cell.type]][, c(1:3, 11)])))
    promoter.interaction.results <- bedTools.2in(bed1=promoters[, c(1:3, 5)], bed2=interaction.loci, opt.string="-wb")
    for (i in 1:length(gene.ids)) {
      
      if (i %% 5000 == 0)
        print(i)
      
      #
      subset.interactions <- unique(as.character(promoter.interaction.results[as.character(promoter.interaction.results[, 4]) == gene.ids[i], 8]))
      hub.results[[cell.type]][i, 2] <- unique(promoters$gene_type[promoters$gene_id == gene.ids[i]])
      if (gene.ids[i] %in% rownames(expression.summary)) {
        hub.results[[cell.type]][i, 3] <- expression.summary[gene.ids[i], cell.type]
      } else {
        hub.results[[cell.type]][i, 3] <- NA
      }
      hub.results[[cell.type]][i, 4] <- length(unique(subset.interactions))
      
      # Annotate interactions as needed for special count categories.
      if (length(subset.interactions) > 0) {
        
        # Count promoter to distal ATAC-seq peak interactions.
        annotation.indices <- which(interaction.data[[cell.type]]$ID %in% subset.interactions)
        annotation.lhs.promoter <- interaction.annotations[[cell.type]]$lhs_promoter[annotation.indices]
        annotation.lhs.promoter.ATAC.seq <- interaction.annotations[[cell.type]]$'lhs_promoter_ATAC-seq'[annotation.indices]
        annotation.lhs.distal.ATAC.seq <- interaction.annotations[[cell.type]]$'lhs_distal_ATAC-seq'[annotation.indices]
        annotation.rhs.promoter <- interaction.annotations[[cell.type]]$rhs_promoter[annotation.indices]
        annotation.rhs.promoter.ATAC.seq <- interaction.annotations[[cell.type]]$'rhs_promoter_ATAC-seq'[annotation.indices]
        annotation.rhs.distal.ATAC.seq <- interaction.annotations[[cell.type]]$'rhs_distal_ATAC-seq'[annotation.indices]
        subset.Pr_DA <- subset.interactions[((((annotation.lhs.promoter > 0) | (annotation.lhs.promoter.ATAC.seq > 0)) & (annotation.rhs.promoter == 0) & (annotation.rhs.promoter.ATAC.seq == 0) & (annotation.rhs.distal.ATAC.seq > 0)) |
                                               (((annotation.rhs.promoter > 0) | (annotation.rhs.promoter.ATAC.seq > 0)) & (annotation.lhs.promoter == 0) & (annotation.lhs.promoter.ATAC.seq == 0) & (annotation.lhs.distal.ATAC.seq > 0)))]
        hub.results[[cell.type]][i, 5] <- length(subset.Pr_DA)
        
        # Count promoter to promoter and promoter to other interactions.
        hub.results[[cell.type]][i, 6] <- length(subset.interactions[subset.interactions %in% Pr_Pr])
        hub.results[[cell.type]][i, 7] <- length(subset.interactions[subset.interactions %in% Pr_NP])
        
        # Count promoter to enhancer and promoter to repressor interactions.
        annotation.lhs.enhancer <- interaction.annotations[[cell.type]]$lhs_active_enhancer[annotation.indices] +
          interaction.annotations[[cell.type]]$lhs_other_enhancer[annotation.indices] +
          interaction.annotations[[cell.type]]$lhs_genic_enhancer[annotation.indices]
        annotation.lhs.repressive <- interaction.annotations[[cell.type]]$lhs_repressive_heterochromatin[annotation.indices] +
          interaction.annotations[[cell.type]]$lhs_repressive_polycomb[annotation.indices]
        annotation.rhs.enhancer <- interaction.annotations[[cell.type]]$rhs_active_enhancer[annotation.indices] +
          interaction.annotations[[cell.type]]$rhs_other_enhancer[annotation.indices] +
          interaction.annotations[[cell.type]]$rhs_genic_enhancer[annotation.indices]
        annotation.rhs.repressive <- interaction.annotations[[cell.type]]$rhs_repressive_heterochromatin[annotation.indices] +
          interaction.annotations[[cell.type]]$rhs_repressive_polycomb[annotation.indices]
        subset.Pr_PE <- subset.interactions[((((annotation.lhs.promoter > 0) | (annotation.lhs.promoter.ATAC.seq > 0)) & (annotation.rhs.promoter == 0) & (annotation.rhs.promoter.ATAC.seq == 0) & (annotation.rhs.enhancer > 0)) |
                                               (((annotation.rhs.promoter > 0) | (annotation.rhs.promoter.ATAC.seq > 0)) & (annotation.lhs.promoter == 0) & (annotation.lhs.promoter.ATAC.seq == 0) & (annotation.lhs.enhancer > 0)))]
        subset.Pr_PR <- subset.interactions[((((annotation.lhs.promoter > 0) | (annotation.lhs.promoter.ATAC.seq > 0)) & (annotation.rhs.promoter == 0) & (annotation.rhs.promoter.ATAC.seq == 0) & (annotation.rhs.repressive > 0)) |
                                               (((annotation.rhs.promoter > 0) | (annotation.rhs.promoter.ATAC.seq > 0)) & (annotation.lhs.promoter == 0) & (annotation.lhs.promoter.ATAC.seq == 0) & (annotation.lhs.repressive > 0)))]
        hub.results[[cell.type]][i, 8] <- length(subset.Pr_PE)
        hub.results[[cell.type]][i, 9] <- length(subset.Pr_PR)
        
      } else {
        
        hub.results[[cell.type]][i, 5] <- 0
        hub.results[[cell.type]][i, 6] <- 0
        hub.results[[cell.type]][i, 7] <- 0
        hub.results[[cell.type]][i, 8] <- 0
        hub.results[[cell.type]][i, 9] <- 0
        
      }
      
    }
    
  }
  
  #
  return(hub.results)
  
}

#
# Inputs: Interactions organized by cell type.
# Outputs: Two layer lists for each permutation of cell types, e.g. specificity[[cell.type.1]][[cell.type.2]], containing overlap matrix between all interactions.
getInteractionSpecificity <- function(interaction.data, overlap.threshold) {
  
  # Analyze each combination of cell types.
  cell.types <- names(interaction.data)
  specificity <- list()
  for (cell.type.1 in cell.types) {
    
    specificity[[cell.type.1]] <- list()
    for (cell.type.2 in cell.types) {
      
      # Create a matrix to hold specificity results for each cell type.
      specificity[[cell.type.1]][[cell.type.2]] <- matrix(FALSE, length(interaction.data[[cell.type.1]][, 1]), length(interaction.data[[cell.type.2]][, 1]))
      print(paste0("Overlapping ", length(interaction.data[[cell.type.1]][, 1]), " ", cell.type.1, " and ", 
                   length(interaction.data[[cell.type.2]][, 1]), " ", cell.type.2, " interactions..."))
      
      # Grab IDs and interactions for each cell type.
      ids.1 <- as.character(interaction.data[[cell.type.1]][, 11])
      ids.2 <- as.character(interaction.data[[cell.type.2]][, 11])
      lhs.interactions.1 <- interaction.data[[cell.type.1]][, c(1:3, 11)]
      lhs.interactions.2 <- interaction.data[[cell.type.2]][, c(1:3, 11)]
      rhs.interactions.1 <- interaction.data[[cell.type.1]][, c(5:7, 11)]
      rhs.interactions.2 <- interaction.data[[cell.type.2]][, c(5:7, 11)]
      
      # Intersect all lhs intersections and all rhs intersections between the two cell types.
      lhs.overlap.results <- bedTools.2in(bed1=lhs.interactions.1, bed2=lhs.interactions.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
      rhs.overlap.results <- bedTools.2in(bed1=rhs.interactions.1, bed2=rhs.interactions.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
      
      # Fill in overlap matrix for every interaction.
      for (k in 1:length(ids.1)) {
        
        # Find interactions in other cell type that intersect both ends of current interaction.
        lhs.intersections <- as.character(lhs.overlap.results[lhs.overlap.results[, 4] == ids.1[k], 8])
        rhs.intersections <- as.character(rhs.overlap.results[rhs.overlap.results[, 4] == ids.1[k], 8])
        both.intersections <- lhs.intersections[lhs.intersections %in% rhs.intersections]
        
        # Fill out table of overlaps.
        specificity[[cell.type.1]][[cell.type.2]][k, ids.2 %in% both.intersections] <- TRUE
        
      }
      
      # Add row and column names to make it easier to query table.
      rownames(specificity[[cell.type.1]][[cell.type.2]]) <- ids.1
      colnames(specificity[[cell.type.1]][[cell.type.2]]) <- ids.2
      
    }
    
  }
  
  # Return specificity results.
  return(specificity)
  
}

# Gets all combinations of cell types and creates a matrix for each cell type where rows are interactions and columns represent overlap with interactions in each cell type.
# Inputs:
# Outputs: 
getSpecificityCategories <- function(specificity.matrix) {
  
  # Process specificity categories for each cell type using the specificity matrix.
  cell.types <- names(specificity.matrix)
  specificity.categories <- list()
  for (cell.type in cell.types) {
    
    # Create results matrix for current cell type, and of course the interactions in the current cell type overlap themselves so set that column as TRUE.
    specificity.categories[[cell.type]] <- data.frame(matrix(FALSE, length(specificity.matrix[[cell.type]][[cell.type]][, 1]), length(cell.types)), stringsAsFactors=F)
    rownames(specificity.categories[[cell.type]]) <- rownames(specificity.matrix[[cell.type]][[cell.type]])
    colnames(specificity.categories[[cell.type]]) <- cell.types
    specificity.categories[[cell.type]][, cell.type] <- TRUE
    
    # Iterate over other cell types.
    other.cell.types <- cell.types[!(cell.types %in% cell.type)]
    for (other.cell.type in other.cell.types) {
      
      # For each interaction in the current cell type, see if it overlaps any interactions in the other cell type, then set the corresponding results value to TRUE/FALSE.
      current.matrix <- specificity.matrix[[cell.type]][[other.cell.type]]
      for (i in 1:length(current.matrix[, 1])) {
        specificity.categories[[cell.type]][i, other.cell.type] <- any(current.matrix[i, ])
      }
      
    }
    
  }
  
  # Return specificity categories.
  return(specificity.categories)
  
}

#
# Inputs:
# Outputs: Returns a character vector of interaction IDs for interactions that are specific to the cell types in specific.cell.types.
getSpecificityCategoryInteractionIDs <- function(specificity.categories, specific.cell.types) {
  
  #
  cell.types <- names(specificity.categories)
  pattern <- rep(FALSE, length(cell.types))
  pattern[cell.types %in% specific.cell.types] <- TRUE
  specific.cell.type.ids <- c()
  for (cell.type in cell.types) {
    
    if (!all(colnames(specificity.categories[[cell.type]]) == cell.types))
      print("Column mismatch error")
    match <- rep(TRUE, length(specificity.categories[[cell.type]][, 1]))
    for (i in 1:length(cell.types)) {
      match <- (match & specificity.categories[[cell.type]][, i] == pattern[i])
    }
    specific.cell.type.ids <- c(specific.cell.type.ids, rownames(specificity.categories[[cell.type]])[match])
    
  }
  
  #
  return(specific.cell.type.ids)
  
}


#
# Inputs: Interaction sets for each cell type (both significant and all), specificity category results, and plotting directory.
# Outputs: 
writeHeatmapFiles <- function(interaction.data, interaction.data.all, chunk.size, output.folder) {
  
  # For each cell type, build lists consisting of significant interactions for one cell type and all interactions for other two cell types to be intersected.
  cell.types <- names(interaction.data)
  for (cell.type in cell.types) {
    
    # 
    print(paste0("Processing cell type: ", cell.type))
    
    # Filter unnecessary columns and filter out interactions with score < 1.
    interaction.data[[cell.type]] <- interaction.data[[cell.type]][, c(1:3, 5:7, 10:11)]
    interaction.data.all[[cell.type]] <- interaction.data.all[[cell.type]][, c(1:3, 5:7, 10:11)]
    interaction.data.all[[cell.type]] <- interaction.data.all[[cell.type]][interaction.data.all[[cell.type]]$score > 0, ]
    
    # Create files for each chunk.
    num.chunks <- ceiling(length(interaction.data[[cell.type]][, 1])/chunk.size)
    print(paste0("Splitting ", length(interaction.data[[cell.type]][, 1]), " interactions into ", num.chunks, " chunks of chunk size ", chunk.size))
    for (i in 1:num.chunks) {
      
      # Get start and end indices
      start.index <- chunk.size*(i-1) + 1
      end.index <- chunk.size*i
      if (end.index > length(interaction.data[[cell.type]][, 1]))
        end.index <- length(interaction.data[[cell.type]][, 1])
      
      #
      current.subset <- interaction.data[[cell.type]][start.index:end.index, ]
      
      #
      print(paste0("Writing chunk from ", start.index, " to ", end.index))
      save(current.subset, file=paste0(output.folder, "/current.subset.", cell.type, ".", start.index, "-", end.index, ".Rdata"))
      
    }
    
  }
  
  #
  save(interaction.data, file=paste0(output.folder, "/interaction.data.Rdata"))
  save(interaction.data.all, file=paste0(output.folder, "/interaction.data.all.Rdata"))
  
}

#
# Inputs:
# Outputs:
processHeatmapFiles <- function(cell.type, start.index, end.index, overlap.threshold=10^-9, input.folder, output.folder) {
  
  #
  original.cell.type <- cell.type
  load(file=paste0(input.folder, "/current.subset.", cell.type, ".", start.index, "-", end.index, ".Rdata"))
  load(file=paste0(input.folder, "/interaction.data.Rdata"))
  load(file=paste0(input.folder, "/interaction.data.all.Rdata"))
  
  #
  print(paste0("Currently processing cell type: ", cell.type))
  print(paste0("Analyzing interactions: ", start.index, " to ", end.index))
  num.interactions <- length(current.subset[, 1])
  print(paste0(num.interactions, " interactions total to be compared"))
  
  #
  cell.types <- names(interaction.data)
  print(paste0("Significant interactions:"))
  for (cell.type in cell.types) {
    print(paste0("Comparing against ", length(interaction.data[[cell.type]][, 1]), " ", cell.type, " interactions"))
  }
  print(paste0("All interactions:"))
  for (cell.type in cell.types) {
    print(paste0("Comparing against ", length(interaction.data.all[[cell.type]][, 1]), " ", cell.type, " interactions"))
  }
  
  #
  results.matrix <- data.frame(matrix(NA, num.interactions, 4), stringsAsFactors=F) # This stores the mean score across all significant interactions (or all interactions if no significant interactions).
  colnames(results.matrix) <- c("ID", cell.types)
  results.matrix[, 1] <- as.character(results.matrix[, 1])
  results.matrix[, 2] <- as.numeric(results.matrix[, 2])
  results.matrix[, 3] <- as.numeric(results.matrix[, 3])
  results.matrix[, 4] <- as.numeric(results.matrix[, 4])
  results.matrix[, 1] <- current.subset[, 8]
  specificity.matrix <- data.frame(matrix(NA, num.interactions, 4), stringsAsFactors=F) # This stores cell type specificity of interaction, e.g. C, HA, or CHA.
  colnames(specificity.matrix) <- c("ID", cell.types)
  specificity.matrix[, 1] <- as.character(results.matrix[, 1])
  specificity.matrix[, 2] <- as.logical(results.matrix[, 2])
  specificity.matrix[, 3] <- as.logical(results.matrix[, 3])
  specificity.matrix[, 4] <- as.logical(results.matrix[, 4])
  specificity.matrix[, 1] <- current.subset[, 8]
  for (cell.type in cell.types) {
    
    #
    print(cell.type)
    
    # Get IDs and coordinates of lhs and rhs loci for each set of interactions.
    subset.ids <- as.character(current.subset[, 8])
    subset.lhs <- current.subset[, c(1:3, 8)]
    subset.rhs <- current.subset[, c(4:6, 8)]
    compare.sig.ids <- as.character(interaction.data[[cell.type]][, 8])
    compare.sig.lhs <- interaction.data[[cell.type]][, c(1:3, 8)]
    compare.sig.rhs <- interaction.data[[cell.type]][, c(4:6, 8)]
    compare.all.ids <- as.character(interaction.data.all[[cell.type]][, 8])
    compare.all.lhs <- interaction.data.all[[cell.type]][, c(1:3, 8)]
    compare.all.rhs <- interaction.data.all[[cell.type]][, c(4:6, 8)]
    
    # Process each interaction in the current chunk one by one.
    for (i in 1:num.interactions) {
      
      #
      print(i)
      
      # Intersect current interaction with significant interactions from this cell type.
      lhs.results <- bedTools.2in.compact(bed1=subset.lhs[i, ], bed2=compare.sig.lhs, opt.string=paste0('-wa -wb -f ', overlap.threshold))
      rhs.results <- bedTools.2in.compact(bed1=subset.rhs[i, ], bed2=compare.sig.rhs, opt.string=paste0('-wa -wb -f ', overlap.threshold))
      lhs.sig.intersections <- as.character(lhs.results[lhs.results[, 4] == subset.ids[i], 8])
      rhs.sig.intersections <- as.character(rhs.results[rhs.results[, 4] == subset.ids[i], 8])
      both.sig.intersections <- lhs.sig.intersections[lhs.sig.intersections %in% rhs.sig.intersections]
      
      # If there are any interactions, get the mean score and set specificity matrix entry to TRUE.
      if (length(both.sig.intersections) > 0) {
        
        specificity.matrix[[cell.type]][i] <- TRUE
        results.matrix[[cell.type]][i] <- mean(interaction.data[[cell.type]][interaction.data[[cell.type]][, 8] %in% both.sig.intersections, 7])
        
      } else if (length(both.sig.intersections) == 0) {
        
        # Set specificity matrix entry to FALSE.
        specificity.matrix[[cell.type]][i] <- FALSE
        
        # Intersect current interaction with all interactions from this cell type.
        lhs.results <- bedTools.2in.compact(bed1=subset.lhs[i, ], bed2=compare.all.lhs, opt.string=paste0('-wa -wb -f ', overlap.threshold))
        rhs.results <- bedTools.2in.compact(bed1=subset.rhs[i, ], bed2=compare.all.rhs, opt.string=paste0('-wa -wb -f ', overlap.threshold))
        lhs.all.intersections <- as.character(lhs.results[lhs.results[, 4] == subset.ids[i], 8])
        rhs.all.intersections <- as.character(rhs.results[rhs.results[, 4] == subset.ids[i], 8])
        both.all.intersections <- lhs.all.intersections[lhs.all.intersections %in% rhs.all.intersections]
        
        # If there are any interactions, get the mean score and fill in results table.
        if (length(both.all.intersections) > 0) {
          
          results.matrix[[cell.type]][i] <- mean(interaction.data.all[[cell.type]][interaction.data.all[[cell.type]][, 8] %in% both.all.intersections, 7])
          
        } else if (length(both.all.intersections) == 0) {
          
          results.matrix[[cell.type]][i] <- 0
          
        } else {
          print("Error")
        }
        
      } else {
        print("Error")
      }
      
    }
    
  }
  
  # Save results to file for loading later on and stitching back together.
  save(results.matrix, file=paste0(output.folder, "/results.matrix.", original.cell.type, ".", start.index, "-", end.index, ".Rdata"))
  save(specificity.matrix, file=paste0(output.folder, "/specificity.matrix.", original.cell.type, ".", start.index, "-", end.index, ".Rdata"))
  
}

#
# Inputs:
# Outputs:
plotHeatmapResults <- function(interactions.res, input.folder, output.prefix) {
  
  #
  setwd(input.folder)
  
  # 
  chunk.size <- 1000
  num.chunks <- c(28, 40, 56)
  start.lasts <- c(28001, 40001, 56001)
  end.lasts <- c(28362, 40445, 56811)
  
  #
  cell.types <- names(interactions.res)
  
  # Stitch together results for each cell type.
  results.matrix.all <- list()
  specificity.matrix.all <- list()
  for (i in 1:length(cell.types)) {

    #
    cell.type <- cell.types[i]
    print(paste0("Processing cell type: ", cell.type))

    #
    results.matrix.all[[cell.type]] <- data.frame(matrix(NA, end.lasts[i], 4), stringsAsFactors=F)
    specificity.matrix.all[[cell.type]] <- data.frame(matrix(NA, end.lasts[i], 4), stringsAsFactors=F)
    results.matrix.all[[cell.type]][, 1] <- as.character(results.matrix.all[[cell.type]][, 1])
    results.matrix.all[[cell.type]][, 2] <- as.numeric(results.matrix.all[[cell.type]][, 2])
    results.matrix.all[[cell.type]][, 3] <- as.numeric(results.matrix.all[[cell.type]][, 3])
    results.matrix.all[[cell.type]][, 4] <- as.numeric(results.matrix.all[[cell.type]][, 4])
    specificity.matrix.all[[cell.type]][, 1] <- as.character(specificity.matrix.all[[cell.type]][, 1])
    specificity.matrix.all[[cell.type]][, 2] <- as.numeric(specificity.matrix.all[[cell.type]][, 2])
    specificity.matrix.all[[cell.type]][, 3] <- as.numeric(specificity.matrix.all[[cell.type]][, 3])
    specificity.matrix.all[[cell.type]][, 4] <- as.numeric(specificity.matrix.all[[cell.type]][, 4])

    #
    for (j in 1:num.chunks[i]) {

      # Read in results for current chunk.
      start.index <- ((j-1)*chunk.size+1)
      end.index <- (j*chunk.size)
      print(paste0("Reading chunk from ", start.index, " to ", end.index, "..."))
      load(paste0("results.matrix.", cell.type, ".", start.index, "-", end.index, ".Rdata"))
      load(paste0("specificity.matrix.", cell.type, ".", start.index, "-", end.index, ".Rdata"))
      results.matrix[, 1] <- as.character(results.matrix[, 1])
      specificity.matrix[, 1] <- as.character(specificity.matrix[, 1])
      results.matrix.all[[cell.type]][start.index:end.index, ] <- results.matrix
      specificity.matrix.all[[cell.type]][start.index:end.index, ] <- specificity.matrix

    }

    #
    print(paste0("Reading last chunk..."))
    load(paste0("results.matrix.", cell.type, ".", start.lasts[i], "-", end.lasts[i], ".Rdata"))
    load(paste0("specificity.matrix.", cell.type, ".", start.lasts[i], "-", end.lasts[i], ".Rdata"))
    results.matrix[, 1] <- as.character(results.matrix[, 1])
    specificity.matrix[, 1] <- as.character(specificity.matrix[, 1])
    results.matrix.all[[cell.type]][start.lasts[i]:end.lasts[i], ] <- results.matrix
    specificity.matrix.all[[cell.type]][start.lasts[i]:end.lasts[i], ] <- specificity.matrix

    #
    print(table(is.na(as.vector(results.matrix.all[[cell.type]]))))
    print(table(is.na(as.vector(specificity.matrix.all[[cell.type]]))))

  }
  
  # Checkpoint.
  save(results.matrix.all, file="results.matrix.all.Rdata")
  save(specificity.matrix.all, file="specificity.matrix.all.Rdata")
  load(file="results.matrix.all.Rdata")
  load(file="specificity.matrix.all.Rdata")
  
  #
  for (cell.type in cell.types) {
    rownames(specificity.matrix.all[[cell.type]]) <- specificity.matrix.all[[cell.type]][, 1]
    specificity.matrix.all[[cell.type]] <- specificity.matrix.all[[cell.type]][, 2:4]
    colnames(specificity.matrix.all[[cell.type]]) <- cell.types
    print(table(specificity.matrix.all[[cell.type]] == specificity.categories[[cell.type]]))
  }
  specificity.categories <- specificity.matrix.all
  
  # Get combinations of cell types.
  combinations <- c()
  for (i in 1:length(cell.types)) {
    temp <- combn(cell.types, i)
    for (j in 1:length(temp[1, ])) {
      combinations <- c(combinations, paste0(temp[, j], collapse="&"))
    }
  }
  
  # Sort interactions from each cell type into each combination.
  counts <- rep(0, length(combinations))
  combination.loci <- list()
  combination.ids <- list()
  for (i in 1:length(combinations)) {
    
    # Get boolean pattern for current combination.
    current.types <- unlist(strsplit(combinations[i], split="&"))
    pattern <- rep(FALSE, length(cell.types))
    pattern[cell.types %in% current.types] <- TRUE
    
    #
    combination.loci[[combinations[i]]] <- c(-1)
    combination.ids[[combinations[i]]] <- c()
    for (cell.type in cell.types) {
      
      # Get the interactions for each cell type that match the current combination/pattern.
      if (!all(colnames(specificity.categories[[cell.type]]) == cell.types))
        print("Column mismatch error")
      match <- rep(TRUE, length(specificity.categories[[cell.type]][, 1]))
      for (j in 1:length(cell.types)) {
        match <- (match & specificity.categories[[cell.type]][, j] == pattern[j])
      }
      
      #
      interaction.ids <- rownames(specificity.categories[[cell.type]])[which(match)]
      length(interaction.ids)
      length(which(match))
      if (length(interaction.data[[cell.type]][interaction.data[[cell.type]]$ID %in% interaction.ids, c(1:3, 5:7)][, 1]) > 0) {
        combination.loci[[combinations[i]]] <- rbind(combination.loci[[combinations[i]]], interaction.data[[cell.type]][interaction.data[[cell.type]]$ID %in% interaction.ids, c(1:3, 5:7)])
      }
      combination.ids[[combinations[i]]] <- c(combination.ids[[combinations[i]]], interaction.ids)
      
    }
    
    #
    combination.loci[[combinations[i]]] <- combination.loci[[combinations[i]]][-1, ]
    counts[i] <- length(unique(combination.loci[[combinations[i]]])[, 1])
    combination.ids[[combinations[i]]] <- combination.ids[[combinations[i]]][!duplicated(combination.loci[[combinations[i]]])]
    
  }
  
  #
  my_palette <- colorRampPalette(c("#6495ED", "#FFFFFF", "#FFFFFF", "#FF3030", "#EE0000", "#B22222"))(n=599)
  col_breaks = c(seq(0, 2.5, length=100),
                 seq(2.51, 4.75, length=100),
                 seq(4.76, 5, length=100),
                 seq(5.01, 7.5, length=100),
                 seq(7.51, 15, length=100),
                 seq(15.01, 100, length=100))
  
  #
  pdf(paste0(output.prefix, "/heatmap_color_bar.pdf"), width=2, height=6) 
  color.bar(my_palette, min=0, max=100, nticks=6, labels=c(0, 2.5, 4.5, 5, 15, 100), title="CHiCAGO score")
  dev.off()
  
  # Now assemble separate matrices for each unique combination of specificities.
  master <- c()
  sections <- list()
  results.merged <- rbind(results.matrix.all[["cortical"]], results.matrix.all[["hippocampal"]], results.matrix.all[["astrocyte"]])
  colnames(results.merged) <- c("ID", cell.types)
  for (i in 1:length(combinations)) {
    
    #
    current.types <- unlist(strsplit(combinations[i], split="&"))
    pattern <- rep(FALSE, length(cell.types))
    pattern[cell.types %in% current.types] <- TRUE
    column.order <- c((1:length(cell.types))[which(pattern)], (1:length(cell.types))[which(!pattern)])
    
    #
    sections[[combinations[i]]] <- results.merged[results.merged[, 1] %in% combination.ids[[combinations[i]]], 2:4]
    head(sections[[combinations[i]]])
    colnames(sections[[combinations[i]]]) <- cell.types
    sections[[combinations[i]]] <- data.matrix(sections[[combinations[i]]])
    master <- rbind(master, sections[[combinations[i]]])
    print(dim(sections[[combinations[i]]]))
    
    #
    #downsampled.indices <- (1:floor(length(sections[[combinations[i]]][, 1])/2))*2-1
    #sections[[combinations[i]]] <- sections[[combinations[i]]][downsampled.indices, ]
    height.factor <- length(sections[[combinations[i]]][, 1])/length(results.merged[, 1])
    print(height.factor)
    print(i)
    print(paste0("entries in section: ", length(sections[[combinations[i]]][, 1])))
    jpeg(paste0("results/", combinations[i],".jpg"), width=6*300, height=6*300*height.factor*10, res=300, pointsize=8)
    heatmap.2(sections[[combinations[i]]],
              main=NULL,
              density.info="none",
              trace="none",
              col=my_palette,
              breaks=col_breaks,
              margins=c(2,2),
              dendrogram="none",
              Colv=NA,
              Rowv=NA,
              key=FALSE,
              labRow=FALSE,
              labCol=FALSE)
    dev.off()
    
  }
  
  # Prepare dendorgram of hclust results.
  pdf(paste0("interaction_specificity_dendrogram.pdf"), width=6, height=6)     
  scores <- t(results.merged[, 2:4])
  rownames(scores) <- c("hippocampal", "cortical", "astrocyte")
  h <- hclust(dist(scores))
  dd <- as.dendrogram(h)
  plot(rev(dd), main="interaction specificity dendrogram", axes=T, cex=1, cex.main=1)
  dev.off()
  
}

#
# Inputs: List with two cell types.
# Outputs: List of interactions from first cell type meeting specificity criteria.
getSpecificInteractions <- function(interaction.data, overlap.threshold, main.type) {
  
  # Specific implementation of bedtools in R.
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
    result = read.table(out.file, header=F, sep="\t", col.names=paste0("V", 1:8))
    unlink(a.file)
    unlink(b.file)
    unlink(out.file)
    
    # Return results.
    return(result)
    
  }
  
  # Get sets of interactions to be compared.  We will return main specific interactions.
  cell.types <- names(interaction.data)
  if (length(cell.types) !=  2) {
    print("Error, should only be two cell types for this function")
  }
  main.interactions <- interaction.data[[cell.types[cell.types == main.type]]]
  compare.interactions <- interaction.data[[cell.types[cell.types != main.type]]]
  
  # Get IDs and coordinates of lhs and rhs loci for each set of interactions.
  main.ids <- as.character(main.interactions[, 11])
  compare.ids <- as.character(compare.interactions[, 11])
  main.lhs.loci <- main.interactions[, c(1:3, 11)]
  compare.lhs.loci <- compare.interactions[, c(1:3, 11)]
  main.rhs.loci <- main.interactions[, c(5:7, 11)]
  compare.rhs.loci <- compare.interactions[, c(5:7, 11)]
  
  #
  main.specific <- c()
  main.common <- c()
  compare.common <- c()
  for (i in 1:length(main.interactions[, 1])) {
    
    #
    print(i)
    
    # Intersect all lhs intersections and all rhs intersections between the two cell types.
    lhs.overlap.results <- bedTools.2in(bed1=main.lhs.loci[i, ], bed2=compare.lhs.loci, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    rhs.overlap.results <- bedTools.2in(bed1=main.rhs.loci[i, ], bed2=compare.rhs.loci, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    lhs.intersections <- as.character(lhs.overlap.results[lhs.overlap.results[, 4] == main.ids[i], 8])
    rhs.intersections <- as.character(rhs.overlap.results[rhs.overlap.results[, 4] == main.ids[i], 8])
    both.intersections <- lhs.intersections[lhs.intersections %in% rhs.intersections]
    if (length(both.intersections) == 0) {
      main.specific <- c(main.specific, main.ids[i])
    } else if (length(both.intersections > 0)) {
      main.common <- c(main.common, main.ids[i])
      compare.common <- c(compare.common, both.intersections)
    }
    
  }
  
  #
  print(paste0(length(main.interactions[, 1]), " interactions for cell type: ", cell.types[cell.types == main.type]))
  print(paste0(length(compare.interactions[, 1]), " interactions for cell type: ", cell.types[cell.types != main.type]))
  print(paste0(length(main.specific), " interactions are specific to cell type: ", cell.types[cell.types == main.type]))
  print(paste0(length(main.common), " interactions are shared for cell type: ", cell.types[cell.types == main.type]))
  print(paste0(length(unique(compare.common)), " interactions are shared for cell type: ", cell.types[cell.types != main.type]))
  print(paste0(length(compare.interactions[, 1]) - length(unique(compare.common)), " interactions are specific to cell type: ", cell.types[cell.types != main.type]))
  
  #
  return(main.interactions[main.interactions$ID %in% main.common, ])
  
}

# Similar to analyzePromoterHubs() but focuses on individual interactions versus being gene promoter-centric.
# Inputs:
# Outputs:
analyzeInteractionHubs <- function(interaction.data, output.prefix) {
  
  # For each interaction, intersect its lhs and rhs with all other coordinates lhs and rhs and see how many other interactions come up besides itself.
  cell.types <- names(interaction.data)
  counts <- matrix(0, length(cell.types), 2)
  rownames(counts) <- cell.types
  colnames(counts) <- c("one to one", "one to many")
  for (i in 1:length(cell.types)) {
    
    # Intersect all lhs intersections and all rhs intersections between the two cell types.
    interaction.ids <- as.character(interaction.data[[cell.types[i]]][, 11])
    lhs.interactions <- interaction.data[[cell.types[i]]][, c(1:3, 11)]
    rhs.interactions <- interaction.data[[cell.types[i]]][, c(5:7, 11)]
    lhs.lhs.overlap.results <- bedTools.2in(bed1=lhs.interactions, bed2=lhs.interactions, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    lhs.rhs.overlap.results <- bedTools.2in(bed1=lhs.interactions, bed2=rhs.interactions, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    rhs.lhs.overlap.results <- bedTools.2in(bed1=rhs.interactions, bed2=lhs.interactions, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    rhs.rhs.overlap.results <- bedTools.2in(bed1=rhs.interactions, bed2=rhs.interactions, opt.string=paste0('-wa -wb -f ', overlap.threshold))
    
    # Process intersection results.
    for (j in 1:length(interaction.ids)) {
      
      lhs.lhs.intersections <- as.character(lhs.lhs.overlap.results[lhs.lhs.overlap.results[, 4] == interaction.ids[j], 8])
      lhs.rhs.intersections <- as.character(lhs.rhs.overlap.results[lhs.rhs.overlap.results[, 4] == interaction.ids[j], 8])
      rhs.lhs.intersections <- as.character(rhs.lhs.overlap.results[rhs.lhs.overlap.results[, 4] == interaction.ids[j], 8])
      rhs.rhs.intersections <- as.character(rhs.rhs.overlap.results[rhs.rhs.overlap.results[, 4] == interaction.ids[j], 8])
      other.ids <- unique(c(lhs.lhs.intersections, lhs.rhs.intersections, rhs.lhs.intersections, rhs.rhs.intersections))
      
      # Determine whether or not it just intersects with itself or with other intersections as well.
      if ((length(other.ids) == 1) && (other.ids == interaction.ids[j])) {
        counts[i, 1] <- counts[i, 1] + 1
      } else if (length(other.ids) > 1) {
        counts[i, 2] <- counts[i, 2] + 1
      } else {
        print(other.ids[j])
      }
      
    }
    
  }
  
  # Plot bar graph of results.
  png(filename=paste0(output.prefix, "/interaction_hubs_analysis_one_many_counts.png"), res=300, width=6, height=4, units="in")
  par(xpd=T, mar=par()$mar + c(0, 0, 0, 12))
  percentages <- apply(counts, 1, function(x) {x*100/sum(x)})
  rownames(percentages) <- c("one to one", "one to many")
  colnames(percentages) <- cell.types
  colors <- c('#FFEFDB', '#FFD39B')
  p <- barplot(percentages, col=colors, border="white", main=paste0("one to one versus one to many interactions"),
               space=rep(1, length(cell.types)), cex.main=0.6, cex.axis=0.6, yaxt='n', cex.names=0.6, width=1)
  legend(2.5*length(cell.types), 100, c("# one to one interactions", "# interactions in hubs"), col=colors, cex=0.6, pch=15)
  x <- c()
  y <- c()
  numbers <- c()
  totals <- c()
  for (j in 1:length(cell.types)) {
    x <- c(x, rep(p[j], 2))
    y <- c(y, percentages[1, j], sum(percentages[1:2, j]))
    numbers <- c(numbers, counts[j, ])
    totals <- c(totals, sum(counts[j, ]))
  }
  text(x, y-6, labels=numbers, col="black", cex=0.5)
  text(p, -7.5, labels=totals, col="black", cex=0.6, font=2 )
  par(mar=c(3, 3, 3, 3) + 0.1)
  dev.off()
  
}

#
# Inputs:
# Outputs:
analyzeTADs <- function(interaction.data, TADs) {
  
  #
  TADs[TADs[, 1] == "chr23", 1] <- "chrX"
  TADs[TADs[, 1] == "chr24", 1] <- "chrY"
  
  # For each cell type, find out how many are inside TADs.
  for (cell.type in cell.types)  {
    
    #
    print(paste0("Analyzing TAD overlaps for cell type: ", cell.type))
    
    #
    lhs.mid <- (interaction.data[[cell.type]]$bait_start + interaction.data[[cell.type]]$bait_end)/2
    rhs.mid <- (interaction.data[[cell.type]]$otherEnd_start + interaction.data[[cell.type]]$otherEnd_end)/2
    lr <- interaction.data[[cell.type]][lhs.mid > rhs.mid, ]
    rl <- interaction.data[[cell.type]][rhs.mid > lhs.mid, ]
    stringent.bed <- rbind(cbind(rl$bait_chr, rl$bait_start, rl$otherEnd_end), 
                           cbind(lr$bait_chr, lr$otherEnd_start, lr$bait_end))
    stringent.bed <- data.frame(stringent.bed, stringsAsFactors=F)
    stringent.bed[, 1] <- as.character(stringent.bed[, 1])
    stringent.bed[, 2] <- as.numeric(stringent.bed[, 2])
    stringent.bed[, 3] <- as.numeric(stringent.bed[, 3])
    stringent.bed <- cbind(stringent.bed, paste0("id", 1:length(stringent.bed[, 1])))
    normal.bed <- rbind(cbind(rl$bait_chr, round(as.numeric((rl$bait_start + rl$bait_end)/2)), round(as.numeric((rl$otherEnd_start + rl$otherEnd_end)/2))), 
                        cbind(lr$bait_chr, round(as.numeric((lr$otherEnd_start + lr$otherEnd_end)/2)), round(as.numeric((lr$bait_start + lr$bait_end)/2))))
    normal.bed <- data.frame(normal.bed, stringsAsFactors=F)
    normal.bed[, 1] <- as.character(normal.bed[, 1])
    normal.bed[, 2] <- as.numeric(normal.bed[, 2])
    normal.bed[, 3] <- as.numeric(normal.bed[, 3])
    normal.bed <- cbind(normal.bed, paste0("id", 1:length(normal.bed[, 1])))
    
    #
    stringent.intersection <- bedTools.2in(bed1=stringent.bed, bed2=TADs, opt.string=" -wb")
    normal.intersection <- bedTools.2in(bed1=normal.bed, bed2=TADs, opt.string=" -wb")
    
    #
    num.stringent.overlaps <- c()
    for (i in 1:length(stringent.bed[, 1])) {
      num.stringent.overlaps <- c(num.stringent.overlaps, length(which(as.character(stringent.intersection[, 4]) == stringent.bed[i, 4])))
    }
    num.normal.overlaps <- c()
    for (i in 1:length(normal.bed[, 1])) {
      num.normal.overlaps <- c(num.normal.overlaps, length(which(as.character(normal.intersection[, 4]) == normal.bed[i, 4])))
    }
    
    #
    print(length(num.stringent.overlaps))
    print(table(num.stringent.overlaps))
    print(table(num.stringent.overlaps > 1))
    print(length(num.normal.overlaps))
    print(table(num.normal.overlaps))
    print(table(num.normal.overlaps > 1))
    
  }
  
}

