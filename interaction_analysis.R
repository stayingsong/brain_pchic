# Annotate ATAC-seq peaks according to their overlap with and distance from interactions, promoters, exons, and introns.
# Notes: Requires unexpanded peaks for overlap with exons and introns.
#
annotateATACSeqPeaks <- function(peaks.expanded, peaks.unexpanded, interactions, promoters, promoters.other, exons, introns) {
  
  # Create a data frame to hold annotation results.
  num.peaks <- length(peaks.expanded[, 1])
  peak.annotation <- data.frame(matrix("", num.peaks, 9), stringsAsFactors=F)
  rownames(peak.annotation) <- peaks.expanded$atac_seq_peak_ID
  colnames(peak.annotation) <- c("peak_type", "promoter_genes", "num_interactions", "proximal_genes", "proximal_other_genes", "distal_genes", "distal_other_genes", "exon_overlap", "intron_overlap")
  peak.annotation$peak_type <- as.character(peak.annotation$peak_type)
  peak.annotation$promoter_genes <- as.character(peak.annotation$promoter_genes)
  peak.annotation$num_interactions <- as.integer(peak.annotation$num_interactions)
  peak.annotation$proximal_genes <- as.character(peak.annotation$proximal_genes)
  peak.annotation$proximal_other_genes <- as.character(peak.annotation$proximal_other_genes)
  peak.annotation$distal_genes <- as.character(peak.annotation$distal_genes)
  peak.annotation$distal_other_genes <- as.character(peak.annotation$distal_other_genes)
  peak.annotation$exon_overlap <- as.logical(peak.annotation$exon_overlap)
  peak.annotation$intron_overlap <- as.logical(peak.annotation$intron_overlap)
  
  # Create data frames to hold IDs of genes whose promoters interact with the lhs or rhs of interactions. Consider strict and other sets of genes separately.
  num.interactions <- length(interactions[, 1])
  interaction.genes <- data.frame(matrix("", num.interactions, 2), stringsAsFactors=F)
  rownames(interaction.genes) <- interactions$interaction_ID
  colnames(interaction.genes) <- c("lhs_genes", "rhs_genes")
  interaction.genes$lhs_genes <- as.character(interaction.genes$lhs_genes)
  interaction.genes$rhs_genes <- as.character(interaction.genes$rhs_genes)
  interaction.other.genes <- data.frame(matrix("", num.interactions, 2), stringsAsFactors=F)
  rownames(interaction.other.genes) <- interactions$interaction_ID
  colnames(interaction.other.genes) <- c("lhs_genes", "rhs_genes")
  interaction.other.genes$lhs_genes <- as.character(interaction.other.genes$lhs_genes)
  interaction.other.genes$rhs_genes <- as.character(interaction.other.genes$rhs_genes)
  
  # Intersect the lhs and rhs bins of each interaction with the promoters.
  interactions.promoters.lhs.results <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.promoters.other.lhs.results <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters.other[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.promoters.rhs.results <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.promoters.other.rhs.results <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters.other[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  for (i in 1:num.interactions) {
    
    # Grab the current interaction ID.
    current.interaction.id <- as.character(interactions$interaction_ID[i])
    
    # Grab the gene IDs of all genes whose promoters intersect either the lhs and rhs bins of each interaction.
    gene.ids.lhs <- as.character(interactions.promoters.lhs.results[as.character(interactions.promoters.lhs.results[, 4]) == current.interaction.id, 8])
    gene.ids.rhs <- as.character(interactions.promoters.rhs.results[as.character(interactions.promoters.rhs.results[, 4]) == current.interaction.id, 8])
    if (length(gene.ids.lhs) > 0) {
      interaction.genes$lhs_genes[i] <- paste0(unique(gene.ids.lhs), collapse=",")
    }
    if (length(gene.ids.rhs) > 0) {
      interaction.genes$rhs_genes[i] <- paste0(unique(gene.ids.rhs), collapse=",")
    }
    other.gene.ids.lhs <- as.character(interactions.promoters.other.lhs.results[as.character(interactions.promoters.other.lhs.results[, 4]) == current.interaction.id, 8])
    other.gene.ids.rhs <- as.character(interactions.promoters.other.rhs.results[as.character(interactions.promoters.other.rhs.results[, 4]) == current.interaction.id, 8])
    if (length(other.gene.ids.lhs) > 0) {
      interaction.other.genes$lhs_genes[i] <- paste0(unique(other.gene.ids.lhs), collapse=",")
    }
    if (length(other.gene.ids.rhs) > 0) {
      interaction.other.genes$rhs_genes[i] <- paste0(unique(other.gene.ids.rhs), collapse=",")
    }
    
  }
  
  # Intersect ATAC-seq peaks with interactions, promoters, exons, and introns.
  peaks.promoters.results <- bedTools.2in(bed1=peaks.expanded[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  peaks.exons.results <- bedTools.2in(bed1=peaks.unexpanded[, c(1:3)], bed2=exons[, c(1:3)], opt.string="-c")
  peaks.introns.results <- bedTools.2in(bed1=peaks.unexpanded[, c(1:3)], bed2=introns[, c(1:3)], opt.string="-c")
  peaks.interactions.lhs.results <- bedTools.2in(bed1=peaks.expanded[, c(1:3, 11)], bed2=interactions[, c(1:3, 11)], opt.string="-wb")
  peaks.interactions.rhs.results <- bedTools.2in(bed1=peaks.expanded[, c(1:3, 11)], bed2=interactions[, c(5:7, 11)], opt.string="-wb")
  print(paste0(num.peaks, " peaks to be annotated..."))
  for (i in 1:num.peaks) {
    
    # Track progress.
    if (i %% 10000 == 0)
      print(paste0(i, " peaks annotated"))
    
    # Process intersections with exons/introns first.
    if (peaks.exons.results[i, 4] > 0) {
      peak.annotation$exon_overlap[i] <- TRUE
    } else {
      peak.annotation$exon_overlap[i] <- FALSE
    } 
    if (peaks.introns.results[i, 4] > 0) {
      peak.annotation$intron_overlap[i] <- TRUE
    } else {
      peak.annotation$intron_overlap[i] <- FALSE
    }
    
    # Grab the current ATAC-seq peak ID.
    current.peak.id <- as.character(peaks.expanded$atac_seq_peak_ID[i])
    
    # First determine whether the ATAC-seq peak intersects any promoters, and if so, grab the gene IDs of the genes associated with the promoters. This is using the strict set only.
    gene.ids.promoters <- as.character(peaks.promoters.results[as.character(peaks.promoters.results[, 4]) == current.peak.id, 8])
    if (length(gene.ids.promoters) > 0) {
      peak.annotation$peak_type[i] <- "promoter"
      peak.annotation$promoter_genes[i] <- paste0(unique(gene.ids.promoters), collapse=",")
    } else {
      peak.annotation$peak_type[i] <- "distal"
    }
    
    # Next determine which interactions intersect the ATAC-seq peak, and grab any proximal and distal genes overlapping the interacting bins, if any.
    interaction.ids.lhs <- as.character(peaks.interactions.lhs.results[as.character(peaks.interactions.lhs.results[, 4]) == current.peak.id, 8])
    interaction.ids.rhs <- as.character(peaks.interactions.rhs.results[as.character(peaks.interactions.rhs.results[, 4]) == current.peak.id, 8])
    interaction.ids.unique <- unique(c(interaction.ids.lhs, interaction.ids.rhs))
    peak.annotation$num_interactions[i] <- length(interaction.ids.unique)
    if (peak.annotation$num_interactions[i] > 0) {
      
      # Go through each interaction and fill in the proximal and distal genes overlapping the interacting bins.
      for (j in 1:peak.annotation$num_interactions[i]) {
        
        # Get the gene IDs associated with the lhs and rhs bins of the interaction, and also figure out which ends of the interaction intersect the ATAC-seq peaak.
        gene.ids.retrieved <- interaction.genes[interaction.ids.unique[j], ]
        other.gene.ids.retrieved <- interaction.other.genes[interaction.ids.unique[j], ]
        in.lhs <- interaction.ids.unique[j] %in% interaction.ids.lhs
        in.rhs <- interaction.ids.unique[j] %in% interaction.ids.rhs
        
        # If the lhs bin of the interaction intersects the ATAC-seq peak, add the corresponding proximal (lhs) and distal (rhs) genes, and vice versa for the rhs bin.
        if (in.lhs) {
          
          peak.annotation$proximal_genes[i] <- paste0(c(peak.annotation$proximal_genes[i], gene.ids.retrieved[1]), collapse=",")
          peak.annotation$distal_genes[i] <- paste0(c(peak.annotation$distal_genes[i], gene.ids.retrieved[2]), collapse=",")
          peak.annotation$proximal_other_genes[i] <- paste0(c(peak.annotation$proximal_other_genes[i], other.gene.ids.retrieved[1]), collapse=",")
          peak.annotation$distal_other_genes[i] <- paste0(c(peak.annotation$distal_other_genes[i], other.gene.ids.retrieved[2]), collapse=",")
          
        }
        if (in.rhs) {
          
          peak.annotation$proximal_genes[i] <- paste0(c(peak.annotation$proximal_genes[i], gene.ids.retrieved[2]), collapse=",")
          peak.annotation$distal_genes[i] <- paste0(c(peak.annotation$distal_genes[i], gene.ids.retrieved[1]), collapse=",")
          peak.annotation$proximal_other_genes[i] <- paste0(c(peak.annotation$proximal_other_genes[i], other.gene.ids.retrieved[2]), collapse=",")
          peak.annotation$distal_other_genes[i] <- paste0(c(peak.annotation$distal_other_genes[i], other.gene.ids.retrieved[1]), collapse=",")
          
        }
        
        # Clean up duplicate gene IDs.
        proximal.genes.temp <- unique(unlist(strsplit(peak.annotation$proximal_genes[i], split=",")))
        proximal.genes.temp <- proximal.genes.temp[proximal.genes.temp != ""]
        peak.annotation$proximal_genes[i] <- paste0(proximal.genes.temp, collapse=",")
        distal.genes.temp <- unique(unlist(strsplit(peak.annotation$distal_genes[i], split=",")))
        distal.genes.temp <- distal.genes.temp[distal.genes.temp != ""]
        peak.annotation$distal_genes[i] <- paste0(distal.genes.temp, collapse=",")
        proximal.other.genes.temp <- unique(unlist(strsplit(peak.annotation$proximal_other_genes[i], split=",")))
        proximal.other.genes.temp <- proximal.other.genes.temp[proximal.other.genes.temp != ""]
        peak.annotation$proximal_other_genes[i] <- paste0(proximal.other.genes.temp, collapse=",")
        distal.other.genes.temp <- unique(unlist(strsplit(peak.annotation$distal_other_genes[i], split=",")))
        distal.other.genes.temp <- distal.other.genes.temp[distal.other.genes.temp != ""]
        peak.annotation$distal_other_genes[i] <- paste0(distal.other.genes.temp, collapse=",")
        
      }
      
    }
    
  }
  
  # Return ATAC-seq peak annotation results.
  return(peak.annotation)
  
}

# Annotate interactions according to their overlap with promoters, ATAC-seq peaks, and features.
# Notes:
#
annotateInteractions <- function(interactions, atac.seq.peaks, atac.seq.annotations, features=NULL, promoters, promoters.other) {
  
  # Interaction annotations are stored in a list organized by feature types (each annotation follows the same ordering of interactions in the original set of interactions).
  interaction.annotation <- list()
  num.interactions <- length(interactions[, 1])
  
  # Process features if specified.
  if (!is.null(features)) {
    
    # Get a list of all features to be used for annotating interactions.
    feature.categories <- names(features)
    feature.table <- c()
    for (i in 1:length(feature.categories)) {
      
      feature.names <- names(features[[feature.categories[i]]])
      feature.table <- rbind(feature.table, cbind(rep(feature.categories[i], length(feature.names)), feature.names))
      
    }
    
    # Perform intersections between features and interaction lhs/rhs bins.
    num.features <- length(feature.table[, 1])
    feature.results.lhs <- list()
    feature.results.rhs <- list()
    for (i in 1:num.features) {
      
      feature.results.lhs[[feature.table[i, 2]]] <- bedTools.2in(bed1=interactions[, c(1:3, 11)], 
                                                                 bed2=cbind(features[[feature.table[i, 1]]][[feature.table[i, 2]]][, c(1:3)], features[[feature.table[i, 1]]][[feature.table[i, 2]]]$feature_ID), 
                                                                 opt.string="-wb", num.cols=8)
      feature.results.rhs[[feature.table[i, 2]]] <- bedTools.2in(bed1=interactions[, c(5:7, 11)], 
                                                                 bed2=cbind(features[[feature.table[i, 1]]][[feature.table[i, 2]]][, c(1:3)], features[[feature.table[i, 1]]][[feature.table[i, 2]]]$feature_ID),
                                                                 opt.string="-wb", num.cols=8)
      interaction.annotation[[paste0(feature.table[i, 2], "_lhs")]] <- rep(0, num.interactions)
      interaction.annotation[[paste0(feature.table[i, 2], "_lhs_ids")]] <- rep("", num.interactions)
      interaction.annotation[[paste0(feature.table[i, 2], "_rhs")]] <- rep(0, num.interactions)
      interaction.annotation[[paste0(feature.table[i, 2], "_rhs_ids")]] <- rep("", num.interactions)
      
    }
    
  }
  
  # Perform intersections between promoters/ATAC-seq peaks and interaction lhs/rhs bins.
  interactions.promoters.results.lhs <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.promoters.other.results.lhs <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters.other[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.peaks.results.lhs <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=atac.seq.peaks[, c(1:3, 11)], opt.string="-wb", num.cols=8)
  interactions.promoters.results.rhs <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.promoters.other.results.rhs <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters.other[, c(1:3, 5)], opt.string="-wb", num.cols=8)
  interactions.peaks.results.rhs <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=atac.seq.peaks[, c(1:3, 11)], opt.string="-wb", num.cols=8)
  interaction.annotation[["promoter_lhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_lhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["promoter_other_lhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_other_lhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["promoter_ATAC-seq_lhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_ATAC-seq_lhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["distal_ATAC-seq_lhs"]] <- rep(0, num.interactions)
  interaction.annotation[["distal_ATAC-seq_lhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["promoter_rhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_rhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["promoter_other_rhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_other_rhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["promoter_ATAC-seq_rhs"]] <- rep(0, num.interactions)
  interaction.annotation[["promoter_ATAC-seq_rhs_ids"]] <- rep("", num.interactions)
  interaction.annotation[["distal_ATAC-seq_rhs"]] <- rep(0, num.interactions)
  interaction.annotation[["distal_ATAC-seq_rhs_ids"]] <- rep("", num.interactions)
  
  # Intersect interactions with promoters, ATAC-seq peaks, and features.
  print(paste0(num.interactions, " interactions to be annotated..."))
  for (i in 1:num.interactions) {
    
    # Track progress.
    if (i %% 10000 == 0)
      print(paste0(i, " interactions annotated"))
    
    # Grab current interaction ID.
    current.interaction.ID <- as.character(interactions$interaction_ID[i])
    
    # Check intersection with promoter intervals.
    results <- as.character(interactions.promoters.results.lhs[as.character(interactions.promoters.results.lhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      interaction.annotation[["promoter_lhs"]][i] <- length(unique(results))
      interaction.annotation[["promoter_lhs_ids"]][i] <- paste0(unique(results), collapse=",")
      
    }
    results <- as.character(interactions.promoters.other.results.lhs[as.character(interactions.promoters.other.results.lhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      interaction.annotation[["promoter_other_lhs"]][i] <- length(unique(results))
      interaction.annotation[["promoter_other_lhs_ids"]][i] <- paste0(unique(results), collapse=",")
      
    }
    results <- as.character(interactions.promoters.results.rhs[as.character(interactions.promoters.results.rhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      interaction.annotation[["promoter_rhs"]][i] <- length(unique(results))
      interaction.annotation[["promoter_rhs_ids"]][i] <- paste0(unique(results), collapse=",")
      
    }
    results <- as.character(interactions.promoters.other.results.rhs[as.character(interactions.promoters.other.results.rhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      interaction.annotation[["promoter_other_rhs"]][i] <- length(unique(results))
      interaction.annotation[["promoter_other_rhs_ids"]][i] <- paste0(unique(results), collapse=",")
      
    }
    
    # Check intersection with ATAC-seq peaks. 
    results <- as.character(interactions.peaks.results.lhs[as.character(interactions.peaks.results.lhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      promoter.results <- results[atac.seq.annotations[results, "peak_type"] == "promoter"]
      interaction.annotation[["promoter_ATAC-seq_lhs"]][i] <- length(unique(promoter.results))
      interaction.annotation[["promoter_ATAC-seq_lhs_ids"]][i] <- paste0(unique(promoter.results), collapse=",")
      distal.results <- results[atac.seq.annotations[results, "peak_type"] == "distal"]
      interaction.annotation[["distal_ATAC-seq_lhs"]][i] <- length(unique(distal.results))
      interaction.annotation[["distal_ATAC-seq_lhs_ids"]][i] <- paste0(unique(distal.results), collapse=",")
      
    }
    results <- as.character(interactions.peaks.results.rhs[as.character(interactions.peaks.results.rhs[, 4]) == current.interaction.ID, 8])
    if (length(results) > 0) {
      
      promoter.results <- results[atac.seq.annotations[results, "peak_type"] == "promoter"]
      interaction.annotation[["promoter_ATAC-seq_rhs"]][i] <- length(unique(promoter.results))
      interaction.annotation[["promoter_ATAC-seq_rhs_ids"]][i] <- paste0(unique(promoter.results), collapse=",")
      distal.results <- results[atac.seq.annotations[results, "peak_type"] == "distal"]
      interaction.annotation[["distal_ATAC-seq_rhs"]][i] <- length(unique(distal.results))
      interaction.annotation[["distal_ATAC-seq_rhs_ids"]][i] <- paste0(unique(distal.results), collapse=",")
      
    }
    
    # Check intersection with features.
    if (!is.null(features)) {
      
      # Iterate through each feature in the feature table.
      for (j in 1:num.features) {
        
        results <- as.character(feature.results.lhs[[feature.table[j, 2]]][as.character(feature.results.lhs[[feature.table[j, 2]]][, 4]) == current.interaction.ID, 8])
        if (length(results) > 0) {
          
          interaction.annotation[[paste0(feature.table[j, 2], "_lhs")]][i] <- length(unique(results))
          interaction.annotation[[paste0(feature.table[j, 2], "_lhs_ids")]][i] <- paste0(unique(results), collapse=",")
          
        }
        results <- as.character(feature.results.rhs[[feature.table[j, 2]]][as.character(feature.results.rhs[[feature.table[j, 2]]][, 4]) == current.interaction.ID, 8])
        if (length(results) > 0) {
          
          interaction.annotation[[paste0(feature.table[j, 2], "_rhs")]][i] <- length(unique(results))
          interaction.annotation[[paste0(feature.table[j, 2], "_rhs_ids")]][i] <- paste0(unique(results), collapse=",")
          
        }
        
      }
      
    }
    
  }
  
  # Return interaction annotation results.
  return(interaction.annotation)
  
}

# Grab score of most significant interaction at each cell type overlapping each loci that is significant in at least one cell type.
# Notes: Takes expanded interaction loci at which to make heatmap scores from cell types in expanded complete interaction data.
#
makeLociSpecificityHeatmap <- function(significant.loci, interaction.data, overlap.threshold) {
  
  # First concatenate all interactin loci coordinates into one variable.
  all.loci <- c()
  for (cell.type in names(significant.loci)) {
    
    print(paste0("Inspecting ", length(significant.loci[[cell.type]][, 1]), " significant interacting loci from cell type: ", cell.type))
    all.loci <- rbind(all.loci, significant.loci[[cell.type]][, c(1:8, 11)])
    
  }
  
  # Get list of unique loci. For each unique loci, concatenate all the interaction IDs so we can determine which cell types they came from.
  print(dim(unique(all.loci[, 1:8]))[1] == dim(unique(all.loci[, c(1:3, 5:7)]))[1]) # Should be TRUE.
  unique.loci <- all.loci %>% group_by_(.dots=c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name")) %>% 
    summarize(interaction_ID=paste0(interaction_ID, collapse=","))
  unique.loci <- data.frame(unique.loci)
  num.unique.loci <- length(unique.loci[, 1])
  print(paste0("There are a total of ", num.unique.loci, " unique significant interacting loci to analyze."))
  print(length(unique(unlist(strsplit(paste0(unique.loci$interaction_ID, collapse=","), split=",")))) == length(all.loci[, 1])) # Should be TRUE.
  print(table(unique.loci[,1] == unique.loci[,5])) # Should be TRUE.
  
  # Create matrix for holding results (rows are significant loci, columns are cell types).
  specificity.heatmap <- matrix(NA, num.unique.loci, length(cell.types))
  rownames(specificity.heatmap) <- unique.loci$interaction_ID
  colnames(specificity.heatmap) <- cell.types
  print(dim(specificity.heatmap))
  
  # Break it down by chromosome.
  for (chr in unique(unique.loci$bait_chr)) {
    
    # Get subsets of unique loci by chromosome.
    unique.loci.by.chr <- unique.loci[unique.loci$bait_chr == chr, ]
    unique.loci.by.chr.rows <- which(unique.loci$bait_chr == chr)
    print(paste0("Processing ", chr, " with ", length(unique.loci.by.chr[, 1]), " entries..."))
    
    # Process each cell type separately.
    for (cell.type in cell.types) {
      
      # Get subsets of interaction data for current cell type by chromosome.
      current.interaction.data <- interaction.data[[cell.type]][interaction.data[[cell.type]]$bait_chr == chr, ]
      
      # Intersect lhs coordinates of significant loci with lhs coordinates of all loci for current cell type & chr. Do same with rhs.
      lhs.results <- bedTools.2in(bed1=unique.loci.by.chr[, c(1:3, 9)], bed2=current.interaction.data[, c(1:3, 11)], opt.string=paste0('-wa -wb -f ', overlap.threshold))
      rhs.results <- bedTools.2in(bed1=unique.loci.by.chr[, c(5:7, 9)], bed2=current.interaction.data[, c(5:7, 11)], opt.string=paste0('-wa -wb -f ', overlap.threshold))
      
      # For each significant loci, see what it intersects with and grab the mean score. Then update the table.
      print(paste0("Processing ", length(unique.loci.by.chr[, 1]), " entries for cell type: ", cell.type))
      for (i in 1:length(unique.loci.by.chr[, 1])) {
        
        # Track progress.
        if (i %% 1000 == 0)
          print(paste0(i, " loci analyzed"))
        
        # An interaction overlaps if both the lhs and rhs overlap. Get IDs of interactions overlapping with lhs, then rhs, then see common.
        current.id <- unique.loci.by.chr$interaction_ID[i]
        lhs.all.intersections <- as.character(lhs.results[lhs.results[, 4] == current.id, 8])
        rhs.all.intersections <- as.character(rhs.results[rhs.results[, 4] == current.id, 8])
        both.all.intersections <- lhs.all.intersections[lhs.all.intersections %in% rhs.all.intersections]
        
        # Grab the maximum score from all interactions overlapping the current loci.
        if (length(both.all.intersections) > 0) {
          specificity.heatmap[unique.loci.by.chr.rows[i], cell.type] <- max(current.interaction.data$score[current.interaction.data$interaction_ID %in% both.all.intersections])
        } else if (length(both.all.intersections) == 0) {
          specificity.heatmap[unique.loci.by.chr.rows[i], cell.type] <- 0
        } else {
          print("Error")
        }
        
      }
      
    }
    
  }
  
  # Return specificity heatmap.
  return(specificity.heatmap)
  
}

# Sample all interactions in a distance-matched manner.
# Notes:
#
sampleDistanceMatchedInteractions <- function(interactions.sig, interactions.all, distance.bins, fold, score.bounds) {
  
  # Set distance bin starts and ends.
  interactions.by.distance <- matrix(0, length(distance.bins), 3)
  interactions.by.distance[, 1] <- distance.bins
  interactions.by.distance[, 2] <- c(distance.bins[2:length(distance.bins)], Inf)
  colnames(interactions.by.distance) <- c("lower", "upper", "num_interactions")
  
  # Calculate distances for the significant interactions as well as for all interactions.
  distances.sig <- abs((interactions.sig[, 6] + interactions.sig[, 7])/2 - (interactions.sig[, 2] + interactions.sig[, 3])/2)
  distances.all <- abs((interactions.all[, 6] + interactions.all[, 7])/2 - (interactions.all[, 2] + interactions.all[, 3])/2)

  # For each distance bin, count the number of significant interactions and sample that number of nonsignificant interactions*fold.
  interactions.sampled <- c()
  for (i in 1:length(distance.bins)) {
    
    # Count how many significant interactions fall in the current distance bin.
    interactions.by.distance[i, 3] <- length(which((distances.sig >= interactions.by.distance[i, 1]) & (distances.sig < interactions.by.distance[i, 2])))
    # print(interactions.by.distance[i, ])
    
    # Take the subset of nonsignificant interactions matching the current distance cutoffs.
    subset.all <- interactions.all[which((distances.all >= interactions.by.distance[i, 1]) & (distances.all < interactions.by.distance[i, 2])), ]
    interactions.sampled <- rbind(interactions.sampled, sampleInteractionsByScore(subset.all, score.bounds[1], score.bounds[2], interactions.by.distance[i, 3]*fold))
    
  }

  # Return sampled interactions.
  return(interactions.sampled)
  
}

# Analyze gene promoter hub-related phenomena.
# Notes:
#
analyzeGenePromoterHubs <- function(interaction.data, annotation.data, promoters, expression.data, resolution) {
  
  # Iterate through each cell type.
  hub.results <- list()
  gene.ids <- unique(promoters$gene_id)
  cell.types <- names(interaction.data)
  for (cell.type in cell.types) {
    
    # Create a data frame for each cell type which has a row for each gene ID and columns for gene ID, # of interactions by class, expression value, and gene type.
    print(paste0("Analyzing gene promoter hubs for ", length(gene.ids), " genes for cell type: ", cell.type))
    hub.results[[cell.type]] <- data.frame(matrix("", length(gene.ids), 14), stringsAsFactors=F)
    rownames(hub.results[[cell.type]]) <- gene.ids
    colnames(hub.results[[cell.type]]) <- c("gene_type", "gene_expression", 
                                            "num_interactions_all", "ids_interactions_all", "num_interactions_PP", "ids_interactions_PP", "num_interactions_PO", "ids_interactions_PO", 
                                            "num_interactions_PA", "ids_interactions_PA", "num_interactions_PE", "ids_interactions_PE", "num_interactions_PR", "ids_interactions_PR")
    hub.results[[cell.type]]$gene_type <- as.character(hub.results[[cell.type]]$gene_type)
    hub.results[[cell.type]]$gene_expression <- as.numeric(hub.results[[cell.type]]$gene_expression)
    hub.results[[cell.type]]$num_interactions_all <- as.integer(hub.results[[cell.type]]$num_interactions_all)
    hub.results[[cell.type]]$ids_interactions_all <- as.character(hub.results[[cell.type]]$ids_interactions_all)
    hub.results[[cell.type]]$num_interactions_PP <- as.integer(hub.results[[cell.type]]$num_interactions_PP)
    hub.results[[cell.type]]$ids_interactions_PP <- as.character(hub.results[[cell.type]]$ids_interactions_PP)
    hub.results[[cell.type]]$num_interactions_PO <- as.integer(hub.results[[cell.type]]$num_interactions_PO)
    hub.results[[cell.type]]$ids_interactions_PO <- as.character(hub.results[[cell.type]]$ids_interactions_PO)
    hub.results[[cell.type]]$num_interactions_PA <- as.integer(hub.results[[cell.type]]$num_interactions_PA)
    hub.results[[cell.type]]$ids_interactions_PA <- as.character(hub.results[[cell.type]]$ids_interactions_PA)
    hub.results[[cell.type]]$num_interactions_PE <- as.integer(hub.results[[cell.type]]$num_interactions_PE)
    hub.results[[cell.type]]$ids_interactions_PE <- as.character(hub.results[[cell.type]]$ids_interactions_PE)
    hub.results[[cell.type]]$num_interactions_PR <- as.integer(hub.results[[cell.type]]$num_interactions_PR)
    hub.results[[cell.type]]$ids_interactions_PR <- as.character(hub.results[[cell.type]]$ids_interactions_PR)
    
    # Get IDs for promoter to promoter and promoter to other interactions for the current cell type.
    ann <- annotation.data[[cell.type]]
    Pr_Pr.ids <- as.character(interaction.data[[cell.type]][((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) &
                                                              ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0)), 11])
    Pr_Np.ids <- as.character(interaction.data[[cell.type]][(((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) &
                                                               (ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0)) |
                                                              ((ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0) &
                                                                 ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0))), 11])
    
    # For each gene, determine how many interactions its promoters overlap with, and also fill in the rest of the results table.
    interacting.loci <- rbind(interaction.data[[cell.type]][, c(1:3, 11)], setNames(interaction.data[[cell.type]][, c(5:7, 11)], names(interaction.data[[cell.type]][, c(1:3, 11)])))
    promoters.interactions.results <- bedTools.2in(bed1=promoters[, c(1:3, 5)], bed2=interacting.loci, opt.string="-wb")
    for (i in 1:length(gene.ids)) {
      
      # Track progress.
      if (i %% 10000 == 0)
        print(paste0(i, " genes analyzed"))
      
      # Grab the type and expression of the current gene.
      hub.results[[cell.type]]$gene_type[i] <- unique(promoters$gene_type[promoters$gene_id == gene.ids[i]])
      if (gene.ids[i] %in% rownames(expression.data)) {
        hub.results[[cell.type]]$gene_expression[i] <- expression.data[gene.ids[i], cell.type]
      } else {
        hub.results[[cell.type]]$gene_expression[i] <- NA
      }
      
      # Get # of interactions intersecting with the current gene's promoters.
      unique.interaction.ids <- unique(as.character(promoters.interactions.results[as.character(promoters.interactions.results[, 4]) == gene.ids[i], 8]))
      unique.interaction.ids <- unique.interaction.ids[unique.interaction.ids != ""]
      hub.results[[cell.type]]$num_interactions_all[i] <- length(unique.interaction.ids)
      hub.results[[cell.type]]$ids_interactions_all[i] <- paste0(unique.interaction.ids, collapse=",")
      
      # If there are any intersecting interactions, perform a more detailed breakdown of the interaction types.
      if (length(unique.interaction.ids) > 0) {
        
        # Count # of promoter to distal ATAC-seq peak interactions intersecting with the current gene's promoters.
        intersecting.indices <- which(interaction.data[[cell.type]]$interaction_ID %in% unique.interaction.ids)
        intersecting.promoter.lhs <- annotation.data[[cell.type]]$promoter_lhs[intersecting.indices]
        intersecting.promoter.ATAC.seq.lhs <- annotation.data[[cell.type]]$'promoter_ATAC-seq_lhs'[intersecting.indices]
        intersecting.distal.ATAC.seq.lhs <- annotation.data[[cell.type]]$'distal_ATAC-seq_lhs'[intersecting.indices]
        intersecting.promoter.rhs <- annotation.data[[cell.type]]$promoter_rhs[intersecting.indices]
        intersecting.promoter.ATAC.seq.rhs <- annotation.data[[cell.type]]$'promoter_ATAC-seq_rhs'[intersecting.indices]
        intersecting.distal.ATAC.seq.rhs <- annotation.data[[cell.type]]$'distal_ATAC-seq_rhs'[intersecting.indices]
        Pr_Da.ids <- unique.interaction.ids[(((intersecting.promoter.lhs > 0) | (intersecting.promoter.ATAC.seq.lhs > 0)) & 
                                               (intersecting.promoter.rhs == 0) & (intersecting.promoter.ATAC.seq.rhs == 0) & (intersecting.distal.ATAC.seq.rhs > 0)) |
                                              ((intersecting.promoter.lhs == 0) & (intersecting.promoter.ATAC.seq.lhs == 0) & (intersecting.distal.ATAC.seq.lhs > 0) & 
                                                 ((intersecting.promoter.rhs > 0) | (intersecting.promoter.ATAC.seq.rhs > 0)))]
        Pr_Da.ids <- Pr_Da.ids[Pr_Da.ids != ""]
        hub.results[[cell.type]]$num_interactions_PA[i] <- length(Pr_Da.ids)
        hub.results[[cell.type]]$ids_interactions_PA[i] <- paste0(Pr_Da.ids, collapse=",")
        
        # Count # of promoter to promoter and promoter to other interactions intersecting with the current gene's promoters.
        hub.results[[cell.type]]$num_interactions_PP[i] <- length(unique.interaction.ids[unique.interaction.ids %in% Pr_Pr.ids])
        hub.results[[cell.type]]$ids_interactions_PP[i] <- paste0(unique.interaction.ids[unique.interaction.ids %in% Pr_Pr.ids], collapse=",")
        hub.results[[cell.type]]$num_interactions_PO[i] <- length(unique.interaction.ids[unique.interaction.ids %in% Pr_Np.ids])
        hub.results[[cell.type]]$ids_interactions_PO[i] <- paste0(unique.interaction.ids[unique.interaction.ids %in% Pr_Np.ids], collapse=",")
        
        # Count # of promoter to enhancer and promoter to repressive interactions intersecting with the current gene's promoters.
        intersecting.enhancer.lhs <- annotation.data[[cell.type]]$Enh_lhs[intersecting.indices] +
          annotation.data[[cell.type]]$EnhG_lhs[intersecting.indices] + annotation.data[[cell.type]]$EnhBiv_lhs[intersecting.indices]
        intersecting.repressive.lhs <- annotation.data[[cell.type]]$Het_lhs[intersecting.indices] +
          annotation.data[[cell.type]]$ReprPC_lhs[intersecting.indices]
        intersecting.enhancer.rhs <- annotation.data[[cell.type]]$Enh_rhs[intersecting.indices] +
          annotation.data[[cell.type]]$EnhG_rhs[intersecting.indices] + annotation.data[[cell.type]]$EnhBiv_rhs[intersecting.indices]
        intersecting.repressive.rhs <- annotation.data[[cell.type]]$Het_rhs[intersecting.indices] +
          annotation.data[[cell.type]]$ReprPC_rhs[intersecting.indices]
        Pr_Pe.ids <- unique.interaction.ids[(((intersecting.promoter.lhs > 0) | (intersecting.promoter.ATAC.seq.lhs > 0)) & 
                                               (intersecting.promoter.rhs == 0) & (intersecting.promoter.ATAC.seq.rhs == 0) & (intersecting.enhancer.rhs > 0)) |
                                              ((intersecting.promoter.lhs == 0) & (intersecting.promoter.ATAC.seq.lhs == 0) & (intersecting.enhancer.lhs > 0) & 
                                                 ((intersecting.promoter.rhs > 0) | (intersecting.promoter.ATAC.seq.rhs > 0)))]
        Pr_PR.ids <- unique.interaction.ids[(((intersecting.promoter.lhs > 0) | (intersecting.promoter.ATAC.seq.lhs > 0)) & 
                                               (intersecting.promoter.rhs == 0) & (intersecting.promoter.ATAC.seq.rhs == 0) & (intersecting.repressive.rhs > 0)) |
                                              ((intersecting.promoter.lhs == 0) & (intersecting.promoter.ATAC.seq.lhs == 0) & (intersecting.repressive.lhs > 0) & 
                                                 ((intersecting.promoter.rhs > 0) | (intersecting.promoter.ATAC.seq.rhs > 0)))]
        Pr_Pe.ids <- Pr_Pe.ids[Pr_Pe.ids != ""]
        Pr_PR.ids <- Pr_PR.ids[Pr_PR.ids != ""]
        hub.results[[cell.type]]$num_interactions_PE[i] <- length(Pr_Pe.ids)
        hub.results[[cell.type]]$ids_interactions_PE[i] <- paste0(Pr_Pe.ids, collapse=",")
        hub.results[[cell.type]]$num_interactions_PR[i] <- length(Pr_PR.ids)
        hub.results[[cell.type]]$ids_interactions_PR[i] <- paste0(Pr_PR.ids, collapse=",")
        
      } else {
        
        # If there aren't any intersecting interactions, set all the counts to 0.
        hub.results[[cell.type]]$num_interactions_PP[i] <- 0
        hub.results[[cell.type]]$ids_interactions_PP[i] <- ""
        hub.results[[cell.type]]$num_interactions_PO[i] <- 0
        hub.results[[cell.type]]$ids_interactions_PO[i] <- ""
        hub.results[[cell.type]]$num_interactions_PA[i] <- 0
        hub.results[[cell.type]]$ids_interactions_PA[i] <- ""
        hub.results[[cell.type]]$num_interactions_PE[i] <- 0
        hub.results[[cell.type]]$ids_interactions_PE[i] <- ""
        hub.results[[cell.type]]$num_interactions_PR[i] <- 0
        hub.results[[cell.type]]$ids_interactions_PR[i] <- ""
        
      }
      
    }
    
  }
  
  # Return gene promoter hub results.
  return(hub.results)
  
}

# Filters hub results based on genes meeting certain criteria.
# Notes:
#
filterHubResults <- function(hub.results, gene.ids, min.expression, min.num.interactions) {
  
  # Filter list of gene promoters to be analyzed for each cell type.
  cell.types <- names(hub.results)
  hub.results.filtered <- list()
  for (cell.type in cell.types) {
    
    hub.results.filtered[[cell.type]] <- hub.results[[cell.type]][hub.results[[cell.type]]$gene_expression > min.expression, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$num_interactions_all >= min.num.interactions, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][rownames(hub.results.filtered[[cell.type]]) %in% gene.ids, ]
    
  }
  
  # Return filtered hub results.
  return(hub.results.filtered)
  
}

# Analyze overlap between Vista elements, interactions, and gene promoters.
# Notes: Uses first promoter set (strict) for nearest/same fragment gene analysis, but uses second promoter set (all) for annotating target genes.
#
analyzeVistaOverlap <- function(vista.data, interaction.data, annotation.data, atac.seq.data, promoters, promoters.all, rmap, resolution) {
  
  # Create results table.
  num.elements <- length(vista.data[, 1])
  cell.types <- names(interaction.data)
  vista.results <- data.frame(matrix(NA, num.elements, 5 + 3*length(cell.types)), stringsAsFactors=F)
  rownames(vista.results) <- vista.data$vista_ID
  colnames(vista.results) <- c("nearest_gene", "nearest_gene_distance", "same_fragment_genes", "same_fragment_other_genes", "same_fragment_ambiguity", 
                               c(rbind(paste0(cell.types, "_atac_seq_peak"), paste0(cell.types, "_interactions"), paste0(cell.types, "_target_genes"))))
  vista.results$nearest_gene <- as.character(vista.results$nearest_gene)
  vista.results$nearest_gene_distance <- as.numeric(vista.results$nearest_gene_distance)
  vista.results$same_fragment_genes <- as.character(vista.results$same_fragment_genes)
  vista.results$same_fragment_other_genes <- as.character(vista.results$same_fragment_other_genes)
  vista.results$same_fragment_ambiguity <- as.logical(vista.results$same_fragment_ambiguity)
  for (cell.type in cell.types) {
    
    vista.results[, paste0(cell.type, "_atac_seq_peak")] <- as.logical(vista.results[, paste0(cell.type, "_atac_seq_peak")])
    vista.results[, paste0(cell.type, "_interactions")] <- as.integer(vista.results[, paste0(cell.type, "_interactions")])
    vista.results[, paste0(cell.type, "_target_genes")] <- as.character(vista.results[, paste0(cell.type, "_target_genes")])
    
  }
  
  # Intersect Vista elements and promoters for reference later (to identify nearest gene).
  vista.sorted <- sortBed(unique(vista.data[, 1:4]))
  promoters.sorted <- sortBed(promoters)
  promoters.sorted.filtered <- promoters.sorted[promoters.sorted[, 1] %in% unique(vista.sorted[, 1]), ]
  vista.promoters.results <- bedTools.2in(functionstring="bedtools closest", bed1=vista.sorted, bed2=promoters.sorted.filtered, opt.string="-t first -d")
  
  # Assign a promoter ID to each promoter for reference later. Also separate promoter sets.
  promoters.with.id <- cbind(promoters, paste0("promoter", 1:length(promoters[, 1])))
  promoters.other <- promoters.all[!(promoters.all$gene_type %in% promoters$gene_type), ]
  promoters.other.with.id <- cbind(promoters.other, paste0("promoter_other", 1:length(promoters.other[, 1])))
  promoters.all.with.id <- rbind(promoters.with.id, setNames(promoters.other.with.id, names(promoters.with.id)))
  
  # Expand Vista elements to minimum resolution and intersect with interactions.
  vista.data.res <- expandFeatures(vista.data, resolution)
  interactions.vista.results.lhs <- list()
  interactions.vista.results.rhs <- list()
  for (cell.type in cell.types) {
    
    interactions.vista.results.lhs[[cell.type]] <- bedTools.2in(bed1=interaction.data[[cell.type]][, c(1:3, 11)], bed2=vista.data.res[, 1:4], opt.string="-wb", num.cols=8)
    interactions.vista.results.rhs[[cell.type]] <- bedTools.2in(bed1=interaction.data[[cell.type]][, c(5:7, 11)], bed2=vista.data.res[, 1:4], opt.string="-wb", num.cols=8)
    
  }
  
  # Intersect Vista elements and ATAC-seq peaks.
  for (cell.type in cell.types) {
    
    vista.peaks.results <- bedTools.2in(bed1=vista.data.res[, 1:4], bed2=atac.seq.data[[cell.type]][, 1:4], opt.string="-c")
    vista.results[, paste0(cell.type, "_atac_seq_peak")] <- vista.peaks.results[, 5] > 0
    
  }
  
  # Intersect expanded Vista elements with restriction fragments.
  rmap.vista.results <- bedTools.2in(bed1=rmap[, 1:4], bed2=vista.data.res[, 1:4], opt.string="-wb", num.cols=8)
  
  # Process Vista elements one by one.
  print(paste0(num.elements, " Vista elements to process..."))
  for (i in 1:length(vista.results[, 1])) {
    
    # Track progress.
    if (i %% 100 == 0)
      print(paste0(i, " elements processed"))
    
    # Get current Vista element ID.
    vista.id <- rownames(vista.results)[i]
    
    # Get nearest gene results.
    vista.results$nearest_gene[i] <- as.character(vista.promoters.results[as.character(vista.promoters.results[, 4]) == vista.id, 11])
    vista.results$nearest_gene_distance[i] <- as.numeric(vista.promoters.results[as.character(vista.promoters.results[, 4]) == vista.id, 12])
    
    # For each Vista element, grab which restriction fragments it overlaps and which gene promoters (strict set) overlap those restriction fragments. This is to resolve same fragment ambiguity later.
    overlapping.bins <- rmap[rmap$rmap_ID %in% rmap.vista.results[as.character(rmap.vista.results[, 8]) == vista.id, 4], ]
    overlapping.bins.promoters.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.with.id, opt.string="-wb", num.cols=12)
    vista.results$same_fragment_genes[i] <- paste0(unique(as.character(overlapping.bins.promoters.results[, 11])), collapse=",")
    
    # Use list of ALL gene promoters for determining same fragment genes.
    overlapping.bins.promoters.other.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.other, opt.string="-wb", num.cols=11)
    vista.results$same_fragment_other_genes[i] <- paste0(unique(as.character(overlapping.bins.promoters.other.results[, 11])), collapse=",")

    # For storing lists of target promoters across all cell types.
    target.promoter.ids.cell.type <- list()
    target.promoter.ids.all <- c()
    
    # Analyze interactions for each cell type.
    for (cell.type in cell.types) {
      
      # For storing lists of target promoters for each cell type.
      target.promoter.ids.cell.type[[cell.type]] <- c()
      
      # Identify the interactions whose lhs intersects the expanded Vista element, then identify their target genes and gene promoters on the rhs.
      target.genes.rhs <- c()
      interaction.ids.lhs <- as.character(interactions.vista.results.lhs[[cell.type]][as.character(interactions.vista.results.lhs[[cell.type]][, 8]) == vista.id, 4])
      num.interactions.lhs <- length(unique(interaction.ids.lhs))
      if (num.interactions.lhs > 0) {
        
        # Find which genes are on the other ends of the overlapping interactions.
        selector <- interaction.data[[cell.type]]$interaction_ID %in% interaction.ids.lhs
        target.genes.rhs <- c(as.character(annotation.data[[cell.type]]$promoter_rhs_ids[selector]), 
                              as.character(annotation.data[[cell.type]]$promoter_other_rhs_ids[selector]))

        # Determine the exact promoters interacting with the other ends of the overlapping interactions.
        interacting.ends.rhs <- interaction.data[[cell.type]][selector, 5:7]
        promoters.interacting.ends.results <- bedTools.2in(bed1=promoters.with.id, bed2=interacting.ends.rhs, opt.string="-c")
        promoters.other.interacting.ends.results <- bedTools.2in(bed1=promoters.other.with.id, bed2=interacting.ends.rhs, opt.string="-c")
        target.promoter.ids.rhs <- as.character(promoters.interacting.ends.results[promoters.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.all <- c(target.promoter.ids.all, target.promoter.ids.rhs)
        
        # Get specific promoters of both sets to determine what target genes should be removed later.
        target.promoter.other.ids.rhs <- as.character(promoters.other.interacting.ends.results[promoters.other.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.cell.type[[cell.type]] <- c(target.promoter.ids.cell.type[[cell.type]], target.promoter.ids.rhs, target.promoter.other.ids.rhs)
        
      }
      
      # Identify the interactions whose rhs intersects the expanded Vista element, then identify their target genes and gene promoters on the lhs.
      target.genes.lhs <- c()
      interaction.ids.rhs <- as.character(interactions.vista.results.rhs[[cell.type]][as.character(interactions.vista.results.rhs[[cell.type]][, 8]) == vista.id, 4])
      num.interactions.rhs <- length(unique(interaction.ids.rhs))
      if (num.interactions.rhs > 0) {
        
        # Find which genes are on the other ends of the overlapping interactions.
        selector <- interaction.data[[cell.type]]$interaction_ID %in% interaction.ids.rhs
        target.genes.lhs <- c(as.character(annotation.data[[cell.type]]$promoter_lhs_ids[selector]),
                              as.character(annotation.data[[cell.type]]$promoter_other_lhs_ids[selector]))

        # Determine the exact promoters interacting with the other ends of the overlapping interactions.
        interacting.ends.lhs <- interaction.data[[cell.type]][selector, 1:3]
        promoters.interacting.ends.results <- bedTools.2in(bed1=promoters.with.id, bed2=interacting.ends.lhs, opt.string="-c")
        promoters.other.interacting.ends.results <- bedTools.2in(bed1=promoters.other.with.id, bed2=interacting.ends.lhs, opt.string="-c")
        target.promoter.ids.lhs <- as.character(promoters.interacting.ends.results[promoters.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.all <- c(target.promoter.ids.all, target.promoter.ids.lhs)
        
        # Get specific promoters of both sets to determine what target genes should be removed later.
        target.promoter.other.ids.lhs <- as.character(promoters.other.interacting.ends.results[promoters.other.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.cell.type[[cell.type]] <- c(target.promoter.ids.cell.type[[cell.type]], target.promoter.ids.lhs, target.promoter.other.ids.lhs)
        
      }
      
      # Store how many interactions overlap the Vista element in total (count unique interaction IDs).
      vista.results[i, paste0(cell.type, "_interactions")] <- length(unique(c(interaction.ids.lhs, interaction.ids.rhs)))
      
      # Translate target genes to gene name and store in results table.
      target.genes.all <- sort(unique(unlist(strsplit(paste0(c(target.genes.lhs, target.genes.rhs), collapse=","), split=","))))
      target.genes.all <- target.genes.all[target.genes.all != ""]
      target.genes.translated <- unique(promoters.all$gene_name[promoters.all$gene_id %in% target.genes.all])
      target.genes.translated <- target.genes.translated[target.genes.translated != ""]
      if (length(target.genes.translated) > 0) {
        vista.results[i, paste0(cell.type, "_target_genes")] <- paste0(sort(target.genes.translated), collapse=",")
      } else {
        vista.results[i, paste0(cell.type, "_target_genes")] <- ""
      }
      
    }
    
    # Need to resolve cases where the promoters of target genes are falling on the exact same fragment(s) as the Vista element.
    updated <- FALSE
    vista.entry.old <- vista.results[i, ]
    for (cell.type in cell.types) {
      
      # First verify the target genes for the current cell type agree with the promoters info.
      target.genes <- unlist(strsplit(vista.results[i, paste0(cell.type, "_target_genes")], split=","))
      promoter.id.genes <- unique(promoters.all.with.id$gene_name[promoters.all.with.id[, 8] %in% target.promoter.ids.cell.type[[cell.type]]])
      if (!((length(target.genes) == length(promoter.id.genes)) && all(target.genes %in% promoter.id.genes)))
        print(paste0("Inconsistency with entry: ", i, " in cell type: ", cell.type))
      
      # Now get all promoters in same fragment(s) as Vista element and filter out all target gene promoters that coincide.
      overlapping.bins.promoters.all.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.all.with.id, opt.string="-wb", num.cols=12)
      filtered <- target.promoter.ids.cell.type[[cell.type]][!(target.promoter.ids.cell.type[[cell.type]] %in% unique(as.character(overlapping.bins.promoters.all.results[, 12])))]
      new.target.genes <- unique(promoters.all.with.id$gene_name[promoters.all.with.id[, 8] %in% filtered])
      new.target.genes <- new.target.genes[new.target.genes != ""]
      if (!((length(target.genes) == length(new.target.genes)) && all(target.genes %in% new.target.genes))) {
        
        print(paste0("Entry ", i, " has been updated in cell type: ", cell.type))
        updated <- TRUE
        
      }
        
      # Update new entry.
      vista.results[i, paste0(cell.type, "_target_genes")] <- paste0(new.target.genes, collapse=",")
      
    }
    
    # Print diagnostic message.
    if (updated) {
      
      print(vista.entry.old)
      print(vista.results[i, ])
      
    }
    
    # Evaluate same fragment ambiguity. This is TRUE only if the nearest gene is also a same fragment gene AND ALL promoters of the nearest gene fall on the same fragment(s) as the Vista element.
    if (vista.results$nearest_gene[i] %in% unique(as.character(overlapping.bins.promoters.results[, 11]))) {
      
      # Check overlap between bins intersecting nearest gene promoter versus bins intersecting (expanded) Vista element.
      nearest.gene.promoters <- promoters.with.id[promoters.with.id$gene_name == vista.results$nearest_gene[i], ]
      nearest.gene.overlapping.bins.results <- bedTools.2in(bed1=nearest.gene.promoters, bed2=rmap, opt.string="-wb", num.cols=12)
      non.overlapping <- nearest.gene.overlapping.bins.results[, 12][!(nearest.gene.overlapping.bins.results[, 12] %in% overlapping.bins.promoters.results[, 4])]
      if (length(non.overlapping) > 0) {
        
        # This implies that an interaction is possible between a bin overlapping the Vista element and a non-overlapping bin coincident with a promoter of the nearest gene.
        vista.results$same_fragment_ambiguity[i] <- FALSE
        
      } else {
        
        # In some cases the nearest gene will still be reported as an interaction target due to the (expanded) Vista element covering many smaller bins. In this case, filter out the target genes.
        vista.results$same_fragment_ambiguity[i] <- TRUE
        for (cell.type in cell.types) {
          
          # May no longer be needed...
          temp <- unlist(strsplit(vista.results[i, paste0(cell.type, "_target_genes")], split=","))
          temp <- temp[temp != vista.results$nearest_gene[i]]
          vista.results[i, paste0(cell.type, "_target_genes")] <- paste0(temp, collapse=",")
          
        }
        
      }

    } else {
      
      # If the nearest gene is not on the same fragment(s) as the Vista element, then there is no same gene ambiguity (can't determine if interaction exists between Vista element and nearest gene).
      vista.results$same_fragment_ambiguity[i] <- FALSE
      
    }
    
  }
  
  # Return Vista element overlap results.
  return(vista.results)
  
}

# Analyze overlap between all GWAS SNPs, interactions, and gene promoters.
# Notes: Uses first promoter set (strict) for nearest/same fragment gene analysis, but uses second promoter set (all) for annotating target genes.
#
analyzeAllSNPOverlap <- function(disease.data, interaction.data, interaction.annotation, atac.seq.data, atac.seq.annotation, promoters, promoters.all, rmap, resolution) {
  
  # Print out information about disease.
  num.all.snps <- length(disease.data[, 1])
  print(paste0("Analyzing ", num.all.snps, " total SNPs."))
  cell.types <- names(interaction.data)
  tag.snp.rsids <- unlist(strsplit(paste0(disease.data$query_snp_rsid[disease.data$is_query_snp == TRUE], collapse=","), split=","))
  
  # Create results table for all SNPs.
  all.snp.results <- data.frame(matrix(NA, num.all.snps, 6 + length(cell.types)*4), stringsAsFactors=F)
  rownames(all.snp.results) <- disease.data$rsid
  colnames(all.snp.results) <- c("nearest_gene", "nearest_gene_distance", "same_fragment_genes", "same_fragment_other_genes", "same_fragment_ambiguity", "is_tag_snp",
                                 c(rbind(paste0(cell.types, "_promoter_atac_seq_peak"), paste0(cell.types, "_distal_atac_seq_peak"), paste0(cell.types, "_interactions"), paste0(cell.types, "_target_genes"))))
  all.snp.results$is_tag_snp <- disease.data$is_query_snp
  
  # Intersect all SNPs and promoters for reference later (to identify nearest gene).
  snps.sorted <- sortBed(unique(disease.data[, 1:4]))
  print(length(unique(snps.sorted$rsid)) == num.all.snps)
  promoters.sorted <- sortBed(promoters)
  promoters.sorted.filtered <- promoters.sorted[promoters.sorted[, 1] %in% unique(snps.sorted[, 1]), ]
  snps.promoters.results <- bedTools.2in(functionstring="bedtools closest", bed1=snps.sorted, bed2=promoters.sorted.filtered, opt.string="-t first -d")
  
  # Assign a promoter ID to each promoter for reference later. Also separate promoter sets.
  promoters.with.id <- cbind(promoters, paste0("promoter", 1:length(promoters[, 1])))
  promoters.other <- promoters.all[!(promoters.all$gene_type %in% promoters$gene_type), ]
  promoters.other.with.id <- cbind(promoters.other, paste0("promoter_other", 1:length(promoters.other[, 1])))
  promoters.all.with.id <- rbind(promoters.with.id, setNames(promoters.other.with.id, names(promoters.with.id)))
  
  # Expand SNPs to minimum resolution and intersect with interactions.
  snps.res <- expandFeatures(disease.data, resolution)
  interactions.snps.results.lhs <- list()
  interactions.snps.results.rhs <- list()
  for (cell.type in cell.types) {
    
    interactions.snps.results.lhs[[cell.type]] <- bedTools.2in(bed1=interaction.data[[cell.type]][, c(1:3, 11)], bed2=snps.res[, 1:4], opt.string="-wb", num.cols=8)
    interactions.snps.results.rhs[[cell.type]] <- bedTools.2in(bed1=interaction.data[[cell.type]][, c(5:7, 11)], bed2=snps.res[, 1:4], opt.string="-wb", num.cols=8)
    
  }
  
  # Intersect SNPs and ATAC-seq peaks.
  for (cell.type in cell.types) {
    
    # Split ATAC-seq peaks into promoter and distal peaks and expand the intervals by the SNP resolution.
    promoter.atac.seq.peaks <- atac.seq.data[[cell.type]][atac.seq.annotation[[cell.type]]$peak_type == "promoter", ]
    promoter.atac.seq.peaks <- expandFeaturesByWidth(promoter.atac.seq.peaks, resolution)
    distal.atac.seq.peaks <- atac.seq.data[[cell.type]][atac.seq.annotation[[cell.type]]$peak_type == "distal", ]
    distal.atac.seq.peaks <- expandFeaturesByWidth(distal.atac.seq.peaks, resolution)
    
    # Intersect ATAC-seq peaks with SNPs and fill in the results table.
    snps.promoter.peaks.results <- bedTools.2in(bed1=snps.res[, 1:4], bed2=promoter.atac.seq.peaks[, 1:4], opt.string="-c")
    snps.distal.peaks.results <- bedTools.2in(bed1=snps.res[, 1:4], bed2=distal.atac.seq.peaks[, 1:4], opt.string="-c")
    all.snp.results[, paste0(cell.type, "_promoter_atac_seq_peak")] <- snps.promoter.peaks.results[, 5] > 0
    all.snp.results[, paste0(cell.type, "_distal_atac_seq_peak")] <- snps.distal.peaks.results[, 5] > 0
    
  }
  
  # Intersect expanded SNPs with restriction fragments.
  rmap.snps.results <- bedTools.2in(bed1=rmap[, 1:4], bed2=snps.res[, 1:4], opt.string="-wb", num.cols=8)
  
  # Process all SNPs one by one.
  print(paste0(num.all.snps, " SNPs to process..."))
  for (i in 1:num.all.snps) {
    
    # Track progress.
    if (i %% 100 == 0)
      print(paste0(i, " SNPs processed"))
    
    # Get current rsid.
    current.rsid <- rownames(all.snp.results)[i]
    
    # Get nearest gene results.
    all.snp.results$nearest_gene[i] <- as.character(snps.promoters.results[as.character(snps.promoters.results[, 4]) == current.rsid, 11])
    all.snp.results$nearest_gene_distance[i] <- as.numeric(snps.promoters.results[as.character(snps.promoters.results[, 4]) == current.rsid, 12])
    
    # For each SNP, grab which restriction fragments it overlaps and which gene promoters (strict set) overlap those restriction fragments. This is to resolve same fragment ambiguity later.
    overlapping.bins <- rmap[rmap$rmap_ID %in% rmap.snps.results[as.character(rmap.snps.results[, 8]) == current.rsid, 4], ]
    overlapping.bins.promoters.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.with.id, opt.string="-wb", num.cols=12)
    all.snp.results$same_fragment_genes[i] <- paste0(unique(as.character(overlapping.bins.promoters.results[, 11])), collapse=",")
    
    # Use list of ALL gene promoters for determining same fragment genes.
    overlapping.bins.promoters.other.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.other, opt.string="-wb", num.cols=11)
    all.snp.results$same_fragment_other_genes[i] <- paste0(unique(as.character(overlapping.bins.promoters.other.results[, 11])), collapse=",")
    
    # For storing lists of target promoters across all cell types.
    target.promoter.ids.cell.type <- list()
    target.promoter.ids.all <- c()
    
    # Analyze interactions for each cell type.
    for (cell.type in cell.types) {
      
      # For storing lists of target promoters for each cell type.
      target.promoter.ids.cell.type[[cell.type]] <- c()
      
      # Identify the interactions whose lhs intersects the expanded SNP, then identify their target genes and gene promoters on the rhs.
      target.genes.rhs <- c()
      interaction.ids.lhs <- as.character(interactions.snps.results.lhs[[cell.type]][as.character(interactions.snps.results.lhs[[cell.type]][, 8]) == current.rsid, 4])
      num.interactions.lhs <- length(unique(interaction.ids.lhs))
      if (num.interactions.lhs > 0) {
        
        # Find which genes are on the other ends of the overlapping interactions.
        selector <- interaction.data[[cell.type]]$interaction_ID %in% interaction.ids.lhs
        target.genes.rhs <- c(as.character(interaction.annotation[[cell.type]]$promoter_rhs_ids[selector]), 
                              as.character(interaction.annotation[[cell.type]]$promoter_other_rhs_ids[selector]))
        
        # Determine the exact promoters interacting with the other ends of the overlapping interactions.
        interacting.ends.rhs <- interaction.data[[cell.type]][selector, 5:7]
        promoters.interacting.ends.results <- bedTools.2in(bed1=promoters.with.id, bed2=interacting.ends.rhs, opt.string="-c")
        promoters.other.interacting.ends.results <- bedTools.2in(bed1=promoters.other.with.id, bed2=interacting.ends.rhs, opt.string="-c")
        target.promoter.ids.rhs <- as.character(promoters.interacting.ends.results[promoters.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.all <- c(target.promoter.ids.all, target.promoter.ids.rhs)
        
        # Get specific promoters of both sets to determine what target genes should be removed later.
        target.promoter.other.ids.rhs <- as.character(promoters.other.interacting.ends.results[promoters.other.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.cell.type[[cell.type]] <- c(target.promoter.ids.cell.type[[cell.type]], target.promoter.ids.rhs, target.promoter.other.ids.rhs)
        
      }
      
      # Identify the interactions whose rhs intersects the expanded SNP, then identify their target genes and gene promoters on the lhs.
      target.genes.lhs <- c()
      interaction.ids.rhs <- as.character(interactions.snps.results.rhs[[cell.type]][as.character(interactions.snps.results.rhs[[cell.type]][, 8]) == current.rsid, 4])
      num.interactions.rhs <- length(unique(interaction.ids.rhs))
      if (num.interactions.rhs > 0) {
        
        # Find which genes are on the other ends of the overlapping interactions.
        selector <- interaction.data[[cell.type]]$interaction_ID %in% interaction.ids.rhs
        target.genes.lhs <- c(as.character(interaction.annotation[[cell.type]]$promoter_lhs_ids[selector]),
                              as.character(interaction.annotation[[cell.type]]$promoter_other_lhs_ids[selector]))
        
        # Determine the exact promoters interacting with the other ends of the overlapping interactions.
        interacting.ends.lhs <- interaction.data[[cell.type]][selector, 1:3]
        promoters.interacting.ends.results <- bedTools.2in(bed1=promoters.with.id, bed2=interacting.ends.lhs, opt.string="-c")
        promoters.other.interacting.ends.results <- bedTools.2in(bed1=promoters.other.with.id, bed2=interacting.ends.lhs, opt.string="-c")
        target.promoter.ids.lhs <- as.character(promoters.interacting.ends.results[promoters.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.all <- c(target.promoter.ids.all, target.promoter.ids.lhs)
        
        # Get specific promoters of both sets to determine what target genes should be removed later.
        target.promoter.other.ids.lhs <- as.character(promoters.other.interacting.ends.results[promoters.other.interacting.ends.results[, 9] > 0, 8])
        target.promoter.ids.cell.type[[cell.type]] <- c(target.promoter.ids.cell.type[[cell.type]], target.promoter.ids.lhs, target.promoter.other.ids.lhs)
        
      }
      
      # Store how many interactions overlap the Vista element in total (count unique interaction IDs).
      all.snp.results[i, paste0(cell.type, "_interactions")] <- length(unique(c(interaction.ids.lhs, interaction.ids.rhs)))
      
      # Translate target genes to gene name and store in results table.
      target.genes.all <- sort(unique(unlist(strsplit(paste0(c(target.genes.lhs, target.genes.rhs), collapse=","), split=","))))
      target.genes.all <- target.genes.all[target.genes.all != ""]
      target.genes.translated <- unique(promoters.all$gene_name[promoters.all$gene_id %in% target.genes.all])
      target.genes.translated <- target.genes.translated[target.genes.translated != ""]
      if (length(target.genes.translated) > 0) {
        all.snp.results[i, paste0(cell.type, "_target_genes")] <- paste0(sort(target.genes.translated), collapse=",")
      } else {
        all.snp.results[i, paste0(cell.type, "_target_genes")] <- ""
      }
      
    }
    
    # Need to resolve cases where the promoters of target genes are falling on the exact same fragment(s) as the SNP.
    updated <- FALSE
    snps.entry.old <- all.snp.results[i, ]
    for (cell.type in cell.types) {
      
      # First verify the target genes for the current cell type agree with the promoters info.
      target.genes <- unlist(strsplit(all.snp.results[i, paste0(cell.type, "_target_genes")], split=","))
      promoter.id.genes <- unique(promoters.all.with.id$gene_name[promoters.all.with.id[, 8] %in% target.promoter.ids.cell.type[[cell.type]]])
      if (!((length(target.genes) == length(promoter.id.genes)) && all(target.genes %in% promoter.id.genes)))
        print(paste0("Inconsistency with entry: ", i, " in cell type: ", cell.type))
      
      # Now get all promoters in same fragment(s) as SNP and filter out all target gene promoters that coincide.
      overlapping.bins.promoters.all.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.all.with.id, opt.string="-wb", num.cols=12)
      filtered <- target.promoter.ids.cell.type[[cell.type]][!(target.promoter.ids.cell.type[[cell.type]] %in% unique(as.character(overlapping.bins.promoters.all.results[, 12])))]
      new.target.genes <- unique(promoters.all.with.id$gene_name[promoters.all.with.id[, 8] %in% filtered])
      new.target.genes <- new.target.genes[new.target.genes != ""]
      if (!((length(target.genes) == length(new.target.genes)) && all(target.genes %in% new.target.genes))) {
        
        print(paste0("Entry ", i, " has been updated in cell type: ", cell.type))
        updated <- TRUE
        
      }
      
      # Update new entry.
      all.snp.results[i, paste0(cell.type, "_target_genes")] <- paste0(new.target.genes, collapse=",")
      
    }
    
    # Print diagnostic message.
    if (updated) {
      
      print(snps.entry.old)
      print(all.snp.results[i, ])
      
    }
    
    # Evaluate same fragment ambiguity. This is TRUE only if the nearest gene is also a same fragment gene AND ALL promoters of the nearest gene fall on the same fragment(s) as the SNP.
    if (all.snp.results$nearest_gene[i] %in% unique(as.character(overlapping.bins.promoters.results[, 11]))) {
      
      # Check overlap between bins intersecting nearest gene promoter versus bins intersecting (expanded) SNP.
      nearest.gene.promoters <- promoters.with.id[promoters.with.id$gene_name == all.snp.results$nearest_gene[i], ]
      nearest.gene.overlapping.bins.results <- bedTools.2in(bed1=nearest.gene.promoters, bed2=rmap, opt.string="-wb", num.cols=12)
      non.overlapping <- nearest.gene.overlapping.bins.results[, 12][!(nearest.gene.overlapping.bins.results[, 12] %in% overlapping.bins.promoters.results[, 4])]
      if (length(non.overlapping) > 0) {
        
        # This implies that an interaction is possible between a bin overlapping the Vista element and a non-overlapping bin coincident with a promoter of the nearest gene.
        all.snp.results$same_fragment_ambiguity[i] <- FALSE
        
      } else {
        
        # In some cases the nearest gene will still be reported as an interaction target due to the (expanded) SNP covering many smaller bins. In this case, filter out the target genes.
        all.snp.results$same_fragment_ambiguity[i] <- TRUE
        for (cell.type in cell.types) {
          
          # May no longer be needed...
          temp <- unlist(strsplit(all.snp.results[i, paste0(cell.type, "_target_genes")], split=","))
          temp <- temp[temp != all.snp.results$nearest_gene[i]]
          all.snp.results[i, paste0(cell.type, "_target_genes")] <- paste0(temp, collapse=",")
          
        }
        
      }
      
    } else {
      
      # If the nearest gene is not on the same fragment(s) as the SNP, then there is no same gene ambiguity (can't determine if interaction exists between SNP and nearest gene).
      all.snp.results$same_fragment_ambiguity[i] <- FALSE
      
    }
    
  }
  
  # Return SNP annotations.
  return(all.snp.results)
  
}

# Analyze overlap between tag GWAS SNPs, interactions, and gene promoters.
# Notes: Uses first promoter set (strict) for nearest/same fragment gene analysis, but uses second promoter set (all) for annotating target genes.
#
analyzeTagSNPOverlap <- function(disease.data, all.snp.results, promoters, promoters.all, rmap, resolution) {
  
  # Print out information about disease.
  num.all.snps <- length(disease.data[, 1])
  tag.snp.data <- disease.data[disease.data$is_query_snp == TRUE, ]
  num.tag.snps <- length(tag.snp.data[, 1])
  print(paste0("Analyzing ", num.tag.snps, " tag SNPs from ", num.all.snps, " total SNPs."))
  tag.snps.ids <- unique(unlist(strsplit(paste0(disease.data$query_snp_rsid, collapse=","), split=",")))
  all.snps.ids <- unique(unlist(strsplit(paste0(disease.data$rsid, collapse=","), split=",")))
  print(table(tag.snps.ids %in% all.snps.ids)) # Should be TRUE.
  cell.types <- gsub("_interactions", "", grep("_interactions", colnames(all.snp.results), value=TRUE))

  # Should be the same (no duplicate rsids in entries). CHECK
  print(length(unlist(strsplit(paste0(disease.data$rsid, collapse=","), split=","))))
  print(length(unique(unlist(strsplit(paste0(disease.data$rsid, collapse=","), split=",")))))
  
  # Each position/entry might have multiple rsids/query rsids. Split these and create a mapping of individual rsids and indices.
  indices.by.rsid <- c()
  # indices.by.query.rsid <- c()
  for (i in 1:length(disease.data[, 1])) {
    
    # Process rsids.
    individual.rsids <- unlist(strsplit(disease.data$rsid[i], split=","))
    # if (length(individual.rsids) > 1)
    #   print(i)
    indices.by.rsid <- rbind(indices.by.rsid, cbind(individual.rsids, rep(i, length(individual.rsids))))
    
    # Process query rsids.
    # individual.query.rsids <- unlist(strsplit(disease.data$query_snp_rsid[i], split=","))
    # indices.by.query.rsid <- rbind(indices.by.query.rsid, cbind(individual.query.rsids, rep(i, length(individual.query.rsids))))
    
  }
  
  # Create a translation table between all SNPs and query SNPs. One entry per SNP to query SNP mapping.
  query.translation <- c()
  for (i in 1:num.all.snps) {
    
    current.rsids <- unlist(strsplit(disease.data$rsid[i], split=","))
    current.query.rsids <- unlist(strsplit(disease.data$query_snp_rsid[i], split=","))
    query.translation <- rbind(query.translation, expand.grid(current.rsids, current.query.rsids))
    
  }
  query.translation <- unique(query.translation) # Not technically needed.
  colnames(query.translation) <- c("rsid", "query_snp_rsid")
  
  # Create results table for tag SNPs.
  tag.snp.results <- data.frame(matrix(NA, num.tag.snps, 6 + length(cell.types)*4), stringsAsFactors=F)
  rownames(tag.snp.results) <- tag.snp.data$rsid
  colnames(tag.snp.results) <- c("nearest_gene", "nearest_gene_distance", "same_fragment_genes", "same_fragment_other_genes", "same_fragment_ambiguity", "num_linked_snps",
                                 c(rbind(paste0(cell.types, "_interactions"), paste0(cell.types, "_target_genes"), paste0(cell.types, "_interactions_imputed"), paste0(cell.types, "_target_genes_imputed"))))
  
  # Assign a promoter ID to each promoter for reference later. Also separate promoter sets.
  promoters.with.id <- cbind(promoters, paste0("promoter", 1:length(promoters[, 1])))
  promoters.other <- promoters.all[!(promoters.all$gene_type %in% promoters$gene_type), ]
  promoters.other.with.id <- cbind(promoters.other, paste0("promoter_other", 1:length(promoters.other[, 1])))
  promoters.all.with.id <- rbind(promoters.with.id, setNames(promoters.other.with.id, names(promoters.with.id)))
  
  # Expand SNPs to minimum resolution and intersect expanded SNPs with restriction fragments.
  snps.res <- expandFeatures(disease.data, resolution)
  rmap.snps.results <- bedTools.2in(bed1=rmap[, 1:4], bed2=snps.res[, 1:4], opt.string="-wb", num.cols=8)
  
  # Process tag SNPs one by one.
  print(paste0(num.tag.snps, " tag SNPs to process..."))
  for (i in 1:num.tag.snps) {
    
    # Track progress.
    if (i %% 100 == 0)
      print(paste0(i, " tag SNPs processed"))
    
    # Get current rsid.
    current.query.rsid <- as.character(unlist(strsplit(rownames(tag.snp.results)[i], split=",")))
    if (length(current.query.rsid) > 1)
      print(paste0("Detected entry with multiple rsids, at least one of which is a query SNP: ", i))
    
    # Find all rsids for which the current rsid is the query SNP. Count unique positions as opposed to unique rsids as multiple rsids may be at the same position.
    imputed.rsids <- as.character(query.translation[query.translation[, 2] %in% current.query.rsid, 1])
    # tag.snp.results$num_linked_snps[i] <- length(imputed.rsids)
    tag.snp.results$num_linked_snps[i] <- length(unique(indices.by.rsid[indices.by.rsid[, 1] %in% imputed.rsids, 2]))

    # Fill in information for the query SNP.
    query.snp.index <- unique(as.integer(indices.by.rsid[indices.by.rsid[, 1] %in% current.query.rsid, 2]))
    if (length(query.snp.index) != 1)
      print(paste0("Cannot find unambiguous information for query SNP: ", i))
    tag.snp.results[i, 1:4] <- all.snp.results[query.snp.index, 1:4]
    if (!all.snp.results$is_tag_snp[query.snp.index])
      print("Inconsistency in identification of tag SNP")
    for (cell.type in cell.types) {
      tag.snp.results[i, paste0(cell.type, "_interactions")] <- all.snp.results[query.snp.index, paste0(cell.type, "_interactions")]
      tag.snp.results[i, paste0(cell.type, "_target_genes")] <- all.snp.results[query.snp.index, paste0(cell.type, "_target_genes")]
    }
    
    # Fill in information for all imputed SNPs.
    imputed.snp.indices <- as.integer(indices.by.rsid[indices.by.rsid[, 1] %in% imputed.rsids, 2])
    if (!(query.snp.index %in% imputed.snp.indices))
      print("Query SNP index is not in list of linked SNP indices")
    imputed.records <- all.snp.results[imputed.snp.indices, ]
    for (cell.type in cell.types) {
      
      tag.snp.results[i, paste0(cell.type, "_interactions_imputed")] <- any(all.snp.results[imputed.snp.indices, paste0(cell.type, "_interactions")] > 0)
      target.genes.all <- unique(unlist(strsplit(paste0(all.snp.results[imputed.snp.indices, paste0(cell.type, "_target_genes")], collapse=","), split=",")))
      target.genes.all <- target.genes.all[target.genes.all != ""]
      if (length(target.genes.all) > 0) {
        tag.snp.results[i, paste0(cell.type, "_target_genes_imputed")] <- paste0(target.genes.all, collapse=",")
      } else {
        tag.snp.results[i, paste0(cell.type, "_target_genes_imputed")] <- ""
      }
      
    }
    
    # Calculate same fragment ambiguity. First, if the tag SNP is only linked to itself and no other SNPs, then just use the individual SNP results.
    if (tag.snp.results$num_linked_snps[i] == 1) {
      
      tag.snp.results$same_fragment_ambiguity[i] <- all.snp.results$same_fragment_ambiguity[query.snp.index]

    } else {
      
      # Second, if any of the linked SNPs do not have the nearest gene as a same fragment gene, then there is no ambiguity. Set 'match' to FALSE if nearest gene is not a same fragment gene for any linked SNP.
      same.fragment.genes.all <- all.snp.results$same_fragment_genes[imputed.snp.indices]
      match <- TRUE
      for (j in 1:length(same.fragment.genes.all)) {
        
        if (!(tag.snp.results$nearest_gene[i] %in% unlist(strsplit(same.fragment.genes.all[j], split=","))))
          match <- FALSE
        
      }
      
      # Final case if 'match' is FALSE is when there is >1 linked SNP and all SNPs are on the same fragment(s) as at least one of the nearest gene's promoters. Then we have to evaluate manually.
      if (match == FALSE ) {
        
        tag.snp.results$same_fragment_ambiguity[i] <- FALSE
        
      } else {
        
        # Diagnostic message.
        # print(paste0("Manually evaluating tag SNP entry: ", i))
        
        # Get all the bins overlapping promoters of the nearest gene.
        nearest.gene.promoters <- promoters.with.id[promoters.with.id$gene_name == tag.snp.results$nearest_gene[i], ]
        nearest.gene.overlapping.bins.results <- bedTools.2in(bed1=nearest.gene.promoters, bed2=rmap, opt.string="-wb", num.cols=12)
        
        # Now process the linked SNPs one by one--if for any of them, they do not overlap one bin overlapped by a nearest gene promoter, then there is no same fragment ambiguity.
        possible.interaction <- FALSE
        for (j in 1:length(imputed.snp.indices)) {
          
          # Now get all bins overlapping any of the linked SNPs.
          overlapping.bins <- rmap[rmap$rmap_ID %in% rmap.snps.results[as.character(rmap.snps.results[, 8]) %in% disease.data$rsid[imputed.snp.indices[j]], 4], ]
          overlapping.bins.promoters.results <- bedTools.2in(bed1=overlapping.bins, bed2=promoters.with.id, opt.string="-wb", num.cols=12)
          
          # Now see if there are any bins overlapping promoters of the nearest gene which are also NOT overlapping any linked SNPs. If this is true, then there is no ssame fragment ambiguity.
          non.overlapping <- nearest.gene.overlapping.bins.results[, 12][!(nearest.gene.overlapping.bins.results[, 12] %in% overlapping.bins.promoters.results[, 4])]
          if (length(non.overlapping) > 0)
            possible.interaction <- TRUE

        }
        
        # Set the evaluated same fragment ambiguity.
        if (possible.interaction) {
          
          tag.snp.results$same_fragment_ambiguity[i] <- FALSE
          
        } else {
          
          # Quickly check to make sure the nearest gene is not a target gene, since that would be impossible given there is same fragment ambiguity.
          # print(paste0(i, ": TRUE"))
          tag.snp.results$same_fragment_ambiguity[i] <- TRUE
          if (tag.snp.results$nearest_gene[i] %in% unlist(strsplit(paste0(tag.snp.results[i, paste0(cell.types, "_target_genes_imputed")], collapse=","), split=",")))
            print(paste0(i, ": ERROR"))
          
        }
        
      }
      
    }

  }
  
  # Return tag SNP annotations.
  return(tag.snp.results)
  
}

# Combine tag SNP results across all diseases.
# Notes:
#
combineTagSNPResults <- function(tag.snp.results, features, diseases) {
  
  # Merge the data and annotation for each tag SNP across all diseases.
  cell.types <- names(features)
  tag.snp.data.merged <- c()
  tag.snp.results.merged <- c()
  for (current.disease in diseases) {
    
    # Grab the data and annotation for the current disease and concatenate both into their merged records.
    disease.data <- features[[cell.types[1]]][["snp"]][[current.disease]]
    tag.snp.data <- disease.data[disease.data$is_query_snp == TRUE, ]
    tag.snp.data.merged <- rbind(tag.snp.data.merged, tag.snp.data)
    tag.snp.results.merged <- rbind(tag.snp.results.merged, cbind(tag.snp.results[[current.disease]], tag.snp.data$feature_ID))
    
    # Print checks (should all be TRUE).
    # print(length(tag.snp.data[, 1]) == length(tag.snp.results[[current.disease]][, 1]))
    # print(table(rownames(tag.snp.results[[current.disease]]) == tag.snp.data$rsid))
    # print(length(unique(tag.snp.data[, 1:3])[, 1]) == length(tag.snp.results[[current.disease]][, 1]))
    
  }
  
  # Split the duplicated and non-duplicated tag SNP records.
  if (length(unique(tag.snp.data.merged[, 1:3])[, 1]) != length(unique(tag.snp.data.merged[, 1:4])[, 1]))
    print("Tag SNPs detected with same positions but different rsids...")
  print(table(duplicated(tag.snp.data.merged[, 1:3]))[["TRUE"]] == table(duplicated(tag.snp.data.merged[, 4]))[["TRUE"]])
  print(paste0("There are a total of ", length(unique(tag.snp.data.merged[, 4])), " unique tag SNP rsids."))
  duplicated.indices <- duplicated(tag.snp.data.merged[, 4], fromLast=F) | duplicated(tag.snp.data.merged[, 4], fromLast=T)
  print(table(duplicated.indices))
  non.duplicated.tag.snp.data.merged <- tag.snp.data.merged[!duplicated.indices, ]
  non.duplicated.tag.snp.results.merged <- tag.snp.results.merged[!duplicated.indices, ]
  duplicated.tag.snp.data.merged <- tag.snp.data.merged[duplicated.indices, ]
  duplicated.tag.snp.results.merged <- tag.snp.results.merged[duplicated.indices, ]
  print(paste0("There are a total of ", length(unique(duplicated.tag.snp.data.merged[, 4])), " duplicated tag SNP rsids."))
  print(paste0("There are a total of ", length(duplicated.tag.snp.data.merged[, 4]), " records associated with duplicated tag SNP rsids."))
  print(table(duplicated.tag.snp.data.merged[, 4] %in% non.duplicated.tag.snp.data.merged[, 4]))
  
  # Reconcile the duplicated records.
  duplicated.rsids <- unique(duplicated.tag.snp.data.merged[, 4])
  num.duplicated.rsids <- length(duplicated.rsids)
  reconciled.data <- c()
  reconciled.results <- c()
  for (i in 1:num.duplicated.rsids) {
    
    # Get current duplicated rsid.
    current.indices <- which(duplicated.tag.snp.data.merged[, 4] == duplicated.rsids[i])
    if (length(current.indices) < 2)
      print(paste0("Error 1: ", i))

    # Reconcile tag SNP data.
    # duplicated.tag.snp.data.merged[current.indices, ]
    temp <- duplicated.tag.snp.data.merged[current.indices[1], ]
    temp$query_snp_rsid <- paste0(sort(unique(unlist(strsplit(paste0(duplicated.tag.snp.data.merged$query_snp_rsid[current.indices], collapse=","), split=",")))), collapse=",")
    temp$feature_ID <- paste0(sort(unique(unlist(strsplit(paste0(duplicated.tag.snp.data.merged$feature_ID[current.indices], collapse=","), split=",")))), collapse=",")
    reconciled.data <- rbind(reconciled.data, temp)
    
    # Reconcile tag SNP results.
    # duplicated.tag.snp.results.merged[current.indices, ]
    if (dim(unique(duplicated.tag.snp.results.merged[current.indices, c(1:5, 7:8, 11:12, 15:16, 19:20)]))[1] != 1)
      print(paste0("Error 2: ", i))
    temp <- duplicated.tag.snp.results.merged[current.indices[1], ]
    temp$num_linked_snps <- max(duplicated.tag.snp.results.merged$num_linked_snps[current.indices])
    temp$`tag.snp.data$feature_ID` <- paste0(sort(unique(unlist(strsplit(paste0(duplicated.tag.snp.results.merged$`tag.snp.data$feature_ID`[current.indices], collapse=","), split=",")))), collapse=",")
    for (cell.type in cell.types) {
      
      temp[paste0(cell.type, "_interactions_imputed")] <- any(duplicated.tag.snp.results.merged[current.indices, paste0(cell.type, "_interactions_imputed")])
      temp[paste0(cell.type, "_target_genes_imputed")] <- paste0(sort(unique(unlist(strsplit(paste0(duplicated.tag.snp.results.merged[current.indices, paste0(cell.type, "_target_genes_imputed")], collapse=","), split=",")))), collapse=",")
      
    }
    reconciled.results <- rbind(reconciled.results, temp)
    
  }
  
  # Merge unique and reconciled results.
  print(table(reconciled.data$feature_ID == reconciled.results$`tag.snp.data$feature_ID`))
  merged.results <- rbind(non.duplicated.tag.snp.results.merged, reconciled.results)
  print(dim(merged.results))
  
  # Return combined tag SNP annotations.
  return(merged.results)
  
}

# Return a list of interaction IDs matching a list of specificity patterns.
# Notes:
#
getCellTypeSpecificInteractions <- function(specificity.entries, specificity.patterns, cell.types.filter, cell.types) {
  
  # Keep only interactions in the specified cell types.
  keep <- c()
  for (cell.type.filter in cell.types.filter)
    keep <- rbind(keep, specificity.entries[grepl(cell.type.filter, specificity.entries[, 1]), ])
  specificity.entries.filter <- keep
  
  # Keep only interaction IDs matching the specificity patterns.
  keep <- c()
  for (specificity.pattern in specificity.patterns)
    keep <- c(keep, specificity.entries.filter[specificity.entries.filter[, 2] == specificity.pattern, 1])
  
  # Return filtered interaction IDs.
  return(keep)
  
}

# Retrieves genes participating in interactions across all cell types with the provided specificity pattern(s).
# Notes: 'mode' can be "all", "other", or "distal_atac". 'specificity.patterns' should be sorted in the same order as the 'cell.types' variable.
#
getSpecificInteractingGenes <- function(interactions.sig, interactions.ann, atac.seq.peaks.ann, specificity.entries, specificity.patterns, promoters, expression.data, expression.threshold, mode="distal_atac") {
  
  # Get the interacting genes matching interactions with the provided specificity pattern(s) in each cell type.
  interacting.genes <- c()
  cell.types <- names(interactions.sig)
  for (cell.type in cell.types) {
    
    # Get the IDs of the interactions matching the provided specificity pattern(s).
    specificity.entries.cell.type <- specificity.entries[grepl(cell.type, specificity.entries[, 1]), ]
    matching.ids <- specificity.entries.cell.type[specificity.entries.cell.type[, 2] %in% specificity.patterns, 1]
    selector <- interactions.sig[[cell.type]]$interaction_ID %in% matching.ids
    
    # Get the interacting genes depending on the type of interactions we are selecting for (all interacting genes, only promoter to other interactions, or only promoter to distal ATAC-seq peak interactions).
    if (mode == "all") {
      
      current.genes <- c(interactions.ann[[cell.type]]$promoter_lhs_ids, interactions.ann[[cell.type]]$promoter_rhs_ids)

    } else if (mode == "other") {
      
      promoter.lhs.selector <- ((interactions.ann[[cell.type]]$promoter_lhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_rhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] == 0)
      promoter.rhs.selector <- ((interactions.ann[[cell.type]]$promoter_rhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_lhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] == 0)
      current.genes <- c(interactions.ann[[cell.type]]$promoter_lhs_ids[selector][promoter.lhs.selector], interactions.ann[[cell.type]]$promoter_rhs_ids[selector][promoter.rhs.selector])
      
    } else if (mode == "distal_atac") {
      
      promoter.lhs.distal.rhs.selector <- ((interactions.ann[[cell.type]]$promoter_lhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_rhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] == 0) & 
        (interactions.ann[[cell.type]]$`distal_ATAC-seq_rhs`[selector] > 0)
      promoter.rhs.distal.lhs.selector <- ((interactions.ann[[cell.type]]$promoter_rhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_lhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] == 0) & 
        (interactions.ann[[cell.type]]$`distal_ATAC-seq_lhs`[selector] > 0)
      current.genes <- c(interactions.ann[[cell.type]]$promoter_lhs_ids[selector][promoter.lhs.distal.rhs.selector], interactions.ann[[cell.type]]$promoter_rhs_ids[selector][promoter.rhs.distal.lhs.selector])
      
    } else if (mode == "both_atac") {
      
      promoter.lhs.distal.rhs.selector <- ((interactions.ann[[cell.type]]$promoter_lhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_rhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] == 0) & 
        (interactions.ann[[cell.type]]$`distal_ATAC-seq_rhs`[selector] > 0)
      promoter.rhs.distal.lhs.selector <- ((interactions.ann[[cell.type]]$promoter_rhs[selector] > 0) | (interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs`[selector] > 0)) & 
        (interactions.ann[[cell.type]]$promoter_lhs[selector] == 0) & (interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs`[selector] == 0) & 
        (interactions.ann[[cell.type]]$`distal_ATAC-seq_lhs`[selector] > 0)
      current.peaks <- c(interactions.ann[[cell.type]]$`promoter_ATAC-seq_lhs_ids`[selector][promoter.lhs.distal.rhs.selector],
                         interactions.ann[[cell.type]]$`promoter_ATAC-seq_rhs_ids`[selector][promoter.rhs.distal.lhs.selector])
      current.genes <- atac.seq.peaks.ann[[cell.type]][current.peaks, "promoter_genes"]
      
    }
    
    #
    current.genes <- unique(unlist(strsplit(paste0(current.genes, collapse=","), split=",")))
    current.genes <- current.genes[current.genes != ""]
    current.genes <- current.genes[expression.data[current.genes, cell.type] > expression.threshold]
    interacting.genes <- c(interacting.genes, current.genes)
    
  }
  
  #
  interacting.genes <- unique(promoters$gene_name[promoters$gene_id %in% interacting.genes])
  print(paste0(length(interacting.genes), " identified."))
  
  #
  return(interacting.genes)
  
}


# End ---------------------------------------------------------------------

