#
# Inputs: Peaks, interactions, and promoter region coordinates.
# Outputs: Data frame with annotation of peak as promoter or distal, gene names associated with promoters, # of interactions intersecting with peak, and gene names associated with distal promoters.
annotateATACPeaks <- function(peaks, peaks.unexpanded, interactions, promoters, introns, exons) {
  
  # Create data frame to hold annotation results.
  annotation <- data.frame(matrix("", length(peaks[, 1]), 7), stringsAsFactors=F)
  rownames(annotation) <- peaks$ID
  colnames(annotation) <- c("annotation", "promoter_genes", "num_interactions", "proximal_genes", "distal_genes", "introns", "exons")
  annotation[, 1] <- as.character(annotation[, 1]) # Promoter or distal.
  annotation[, 2] <- as.character(annotation[, 2]) # Genes whose promoter regions intersect peak.
  annotation[, 3] <- as.integer(annotation[, 3]) # Number of interactions that intersect peak.
  annotation[, 4] <- as.character(annotation[, 4]) # Genes in interaction end intersecting peak.
  annotation[, 5] <- as.character(annotation[, 5]) # Genes in interaction other end not intersecting peak.
  annotation[, 6] <- as.character(annotation[, 6])
  annotation[, 7] <- as.character(annotation[, 7])
  
  # Intersect interaction lhs/rhs coordinates with promoters.
  interaction.genes <- data.frame(matrix("", length(interactions[, 1]), 2), stringsAsFactors=F)
  rownames(interaction.genes) <- interactions[, 11]
  colnames(interaction.genes) <- c("lhs_genes", "rhs_genes")
  interaction.genes[, 1] <- as.character(interaction.genes[, 1])
  interaction.genes[, 2] <- as.character(interaction.genes[, 2])
  lhs.interaction.promoters <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  rhs.interaction.promoters <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  for (i in 1:length(interactions[, 1])) {
    
    # Grab current interaction ID.
    current.interaction.id <- as.character(interactions[i, 11])
    
    # Get gene IDs of all genes whose promoters intersect lhs and rhs of interaction and fill into table.
    lhs.subset.gene.ids <- as.character(lhs.interaction.promoters[as.character(lhs.interaction.promoters[, 4]) == current.interaction.id, 8])
    rhs.subset.gene.ids <- as.character(rhs.interaction.promoters[as.character(rhs.interaction.promoters[, 4]) == current.interaction.id, 8])
    if (length(lhs.subset.gene.ids) > 0) {
      interaction.genes[i, 1] <- paste0(unique(lhs.subset.gene.ids), collapse=",")
    }
    if (length(rhs.subset.gene.ids) > 0) {
      interaction.genes[i, 2] <- paste0(unique(rhs.subset.gene.ids), collapse=",")
    }
    
  }
  
  # Intersect peak coordinates with promoter/interaction coordinates.
  promoter.intersection <- bedTools.2in(bed1=peaks[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  intron.intersection <- bedTools.2in(bed1=peaks.unexpanded[, c(1:3)], bed2=introns[, c(1:3)], opt.string="-c")
  exon.intersection <- bedTools.2in(bed1=peaks.unexpanded[, c(1:3)], bed2=exons[, c(1:3)], opt.string="-c")
  lhs.interaction.peaks <- bedTools.2in(bed1=peaks[, c(1:3, 11)], bed2=interactions[, c(1:3, 11)], opt.string="-wb")
  rhs.interaction.peaks <- bedTools.2in(bed1=peaks[, c(1:3, 11)], bed2=interactions[, c(5:7, 11)], opt.string="-wb")
  print(paste0(length(peaks[, 1]), " peaks to annotate..."))
  for (i in 1:length(peaks[, 1])) {
    
    if (i %% 10000 == 0)
      print(i)
    
    # Process intron/exon overlaps first.
    if (intron.intersection[i, 4] > 0) {
      annotation[i, 6] <- TRUE
    } else {
      annotation[i, 6] <- FALSE
    }
    if (exon.intersection[i, 4] > 0) {
      annotation[i, 7] <- TRUE
    } else {
      annotation[i, 7] <- FALSE
    } 
    
    # Grab current peak ID.
    current.peak.id <- as.character(peaks[i, 11])
    
    # First determine whether or not the peak intersects promoter region(s), and if so grab the gene IDs of all associated genes.
    promoter.subset.gene.ids <- as.character(promoter.intersection[as.character(promoter.intersection[, 4]) == current.peak.id, 8])
    if (length(promoter.subset.gene.ids) > 0) {
      annotation[i, 1] <- "promoter"
      annotation[i, 2] <- paste0(unique(promoter.subset.gene.ids), collapse=",")
    } else {
      annotation[i, 1] <- "distal"
    }
    
    # Next determine how many interactions intersect peak, and grab proximal and distal genes associated with interaction (if any).
    lhs.subset.interaction.ids <- as.character(lhs.interaction.peaks[as.character(lhs.interaction.peaks[, 4]) == current.peak.id, 8])
    rhs.subset.interaction.ids <- as.character(rhs.interaction.peaks[as.character(rhs.interaction.peaks[, 4]) == current.peak.id, 8])
    unique.interaction.ids <- unique(c(lhs.subset.interaction.ids, rhs.subset.interaction.ids))
    annotation[i, 3] <- length(unique.interaction.ids)
    if (annotation[i, 3] > 0) {
      
      # Go through each interaction intersecting peak and fill in the proximal and distal genes according to the table we made earlier.
      for (j in 1:length(unique.interaction.ids)) {
        
        retrieved.gene.ids <- interaction.genes[unique.interaction.ids[j], ]
        in.lhs <- unique.interaction.ids[j] %in% lhs.subset.interaction.ids
        in.rhs <- unique.interaction.ids[j] %in% rhs.subset.interaction.ids
        if (in.lhs & in.rhs) {
          
          annotation[i, 4] <- paste0(c(annotation[i, 4], "CO"), collapse=",")
          annotation[i, 5] <- paste0(c(annotation[i, 5], "CO"), collapse=",")
          
        } else if (in.lhs) {
          
          annotation[i, 4] <- paste0(c(annotation[i, 4], retrieved.gene.ids[1]), collapse=",")
          annotation[i, 5] <- paste0(c(annotation[i, 5], retrieved.gene.ids[2]), collapse=",")
          
        } else if (in.rhs) {
          
          annotation[i, 4] <- paste0(c(annotation[i, 4], retrieved.gene.ids[2]), collapse=",")
          annotation[i, 5] <- paste0(c(annotation[i, 5], retrieved.gene.ids[1]), collapse=",")
          
        }
        
        # Clean up duplicate gene IDs.
        temp.proximal <- unique(unlist(strsplit(annotation[i, 4], split=",")))
        temp.proximal <- temp.proximal[temp.proximal != ""]
        annotation[i, 4] <- paste0(temp.proximal, collapse=",")
        temp.distal <- unique(unlist(strsplit(annotation[i, 5], split=",")))
        temp.distal <- temp.distal[temp.distal != ""]
        annotation[i, 5] <- paste0(temp.distal, collapse=",")
        
      }
      
    }
    
  }
  
  # Return annotation results.
  return(annotation)
  
}

#
# Inputs:
# Outputs:
annotateInteractions <- function(interactions, promoters, atac.peaks, atac.annotations, features=NULL) {
  
  #
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
  
  # Annotations will be stored in a list organized by feature names.
  annotation <- list()
  num.interactions <- length(interactions[, 1])
  
  # Process features if they are specified.
  if (!is.null(features)) {
    
    # Get list of all features to be used for annotation.  Table keeps track of unique pairs of feature types + feature names to be iterated over.
    feature.classes <- names(features)
    feature.pairs <- c()
    for (i in 1:length(feature.classes)) {
      
      feature.names <- names(features[[feature.classes[i]]])
      feature.pairs <- rbind(feature.pairs, cbind(rep(feature.classes[i], length(feature.names)), feature.names))
      
    }
    
    # Perform intersections to determine whether the lhs and rhs coordinates of interactions intersect each feature.
    lhs.feature.results <- list()
    rhs.feature.results <- list()
    for (i in 1:length(feature.pairs[, 1])) {
      
      lhs.feature.results[[feature.pairs[i, 2]]] <- bedTools.2in(bed1=interactions[, c(1:3, 11)], 
                                                                 bed2=cbind(features[[feature.pairs[i, 1]]][[feature.pairs[i, 2]]][, c(1:3)], features[[feature.pairs[i, 1]]][[feature.pairs[i, 2]]]$ID), 
                                                                 opt.string="-wb")
      rhs.feature.results[[feature.pairs[i, 2]]] <- bedTools.2in(bed1=interactions[, c(5:7, 11)], 
                                                                 bed2=cbind(features[[feature.pairs[i, 1]]][[feature.pairs[i, 2]]][, c(1:3)], features[[feature.pairs[i, 1]]][[feature.pairs[i, 2]]]$ID),
                                                                 opt.string="-wb")
      annotation[[paste0("lhs_", feature.pairs[i, 2])]] <- rep(0, num.interactions)
      annotation[[paste0("lhs_", feature.pairs[i, 2], "_ids")]] <- rep("", num.interactions)
      annotation[[paste0("rhs_", feature.pairs[i, 2])]] <- rep(0, num.interactions)
      annotation[[paste0("rhs_", feature.pairs[i, 2], "_ids")]] <- rep("", num.interactions)
      
    }
    
  }
  
  # Perform intersections to determine whether the lhs and rhs coordinates of interactions intersect promoters regions/ATAC-seq peaks.
  lhs.promoter.results <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  lhs.atac.results <- bedTools.2in(bed1=interactions[, c(1:3, 11)], bed2=atac.peaks[, c(1:3, 11)], opt.string="-wb")
  rhs.promoter.results <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=promoters[, c(1:3, 5)], opt.string="-wb")
  rhs.atac.results <- bedTools.2in(bed1=interactions[, c(5:7, 11)], bed2=atac.peaks[, c(1:3, 11)], opt.string="-wb")
  annotation[["lhs_promoter"]] <- rep(0, num.interactions)
  annotation[["lhs_promoter_ids"]] <- rep("", num.interactions)
  annotation[["lhs_promoter_ATAC-seq"]] <- rep(0, num.interactions)
  annotation[["lhs_promoter_ATAC-seq_ids"]] <- rep("", num.interactions)
  annotation[["lhs_distal_ATAC-seq"]] <- rep(0, num.interactions)
  annotation[["lhs_distal_ATAC-seq_ids"]] <- rep("", num.interactions)
  annotation[["rhs_promoter"]] <- rep(0, num.interactions)
  annotation[["rhs_promoter_ids"]] <- rep("", num.interactions)
  annotation[["rhs_promoter_ATAC-seq"]] <- rep(0, num.interactions)
  annotation[["rhs_promoter_ATAC-seq_ids"]] <- rep("", num.interactions)
  annotation[["rhs_distal_ATAC-seq"]] <- rep(0, num.interactions)
  annotation[["rhs_distal_ATAC-seq_ids"]] <- rep("", num.interactions)
  
  # Annotate interactions with intersection results.
  print(paste0(num.interactions, " interactions to annotate..."))
  for (i in 1:num.interactions) {
    
    if (i %% 5000 == 0)
      print(i)
    
    # Grab current interaction ID to search intersection results.
    current.interaction.ID <- as.character(interactions[i, 11])
    
    # Check intersection with promoter regions.
    res <- as.character(lhs.promoter.results[as.character(lhs.promoter.results[, 4]) == current.interaction.ID, 8])
    if (length(res) > 0) {
      annotation[["lhs_promoter"]][i] <- length(res)
      annotation[["lhs_promoter_ids"]][i] <- paste0(unique(res), collapse=",")
    }
    res <- as.character(rhs.promoter.results[as.character(rhs.promoter.results[, 4]) == current.interaction.ID, 8])
    if (length(res) > 0) {
      annotation[["rhs_promoter"]][i] <- length(res)
      annotation[["rhs_promoter_ids"]][i] <- paste0(unique(res), collapse=",")
    }
    
    # Check intersection with ATAC-seq peaks. 
    res <- as.character(lhs.atac.results[as.character(lhs.atac.results[, 4]) == current.interaction.ID, 8])
    if (length(res) > 0) {
      
      promoter.res <- res[atac.annotations[res, "annotation"] == "promoter"]
      annotation[["lhs_promoter_ATAC-seq"]][i] <- length(promoter.res)
      annotation[["lhs_promoter_ATAC-seq_ids"]][i] <- paste0(unique(promoter.res), collapse=",")
      distal.res <- res[atac.annotations[res, "annotation"] == "distal"]
      annotation[["lhs_distal_ATAC-seq"]][i] <- length(distal.res)
      annotation[["lhs_distal_ATAC-seq_ids"]][i] <- paste0(unique(distal.res), collapse=",")
      
    }
    res <- as.character(rhs.atac.results[as.character(rhs.atac.results[, 4]) == current.interaction.ID, 8])
    if (length(res) > 0) {
      
      promoter.res <- res[atac.annotations[res, "annotation"] == "promoter"]
      annotation[["rhs_promoter_ATAC-seq"]][i] <- length(promoter.res)
      annotation[["rhs_promoter_ATAC-seq_ids"]][i] <- paste0(unique(promoter.res), collapse=",")
      distal.res <- res[atac.annotations[res, "annotation"] == "distal"]
      annotation[["rhs_distal_ATAC-seq"]][i] <- length(distal.res)
      annotation[["rhs_distal_ATAC-seq_ids"]][i] <- paste0(unique(distal.res), collapse=",")
      
    }
    
    # Check intersections with features.
    if (!is.null(features)) {
      
      for (j in 1:length(feature.pairs[, 1])) {
        
        res <- as.character(lhs.feature.results[[feature.pairs[j, 2]]][as.character(lhs.feature.results[[feature.pairs[j, 2]]][, 4]) == current.interaction.ID, 8])
        if (length(res) > 0) {
          annotation[[paste0("lhs_", feature.pairs[j, 2])]][i] <- length(res)
          annotation[[paste0("lhs_", feature.pairs[j, 2], "_ids")]][i] <- paste0(res, collapse=",")
        }
        res <- as.character(rhs.feature.results[[feature.pairs[j, 2]]][as.character(rhs.feature.results[[feature.pairs[j, 2]]][, 4]) == current.interaction.ID, 8])
        if (length(res) > 0) {
          annotation[[paste0("rhs_", feature.pairs[j, 2])]][i] <- length(res)
          annotation[[paste0("rhs_", feature.pairs[j, 2], "_ids")]][i] <- paste0(res, collapse=",")
        }
        
      }
      
    }
    
  }
  
  # Return results.
  return(annotation)
  
}

