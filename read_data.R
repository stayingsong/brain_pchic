# Read in interactions by cell type.
# Notes: Interaction IDs are only added after interaction deduplication.
#
read.interactions <- function(interactions.list, pipeline.dir) {
  
  # Read in file containing interactions list.
  interaction.files <- read.table(paste0(pipeline.dir, "/", interactions.list), sep="\t", header=F, stringsAsFactors=F)
  
  # Read in interactions for each cell type.
  interaction.data <- list()
  for (i in 1:length(interaction.files[, 1])) {

    cell.type <- interaction.files[i, 2]
    interaction.data[[cell.type]] <- read.table(paste0(pipeline.dir, "/", interaction.files[i, 1]), sep="\t", stringsAsFactors=F, header=T,
                                                colClasses=c("character", "integer", "integer", "character", "character", "integer", "integer", "character", "integer", "numeric"))
    colnames(interaction.data[[cell.type]]) <- c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score")
    print(paste0("Read in ", length(interaction.data[[cell.type]][, 1]), " interactions for: ", cell.type))
    
  }
  
  # Return interactions.
  return(interaction.data)
  
}

# Read in ATAC-seq peaks by cell type, removing peaks with duplicate BED coordinates.
# Notes: 
#
read.atac.seq.data <- function(atac.seq.peak.list, pipeline.dir, allowed.chrs, max.atac.seq.peaks) {
  
  # Read in file containing ATAC-seq peaks list.
  atac.seq.peak.files <- read.table(paste0(pipeline.dir, "/", atac.seq.peak.list), sep="\t", header=F, stringsAsFactors=F)
  
  # Read in ATAC-seq peaks for each cell type.
  atac.seq.peaks <- list()
  for (i in 1:length(atac.seq.peak.files[, 1])) {
    
    cell.type <- atac.seq.peak.files[i, 2]
    atac.seq.peaks[[cell.type]] <- read.table(paste0(pipeline.dir, "/", atac.seq.peak.files[i, 1]), sep="\t", stringsAsFactors=F, header=F,
                                              colClasses=c("character", "integer", "integer", "character", "numeric", "character", "numeric", "numeric", "numeric", "numeric"))
    colnames(atac.seq.peaks[[cell.type]]) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
    atac.seq.peaks[[cell.type]] <- atac.seq.peaks[[cell.type]][atac.seq.peaks[[cell.type]]$chrom %in% allowed.chrs, ]
    print(paste0("Read in ", length(atac.seq.peaks[[cell.type]]$chrom), " raw peaks on main chromosomes for: ", cell.type))
    
    # Remove peaks with duplicate BED coordinates, keeping only the peaks with the highest q-values.
    atac.seq.peaks[[cell.type]] <- atac.seq.peaks[[cell.type]][order(-atac.seq.peaks[[cell.type]]$qValue), ]
    is.duplicated <- duplicated(atac.seq.peaks[[cell.type]][, 1:3])
    duplicated.peaks <- atac.seq.peaks[[cell.type]][is.duplicated, ]
    atac.seq.peaks[[cell.type]] <- atac.seq.peaks[[cell.type]][!is.duplicated, ]
    print(paste0(length(atac.seq.peaks[[cell.type]]$chrom), " peaks retained after deduplication"))
    
    # Downsample ATAC-seq peaks based on score (if a maximum # of peaks is specified).
    if (!is.null(max.atac.seq.peaks) && (length(atac.seq.peaks[[cell.type]]$chrom) > max.atac.seq.peaks)) {
      
      atac.seq.peaks[[cell.type]] <- atac.seq.peaks[[cell.type]][1:max.atac.seq.peaks, ]
      print(paste0(length(atac.seq.peaks[[cell.type]]$chrom), " peaks retained after downsampling"))
      
    }
    
    # Add ATAC-seq peak IDs.
    atac.seq.peaks[[cell.type]] <- cbind(atac.seq.peaks[[cell.type]], paste0(cell.type, "_", 1:length(atac.seq.peaks[[cell.type]]$chrom)))
    colnames(atac.seq.peaks[[cell.type]])[11] <- "atac_seq_peak_ID"

  }
  
  # Return ATAC-seq peaks.
  return(atac.seq.peaks)
  
}

# Parse TSS annotation file to obtain a list of promoters.
# Notes: The TSS annotation file can be prepared from the GENCODE GTF release file using 'gtf_processing.R'.
#
read.promoters <- function(annotation.file, allowed.chrs, tss.upstream, tss.downstream, allowed.types) {
  
  # Read in TSS annotation file. 
  annotation <- read.table(annotation.file, sep="\t", header=F, stringsAsFactors=F, 
                           colClasses=c("character", "integer", "integer", "character", "character", "character", "character"))
  colnames(annotation) <- c("chrom", "chromStart", "chromEnd", "strand", "gene_id", "gene_type", "gene_name")
  annotation <- annotation[annotation$chrom %in% allowed.chrs, ]
  print(paste0("Read in ", length(unique(annotation$gene_id)), " unique gene IDs on main chromosomes."))
  
  # Expand intervals based on strand information.
  tss.plus <- annotation[annotation[, 4] == "+", ]
  tss.plus[, 2] <- tss.plus[, 2] - tss.upstream
  tss.plus[tss.plus[, 2] < 1, 2] <- 1
  tss.plus[, 3] <- tss.plus[, 3] + tss.downstream
  tss.minus <- annotation[annotation[, 4] == "-", ]
  tss.minus[, 2] <- tss.minus[, 2] - tss.downstream
  tss.minus[tss.minus[, 2] < 1, 2] <- 1
  tss.minus[, 3] <- tss.minus[, 3] + tss.upstream
  promoters <- rbind(tss.plus, tss.minus)

  # Filter promoters by gene type if list of allowed types is specified.
  if (!is.null(allowed.types))
    promoters <- promoters[promoters$gene_type %in% allowed.types, ]
  print(paste0("After filtering for allowed gene types, retained ", length(unique(promoters$gene_id)), " unique gene IDs."))
  
  # Return list of promoters.
  return(promoters[order(promoters$gene_name), ])
  
}

# Read in RNA-seq transcript quantification results from RSEM by replicate/cell type.
# Notes: 'field' can be one of the following: length, effective_length, expected_count, TPM, FPKM.  Rows of the RSEM results files should follow the same order of gene IDs.
#
read.rna.seq.data <- function(rna.seq.list, pipeline.dir, field) {
  
  # Read in file containing RNA-seq results list.
  rna.seq.files <- read.table(paste0(pipeline.dir, "/", rna.seq.list), sep="\t", stringsAsFactors=F, header=F)
  
  # Read in RNA-seq results by cell type.
  cell.types <- unique(rna.seq.files[, 2])
  rna.seq.summary <- list()
  for (cell.type in cell.types) {
    
    # Read in raw data for all replicates corresponding to the current cell type.
    subset.files <- rna.seq.files[rna.seq.files[, 2] == cell.type, ]
    raw.data <- list()
    num.replicates <- length(subset.files[, 1])
    for (i in 1:num.replicates) {
      raw.data[[i]] <- read.table(paste0(pipeline.dir, "/", subset.files[i, 1]), sep="\t", stringsAsFactors=F, header=T, 
                                  colClasses=c("character", "character", rep("numeric", 5)))
    }
    
    # Get list of gene IDs. Rows of the RSEM results files should follow the same order of gene IDs.
    gene.ids <- unlist(strsplit(raw.data[[1]]$gene_id, split="\\."))[(1:length(raw.data[[1]]$gene_id))*2-1]
    num.genes <- length(gene.ids)
    
    # Create and populate results table for each replicate. Rows are gene IDs, columns are replicates.
    rna.seq.summary[[cell.type]] <- data.frame(matrix(NA, num.genes, num.replicates))
    rownames(rna.seq.summary[[cell.type]]) <- gene.ids
    colnames(rna.seq.summary[[cell.type]]) <- c(paste0(cell.type, "_", 1:num.replicates))
    for (i in 1:num.replicates) {
      
      if (all(unlist(strsplit(raw.data[[i]]$gene_id, split="\\."))[(1:length(raw.data[[i]]$gene_id))*2-1] == gene.ids)) {
        rna.seq.summary[[cell.type]][, i] <- raw.data[[i]][, field]
      } else {
        print("ERROR: Rows of the RSEM results files should follow the same order of gene IDs.")
      }
      
    }
    
  }
  
  # Return RNA-seq results.
  return(rna.seq.summary)
  
}

# Read in features by cell type.
# Notes: Returns nested list with the following levels, in order: cell type, feature type (chromhmm, har, snp, vista), feature name.
#
read.features <- function(features.list, pipeline.dir) {
  
  # Read in file containing features list.
  features.files <- read.table(paste0(pipeline.dir, "/", features.list), sep="\t", stringsAsFactors=F, header=F)
  
  # Read in features by cell type.
  cell.types <- unique(features.files[, 2])
  features.all <- list()
  for (cell.type in cell.types) {
    
    # Read in features by feature type.
    subset.files <- features.files[features.files[, 2] == cell.type, ]
    feature.types <- unique(subset.files[, 3])
    features.all[[cell.type]] <- list()
    for (feature.type in feature.types) {
      
      # Read in features according to feature type.
      feature.type.files <- subset.files[subset.files[, 3] == feature.type, ]
      features.all[[cell.type]][[feature.type]] <- list()
      for (i in 1:length(feature.type.files[, 1])) {
        
        # Files for different feature types will be formatted differently (todo: update to read in formats from config file and avoid hard-coding).
        features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]] <- read.table(paste0(pipeline.dir, "/", feature.type.files[i, 1]), sep="\t", stringsAsFactors=F, header=F)
        num.features <- length(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]][, 1])
        print(paste0("Read in ", num.features, " ", feature.type.files[i, 4], " ", feature.type, " features for cell type: ", cell.type))
        if (feature.type == "chromhmm") {
          colnames(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "state")
        } else if (feature.type == "har") {
          colnames(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "har_ID")
        } else if (feature.type == "snp") {
          colnames(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "rsid", "ref", "alt", "query_snp_rsid", "is_query_snp", "exon_overlap")
        } else if (feature.type == "vista") {
          colnames(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "vista_ID", "annotation")
        } 
        
        # Add feature ID as the last column.
        num.columns <- length(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]])
        features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]] <- 
          cbind(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]], paste0(feature.type.files[i, 4], "_", 1:num.features))
        colnames(features.all[[cell.type]][[feature.type]][[feature.type.files[i, 4]]])[num.columns + 1] <- "feature_ID"
        
      }
      
    }
    
  }
  
  # Return features.
  return(features.all)
  
}


# End ---------------------------------------------------------------------

