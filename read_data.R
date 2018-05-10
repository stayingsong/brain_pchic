# Read in cell type specific interactions.
# Inputs: Tab separated list of filenames containing interactions in the first column and cell types in the second column.  Interactions should be in CHiCAGO IBED format (10 columns).
# Outputs: List where the cell types are the keys and the interactions are the values.  Adds and 11th interaction ID column.
read.interactions <- function(interaction.list) {
  
  # Read in file containing interaction lists.
  interaction.files <- read.table(interaction.list, sep="\t", header=F, stringsAsFactors=F)
  
  # Read in interactions for each cell type.
  interaction.data <- list()
  for (i in 1:length(interaction.files[, 1])) {
    
    interaction.data[[interaction.files[i, 2]]] <- read.table(file=interaction.files[i, 1], sep="\t", stringsAsFactors=F, header=T,
                                                              colClasses=c("character", "integer", "integer", "character", "character", "integer", "integer", "character", "integer", "numeric"))
    interaction.data[[interaction.files[i, 2]]] <- cbind(interaction.data[[interaction.files[i, 2]]], paste0(interaction.files[i, 2], "_", 1:length(interaction.data[[interaction.files[i, 2]]][, 1])))
    colnames(interaction.data[[interaction.files[i, 2]]]) <- c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score", "ID")
    print(paste0("Read in ", length(interaction.data[[interaction.files[i, 2]]][, 1]), " interactions for: ", interaction.files[i, 2]))
    
  }
  
  # Return interactions.
  return(interaction.data)
  
}

# Read in cell type specific ATAC-seq peak lists.
# Inputs: Tab separated list of filenames containing ATAC-seq peaks in the first column and cell types in the second column.
# Outputs: ATAC-seq peaks should be in narrowPeak format (10 columns).  Adds an 11th ATAC-seq peak ID column.
read.atac.seq.data <- function(atac.seq.list) {
  
  # Read in file containing ATAC-seq peaks.
  atac.seq.files <- read.table(atac.seq.list, sep="\t", header=F, stringsAsFactors=F)
  
  # For filtering out peaks on disallowed chromosomes.
  allowed.chrs <- paste0("chr", c(1:22, "X")) # Exclude Y since we don't have any interactions in Y.
  
  # Read in ATAC-seq peaks for each cell type.
  atac.seq.peaks <- list()
  for (i in 1:length(atac.seq.files[, 1])) {
    
    atac.seq.peaks[[atac.seq.files[i, 2]]] <- read.table(file=atac.seq.files[i, 1], sep="\t", stringsAsFactors=F, header=F,
                                                         colClasses=c("character", "integer", "integer", "character", "numeric", "character", "numeric", "numeric", "numeric", "numeric"))
    atac.seq.peaks[[atac.seq.files[i, 2]]] <- cbind(atac.seq.peaks[[atac.seq.files[i, 2]]], paste0(file=atac.seq.files[i, 2], "_", 1:length(atac.seq.peaks[[atac.seq.files[i, 2]]][, 1])))
    colnames(atac.seq.peaks[[atac.seq.files[i, 2]]]) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak", "ID")
    print(paste0("Read in ", length(atac.seq.peaks[[atac.seq.files[i, 2]]][, 1]), " peaks for: ", atac.seq.files[i, 2]))
    
    duplicate.test <- atac.seq.peaks[[atac.seq.files[i, 2]]][, 1:3]
    num.peaks <- length(duplicate.test[, 1])
    num.unique.peaks <- length(unique(duplicate.test)[, 1])
    if (num.peaks != num.unique.peaks) {
      print("Duplicate peaks found")
    }
    
    atac.seq.peaks[[atac.seq.files[i, 2]]] <- atac.seq.peaks[[atac.seq.files[i, 2]]][atac.seq.peaks[[atac.seq.files[i, 2]]][, 1] %in% allowed.chrs, ]
    print(paste0("Retained ", length(atac.seq.peaks[[atac.seq.files[i, 2]]][, 1]), " peaks on main chromosomes for: ", atac.seq.files[i, 2]))
    
  }
  
  # Return peaks.
  return(atac.seq.peaks)
  
}

# Read in cell type specific expression data.
# Inputs: Tab separated list of filenames containing RNA-seq expression results in the first column and cell types in the second column.
# Outputs:
read.expression.tpm.data <- function(rna.seq.list) {
  
  # Read in file containing RNA-seq expression results lists.
  rna.seq.files <- read.table(rna.seq.list, sep="\t", stringsAsFactors=F, header=F)
  
  # Read in RNA-seq expression results by cell type.
  cell.types <- unique(rna.seq.files[, 2])
  expression.results <- list()
  rna.seq.data <- list()
  for (cell.type in cell.types) {
    
    # Read in data for all replicates corresponding to the current cell type.
    cell.type.files <- rna.seq.files[rna.seq.files[, 2] == cell.type, ]
    rna.seq.data[[cell.type]] <- list()
    for (i in 1:length(cell.type.files[, 1])) {
      rna.seq.data[[cell.type]][[i]] <- read.table(file=cell.type.files[i, 1], sep="\t", stringsAsFactors=F, header=T, 
                                                   colClasses=c("character", "character", rep("numeric", 15)))
    }
    
    # Get list of gene IDs (assumes all results are in same order with respect to gene ID).
    gene.ids <- unlist(strsplit(rna.seq.data[[cell.type]][[1]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[1]]$gene_id))*2-1]
    num.genes <- length(gene.ids)
    
    # Create and populate table containing count results for each replicate.
    expression.results[[cell.type]] <- data.frame(matrix(NA, num.genes, length(rna.seq.data[[cell.type]])))
    colnames(expression.results[[cell.type]]) <- c(paste0(cell.type, 1:length(rna.seq.data[[cell.type]])))
    rownames(expression.results[[cell.type]]) <- gene.ids
    for (i in 1:length(cell.type.files[, 1])) {
      
      if (all(unlist(strsplit(rna.seq.data[[cell.type]][[i]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[i]]$gene_id))*2-1] == gene.ids)) {
        expression.results[[cell.type]][, i] <- rna.seq.data[[cell.type]][[i]]$TPM
      } else {
        print("Error: gene IDs do not correspond for this replicate.")
      }
      
    }
    
  }
  
  # Return expression results.
  return(expression.results)
  
}

# Read in cell type specific expression data.
# Inputs: Tab separated list of filenames containing RNA-seq expression results in the first column and cell types in the second column.
# Outputs:
read.expression.counts.data <- function(rna.seq.list) {
  
  # Read in file containing RNA-seq expression results lists.
  rna.seq.files <- read.table(rna.seq.list, sep="\t", stringsAsFactors=F, header=F)
  
  # Read in RNA-seq expression results by cell type.
  cell.types <- unique(rna.seq.files[, 2])
  expression.results <- list()
  rna.seq.data <- list()
  for (cell.type in cell.types) {
    
    # Read in data for all replicates corresponding to the current cell type.
    cell.type.files <- rna.seq.files[rna.seq.files[, 2] == cell.type, ]
    rna.seq.data[[cell.type]] <- list()
    for (i in 1:length(cell.type.files[, 1])) {
      rna.seq.data[[cell.type]][[i]] <- read.table(file=cell.type.files[i, 1], sep="\t", stringsAsFactors=F, header=T, 
                                                   colClasses=c("character", "character", rep("numeric", 15)))
    }
    
    # Get list of gene IDs (assumes all results are in same order with respect to gene ID).
    gene.ids <- unlist(strsplit(rna.seq.data[[cell.type]][[1]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[1]]$gene_id))*2-1]
    num.genes <- length(gene.ids)
    
    # Create and populate table containing count results for each replicate.
    expression.results[[cell.type]] <- data.frame(matrix(NA, num.genes, length(rna.seq.data[[cell.type]])))
    colnames(expression.results[[cell.type]]) <- c(paste0(cell.type, 1:length(rna.seq.data[[cell.type]])))
    rownames(expression.results[[cell.type]]) <- gene.ids
    for (i in 1:length(cell.type.files[, 1])) {
      
      if (all(unlist(strsplit(rna.seq.data[[cell.type]][[i]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[i]]$gene_id))*2-1] == gene.ids)) {
        expression.results[[cell.type]][, i] <- rna.seq.data[[cell.type]][[i]]$expected_count
      } else {
        print("Error: gene IDs do not correspond for this replicate.")
      }
      
    }
  
  }
  
  # Return expression results.
  return(expression.results)
  
}

#
# Inputs: Tab separated list of filenames containing RNA-seq expression results in the first column and cell types in the second column.
# Outputs:
read.expression.gene.lengths <- function(rna.seq.list) {
  
  # Read in file containing RNA-seq expression results lists.
  rna.seq.files <- read.table(rna.seq.list, sep="\t", stringsAsFactors=F, header=F)
  
  #
  cell.types <- unique(rna.seq.files[, 2])
  expression.results <- list()
  rna.seq.data <- list()
  for (cell.type in cell.types) {
    
    #
    cell.type.files <- rna.seq.files[rna.seq.files[, 2] == cell.type, ]
    rna.seq.data[[cell.type]] <- list()
    for (i in 1:length(cell.type.files[, 1])) {
      rna.seq.data[[cell.type]][[i]] <- read.table(file=cell.type.files[i, 1], sep="\t", stringsAsFactors=F, header=T, 
                                                   colClasses=c("character", "character", rep("numeric", 15)))
    }
    
    # Get list of gene IDs (assumes all results are in same order with respect to gene ID).
    gene.ids <- unlist(strsplit(rna.seq.data[[cell.type]][[1]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[1]]$gene_id))*2-1]
    num.genes <- length(gene.ids)
    
    # Create and populate table containing count results for each replicate.
    expression.results[[cell.type]] <- data.frame(matrix(NA, num.genes, length(rna.seq.data[[cell.type]])))
    colnames(expression.results[[cell.type]]) <- c(paste0(cell.type, 1:length(rna.seq.data[[cell.type]])))
    rownames(expression.results[[cell.type]]) <- gene.ids
    for (i in 1:length(cell.type.files[, 1])) {
      
      if (all(unlist(strsplit(rna.seq.data[[cell.type]][[i]]$gene_id, split="\\."))[(1:length(rna.seq.data[[cell.type]][[i]]$gene_id))*2-1] == gene.ids)) {
        expression.results[[cell.type]][, i] <- rna.seq.data[[cell.type]][[i]]$length
      } else {
        print("Error: gene IDs do not correspond for this replicate.")
      }
      
    }
    
  }
  
  # Return expression results.
  return(expression.results)
  
}

# Read in cell type specific features.
# Inputs:
# Outputs:
read.features <- function(features.list) {
  
  # Read in file containing features lists.
  features.files <- read.table(features.list, sep="\t", stringsAsFactors=F, header=F)
  
  # First sort by cell type.
  cell.types <- unique(features.files[, 2])
  features <- list()
  for (cell.type in cell.types) {
    
    #
    print(paste0("Reading in features for cell type: ", cell.type))
    
    # Then sort by feature type.
    cell.type.files <- features.files[features.files[, 2] == cell.type, ]
    feature.types <- unique(cell.type.files[, 3])
    features[[cell.type]] <- list()
    for (feature.type in feature.types) {
      
      # Read in features.
      feature.type.files <- cell.type.files[cell.type.files[, 3] == feature.type, ]
      features[[cell.type]][[feature.type]] <- list()
      for (i in 1:length(feature.type.files[, 1])) {
        
        # Different types of features will have differing #s of columns.
        features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]] <- read.table(feature.type.files[i, 1], sep="\t", stringsAsFactors=F, header=F)
        #print(paste0(feature.type.files[i, 4], ": ", length(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]][, 1])))
        
        #
        if (feature.type == "ChIP-seq") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue")
        } else if (feature.type == "DHS-seq") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue")
        } else if (feature.type == "SNP") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "rsID", "ref", "alt", "query_snp_rsid", "is_query_snp")
        } else if (feature.type == "eQTL") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "rsID", "q-value", "gene_id", "gene_name")
        } else if (feature.type == "vista") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "vistaID", "annotation")
        } else if (feature.type == "ChromHMM") {
          colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]]) <- c("chrom", "chromStart", "chromEnd", "state")
        }
        
        # Add feature ID as the last column.
        num.columns <- length(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]])
        num.features <- length(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]][, 1])
        features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]] <- 
          cbind(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]], paste0(feature.type.files[i, 4], "_", 1:num.features))
        colnames(features[[cell.type]][[feature.type]][[feature.type.files[i, 4]]])[num.columns+1] <- "ID"
        
      }
    
    }
    
  }
  
  # Return features.
  return(features)
  
}

