# Makes plots of # of interactions by score.
# Inputs: Interaction set(s) organized by cell type, vector of scores, and output prefix.
# Outputs: Plots of # of interactions by score.
plotNumInteractionsBySignificance <- function(interaction.data, p.values, xlim, ylim, cutoff, output.prefix) {
  
  # Analyze # of interactions as a function of p-value cutoff.
  cell.types <- names(interaction.data)
  nums <- matrix(NA, length(cell.types), length(p.values))
  for (i in 1:length(cell.types)) {
    
    for (j in 1:length(p.values)) {
      
      nums[i, j] <- length(interaction.data[[cell.types[i]]][interaction.data[[cell.types[i]]]$score > p.values[j], 1])
      
    }
    
  }
  
  # Prepare plot.
  pdf(paste0(output.prefix, "/interaction_numbers_by_significance_", paste0(xlim, collapse="-"), "_", paste0(ylim, collapse="-"), ".pdf"), width=6, height=6)
  colors <- c("#3399FF", "#FF9966", "#66FF99", "magenta")
  plot(1, type='n', xlim=xlim, ylim=ylim, xlab='CHiCAGO score cutoff', ylab='# of interactions', main="# of interactions as a function of CHiCAGO score cutoff")
  for (i in 1:length(cell.types)) {
    lines(p.values, nums[i, ], type='l', lwd=2, col=colors[i])
  }
  axis(1, at=seq(xlim[1], xlim[2], 1), labels=seq(xlim[1], xlim[2], 1))
  if (is.null(cutoff) == FALSE)
    abline(v=cutoff, lty=2, col="red")
  legend("topright", legend=cell.types, lty=rep(1, length(cell.types)), col=colors[1:length(cell.types)])
  dev.off()
  
}

# Plots distribution of bin sizes for interactions.
# Inputs: Interaction set(s) organized by cell type, minimum resolution, and output prefix.
# Outputs:
plotBinWidths <- function(interactions, max.width, output.prefix, file.suffix) {
  
  # Plot distribution of expanded bin sizes.
  res <- 2000
  pdf(paste0(output.prefix, "/interaction_bin_size_distribution_", file.suffix, ".pdf"), width=4*length(names(interactions)), height=4)
  par(mfrow=c(1, length(names(interactions))))
  for (cell.type in names(interactions)) {
    
    # Plot distribution of fragment lengths to get an estimate of mean resolution.
    bin.sizes <- c(abs(interactions[[cell.type]][, 2] - interactions[[cell.type]][, 3]), 
                   abs(interactions[[cell.type]][, 6] - interactions[[cell.type]][, 7]))
    hist(bin.sizes, main=paste0(cell.type, " interactions bin size distribution"), sub=paste0("mean bin size: ", signif(mean(bin.sizes), 5)),
         xlab="bin size (bp)", xlim=c(0, max.width), breaks=seq(0, 1000000, res), ylab=NULL, xaxt="n")
    axis(1, at=seq(0, max.width, res*2), labels=seq(0, max.width, res*2))
    abline(v=mean(bin.sizes), col="red")
    
  }
  dev.off()
  
}

#
# Inputs:
# Outputs:
plotDistances <- function(interaction.data, type, output.prefix) {
  
  #
  cell.types <- names(interaction.data)
  
  # 
  if (type == "hist") {
    
    colors=c("#3399FF", "#FF9966", "#66FF99")
    pdf(paste0(output.prefix, "/interaction_distance_distribution_histogram.pdf"), width=4*length(cell.types), height=4)
    par(mfrow=c(1, length(cell.types)))
    i=0
    for (cell.type in cell.types) {
      
      #
      distances <- abs(interaction.data[[cell.type]][, 2] + interaction.data[[cell.type]][, 3] - interaction.data[[cell.type]][, 6] - interaction.data[[cell.type]][, 7])/2000
      bin.size <- 10
      max.dist <- ceiling(max(distances)/bin.size)*bin.size
      i <- i + 1
      ceiling(max(distances)/1000/bin.size)
      hist(distances, main=paste0(cell.type, " interaction distance distribution"),
           xlab="interaction distance (kb)", xlim=c(0, 250), ylab="# interactions", breaks=seq(0, max.dist, bin.size), col=colors[i])
      abline(v=mean(distances), col="red")
      
    }
    dev.off()
    
  } else if (type == "cdf") {
    
    # Get distances for each cell type.
    distances <- list()
    for (cell.type in cell.types) {
      distances[[cell.type]] <- abs(interaction.data[[cell.type]][, 2] + interaction.data[[cell.type]][, 3] - interaction.data[[cell.type]][, 6] - interaction.data[[cell.type]][, 7])/2000
    }
    
    # Create a single chart with all 3 CDF plots.
    pdf(paste0(output.prefix, "/interaction_distance_distribution_cdf_individual.pdf"), width=6, height=6)
    plot(ecdf(distances[["cortical"]]), col="#3399FF", main=NA,  xlim=c(0, 500), ylim=c(0, 1), xlab="interaction distance (kb)", ylab="interaction distance CDF")
    plot(ecdf(distances[["hippocampal"]]), col="#FF9966", add=T)
    plot(ecdf(distances[["astrocyte"]]), col="#66FF99", add=T)
    plot(ecdf(c(distances[["cortical"]], distances[["hippocampal"]], distances[["astrocyte"]])), col="darkgrey", lty=2, add=T)
    legend('right', c(cell.types, "merged"), fill=c(c("#3399FF", "#FF9966", "#66FF99"), "darkgrey"), border=NA)
    dev.off()
    
    # Create a single chart with all 3 CDF plots.
    pdf(paste0(output.prefix, "/interaction_distance_distribution_cdf_all.pdf"), width=6, height=6)
    plot(ecdf(c(distances[["cortical"]], distances[["hippocampal"]], distances[["astrocyte"]])), col="darkgrey", lty=2, main=NA,  xlim=c(0, 500), ylim=c(0, 1), 
         xlab="interaction distance (kb)", ylab="interaction distance CDF")
    dev.off()
    
  } else if (type == "hybrid") {
    
    # Get distances for each cell type.
    distances <- list()
    all.distances <- c()
    for (cell.type in cell.types) {
      distances[[cell.type]] <- abs(interaction.data[[cell.type]][, 2] + interaction.data[[cell.type]][, 3] - interaction.data[[cell.type]][, 6] - interaction.data[[cell.type]][, 7])/2000
      all.distances <- c(all.distances, distances[[cell.type]])
    }
   
    # Get an idea of the y-axis limits for the combined plot.
    max.display <- 500
    bin.size <- 10
    max.dist <- ceiling(max(all.distances)/bin.size)*bin.size
    h_test <- hist(all.distances[all.distances < max.display], breaks=seq(0, max.dist, bin.size))
    top <- max(h_test$counts)*1.1
    top <- 25000
    
    # Prepare histogram and CDF.
    pdf(paste0(output.prefix, "/interaction_distance_distribution_hybrid.pdf"), width=6, height=6)
    h  <- hist(all.distances[all.distances < max.display], breaks=seq(0, max.dist, bin.size), xlim=c(0, max.display-bin.size*1.9), 
               ylim=c(-1, top), col="#F0F0F0", border="#C1C1C1", yaxt="n", xlab="interaction distance (kb)", ylab=NA, main=NA)
    ec <- ecdf(all.distances)
    ec_1 <- ecdf(distances[["cortical"]])
    ec_2 <- ecdf(distances[["hippocampal"]])
    ec_3 <- ecdf(distances[["astrocyte"]])
    axis(side=2, at=seq(0, 25000, 5000), labels=seq(0, 25000, 5000), las=2)
    axis(side=4, at=seq(0, top, top/5), labels=seq(0, 1, 0.2), las=2)
    lines(x = c(0, h$mids), y=c(0, ec(h$mids)*top), col ="darkgrey", lwd=2)
    lines(x = c(0, h$mids), y=c(0, ec_1(h$mids)*top), col ="#3399FF", lwd=2)
    lines(x = c(0, h$mids), y=c(0, ec_2(h$mids)*top), col ="#FF9966", lwd=2)
    lines(x = c(0, h$mids), y=c(0, ec_3(h$mids)*top), col ="#66FF99", lwd=2)
    abline(h=top, col="#5B5B5B", lty=2, lwd=1)
    abline(h=0, col="#5B5B5B", lty=2, lwd=1)
    legend('right', c(cell.types, "combined"), fill=c(c("#3399FF", "#FF9966", "#66FF99"), "darkgrey"), border=NA)
    dev.off()
    
  }
  
}

#
# Inputs:
# Outputs:
printATACPeakAnnotations <- function(peak.annotations, output.prefix) {
  
  # Compile results table for storing location of ATAC-seq peaks and whether or not they are interacting.
  results <- matrix(NA, length(names(peak.annotations)), 16)
  rownames(results) <- names(peak.annotations)
  colnames(results) <- c("distal_none", "distal_intron", "distal_exon", "distal_both", "promoter_none", "promoter_intron", "promoter_exon", "promoter_both",
                         "distal_none_int", "distal_intron_int", "distal_exon_int", "distal_both_int", 
                         "promoter_none_int", "promoter_intron_int", "promoter_exon_int", "promoter_both_int")
  compiled <- list()
  for (cell.type in names(peak.annotations)) {
    
    # Fill in table with the appropriate ATAC-seq annotation categories.
    current <- peak.annotations[[cell.type]]
    distal.none <- table((current$annotation == "distal") & (current$introns == "FALSE") & (current$exons == "FALSE"))[["TRUE"]]
    distal.intron <- table((current$annotation == "distal") & (current$introns == "TRUE") & (current$exons == "FALSE"))[["TRUE"]]
    distal.exon <- table((current$annotation == "distal") & (current$introns == "FALSE") & (current$exons == "TRUE"))[["TRUE"]]
    distal.both <- table((current$annotation == "distal") & (current$introns == "TRUE") & (current$exons == "TRUE"))[["TRUE"]]
    promoter.none <- table((current$annotation == "promoter") & (current$introns == "FALSE") & (current$exons == "FALSE"))[["TRUE"]]
    promoter.intron <- table((current$annotation == "promoter") & (current$introns == "TRUE") & (current$exons == "FALSE"))[["TRUE"]]
    promoter.exon <- table((current$annotation == "promoter") & (current$introns == "FALSE") & (current$exons == "TRUE"))[["TRUE"]]
    promoter.both <- table((current$annotation == "promoter") & (current$introns == "TRUE") & (current$exons == "TRUE"))[["TRUE"]]
    distal.none.interacting <- table((current$annotation == "distal") & (current$introns == "FALSE") & (current$exons == "FALSE") & (current$num_interactions > 0))[["TRUE"]]
    distal.intron.interacting <- table((current$annotation == "distal") & (current$introns == "TRUE") & (current$exons == "FALSE") & (current$num_interactions > 0))[["TRUE"]]
    distal.exon.interacting <- table((current$annotation == "distal") & (current$introns == "FALSE") & (current$exons == "TRUE") & (current$num_interactions > 0))[["TRUE"]]
    distal.both.interacting <- table((current$annotation == "distal") & (current$introns == "TRUE") & (current$exons == "TRUE") & (current$num_interactions > 0))[["TRUE"]]
    promoter.none.interacting <- table((current$annotation == "promoter") & (current$introns == "FALSE") & (current$exons == "FALSE") & (current$num_interactions > 0))[["TRUE"]]
    promoter.intron.interacting <- table((current$annotation == "promoter") & (current$introns == "TRUE") & (current$exons == "FALSE") & (current$num_interactions > 0))[["TRUE"]]
    promoter.exon.interacting <- table((current$annotation == "promoter") & (current$introns == "FALSE") & (current$exons == "TRUE") & (current$num_interactions > 0))[["TRUE"]]
    promoter.both.interacting <- table((current$annotation == "promoter") & (current$introns == "TRUE") & (current$exons == "TRUE") & (current$num_interactions > 0))[["TRUE"]]
    
    # Fill in results table.
    compiled[[cell.type]] <- c(distal.none, distal.intron, distal.exon, distal.both, promoter.none, promoter.intron, promoter.exon, promoter.both,
                               distal.none.interacting, distal.intron.interacting, distal.exon.interacting, distal.both.interacting,
                               promoter.none.interacting, promoter.intron.interacting, promoter.exon.interacting, promoter.both.interacting)
    #print(compiled[[cell.type]][1:8])
    #print(compiled[[cell.type]][9:16])
    results[cell.type, ] <- compiled[[cell.type]]
    
  }
  
  # Write results to file.
  write.table(results, file=paste0(output.prefix, "/ATAC-seq_peaks_annotation.txt"), quote=F, row.names=T, col.names=T, sep="\t")
  
}

# Annotates ATAC-seq peaks.
# Inputs: A set of annotations (annotation, promoter_genes, num_interactions, proximal_genes, distal_genes) organized by cell type.
# Outputs: 
plotATACPeakAnnotations <- function(peak.annotations, output.prefix) {
  
  # Plot results in a series of pie charts.
  cell.types <- names(peak.annotations)
  pdf(paste0(output.prefix, "/ATAC-seq_peaks_annotation.pdf"), width=4*length(names(peak.annotations)), height=4)
  par(mfrow=c(1, length(names(peak.annotations))))
  for (cell.type in cell.types) {
    
    # Count peaks in different categories.
    P_I <- table((peak.annotations[[cell.type]]$annotation == "promoter") & (peak.annotations[[cell.type]]$num_interactions > 0))["TRUE"]
    D_I <- table((peak.annotations[[cell.type]]$annotation == "distal") & (peak.annotations[[cell.type]]$num_interactions > 0))["TRUE"]
    P_NI <- table((peak.annotations[[cell.type]]$annotation == "promoter") & (peak.annotations[[cell.type]]$num_interactions == 0))["TRUE"]
    D_NI <- table((peak.annotations[[cell.type]]$annotation == "distal") & (peak.annotations[[cell.type]]$num_interactions == 0))["TRUE"]
    
    # Settings for pie chart.
    colors=list(promoter='#87CEFA',
                distal='#FCE6C9',
                interaction='#DC143C',
                none='#FFFFFF')
    iniR=0.8
    offset.rad=pi/2
    offset.angle=90
    
    # Start plotting pie chart.
    label.slices <- rev(c(P_NI+P_I, D_I+D_NI))
    labels <- rev(c('promoter', 'distal'))
    labels <- paste0(labels, ", ", format(round(label.slices/sum(label.slices)*100, 1), nsmall=1), "%")
    pie(label.slices, radius=1.05*iniR, init.angle=offset.angle, col=rep('white', length(label.slices)), border=TRUE, labels=labels, cex=1, cex.main=1, main=paste0(cell.type, " atac-seq peak annotation"))
    floating.pie(0, 0, 1, radius=1*iniR, startpos=offset.rad, col='white', border=FALSE)
    
    # Plot interaction slices.
    interaction.slices <- rev(c(P_NI, P_I+D_I, D_NI))
    interaction.colors <- as.character(colors[rev(c('none', 'interaction', 'none'))])
    floating.pie(0, 0, interaction.slices, radius=0.95*iniR, startpos=offset.rad, col=interaction.colors, border=FALSE)
    
    # Plot annotation slices.
    annotation.slices <- rev(c(P_NI+P_I, D_I+D_NI))
    annotation.colors <- as.character(colors[rev(c('promoter', 'distal'))])
    floating.pie(0, 0, annotation.slices, radius=0.85*iniR, startpos=offset.rad, col=annotation.colors, border=TRUE)
    
    # Plot blank interior.
    floating.pie(0, 0, 1, radius=0.45*iniR, startpos=offset.rad, col='black', border=TRUE)
    floating.pie(0, 0, 1, radius=0.446*iniR, startpos=offset.rad, col='white', border=FALSE)
    
  }
  dev.off()
  
}

# Plots distance between different classes of ATAC-seq peaks (promoter, distal) and nearest interaction loci.
# Inputs:
# Outputs:
distanceATACPeaksInteractionLoci <- function(peaks, annotation, interaction.data.all, output.prefix) {
  
  # Analyze multiple score cutoffs for each cell type.
  score.cutoffs <- c(3, 4, 5)
  cell.types <- names(peaks)
  pdf(paste0(output.prefix, paste0("/ATAC-seq_peaks_closest_interaction_loci_cdf_individual.pdf")), width=4*length(cell.types), height=4)
  par(mfrow=c(1, 3))
  for (cell.type in cell.types) {
    
    # Print current cell type.
    print(cell.type)
    
    # Plot CDF for first score cutoff.
    cutoff <- interaction.data.all[[cell.type]][interaction.data.all[[cell.type]]$score > score.cutoffs[1], ]
    interaction.loci <- rbind(cutoff[, c(1:3, 11)], setNames(cutoff[, c(5:7, 11)], colnames(cutoff[, c(1:3, 11)])))
    interaction.loci <- sortBed(interaction.loci)
    sorted.peaks <- sortBed(peaks[[cell.type]][, c(1:3, 11)])
    closest.results <- bedTools.2in(functionstring="bedtools closest", bed1=sorted.peaks, bed2=interaction.loci, opt.string="-t first -d")
    all.distances <- closest.results$V9/1000
    plot(ecdf(all.distances), main=paste0("peak distance to closest interaction loci CDF (", cell.type, ")"), col="#71C671",
         xlab="distance to closest interaction (kb)", xlim=c(0, 100), ylim=c(0, 0.9), ylab="% peaks within distance to closest interaction", cex.main=1.25)
    
    # Plot percentiles for first score cutoff.
    print(table(all.distances < 10)["TRUE"]/length(all.distances))
    print(table(all.distances < 20)["TRUE"]/length(all.distances))
    print(table(all.distances < 50)["TRUE"]/length(all.distances))
    print(table(all.distances < 100)["TRUE"]/length(all.distances))
    
    # Plot CDF for second score cutoff.
    cutoff <- interaction.data.all[[cell.type]][interaction.data.all[[cell.type]]$score > score.cutoffs[2], ]
    interaction.loci <- rbind(cutoff[, c(1:3, 11)], setNames(cutoff[, c(5:7, 11)], colnames(cutoff[, c(1:3, 11)])))
    interaction.loci <- sortBed(interaction.loci)
    sorted.peaks <- sortBed(peaks[[cell.type]][, c(1:3, 11)])
    closest.results <- bedTools.2in(functionstring="bedtools closest", bed1=sorted.peaks, bed2=interaction.loci, opt.string="-t first -d")
    all.distances <- closest.results$V9/1000
    plot(ecdf(all.distances), col="#FF7F24", xlab="distance to closest interaction (kb)", add=T)
    
    # Plot percentiles for second score cutoff.
    print(table(all.distances < 10)["TRUE"]/length(all.distances))
    print(table(all.distances < 20)["TRUE"]/length(all.distances))
    print(table(all.distances < 50)["TRUE"]/length(all.distances))
    print(table(all.distances < 100)["TRUE"]/length(all.distances))
    
    # Plot CDF for third score cutoff.
    cutoff <- interaction.data.all[[cell.type]][interaction.data.all[[cell.type]]$score > score.cutoffs[3], ]
    interaction.loci <- rbind(cutoff[, c(1:3, 11)], setNames(cutoff[, c(5:7, 11)], colnames(cutoff[, c(1:3, 11)])))
    interaction.loci <- sortBed(interaction.loci)
    sorted.peaks <- sortBed(peaks[[cell.type]][, c(1:3, 11)])
    closest.results <- bedTools.2in(functionstring="bedtools closest", bed1=sorted.peaks, bed2=interaction.loci, opt.string="-t first -d")
    all.distances <- closest.results$V9/1000
    plot(ecdf(all.distances), col="#CD5555", xlab="distance to closest interaction (kb)", add=T)
  
    # Plot percentiles for third score cutoff.
    print(table(all.distances < 10)["TRUE"]/length(all.distances))
    print(table(all.distances < 20)["TRUE"]/length(all.distances))
    print(table(all.distances < 50)["TRUE"]/length(all.distances))
    print(table(all.distances < 100)["TRUE"]/length(all.distances))
    
    # Plot legend.
    legend('topleft', c("score cutoff=3", "score cutoff=4", "score cutoff=5"), fill=c("#71C671", "#FF7F24", "#CD5555"), border=NA)
    
  }
  dev.off()
  
  # Prepare specific plot for a single score cutoff but with all cell types.
  score.cutoffs <- 5
  cell.types <- names(peaks)
  names <- cell.types
  colors=c("#3399FF", "#FF9966", "#66FF99")
  all.distances <- list()
  for (i in 1:length(cell.types)) {
    
    # Get distances for current cell type.
    cell.type <- cell.types[i]
    cutoff <- interaction.data.all[[cell.type]][interaction.data.all[[cell.type]]$score > score.cutoffs[1], ]
    interaction.loci <- rbind(cutoff[, c(1:3, 11)], setNames(cutoff[, c(5:7, 11)], colnames(cutoff[, c(1:3, 11)])))
    interaction.loci <- sortBed(interaction.loci)
    sorted.peaks <- sortBed(peaks[[cell.type]][, c(1:3, 11)])
    closest.results <- bedTools.2in(functionstring="bedtools closest", bed1=sorted.peaks, bed2=interaction.loci, opt.string="-t first -d")
    all.distances[[i]] <- closest.results$V9/1000

  }
  
  # Plot CDFs for each cell type on one plot.
  pdf(paste0(output.prefix, paste0("/ATAC-seq_peaks_closest_interaction_loci_cdf_all.pdf")), width=6, height=6)
  plot(ecdf(all.distances[[1]]), main=paste0("peak distance to closest interaction loci CDF"), col=colors[1],
       xlab="distance to closest interaction (kb)", xlim=c(0, 100), ylim=c(0, 0.9), ylab="% peaks within distance to closest interaction", cex.main=1.25)
  plot(ecdf(all.distances[[2]]), col=colors[2], xlab="distance to closest interaction (kb)", add=T)
  plot(ecdf(all.distances[[3]]), col=colors[3], xlab="distance to closest interaction (kb)", add=T)
  legend('topleft', cell.types, fill=colors, border=NA)
  dev.off()
  
}

#
# Inputs:
# Outputs:
plotPromoterOtherInteractionTranscriptTypes <- function(annotation.data, promoters, output.prefix) {
  
  # Determine which interactions are promoter to non-promoter interactions for each cell type.
  cell.types <- names(annotation.data)
  Pr_NP.subset.genes <- list()
  for (i in 1:length(cell.types)) {
    
    # Get annotation data for current cell type.
    annotation <- annotation.data[[cell.types[i]]]
    
    # Retrieve only promoter to non-promoter interactions.
    Pr_NP.selector <- (((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) & (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0))
    NP_Pr.selector <- ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) & ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)))
    Pr_NP.subset.genes[[cell.types[i]]] <- c(annotation[["lhs_promoter_ids"]][Pr_NP.selector], annotation[["rhs_promoter_ids"]][NP_Pr.selector])
    
  }
  
  # Make pie chart of gene types for promoter to non-promoter interactions.
  protein_coding_labels <- c("IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_V_gene", "protein_coding", "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene")
  small_noncoding_RNA_labels <- c("miRNA", "miscRNA", "rRNA", "snoRNA", "snRNA")
  lincRNA_labels <- c("lincRNA")
  pseudogene_labels <- c("IG_C_pseudogene", "IG_J_pseudogene", "IG_V_pseudogene", "polymorphic_pseudogene", "pseudogene", "TR_J_pseudogene", "TR_V_pseudogene")
  other_labels <- c("3prime_overlapping_ncrna", "antisense", "processed_transcript", "sense_intronic", "sense_overlapping")
  pdf(paste0(output.prefix, "/promoter_to_other_transcript_types.pdf"), width=4*length(cell.types), height=4)
  par(mfrow=c(1, length(cell.types)))
  for (cell.type in cell.types) {
    
    # Determine counts for each category.
    protein_coding <- 0
    small_noncoding_RNA <- 0
    lincRNA <- 0
    pseudogene <- 0
    other <- 0
    num.genes <- length(Pr_NP.subset.genes[[cell.type]])
    for (i in 1:num.genes) {
      
      # Get genes for each interaction.
      current.genes <- unlist(strsplit(Pr_NP.subset.genes[[cell.type]][i], split=","))
      for (j in 1:length(current.genes)) {
        
        # Determine gene type for each gene.
        type <- unique(promoters$gene_type[promoters$gene_id == current.genes[j]])
        if (length(type) > 1) 
          print("Error")
        else if (length(type) == 1) {
          if (type %in% protein_coding_labels)
            protein_coding <- protein_coding + 1
          else if (type %in% small_noncoding_RNA_labels)
            small_noncoding_RNA <- small_noncoding_RNA + 1
          else if (type %in% lincRNA_labels)
            lincRNA <- lincRNA + 1
          else if (type %in% pseudogene_labels)
            pseudogene <- pseudogene + 1
          else if (type %in% other_labels)
            other <- other + 1
        }
        
      }
      
    }
    
    # Plot pie chart.
    slices <- c(protein_coding, small_noncoding_RNA, lincRNA, pseudogene, other)
    labels <- c("protein coding genes", "small noncoding RNA", "lincRNA", "pseudogenes", "other")
    pie(slices, labels=labels, main=paste0(cell.type, " promoter to other transcript types"))

  }
  dev.off()
  
}

# This is the same as plotInteractionClasses() except it plots all the cell types in one figure.
# Inputs:
# Outputs:
plotAllInteractionClasses <- function(annotation.data, output.prefix) {
  
  # Get counts for each category for each cell type.
  cell.types <- names(annotation.data)
  Pr_Pr <- rep(NA, length(cell.types))
  Pr_NP <- rep(NA, length(cell.types))
  NP_NP <- rep(NA, length(cell.types))
  PA_PA <- rep(NA, length(cell.types))
  PA_Pn <- rep(NA, length(cell.types))
  Pn_Pn <- rep(NA, length(cell.types))
  Pn_DA <- rep(NA, length(cell.types))
  PA_DA <- rep(NA, length(cell.types))
  Pn_Dn <- rep(NA, length(cell.types))
  PA_Dn <- rep(NA, length(cell.types))
  Pn_DE <- rep(NA, length(cell.types))
  PA_DE <- rep(NA, length(cell.types))
  Pn_NE <- rep(NA, length(cell.types))
  PA_NE <- rep(NA, length(cell.types))
  Pn_EA <- rep(NA, length(cell.types))
  PA_EA <- rep(NA, length(cell.types))
  Pn_Ea <- rep(NA, length(cell.types))
  PA_Ea <- rep(NA, length(cell.types))
  Pn_eA <- rep(NA, length(cell.types))
  PA_eA <- rep(NA, length(cell.types))
  Pn_ea <- rep(NA, length(cell.types))
  PA_ea <- rep(NA, length(cell.types))
  for (i in 1:length(cell.types)) {
    
    # Get annotation data for current cell type.
    annotation <- annotation.data[[cell.types[i]]]
    
    # Count # of interactions that are promoter to promoter, promoter to non-promoter, and non-promoter to non-promoter.
    Pr_Pr[i] <- table(((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
                        ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)))["TRUE"]
    Pr_NP[i] <- table((((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
                         (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)) |
                        ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) &
                           ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0))))["TRUE"]
    NP_NP[i] <- table((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) &
                        (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0))["TRUE"]
    print(NP_NP[i])
    
    # Dissect promoter to promoter interactions into finer categories.
    Pr_Pr.subset <- ((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
      ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0))
    PA_PA[i] <- table((annotation[["lhs_promoter_ATAC-seq"]][Pr_Pr.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_Pr.subset] > 0))["TRUE"]
    PA_Pn[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_Pr.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_Pr.subset] == 0)) |
                        ((annotation[["lhs_promoter_ATAC-seq"]][Pr_Pr.subset] == 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_Pr.subset] > 0)))["TRUE"]
    Pn_Pn[i] <- table((annotation[["lhs_promoter_ATAC-seq"]][Pr_Pr.subset] == 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_Pr.subset] == 0))["TRUE"]
    
    # Dissect promoter to non-promoter interactions into finer categories.
    Pr_NP.subset <- (((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) & (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)) |
      ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) & ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)))
    Pn_DA[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    PA_DA[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    Pn_Dn[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    PA_Dn[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    
    # Annotate enhancer states.
    lhs.enhancer <- annotation$lhs_active_enhancer + annotation$lhs_genic_enhancer + annotation$lhs_other_enhancer
    rhs.enhancer <- annotation$rhs_active_enhancer + annotation$rhs_genic_enhancer + annotation$rhs_other_enhancer
    Pr_NP.subset <- (((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) & (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)) |
      ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) & ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)))
    lhs.enhancer <- lhs.enhancer[Pr_NP.subset]
    rhs.enhancer <- rhs.enhancer[Pr_NP.subset]
    Pn_DE[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer > 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer > 0)))["TRUE"]
    PA_DE[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer > 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer > 0)))["TRUE"]
    Pn_NE[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer == 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer == 0)))["TRUE"]
    PA_NE[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer == 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer == 0)))["TRUE"]
    
    # Annotate enhancer + ATAC-seq states in conjunction.
    Pn_EA[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    PA_EA[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    Pn_Ea[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    PA_Ea[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer > 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer > 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    Pn_eA[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    PA_eA[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] > 0)))["TRUE"]
    Pn_ea[i] <- table(((annotation[["lhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (rhs.enhancer == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter"]][Pr_NP.subset] > 0) & (annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] == 0) & (lhs.enhancer == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    PA_ea[i] <- table(((annotation[["lhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (rhs.enhancer == 0) & (annotation[["rhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)) |
                        ((annotation[["rhs_promoter_ATAC-seq"]][Pr_NP.subset] > 0) & (lhs.enhancer == 0) & (annotation[["lhs_distal_ATAC-seq"]][Pr_NP.subset] == 0)))["TRUE"]
    
  }
  
  # Prepare counts table.
  counts.main <- matrix(0, 2, 3) # For each cell type, the categories are from bottom to top: P to P, P to other
  counts.aux <- matrix(0, 4, 3) # Split P to P into w/ or w/o ATAC-seq support, then split P to other into w/ or w/o distal ATAC-seq peak.
  counts.enh <- matrix(0, 6, 3)
  colnames(counts.main) <- cell.types
  colnames(counts.aux) <- cell.types
  colnames(counts.enh) <- cell.types
  rownames(counts.main) <- c("promoter to promoter", "promoter to other")
  rownames(counts.aux) <- c("w/ promoter ATAC-seq peak", "w/o promoter ATAC-seq peak", "w/ distal ATAC-seq peak", "w/o distal ATAC-seq peak")
  rownames(counts.enh) <- c("w/ distal enhancer state", "w/o distal enhancer state", 
                            "ATAC-seq peak & enhancer state", "ATAC-seq peak only", "enhancer state only", "neither ATAC-seq peak or enhancer state")
  
  # Fill in counts table.
  for (i in 1:length(cell.types)) {
    counts.main[1, i] <- Pr_Pr[i]
    counts.main[2, i] <- Pr_NP[i]
    counts.aux[1, i] <- PA_PA[i] + PA_Pn[i]
    counts.aux[2, i] <- Pn_Pn[i]
    counts.aux[3, i] <- Pn_DA[i] + PA_DA[i]
    counts.aux[4, i] <- Pn_Dn[i] + PA_Dn[i]
    counts.enh[1, i] <- Pn_DE[i] + PA_DE[i]
    counts.enh[2, i] <- Pn_NE[i] + PA_NE[i]
    counts.enh[3, i] <- Pn_EA[i] + PA_EA[i]
    counts.enh[4, i] <- Pn_eA[i] + PA_eA[i]
    counts.enh[5, i] <- Pn_Ea[i] + PA_Ea[i]
    counts.enh[6, i] <- Pn_ea[i] + PA_ea[i]
  }
  
  # Calculate percentages.
  percentages.main <- apply(counts.main, 2, function(x) {x*100/sum(x, na.rm=T)})
  percentages.aux <- apply(counts.aux, 2, function(x) {x*100/sum(x, na.rm=T)})
  percentages.enh <- apply(rbind(counts.aux[1, ], counts.aux[2, ], counts.enh[1, ], counts.enh[2, ]), 2, function(x) {x*100/sum(x, na.rm=T)})
  
  # Add colors and positions for each entry.
  colors.main <- c()
  colors.aux <- c()
  colors.enh <- c()
  d <- c(-1, -1, -1, -1, -1)
  spacing <- 1.75
  offset <- 0.25
  for (i in 1:length(cell.types)) {
    
    cell.type <- cell.types[i]
    d <- rbind(d,
               c(paste0(cell.type, "_main"), "Pr_Pr", percentages.main[1, i], 0.9, i*spacing),
               c(paste0(cell.type, "_main"), "Pr_NP", percentages.main[2, i], 0.9, i*spacing),
               c(paste0(cell.type, "_aux"), "Pa_Pa", percentages.aux[1, i], 0.4, i*spacing+offset),
               c(paste0(cell.type, "_aux"), "Pn_Pn", percentages.aux[2, i], 0.4, i*spacing+offset),
               c(paste0(cell.type, "_aux"), "Pr_DA", percentages.aux[3, i], 0.4, i*spacing+offset),
               c(paste0(cell.type, "_aux"), "Pr_Dn", percentages.aux[4, i], 0.4, i*spacing+offset),
               c(paste0(cell.type, "_enh"), "blank1", percentages.enh[1, i], 0.4, i*spacing+offset*2),
               c(paste0(cell.type, "_enh"), "blank2", percentages.enh[2, i], 0.4, i*spacing+offset*2),
               c(paste0(cell.type, "_enh"), "Pr_DE", percentages.enh[3, i], 0.4, i*spacing+offset*2),
               c(paste0(cell.type, "_enh"), "Pr_NE", percentages.enh[4, i], 0.4, i*spacing+offset*2))
    
  }
  
  # Convert to a data frame for plotting.
  d <- d[-1, ]
  d <- data.frame(d, stringsAsFactors=F)
  d[, 1] <- factor(d[, 1], levels=c(rbind(paste0(cell.types, "_main"), paste0(cell.types, "_aux"), paste0(cell.types, "_enh"))), ordered=T)
  d[, 2] <- factor(d[, 2], levels=c("Pr_NP", "Pr_Pr", "Pr_Dn", "Pr_DA", "Pn_Pn", "Pa_Pa", "Pr_NE", "Pr_DE", "blank1", "blank2"), ordered=T)
  d[, 3] <- as.numeric(d[, 3])
  d[, 4] <- as.numeric(d[, 4])
  d[, 5] <- as.numeric(d[, 5])
  colnames(d) <- c("cell_type", "class", "count", "bar_width", "x_pos")
  
  # Plot bar plot.
  ggplot(data=d, aes(x=cell_type, y=count, fill=class, width=bar_width)) + 
    geom_bar(stat="identity") +
    scale_fill_manual("legend", values = c("Pr_Pr"="#63B8FF", "Pr_NP"="#FFEC8B", "Pr_Dn"="#FFFFFF", 
                                           "Pr_DA"="#AB82FF", "Pn_Pn"="#FFFFFF", "Pa_Pa"="#AB82FF",
                                           "Pr_DE"="#DC143C", "Pr_NE"="#FFFFFF", "blank1"="#FFFFFF", "blank2"="#FFFFFF")) + 
    theme_bw() +
    theme(axis.line = element_line(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    labs(y="% of interactions") +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  ggsave(paste0(output.prefix, "/interaction_categories_bar_graph.pdf"), width=6, height=6)
  write.table(rbind(counts.main, counts.aux, counts.enh), file=paste0(output.prefix, "/interaction_categories_counts.txt"), row.names=T, col.names=T, sep="\t", quote=F)
  
}

# Plots a Venn diagram of cell type specificity of all interactions (this does not filter for the same interaction in multiple cell types).
# Inputs:
# Outputs:
plotSpecificityVenn <- function(specificity.categories, interaction.data, output.prefix) {
  
  # Get combinations of cell types.
  cell.types <- names(specificity.categories)
  combinations <- c()
  for (i in 1:length(cell.types)) {
    temp <- combn(cell.types, i)
    for (j in 1:length(temp[1, ])) {
      combinations <- c(combinations, paste0(temp[, j], collapse="&"))
    }
  }
  
  # Get counts for each combination.
  counts.old <- rep(0, length(combinations)) # Significant interactions between same fragments in two cell types are counted separately.
  counts.new <- rep(0, length(combinations)) # If same two fragments are interacting in two cell types, consider it the same interaction (no double counting).
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
      counts.old[i] <- counts.old[i] + length(which(match))
      
    }
    
    #
    combination.loci[[combinations[i]]] <- combination.loci[[combinations[i]]][-1, ]
    counts.new[i] <- length(unique(combination.loci[[combinations[i]]])[, 1])
    
  }
  
  #
  print(counts.old)
  print(counts.new)
  
  # Make plots.
  pdf(paste0(output.prefix, "/interaction_specificity_venn.pdf"), height=6, width=6)
  values <- counts.new
  names(values) <- combinations
  fit <- euler(values)
  colors <- c("#BFEFFF", "#FFD39B", "#BDFCC9")
  plot(fit, fill_opacity = 0.01, shape="circle", title="Cell type specificity of interactions",
       fill=c(colors[1:length(cell.types)]), border="transparent", fontsize=12, quantities=list(fontsize=12))
  dev.off()
  
}

#
# Inputs:
# Outputs:
plotPromoterHubsGeneral <- function(hub.results, subset.gene.ids, min.expression, output.prefix, output.suffix) {
  
  # Only keep records for promoters with at least one interaction and having an expression above the specified cutoff.
  cell.types <- names(hub.results)
  hub.results.filtered <- list()
  for (cell.type in cell.types) {
    
    hub.results.filtered[[cell.type]] <- hub.results[[cell.type]][!is.na(hub.results[[cell.type]]$expression_value), ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$expression_value > min.expression, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$num_interactions > 0, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$gene_id %in% subset.gene.ids, ]
    
  }
  
  # Plot bar graph representing # of genes with only one interaction and # of genes with more than one interaction for each cell type.
  pdf(paste0(output.prefix, "/promoter_hub_one_versus_many_counts.", output.suffix, ".pdf"), width=6, height=4)
  par(xpd=T, mar=par()$mar + c(0, 0, 0, 12))
  counts <- matrix(NA, length(cell.types), 2)
  for (i in 1:length(cell.types)) {
    
    counts[i, 1] <- table(hub.results[[cell.types[i]]]$num_interactions == 1)["TRUE"]
    counts[i, 2] <- table(hub.results[[cell.types[i]]]$num_interactions > 1)["TRUE"]
    
  }
  percentages <- apply(counts, 1, function(x) {x*100/sum(x)})
  rownames(percentages) <- c("one to one", "one to many")
  colnames(percentages) <- cell.types
  colors <- c('#FFEFDB', '#FFD39B')
  p <- barplot(percentages, col=colors, border="white", main=paste0("one to one versus one to many interactions for each gene"),
               space=rep(1, length(cell.types)), cex.main=0.6, cex.axis=0.6, yaxt='n', cex.names=0.6, width=1)
  legend(2.5*length(cell.types), 100, c("# genes with one interaction", "# genes with multiple interactions"), col=colors, cex=0.6, pch=15)
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
  
  # Plot histograms of # of interactions coincident on each gene.
  pdf(paste0(output.prefix, "/promoter_hub_interactions_per_promoter_histogram.all.", output.suffix, ".pdf"), width=4*length(cell.types), height=4)
  par(mfrow=c(1, length(cell.types)))
  for (cell.type in cell.types) {
    
    # Plot distribution of fragment lengths to get an estimate of mean resolution.
    h <- hist(hub.results[[cell.type]]$num_interactions[hub.results[[cell.type]]$num_interactions > 0], main=paste0("# of interactions per gene\nfor genes with > 1 interaction\n", cell.type),
              xlab="# interactions", xlim=c(0, 25), ylab="# genes", breaks=seq(0, 500, 1), cex.main=1)
    abline(v=mean(hub.results[[cell.type]]$num_interactions[hub.results[[cell.type]]$num_interactions > 0]), col="red")
    print(mean(hub.results[[cell.type]]$num_interactions[hub.results[[cell.type]]$num_interactions > 0]))
    write.table(h$counts, file=paste0(output.prefix, "/promoter_hub_interactions_per_promoter_counts.", cell.type, ".all.txt"), row.names=F, col.names=F, quote=F, sep="\t")
    
  }
  dev.off()
  
  # Plot histograms of # of interactions coincident on each gene (CODING GENE PROMOTERS ONLY).
  hub.results.subset <- list()
  pdf(paste0(output.prefix, "/promoter_hub_interactions_per_promoter_histogram.coding.", output.suffix, ".pdf"), width=4*length(cell.types), height=4)
  par(mfrow=c(1, length(cell.types)))
  for (cell.type in cell.types) {
    
    # Plot distribution of fragment lengths to get an estimate of mean resolution.
    hub.results.subset[[cell.type]] <- hub.results[[cell.type]][hub.results[[cell.type]]$gene_type == "protein_coding", ]
    h <- hist(hub.results.subset[[cell.type]]$num_interactions[hub.results.subset[[cell.type]]$num_interactions > 0], main=paste0("# of interactions per gene\nfor genes with > 1 interaction\n", cell.type),
              xlab="# interactions", xlim=c(0, 25), ylab="# genes", breaks=seq(0, 500, 1), cex.main=1)
    abline(v=mean(hub.results.subset[[cell.type]]$num_interactions[hub.results.subset[[cell.type]]$num_interactions > 0]), col="red")
    print(mean(hub.results.subset[[cell.type]]$num_interactions[hub.results.subset[[cell.type]]$num_interactions > 0]))
    write.table(h$counts, file=paste0(output.prefix, "/promoter_hub_interactions_per_promoter_counts.counts.", cell.type, ".coding.txt"), row.names=F, col.names=F, quote=F, sep="\t")
    
  }
  dev.off()
  
  # Analyze gene overlap between cell types.
  gene.id.lists <- list()
  gene.name.lists <- list()
  for (cell.type in cell.types) {

    gene.ids <- hub.results[[cell.type]]$gene_id[hub.results[[cell.type]]$num_interactions >  0]
    gene.id.lists[[cell.type]] <- unique(gene.ids)

  }
  
  # 
  interacting.genes <- unique(c(gene.id.lists[["cortical"]], gene.id.lists[["hippocampal"]], gene.id.lists[["astrocyte"]]))
  print(length(interacting.genes))
  
  #
  venn.diagram(list(cortical=gene.id.lists[["cortical"]], hippocampal=gene.id.lists[["hippocampal"]], astrocyte=gene.id.lists[["astrocyte"]]),
               fill=c("#3399FF", "#FF9966", "#66FF99"), alpha=rep(0.5, 3), cat.cex=1.5, cex=1.5,
               filename=paste0(output.prefix, "/interacting_gene_id_overlaps.all.png"))
  
  #
  coding.gene.id.lists <- list()
  noncoding.gene.id.lists <- list()
  for (cell.type in cell.types) {
    
    coding.gene.id.lists[[cell.type]] <- c()
    noncoding.gene.id.lists[[cell.type]] <- c()
    for (i in 1:length(gene.id.lists[[cell.type]])) {
      if (unique(promoters$gene_type[promoters$gene_id == gene.id.lists[[cell.type]][i]]) == "protein_coding") {
        coding.gene.id.lists[[cell.type]] <- c(coding.gene.id.lists[[cell.type]], gene.id.lists[[cell.type]][i])
      } else {
        noncoding.gene.id.lists[[cell.type]] <- c(noncoding.gene.id.lists[[cell.type]], gene.id.lists[[cell.type]][i])
      }
    }
    
  }
  
  #
  interacting.genes <- unique(c(coding.gene.id.lists[["cortical"]], coding.gene.id.lists[["hippocampal"]], coding.gene.id.lists[["astrocyte"]]))
  print(length(interacting.genes))
  
  #
  venn.diagram(list(cortical=coding.gene.id.lists[["cortical"]], hippocampal=coding.gene.id.lists[["hippocampal"]], astrocyte=coding.gene.id.lists[["astrocyte"]]),
               fill=c("#3399FF", "#FF9966", "#66FF99"), alpha=rep(0.5, 3), cat.cex=1.5, cex=1.5,
               filename=paste0(output.prefix, "/interacting_gene_id_overlaps.coding.png"))
  
  #
  interacting.genes <- unique(c(noncoding.gene.id.lists[["cortical"]], noncoding.gene.id.lists[["hippocampal"]], noncoding.gene.id.lists[["astrocyte"]]))
  print(length(interacting.genes))
  
  #
  venn.diagram(list(cortical=noncoding.gene.id.lists[["cortical"]], hippocampal=noncoding.gene.id.lists[["hippocampal"]], astrocyte=noncoding.gene.id.lists[["astrocyte"]]),
               fill=c("#3399FF", "#FF9966", "#66FF99"), alpha=rep(0.5, 3), cat.cex=1.5, cex=1.5,
               filename=paste0(output.prefix, "/interacting_gene_id_overlaps.noncoding.png"))
  
}

#
# Inputs:
# Outputs:
plotPromoterHubsSpecific <- function(hub.results, subset.gene.ids, min.expression, bin.size, min.points.per.bin, y.max, max.num.bins, category, output.prefix, output.suffix) {
  
  # Utility function for plotting boxplots.
  plotExpressionBoxPlot <- function(hub.results.filtered, bin.size, min.points.per.bin, y.max, max.num.bins, col.id, output.filename) {

    # Plot bar graph representing gene expression as a function of # of interactions (larger bin size).
    cell.types <- names(hub.results.filtered)
    custom.names <- c("cortical", "hippocampal", "astrocyte")
    pdf(output.filename, width=4*length(cell.types), height=4)
    par(mfrow=c(1, length(cell.types)))
    index <- 1
    for (cell.type in cell.types) {
      
      # Sort the data into bins.
      data <- list()
      i=1
      values <- hub.results.filtered[[cell.type]]$expression_value[(hub.results.filtered[[cell.type]][, col.id] > bin.size*(i - 1)) & (hub.results.filtered[[cell.type]][, col.id] <= bin.size*i)]
      while((length(values) >= min.points.per.bin) && (i <= max.num.bins)) {
        
        data[[i]] <- values
        i <- i + 1
        values <- hub.results.filtered[[cell.type]]$expression_value[(hub.results.filtered[[cell.type]][, col.id] > bin.size*(i - 1)) & (hub.results.filtered[[cell.type]][, col.id] <= bin.size*i)]
        
      }
      num.bins <- i - 1
      
      # Transform data into long form.
      transformed.data <- c()
      means <- c()
      for (i in 1:length(data)) {
        
        transformed.data <- rbind(transformed.data, cbind(data[[i]], rep(i, length(data[[i]])), rep(i, length(data[[i]]))))
        means <- c(means, mean(data[[i]], trim=0.25))
        
      }
      colnames(transformed.data) <- c("value", "group", "position")
      transformed.data <- data.frame(transformed.data)
      transformed.data[, 2] <- as.character(transformed.data[, 2])
      
      # Plot boxplot.
      boxplot(value ~ group, transformed.data, outline=F, boxwex=0.4, col="#FFFAF0", yaxt="n", ylim=c(0, y.max[index]), xlim=c(0.5, num.bins+0.5), xaxt="n",
              medlty=1, medlwd=1, medcol="black",
              main=paste0(cell.type, " mean gene expression versus\n", category, " per promoter"), xlab="# interactions per promoter", ylab="mean gene expression (normalized RPKM)")
      axis(side=1, at=(1:(num.bins+1)-0.5), labels=(0:num.bins)*bin.size)
      axis(side=2, at=seq(0, 300, 20), labels=seq(0, 300, 20))
      meds <- data.frame(x=1:length(data), y=sapply(data, median))
      meds$x <- meds$x*bin.size - bin.size/2
      means <- data.frame(x=1:length(data), y=means)
      means$x <- means$x*bin.size - bin.size/2
      statistics <- lm(y ~ x, data=means)
      means$x <- 1:num.bins
      points(1:num.bins, means$y, pch=23, bg="#FFC1C1", xaxt="n", cex=0.85)
      abline(lm(y ~ x, data=means), col="red", lty=5, lwd=2)
      fit = lm(y ~ x, data=means)
      print(summary(fit))
      
      # Increment index.
      index <- index + 1
      
    }
    dev.off()
    
  }
  
  # Apply minimum expression/interaction cutoffs to filter out noise.
  hub.results.filtered <- list()
  cell.types <- names(hub.results)
  for (cell.type in cell.types) {
    
    hub.results.filtered[[cell.type]] <- hub.results[[cell.type]][!is.na(hub.results[[cell.type]]$expression_value), ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$expression_value > min.expression, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]][, category] > 0, ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$gene_id %in% subset.gene.ids, ]
    
  }
  
  # For enhancer interactions, filter out promoters interacting with the same # of repressive interactions as enhancer interactions.
  if (category == "num_interactions_PE") {
    
    for (cell.type in cell.types) {
      
      selector <- hub.results.filtered[[cell.type]]$num_interactions_PR == 0
      print(table(selector))
      hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][selector, ]
      
    }
    
  }
  
  # Plot boxplots representing gene expression as a function of # of interactions.
  plotExpressionBoxPlot(hub.results.filtered, bin.size, min.points.per.bin, y.max, max.num.bins, category, output.filename=paste0(output.prefix, "/promoter_hub_expression_boxplot.", output.suffix, ".pdf"))
  
}

#
# Inputs:
# Outputs:
plotPromoterHubsComparisons <- function(hub.results, hub.results.control, subset.gene.ids, min.expression, filter.percent, y.max, output.prefix) {
  
  # Apply minimum expression/interaction cutoffs to filter out noise.
  cell.types <- names(hub.results)
  hub.results.filtered <- list()
  for (cell.type in cell.types) {
    hub.results.filtered[[cell.type]] <- hub.results[[cell.type]][!is.na(hub.results[[cell.type]]$expression_value), ]
    hub.results.filtered[[cell.type]] <- hub.results.filtered[[cell.type]][hub.results.filtered[[cell.type]]$expression_value > min.expression, ]
  }
  
  # Filter control set as well.
  hub.results.control.filtered <- list()
  for (cell.type in cell.types) {
    hub.results.control.filtered[[cell.type]] <- hub.results.control[[cell.type]][!is.na(hub.results.control[[cell.type]]$expression_value), ]
    hub.results.control.filtered[[cell.type]] <- hub.results.control.filtered[[cell.type]][hub.results.control.filtered[[cell.type]]$expression_value > min.expression, ]
  }
  
  # Compare repressive versus enhancer interactions.
  index <- 1
  enhancer.means.all <- c()
  repressive.means.all <- c()
  enhancer.control.means.all <- c()
  repressive.control.means.all <- c()
  enhancer.values.all <- c()
  repressive.values.all <- c()
  enhancer.control.values.all <- c()
  repressive.control.values.all <- c()
  for (cell.type in cell.types) {
    
    #
    print(cell.type)
    
    #
    enhancer.values <- hub.results.filtered[[cell.type]]$expression_value[(hub.results.filtered[[cell.type]]$num_interactions_PE > 0) & 
                                                                            (hub.results.filtered[[cell.type]]$num_interactions_PR == 0)]
    repressive.values <- hub.results.filtered[[cell.type]]$expression_value[(hub.results.filtered[[cell.type]]$num_interactions_PR > 0)]
    enhancer.control.values <- hub.results.control.filtered[[cell.type]]$expression_value[(hub.results.control.filtered[[cell.type]]$num_interactions_PE > 0) & 
                                                                            (hub.results.control.filtered[[cell.type]]$num_interactions_PR == 0)]
    repressive.control.values <- hub.results.control.filtered[[cell.type]]$expression_value[(hub.results.control.filtered[[cell.type]]$num_interactions_PR > 0)]
    
    #
    enhancer.values <- filterDistribution(enhancer.values, filter.percent)
    repressive.values <- filterDistribution(repressive.values, filter.percent)
    enhancer.control.values <- filterDistribution(enhancer.control.values, filter.percent)
    repressive.control.values <- filterDistribution(repressive.control.values, filter.percent)
    
    #
    enhancer.values.all <- c(enhancer.values.all, enhancer.values)
    repressive.values.all <- c(repressive.values.all, repressive.values)
    enhancer.control.values.all <- c(enhancer.control.values.all, enhancer.control.values)
    repressive.control.values.all <- c(repressive.control.values.all, repressive.control.values)
    enhancer.means.all <- c(enhancer.means.all, mean(enhancer.values))
    repressive.means.all <- c(repressive.means.all, mean(repressive.values))
    enhancer.control.means.all <- c(enhancer.control.means.all, mean(enhancer.control.values))
    repressive.control.means.all <- c(repressive.control.means.all, mean(repressive.control.values))
    
    #
    t <- t.test(enhancer.values, repressive.values)
    print(t$p.value)
    t <- t.test(enhancer.control.values, repressive.control.values)
    print(t$p.value)
    
    #
    all.data <- cbind(enhancer.values, rep("enhancer", length(enhancer.values)))
    all.data <- rbind(all.data, cbind(repressive.values, rep("repressive", length(repressive.values))))
    all.data <- rbind(all.data, cbind(enhancer.control.values, rep("enhancer control", length(enhancer.control.values))))
    all.data <- rbind(all.data, cbind(repressive.control.values, rep("repressive control", length(repressive.control.values))))
    colnames(all.data) <- c("expression", "group")
    all.data <- data.frame(all.data)
    all.data[, 1] <- as.numeric(as.character(all.data[, 1]))
    all.data[, 2] <- factor(all.data[, 2], ordered=T, levels=c("enhancer", "repressive", "enhancer control", "repressive control"))
    
    #
    p <- ggplot(all.data, aes(x=group, y=expression, fill=group)) + scale_y_log10() + geom_violin(trim=T, color="darkred") + 
      theme_bw() +
      theme(axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank()) +
      geom_boxplot(width=0.1, outlier.shape=NA) + ylim(0, y.max[index]) +
      stat_summary(fun.y=mean, geom="point", shape=23, size=2, color="black", fill="red") +
      scale_fill_manual(values=c("#87CEFA", "#FFFACD", "#F0F8FF", "#FFFFF0")) + 
      geom_hline(yintercept = mean(enhancer.values), linetype="dashed", color="black", size=0.25) + 
      geom_hline(yintercept = mean(repressive.values), linetype="dashed", color="black", size=0.25) +
      geom_hline(yintercept = mean(enhancer.control.values), linetype="dashed", color="grey", size=0.25) + 
      geom_hline(yintercept = mean(repressive.control.values), linetype="dashed", color="grey", size=0.25) +
      theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ggsave(paste0(output.prefix, "/promoter_hub_enhancer_versus_repressive_expression.", cell.type, ".pdf"), width=4.5, height=6)
    
    #
    index <- index + 1

  }
  
  #
  t <- t.test(enhancer.values.all, repressive.values.all)
  print(t$p.value)
  t <- t.test(enhancer.control.values.all, repressive.control.values.all)
  print(t$p.value)
  
  #
  all.data <- cbind(enhancer.values.all, rep("enhancer", length(enhancer.values.all)))
  all.data <- rbind(all.data, cbind(repressive.values.all, rep("repressive", length(repressive.values.all))))
  all.data <- rbind(all.data, cbind(enhancer.control.values.all, rep("enhancer control", length(enhancer.control.values.all))))
  all.data <- rbind(all.data, cbind(repressive.control.values.all, rep("repressive control", length(repressive.control.values.all))))
  colnames(all.data) <- c("expression", "group")
  all.data <- data.frame(all.data)
  all.data[, 1] <- as.numeric(as.character(all.data[, 1]))
  all.data[, 2] <- factor(all.data[, 2], ordered=T, levels=c("enhancer", "repressive", "enhancer control", "repressive control"))
  
  #
  means <- data.frame(matrix(NA, 12, 5), stringsAsFactors=F)
  means[, 1] <- as.numeric(means[, 1])
  means[, 2] <- as.numeric(means[, 2])
  means[, 3] <- as.character(means[, 3])
  means[, 4] <- as.character(means[, 4])
  means[, 5] <- as.character(means[, 5])
  colnames(means) <- c("x", "value", "cell_type", "group", "color")
  means[, 3] <- rep(c("cortical", "hippocampal", "astrocyte"), 4)
  means[, 4] <- c(rep("enhancer", 3), rep("repressive", 3), rep("enhancer control", 3), rep("repressive control", 3))
  means[, 1] <- c(1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4)
  means[, 2] <- c(enhancer.means.all, repressive.means.all, enhancer.control.means.all, repressive.control.means.all)
  means[, 5] <- rep(c("#3399FF", "#FF9966", "#66FF99"), 4)
  
  #
  means <- data.frame(matrix(NA, 4, 4), stringsAsFactors=F)
  means[, 1] <- as.numeric(means[, 1])
  means[, 2] <- as.numeric(means[, 2])
  means[, 3] <- as.character(means[, 3])
  means[, 4] <- as.character(means[, 4])
  colnames(means) <- c("x", "value", "group", "color")
  means[, 3] <- c(rep("enhancer", 1), rep("repressive", 1), rep("enhancer control", 1), rep("repressive control", 1))
  means[, 1] <- 1:4
  means[, 2] <- c(mean(enhancer.means.all), mean(repressive.means.all), mean(enhancer.control.means.all), mean(repressive.control.means.all))
  means[, 4] <- rep("white", 4)
  
  #
  p <- ggplot(all.data, aes(x=group, y=expression, fill=group)) + 
    scale_y_sqrt(breaks=c(0.5, seq(0, 200, 20)), labels=c(0, seq(0, 200, 20)), limits=c(0.5, 180)) +
    geom_violin(trim=T, color="black") + 
    theme_bw() +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank()) +
    geom_boxplot(width=0.1, outlier.shape=NA) +
    geom_point(data = means, mapping = aes(x = x, y = value), shape=23, color="black", fill=means$color, size=2.75) +
    geom_hline(yintercept = mean(enhancer.means.all), linetype="dashed", color="red", size=0.25) + 
    geom_hline(yintercept = mean(repressive.means.all), linetype="dashed", color="black", size=0.25) +
    geom_hline(yintercept = mean(enhancer.control.means.all), linetype="dashed", color="grey", size=0.25) + 
    geom_hline(yintercept = mean(repressive.control.means.all), linetype="dashed", color="grey", size=0.25) +
    scale_fill_manual(values=c("#87CEFA", "#FFFACD", "#F0F8FF", "#FFFFF0")) + 
    theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste0(output.prefix, "/promoter_hub_enhancer_versus_repressive_expression.all.pdf"), width=4.5, height=6)

  #
  print(mean(enhancer.means.all))
  print(mean(repressive.means.all))
  print(mean(enhancer.control.means.all))
  print(mean(repressive.control.means.all))
  
}

#
# Inputs:
# Outputs:
plotMotifs <- function(motif.file, output.prefix) {
  
  # Data frame with TF as row, and columns for p-value for each cell type.  Have a similar data frame for TPM instead of p-value.
  motif.pvalues <- read.table(motif.file, sep=",", header=T, stringsAsFactors=F, colClasses=c("character", rep("numeric", 6)))
  motif.pvalues[, 2:4] <- motif.pvalues[, 2:4]
  color.signal <- motif.pvalues[, 2:4]
  sizes.signal <- motif.pvalues[, 5:7]
  
  #
  my_palette <- colorRampPalette(c("#4169E1", "#CAE1FF", "#FFFFE0", "#FFFF00", "#FFFF00", "#FFCC09", "#FF9912"))(n=699)
  col_breaks = c(seq(-0.01, 2, length=100),
                 seq(2.01,  4.8, length=100),
                 seq(4.81, 9.6, length=100),
                 seq(9.61, 20, length=100),
                 seq(20.01, 40, length=100),
                 seq(40.01, 100, length=100),
                 seq(100.01, 384, length=100))
  original.sizes <- c(seq(-0.00001, 2.5, length=100),
                      seq(2.50001, 5, length=100),
                      seq(5.00001, 10, length=100),
                      seq(10.00001, 15, length=100),
                      seq(15.00001, 25, length=100),
                      seq(25.00001, 50, length=100),
                      seq(50.00001, 200, length=100))
  new.sizes <- c(seq(0, 6, length=100),
                 seq(6.00001, 12, length=100),
                 seq(12.00001, 18, length=100),
                 seq(18.00001, 24, length=100),
                 seq(24.00001, 30, length=100),
                 seq(30.00001, 36, length=100),
                 seq(36.00001, 48, length=100))
  
  #
  x.spacing <- 1.5
  y.spacing <- 25
  scale.factor <- 1
  
  #
  x.values <- matrix(NA, length(color.signal[, 1]), length(color.signal)+1)
  x.values[, 1] <- rep(x.spacing*1, length(x.values[, 1]))+0.5
  x.values[, 2] <- rep(x.spacing*2, length(x.values[, 1]))
  x.values[, 3] <- rep(x.spacing*3, length(x.values[, 1]))
  x.values[, 4] <- rep(x.spacing*4, length(x.values[, 1]))
  x.values <- x.values + x.spacing*1.5
  y.values <- matrix(NA, length(color.signal[, 1]), length(color.signal)+1)
  y.values[, 1] <- seq(-y.spacing, -y.spacing*length(y.values[, 1]), -y.spacing)
  y.values[, 2] <- seq(-y.spacing, -y.spacing*length(y.values[, 1]), -y.spacing)
  y.values[, 3] <- seq(-y.spacing, -y.spacing*length(y.values[, 1]), -y.spacing)
  y.values[, 4] <- seq(-y.spacing, -y.spacing*length(y.values[, 1]), -y.spacing)
  labels <- matrix("", length(color.signal[, 1]), length(color.signal)+1)
  labels[, 1] <- motif.pvalues[, 1]
  colors <- matrix(NA, length(color.signal[, 1]), length(color.signal)+1)
  colors[, 1] <- "white"
  colors[, 2] <- my_palette[.bincode(x=c(unlist(color.signal[, 1])), breaks=col_breaks)]
  colors[, 3] <- my_palette[.bincode(x=c(unlist(color.signal[, 2])), breaks=col_breaks)]
  colors[, 4] <- my_palette[.bincode(x=c(unlist(color.signal[, 3])), breaks=col_breaks)]
  sizes <- matrix(NA, length(color.signal[, 1]), length(color.signal)+1)
  sizes[, 1] <- 0
  sizes[, 2] <- new.sizes[.bincode(x=c(unlist(sizes.signal[, 1])), breaks=original.sizes)]*scale.factor
  sizes[, 3] <- new.sizes[.bincode(x=c(unlist(sizes.signal[, 2])), breaks=original.sizes)]*scale.factor
  sizes[, 4] <- new.sizes[.bincode(x=c(unlist(sizes.signal[, 3])), breaks=original.sizes)]*scale.factor
  
  #
  all.points <- data.frame(matrix(NA, length(unlist(x.values)), 5))
  colnames(all.points) <- c("x", "y", "label", "sizes", "colors")
  all.points[, 1] <- as.numeric(all.points[, 1])
  all.points[, 2] <- as.numeric(all.points[, 2])
  all.points[, 3] <- as.character(all.points[, 3])
  all.points[, 4] <- as.numeric(all.points[, 4])
  all.points[, 5] <- as.numeric(all.points[, 5])
  all.points[, 1] <- c(x.values)
  all.points[, 2] <- c(y.values)
  all.points[, 3] <- c(labels)
  all.points[, 4] <- c(sizes)
  all.points[, 5] <- c(colors)
  custom.colors <- as.character(all.points$colors)
  names(custom.colors) <- as.character(all.points$colors)
  ggplot(all.points, aes(x=x, y=y, size=sizes, fill=colors)) + scale_size_area() +
    geom_point(pch=21, color=c(rep("white", length(labels[, 1])), rep("black", 3*length(labels[, 1])))) +
    xlim(0, range(all.points$x)[2]+x.spacing) + ylim(range(all.points$y)[1]-y.spacing, 0) +
    geom_text(aes(label=c(labels)), hjust=1, vjust=0.5, size=4) +
    scale_fill_manual(values=custom.colors) +
    theme_bw() + theme(legend.position="none") +
    theme(axis.line = element_line(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  ggsave(paste0(output.prefix, "/motif_enrichment_scatterplot.pdf"), width=3.5, height=8)
  
  # Make color bar.
  color.bar(my_palette, min=0, max=100, nticks=8, 
            labels=c(expression(10^0), expression(10^-2), expression(10^-5), expression(10^-10), expression(10^-20), expression(10^-40), expression(10^-100), expression(10^-300)), 
            title="Enrichment p-value")
  
  # Make size bar.
  num.marks <- 8
  x.spacing <- 2
  y.spacing <- 25
  x.values <- matrix(NA, num.marks, 2)
  x.values[, 1] <- seq(x.spacing, num.marks*x.spacing, x.spacing)
  x.values[, 2] <- seq(x.spacing, num.marks*x.spacing, x.spacing)
  y.values <- matrix(NA, num.marks, 2)
  y.values[, 1] <- rep(-y.spacing*1, num.marks)
  y.values[, 2] <- rep(-y.spacing*2, num.marks)
  labels <- matrix("", num.marks, 2)
  labels[, 2] <- c(0.1, 1, 2.5, 5, 10, 20, 50, 150)
  colors <- matrix("#FFFFF0", num.marks, 2)
  sizes <- matrix(0, num.marks, 2)
  sizes[, 1] <- new.sizes[.bincode(x=c(unlist(labels[, 2])), breaks=original.sizes)]*scale.factor
  all.points <- data.frame(matrix(NA, length(unlist(x.values)), 5))
  colnames(all.points) <- c("x", "y", "label", "sizes", "colors")
  all.points[, 1] <- as.numeric(all.points[, 1])
  all.points[, 2] <- as.numeric(all.points[, 2])
  all.points[, 3] <- as.character(all.points[, 3])
  all.points[, 4] <- as.numeric(all.points[, 4])
  all.points[, 5] <- as.numeric(all.points[, 5])
  all.points[, 1] <- c(x.values)
  all.points[, 2] <- c(y.values)
  all.points[, 3] <- c(labels)
  all.points[, 4] <- c(sizes)
  all.points[, 5] <- c(colors)
  custom.colors <- as.character(all.points$colors)
  names(custom.colors) <- as.character(all.points$colors)
  ggplot(all.points, aes(x=x, y=y, size=sizes, fill=colors)) + 
    geom_point(pch=21, color=c(rep("black", length(labels[, 1])), rep("white", length(labels[, 1])))) + scale_size_area(max_size=25) +
    xlim(0, range(all.points$x)[2]+x.spacing) + ylim(range(all.points$y)[1]-y.spacing, 0) +
    geom_text(aes(label=c(labels)), hjust=0.5, vjust=0.5, size=4) +
    scale_fill_manual(values=custom.colors) + scale_size_area() +
    theme_bw() + theme(legend.position="none") +
    theme(axis.line = element_line(colour="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank()) +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    ggtitle("TF expression (normalized RPKM)")
  ggsave(paste0(output.prefix, "/motif_enrichment_size_bar.pdf"), width=3, height=1.5)
  
}

#
# Inputs:
# Outputs:
plotDiseaseEnrichhment <- function(enrichment.file, output.prefix) {
  
  #
  enrichment <- read.table(enrichment.file, sep=",", header=T, colClasses=c("character", "numeric", "numeric", "numeric"))
  values <- data.matrix(enrichment[, 2:4])
  rownames(values) <- c("AD", "ADHD", "ASD", "BD", "EP", "MDD", "PD", "SCZ", "T1D")
  
  # 
  colors = c(seq(0.5, 1.0, length=100),
             seq(1.01, 1.25, length=100),
             seq(1.26, 2.5, length=100))
  my_palette <- colorRampPalette(c("#4169E1", "white", "#CD0000"))(n = 299)
  
  # Prepare plot.
  pdf(paste0(output.prefix, "/SNP_enrichment_plain.pdf"), height=4, width=6)
  heatmap <- heatmap.2(values, Rowv=NA, Colv=NA, col = my_palette, breaks=colors, cexRow=1, cexCol=0.9, margins=c(6, 4), 
                       dendrogram="none", trace="none",  density.info="none",
                       keysize=1.25, key.par = list(cex=0.5))
  dev.off()
  # pdf("/Users/michael/Box Sync/MS_analysis/brain_pchic/figures/publication/SNPs/SNP_enrichment_values.pdf", height=4, width=6)
  # heatmap <- heatmap.2(values, Rowv=NA, Colv=NA, col = my_palette, breaks=colors, cexRow=1, cexCol=0.9, margins=c(6, 4), 
  #                      dendrogram="none", trace="none",  density.info="none", cellnote=format(round(values, 2), nsmall = 2), notecol="black", notecex=0.6,
  #                      keysize=1.25, key.par = list(cex=0.5))
  # dev.off()
  
}

#
# Inputs:
# Outputs:
plotAstrocyteMarkers <- function(astrocyte.markers.file, output.prefix) {
  
  #
  markers <- read.table(astrocyte.markers.file, sep=",", header=T, colClasses=c("character", "numeric", "numeric", "numeric"))
  values <- data.matrix(markers[, 4])
  rownames(values) <- as.character(markers[, 1])
  values <- values[order(values), ]
  values <- data.matrix(values)
  
  #
  colors = c(seq(0, 2.5, length=50),
             seq(2.51, 10, length=50),
             seq(10.01, 25, length=50),
             seq(25.01, 250, length=50))
  my_palette <- colorRampPalette(c("#6495ED", "white", "#E68080", "#CD0000"))(n = 199)
  
  #
  pdf(paste0(output.prefix, "/astrocyte_markers.pdf"), height=6, width=3)
  heatmap <- heatmap.2(cbind(values, values), Rowv=NA, Colv=NA, col = my_palette, breaks=colors, cexRow=0.5, cexCol=0.9, margins=c(5, 5), 
                       dendrogram="none", trace="none",  density.info="none",
                       keysize=1.25, key.par = list(cex=0.5))
  dev.off()
  
  #
  color.bar(my_palette, min=0, max=10, nticks=5, 
            labels=c(0, 2.5, 10, 25, 250), 
            title="Gene expression (TMM-normalized RPKM)")
  
}

# 
# Inputs:
# Outputs:
writeFeatureOverlapResults <- function(interactions.res, interactions.ann, output.prefix) {
  
  #
  cell.types <- names(interactions.res)
  for (cell.type in cell.types) {
    
    # Get annotation data for current cell type.
    annotation <- interactions.ann[[cell.type]]
    
    # Count # of interactions that are promoter to promoter, promoter to non-promoter, and non-promoter to non-promoter.
    Pr_Pr <- ((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
      ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0))
    Pr_NP <- (((annotation[["lhs_promoter"]] > 0) | (annotation[["lhs_promoter_ATAC-seq"]] > 0)) &
                (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)) |
      ((annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) &
         ((annotation[["rhs_promoter"]] > 0) | (annotation[["rhs_promoter_ATAC-seq"]] > 0)))
    NP_NP <- (annotation[["lhs_promoter"]] == 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) &
      (annotation[["rhs_promoter"]] == 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0)
    
    # Dissect promoter to non-promoter interactions into finer categories.
    Pn_DA <- Pr_NP & (((annotation[["lhs_promoter"]] > 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) & (annotation[["rhs_distal_ATAC-seq"]] > 0)) |
      ((annotation[["rhs_promoter"]] > 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0) & (annotation[["lhs_distal_ATAC-seq"]] > 0)))
    PA_DA <- Pr_NP & (((annotation[["lhs_promoter_ATAC-seq"]] > 0) & (annotation[["rhs_distal_ATAC-seq"]] > 0)) |
      ((annotation[["rhs_promoter_ATAC-seq"]] > 0) & (annotation[["lhs_distal_ATAC-seq"]] > 0)))
    Pn_Dn <- Pr_NP & (((annotation[["lhs_promoter"]] > 0) & (annotation[["lhs_promoter_ATAC-seq"]] == 0) & (annotation[["rhs_distal_ATAC-seq"]] == 0)) |
      ((annotation[["rhs_promoter"]] > 0) & (annotation[["rhs_promoter_ATAC-seq"]] == 0) & (annotation[["lhs_distal_ATAC-seq"]] == 0)))
    PA_Dn <- Pr_NP & (((annotation[["lhs_promoter_ATAC-seq"]] > 0) & (annotation[["rhs_distal_ATAC-seq"]] == 0)) |
      ((annotation[["rhs_promoter_ATAC-seq"]] > 0) & (annotation[["lhs_distal_ATAC-seq"]] == 0)))
    
    #
    print(paste0("Outputting ", length(interactions.res[[cell.type]][, 1]), " ", cell.type, " features into report file..."))
    features.list <- names(interactions.ann[[cell.type]])
    ids.indices <- grep("ids", features.list)
    features.list <- features.list[-ids.indices]
    frame <- interactions.sig[[cell.type]]
    frame <- cbind(frame, rep("", length(frame[, 1])))
    for (i in 1:length(features.list)) {
      
      #
      frame <- cbind(frame, interactions.ann[[cell.type]][features.list[i]])
      
    }
    
    #
    colnames(frame) <- c(colnames(interactions.res[[cell.type]]), "category", features.list)
    frame$category <- factor(frame$category, ordered=T, levels=c("Pr_Pr", "PA_DA", "Pn_DA", "PA_Dn", "Pn_Dn", "NP_NP"))
    frame$category[Pr_Pr] <- "Pr_Pr"
    frame$category[NP_NP] <- "NP_NP"
    frame$category[Pn_DA] <- "Pn_DA"
    frame$category[PA_DA] <- "PA_DA"
    frame$category[Pn_Dn] <- "Pn_Dn"
    frame$category[PA_Dn] <- "PA_Dn"
    
    #
    frame <- cbind(frame, interactions.ann[[cell.type]]$lhs_promoter_ids, interactions.ann[[cell.type]]$rhs_promoter_ids,
                   interactions.ann[[cell.type]]$lhs_AD_ids, interactions.ann[[cell.type]]$rhs_AD_ids,
                   interactions.ann[[cell.type]]$lhs_ADHD_ids, interactions.ann[[cell.type]]$rhs_ADHD_ids,
                   interactions.ann[[cell.type]]$lhs_ASD_ids, interactions.ann[[cell.type]]$rhs_ASD_ids,
                   interactions.ann[[cell.type]]$lhs_BD_ids, interactions.ann[[cell.type]]$rhs_BD_ids,
                   interactions.ann[[cell.type]]$lhs_EP_ids, interactions.ann[[cell.type]]$rhs_EP_ids,
                   interactions.ann[[cell.type]]$lhs_MDD_ids, interactions.ann[[cell.type]]$rhs_MDD_ids,
                   interactions.ann[[cell.type]]$lhs_PD_ids, interactions.ann[[cell.type]]$rhs_PD_ids,
                   interactions.ann[[cell.type]]$lhs_SCZ_ids, interactions.ann[[cell.type]]$rhs_SCZ_ids)
    
    #
    if (cell.type == "cortical") {
      colnames(frame)[71:88] <- c("lhs_promoter_ids", "rhs_promoter_ids", "lhs_AD_ids", "rhs_AD_ids", "lhs_ADHD_ids", "rhs_ADHD_ids",
                                  "lhs_ASD_ids", "rhs_ASD_ids", "lhs_BD_ids", "rhs_BD_ids", "lhs_EP_ids", "rhs_EP_ids",
                                  "lhs_MDD_ids", "rhs_MDD_ids", "lhs_PD_ids", "rhs_PD_ids", "lhs_SCZ_ids", "rhs_SCZ_ids")
    } else if (cell.type %in% c("hippocampal", "astrocyte")) {
      colnames(frame)[69:86] <- c("lhs_promoter_ids", "rhs_promoter_ids", "lhs_AD_ids", "rhs_AD_ids", "lhs_ADHD_ids", "rhs_ADHD_ids",
                                  "lhs_ASD_ids", "rhs_ASD_ids", "lhs_BD_ids", "rhs_BD_ids", "lhs_EP_ids", "rhs_EP_ids",
                                  "lhs_MDD_ids", "rhs_MDD_ids", "lhs_PD_ids", "rhs_PD_ids", "lhs_SCZ_ids", "rhs_SCZ_ids")
    }
    
    
    #
    frame[, "lhs_promoter_ids"] <- as.character(frame[, "lhs_promoter_ids"])
    frame[, "rhs_promoter_ids"] <- as.character(frame[, "rhs_promoter_ids"])
    
    #
    for (i in 1:length(frame[, 1])) {
      
      lhs.ids <- unlist(strsplit(as.character(frame[i, "lhs_promoter_ids"]), split=","))
      lhs.ids <- unique(promoters$gene_name[promoters$gene_id %in% lhs.ids])
      lhs.ids <- lhs.ids[lhs.ids != ""]
      frame[i, "lhs_promoter_ids"] <- paste0(lhs.ids, collapse=",")
      rhs.ids <- unlist(strsplit(as.character(frame[i, "rhs_promoter_ids"]), split=","))
      rhs.ids <- unique(promoters$gene_name[promoters$gene_id %in% rhs.ids])
      rhs.ids <- rhs.ids[rhs.ids != ""]
      frame[i, "rhs_promoter_ids"] <- paste0(rhs.ids, collapse=",")
      
    }
    
    # Rearrange columns before writing to file.
    if (cell.type == "cortical") {
      frame <- frame[, c(1:10, 66, 67, 65, 71, 69, 70, 68, 72,
                         47, 51, 53, 61, 49, 55, 57, 63,
                         48, 52, 54, 62, 50, 56, 58, 64,
                         27, 73, 28, 74, 25, 75, 26, 76,
                         29, 77, 30, 78, 31, 79, 32, 80,
                         33, 81, 34, 82, 39, 83, 40, 84,
                         35, 85, 36, 86, 37, 87, 38, 88)]
    } else if (cell.type %in% c("hippocampal", "astrocyte")) {
      frame <- frame[, c(1:10, 64, 65, 63, 69, 67, 68, 66, 70,
                         45, 49, 51, 59, 47, 53, 55, 61,
                         46, 50, 52, 60, 48, 54, 56, 62,
                         27, 71, 28, 72, 25, 73, 26, 74,
                         29, 75, 30, 76, 31, 77, 32, 78,
                         33, 79, 34, 80, 39, 81, 40, 82,
                         35, 83, 36, 84, 37, 85, 38, 86)]
    }
    
    #
    write.table(frame, file=paste0(output.prefix, "/interaction_feature_overlap_report.", cell.type, ".txt"), sep="\t", col.names=T, row.names=F, quote=F)
  
  }
  
}

