# Plot interaction annotation results.
# Notes:
# 
plotInteractionAnnotations <- function(annotation.data, output.prefix) {
  
  # Get counts for each category for each cell type.
  cell.types <- names(annotation.data)
  classes <- matrix(NA, length(cell.types), 10)
  rownames(classes) <- cell.types
  colnames(classes) <- c("Pr_Pr", "Pr_Np", "Np_Np", "Pa_Pa", "Pa_Pn", "Pn_Pn", "Pa_Da", "Pa_Dn", "Pn_Da", "Pn_Dn")
  for (cell.type in cell.types) {
    
    # Get annotation data for current cell type.
    ann <- annotation.data[[cell.type]]
    
    # Count # of interactions that are promoter to promoter, promoter to non-promoter, and non-promoter to non-promoter.
    classes[cell.type, "Pr_Pr"] <- table(((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & 
                                           ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0)))["TRUE"]
    classes[cell.type, "Pr_Np"] <- table((((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & 
                                            ((ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0))) |
                                           (((ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0)) & 
                                              ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0))))["TRUE"]
    classes[cell.type, "Np_Np"] <- table((ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0) & 
                                           (ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0))["TRUE"]
    
    # Dissect promoter to promoter interactions into finer categories.
    Pr_Pr.subset <- ((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0))
    classes[cell.type, "Pa_Pa"] <- table((ann[["promoter_ATAC-seq_lhs"]][Pr_Pr.subset] > 0) & 
                                           (ann[["promoter_ATAC-seq_rhs"]][Pr_Pr.subset] > 0))["TRUE"]
    classes[cell.type, "Pa_Pn"] <- table(((ann[["promoter_ATAC-seq_lhs"]][Pr_Pr.subset] > 0) & 
                                            (ann[["promoter_ATAC-seq_rhs"]][Pr_Pr.subset] == 0)) |
                                           ((ann[["promoter_ATAC-seq_lhs"]][Pr_Pr.subset] == 0) & 
                                              (ann[["promoter_ATAC-seq_rhs"]][Pr_Pr.subset] > 0)))["TRUE"]
    classes[cell.type, "Pn_Pn"] <- table((ann[["promoter_ATAC-seq_lhs"]][Pr_Pr.subset] == 0) & 
                                           (ann[["promoter_ATAC-seq_rhs"]][Pr_Pr.subset] == 0))["TRUE"]

    # Dissect promoter to non-promoter interactions into finer categories.
    Pr_Np.subset <- (((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & ((ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0))) |
      (((ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0)) & ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0)))
    classes[cell.type, "Pa_Da"] <- table(((ann[["promoter_ATAC-seq_lhs"]][Pr_Np.subset] > 0) & 
                                            (ann[["distal_ATAC-seq_rhs"]][Pr_Np.subset] > 0)) |
                                           ((ann[["distal_ATAC-seq_lhs"]][Pr_Np.subset] > 0) & 
                                              (ann[["promoter_ATAC-seq_rhs"]][Pr_Np.subset] > 0)))["TRUE"]
    classes[cell.type, "Pa_Dn"] <- table(((ann[["promoter_ATAC-seq_lhs"]][Pr_Np.subset] > 0) & 
                                            (ann[["distal_ATAC-seq_rhs"]][Pr_Np.subset] == 0)) | 
                                           ((ann[["distal_ATAC-seq_lhs"]][Pr_Np.subset] == 0) & 
                                              (ann[["promoter_ATAC-seq_rhs"]][Pr_Np.subset] > 0)))["TRUE"]
    classes[cell.type, "Pn_Da"] <- table(((ann[["promoter_lhs"]][Pr_Np.subset] > 0) & 
                                            (ann[["promoter_ATAC-seq_lhs"]][Pr_Np.subset] == 0) & 
                                            (ann[["distal_ATAC-seq_rhs"]][Pr_Np.subset] > 0)) |
                                           ((ann[["distal_ATAC-seq_lhs"]][Pr_Np.subset] > 0) & 
                                              (ann[["promoter_rhs"]][Pr_Np.subset] > 0) & 
                                              (ann[["promoter_ATAC-seq_rhs"]][Pr_Np.subset] == 0)))["TRUE"]
    classes[cell.type, "Pn_Dn"] <- table(((ann[["promoter_lhs"]][Pr_Np.subset] > 0) & 
                                            (ann[["promoter_ATAC-seq_lhs"]][Pr_Np.subset] == 0) & 
                                            (ann[["distal_ATAC-seq_rhs"]][Pr_Np.subset] == 0)) |
                                           ((ann[["distal_ATAC-seq_lhs"]][Pr_Np.subset] == 0) & 
                                              (ann[["promoter_rhs"]][Pr_Np.subset] > 0) & 
                                              (ann[["promoter_ATAC-seq_rhs"]][Pr_Np.subset] == 0)))["TRUE"]
    
  }
  
  # Write results to file.
  write.table(classes, file=paste0(output.prefix, "/interaction.classes.all.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  
  # Prepare counts tables.
  counts.basic <- classes[cell.types, c("Pr_Pr", "Pr_Np")]
  counts.basic.prop <- prop.table(counts.basic, 1)
  counts.basic.colors <- c("#1E90FF", "#7B68EE")
  counts.detailed <- classes[cell.types, c("Pa_Pa", "Pa_Pn", "Pn_Pn", "Pa_Da", "Pa_Dn", "Pn_Da", "Pn_Dn")]
  counts.detailed.prop <- prop.table(counts.detailed, 1)
  counts.detailed.colors <- c("#1E90FF", "#00BFFF", "#87CEFA", "#7B68EE", "#9370DB", "#BA55D3", "#EE82EE")
  
  # Prepare stacked barplot (basic).
  pdf(paste0(output.prefix, "/interaction.classes.basic.pdf"), width=4.5, height=6)
  par(mar=c(4.5, 4.5, 4.5, 6), xpd=TRUE)
  b <- barplot(t(counts.basic.prop), col=counts.basic.colors, axes=FALSE, ylim=c(0, 1), las=2, width=0.75,
               ylab="% of interactions", main="interaction classes basic")
  y <- apply(t(counts.basic.prop), 2, cumsum) - 0.02
  text(sort(rep(b, length(counts.basic.colors))), c(y), c(t(counts.basic)))
  axis(2, at=seq(0, 1, 0.1), labels=paste0(seq(0, 100, 10), "%"))
  legend("right", legend=colnames(counts.basic.prop), fill=counts.basic.colors, cex=1, inset=c(-0.5, 0))
  dev.off()
  
  # Prepare stacked barplot (detailed).
  pdf(paste0(output.prefix, "/interaction.classes.detailed.pdf"), width=4.5, height=6)
  par(mar=c(4.5, 4.5, 4.5, 6), xpd=TRUE)
  b <- barplot(t(counts.detailed.prop), col=counts.detailed.colors, axes=FALSE, ylim=c(0, 1), las=2, width=0.75,
               ylab="% of interactions", main="interaction classes detailed")
  y <- apply(t(counts.detailed.prop), 2, cumsum) - 0.02
  text(sort(rep(b, length(counts.detailed.colors))), c(y), c(t(counts.detailed)))
  axis(2, at=seq(0, 1, 0.1), labels=paste0(seq(0, 100, 10), "%"))
  legend("right", legend=colnames(counts.detailed.prop), fill=counts.detailed.colors, cex=1, inset=c(-0.5, 0))
  dev.off()

}

# Plot ATAC-seq peak annotation results.
# Notes:
#
plotATACSeqPeakAnnotations <- function(annotation.data, output.prefix) {
  
  # Get counts for each category for each cell type.
  cell.types <- names(annotation.data)
  classes <- matrix(NA, length(cell.types), 4)
  rownames(classes) <- cell.types
  colnames(classes) <- c("Pr_I", "Di_I", "Pr_N", "Di_N")
  
  # Plot ATAC-seq peak annotation results in a series of pie charts.
  pdf(paste0(output.prefix, "/atac.seq.peak.classes.pdf"), width=4*length(names(annotation.data)), height=4)
  par(mfrow=c(1, length(names(annotation.data))))
  for (cell.type in cell.types) {
    
    # Count peaks in different classes.
    classes[cell.type, "Pr_I"] <- table((annotation.data[[cell.type]]$peak_type == "promoter") & (annotation.data[[cell.type]]$num_interactions > 0))["TRUE"]
    classes[cell.type, "Di_I"] <- table((annotation.data[[cell.type]]$peak_type == "distal") & (annotation.data[[cell.type]]$num_interactions > 0))["TRUE"]
    classes[cell.type, "Pr_N"] <- table((annotation.data[[cell.type]]$peak_type == "promoter") & (annotation.data[[cell.type]]$num_interactions == 0))["TRUE"]
    classes[cell.type, "Di_N"] <- table((annotation.data[[cell.type]]$peak_type == "distal") & (annotation.data[[cell.type]]$num_interactions == 0))["TRUE"]
    
    # Settings for pie chart.
    colors=list(promoter='#87CEFA', distal='#FCE6C9', interaction='#DC143C', none='#FFFFFF')
    iniR=0.8
    offset.rad=pi/2
    offset.angle=90
    
    # Start plotting pie chart.
    label.slices <- rev(c(classes[cell.type, "Pr_N"] + classes[cell.type, "Pr_I"], classes[cell.type, "Di_I"] + classes[cell.type, "Di_N"]))
    labels <- rev(c('promoter', 'distal'))
    labels <- paste0(labels, ", ", format(round(label.slices/sum(label.slices)*100, 1), nsmall=1), "%")
    pie(label.slices, radius=1.05*iniR, init.angle=offset.angle, col=rep('white', length(label.slices)), border=TRUE, labels=labels, cex=1, cex.main=1, main=paste0(cell.type, " atac-seq peak annotation"))
    floating.pie(0, 0, 1, radius=1*iniR, startpos=offset.rad, col='white', border=FALSE)
    
    # Plot interaction slices.
    interaction.slices <- rev(c(classes[cell.type, "Pr_N"], classes[cell.type, "Pr_I"] + classes[cell.type, "Di_I"], classes[cell.type, "Di_N"]))
    interaction.colors <- as.character(colors[rev(c('none', 'interaction', 'none'))])
    floating.pie(0, 0, interaction.slices, radius=0.95*iniR, startpos=offset.rad, col=interaction.colors, border=FALSE)
    
    # Plot annotation slices.
    annotation.slices <- rev(c(classes[cell.type, "Pr_N"] + classes[cell.type, "Pr_I"], classes[cell.type, "Di_I"] + classes[cell.type, "Di_N"]))
    annotation.colors <- as.character(colors[rev(c('promoter', 'distal'))])
    floating.pie(0, 0, annotation.slices, radius=0.85*iniR, startpos=offset.rad, col=annotation.colors, border=TRUE)
    
    # Plot blank interior.
    floating.pie(0, 0, 1, radius=0.45*iniR, startpos=offset.rad, col='black', border=TRUE)
    floating.pie(0, 0, 1, radius=0.446*iniR, startpos=offset.rad, col='white', border=FALSE)
    
  }
  dev.off()
  
  # Write results to file.
  write.table(classes, file=paste0(output.prefix, "/atac.seq.peak.classes.txt"), sep="\t", row.names=T, col.names=T, quote=F)
  
}

# Plot # of interactions by score cutoff.
# Notes:
#
plotNumInteractionsBySignificance <- function(interaction.data, p.values, xlim, ylim, cutoff, output.prefix) {
  
  # Count # of interactions as a function of score cutoff.
  cell.types <- names(interaction.data)
  nums <- matrix(NA, length(cell.types), length(p.values))
  for (i in 1:length(cell.types))
    for (j in 1:length(p.values))
      nums[i, j] <- length(interaction.data[[cell.types[i]]][interaction.data[[cell.types[i]]]$score > p.values[j], 1])
      
  # Plot results.
  pdf(paste0(output.prefix, "/number.interactions.by.score.", paste0(xlim, collapse="-"), ".", paste0(ylim, collapse="-"), ".pdf"), width=6, height=6)
  colors <- c("#3399FF", "#FF9966", "#66FF99", "#D15FEE")
  plot(1, type='n', xlim=xlim, ylim=ylim, xlab='score cutoff', ylab='# of interactions', main="# of interactions as a function of  score cutoff")
  for (i in 1:length(cell.types)) {
    lines(p.values, nums[i, ], type='l', lwd=2, col=colors[i])
  }
  axis(1, at=seq(xlim[1], xlim[2], 1), labels=seq(xlim[1], xlim[2], 1))
  if (is.null(cutoff) == FALSE)
    abline(v=cutoff, lty=2, col="red")
  legend("topright", legend=cell.types, lty=rep(1, length(cell.types)), col=colors[1:length(cell.types)])
  dev.off()
  
}

# Plot distance distribution of interactions for each cell type.
# Notes: 'type' can be "hist", "cdf", or "both". Default unit is kb.
#
plotDistances <- function(interaction.data, type, bin.size, dist.lim, colors, output.prefix) {
  
  # Get distances for each cell type.
  dists <- list()
  dists.merged <- c()
  for (cell.type in cell.types) {
    
    dists[[cell.type]] <- abs(interaction.data[[cell.type]][, 2] + interaction.data[[cell.type]][, 3] - interaction.data[[cell.type]][, 6] - interaction.data[[cell.type]][, 7])/2000
    dists.merged <- c(dists.merged, dists[[cell.type]])
    
  }
  
  # Prepare different plots based on the options.
  cell.types <- names(interaction.data)
  if (type == "hist") {
    
    # Plot interaction distance histograms for each cell type.
    pdf(paste0(output.prefix, "/interaction.distance.histogram.pdf"), width=4*length(cell.types), height=4)
    par(mfrow=c(1, length(cell.types)))
    for (i in 1:length(cell.types)) {
      
      # Get interaction distances for the current cell type.
      cell.type <- cell.types[i]
      max.dist <- ceiling(max(dists[[cell.type]])/bin.size)*bin.size
      
      # Plot histogram and mean interaction distance for the current cell type.
      hist(dists[[cell.type]], main=paste0(cell.type, " interaction distance histogram"),
           xlab="interaction distance (kb)", xlim=c(0, dist.lim), ylab="# interactions", breaks=seq(0, max.dist, bin.size), col=colors[i])
      abline(v=mean(dists[[cell.type]]), col="red")
      
    }
    dev.off()
    
  } else if (type == "cdf") {
    
    # Plot interaction distance CDF plot for each cell type (and also the merged set).
    pdf(paste0(output.prefix, "/interaction.distance.cdf.plot.pdf"), width=6, height=6)
    plot(ecdf(dists[[cell.types[1]]]), col=colors[1], main=NA,  xlim=c(0, dist.lim), ylim=c(0, 1), xlab="interaction distance (kb)", ylab="interaction distance CDF plot")
    for (i in 2:length(cell.types))
      plot(ecdf(dists[[cell.types[i]]]), col=colors[i], add=T)
    plot(ecdf(dists.merged), col="#555555", lty=2, add=T)
    legend('right', c(cell.types, "merged"), fill=c(colors, "#555555"), border=NA)
    dev.off()
    
  }
  
}

# Plot histograms of the # of interactions per gene.
# Notes:
#
plotGenePromoterHubsHistogram <- function(hub.results, colors, output.prefix, output.suffix) {
  
  # Plot histograms of # of interactions coincident on promoters for each gene.
  pdf(paste0(output.prefix, "/gene.promoter.hubs.histogram.", output.suffix, ".pdf"), width=4*length(cell.types), height=4)
  par(mfrow=c(1, length(cell.types)))
  for (i in 1:length(cell.types)) {

    cell.type <- cell.types[i]
    h <- hist(hub.results[[cell.type]]$num_interactions_all, main=paste0("# of interactions per gene for cell type: ", cell.type),
              xlab="# interactions", xlim=c(0, 25), ylab="# genes", breaks=seq(0, 1000, 1), col=colors[i], cex.main=1)
    abline(v=mean(hub.results[[cell.type]]$num_interactions_all), col="red")
    print(paste0(mean(hub.results[[cell.type]]$num_interactions_all), " mean interactions per promoter for cell type: ", cell.type))
    print(sum(h$counts))
    write.table(h$counts, file=paste0(output.prefix, "/gene.promoter.hubs.histogram.", output.suffix, ".", cell.type, ".txt"), sep="\t", row.names=T, col.names=T, quote=F)
    
  }
  dev.off()
  
}

# Plot Venn diagrams of cell type specific interactions.
# Notes: 
#
plotSpecificityVenn <- function(specificity.table, colors, output.prefix) {
  
  # Get unique combinations of cell types.
  cell.types <- colnames(specificity.table)
  combinations <- c()
  for (i in 1:length(cell.types)) {
    
    temp <- combn(cell.types, i)
    for (j in 1:length(temp[1, ]))
      combinations <- c(combinations, paste0(temp[, j], collapse="&"))
    
  }
  
  # Get counts for each combination.
  counts <- rep(0, length(combinations))
  specificity.pattern <- matrix(0, dim(specificity.table)[1], dim(specificity.table)[2])
  specificity.pattern[specificity.table] <- 1
  specificity.pattern <- apply(specificity.pattern, 1, paste0, collapse="")
  combination.ids <- list()
  for (i in 1:length(combinations)) {
    
    # Get boolean pattern for current combination, then count the number of matches.
    current.types <- unlist(strsplit(combinations[i], split="&"))
    pattern <- rep(FALSE, length(cell.types))
    pattern[cell.types %in% current.types] <- TRUE
    check <- rep(0, length(cell.types))
    check[pattern] <- 1
    check <- paste0(check, collapse="")
    counts[i] <- length(which(specificity.pattern == check))
    combination.ids[[combinations[i]]] <- rownames(specificity.matrix)[which(specificity.pattern == check)]
    
  }
  
  # Make circle Venn diagram.
  pdf(paste0(output.prefix, "/interaction_specificity_venn_1.pdf"), height=6, width=6)
  values <- counts
  names(values) <- combinations
  fit <- euler(values)
  print(plot(fit, fill_opacity=1, alpha=0.5, shape="circle", title="Cell type specificity of interactions",
             fill=c(colors[1:length(cell.types)]), border="transparent", fontsize=12, quantities=list(fontsize=12)))
  dev.off()
  
  # Transform the data into a list with unique IDs for counting.
  bins <- list()
  counter <- 1
  for (i in 1:length(combinations)) {
    
    current.types <- unlist(strsplit(combinations[i], split="&"))
    for (current.type in current.types)
      bins[[current.type]] <- c(bins[[current.type]], counter:(counter+counts[i]))
    counter <- counter + counts[i]
    
  }
  
  # Make elliptical Venn diagram.
  pdf(paste0(output.prefix, "/interaction_specificity_venn_2.pdf"), height=6, width=6)
  vp <- venn.diagram(bins, alpha=0.33, filename=NULL, fill=c(colors[1:length(cell.types)]), 
                     main.fontfamily="sans", cat.fontfamily="sans", fontfamily="sans")
  grid.draw(vp)
  dev.off()
  
}

# Print out basic results for Vista element analysis.
# Notes:
#
plotVistaResults <- function(vista.data, vista.results, cns.annotations, colors, output.prefix) {
  
  # Get a list of cell types in the results table.
  cell.types <- gsub("_interactions", "", grep("_interactions", colnames(vista.results), value=TRUE))
  
  # Annotate CNS-related Vista elements.
  num.elements <- length(vista.data[, 1])
  vista.data <- cbind(vista.data[, 1:5], rep(FALSE, num.elements))
  colnames(vista.data) <- c(colnames(vista.data)[1:5], "CNS_related")
  vista.data$CNS_related <- as.logical(vista.data$CNS_related)
  for (i in 1:num.elements) {
    
    current.annotation <- vista.data[i, "annotation"]
    current.annotation <- gsub("ganglion, cranial", "ganglion|cranial", current.annotation)
    current.annotation <- unlist(strsplit(current.annotation, split=","))
    for (j in 1:length(current.annotation))
      current.annotation[j] <- unlist(strsplit(current.annotation[j], split="\\["))[1]
    if (any(current.annotation %in% cns.annotations))
      vista.data$CNS_related[i] <- TRUE
    
  }
  
  # Print the # of CNS-related Vista elements.
  num.cns <- length(which(vista.data$CNS_related == TRUE))
  print(paste0(num.cns, " of ", num.elements, " Vista elements have CNS annotations."))
  cns.ids <- vista.data$vista_ID[vista.data$CNS_related == TRUE]
  
  # Count the # of interacting Vista elements per cell type.
  interacting.vista.ids.by.cell.type <- list()
  interacting.vista.ids.all <- c()
  for (cell.type in cell.types) {
    
    interacting.vista.ids.by.cell.type[[cell.type]] <- rownames(vista.results)[which(vista.results[, paste0(cell.type, "_interactions")] > 0)]
    interacting.vista.ids.all <- c(interacting.vista.ids.all, interacting.vista.ids.by.cell.type[[cell.type]])
    
  }
  
  # Print the # of interacting Vista elements across all cell types and plot a Venn diagram showing the intersections between sets.
  print(paste0("There are a total of ", length(unique(interacting.vista.ids.all)), " interacting Vista elements across all cell types."))
  print(paste0("There are a total of ", length(unique(interacting.vista.ids.all[interacting.vista.ids.all %in% cns.ids])), " CNS-related interacting Vista elements across all cell types."))
  pdf(paste0(output.prefix, "/interacting.vista.elements.venn.pdf"), width=6, height=6) 
  vp <- venn.diagram(interacting.vista.ids.by.cell.type, alpha = 0.33, filename = NULL, fill=c(colors[1:length(cell.types)]), 
                     main.fontfamily="sans", cat.fontfamily="sans", fontfamily="sans")
  grid.draw(vp)  
  dev.off()
  
  # Count the # of Vista elements overlapping ATAC-seq peaks per cell type.
  atac.seq.vista.ids.by.cell.type <- list()
  atac.seq.vista.ids.all <- c()
  for (cell.type in cell.types) {
    
    atac.seq.vista.ids.by.cell.type[[cell.type]] <- rownames(vista.results)[which(vista.results[, paste0(cell.type, "_atac_seq_peak")] == TRUE)]
    atac.seq.vista.ids.all <- c(atac.seq.vista.ids.all, atac.seq.vista.ids.by.cell.type[[cell.type]])
    
  }
  
  # Print the # of Vista elements overlapping ATAC-seq peaks across all cell types and plot a Venn diagram showing the intersections between sets.
  print(paste0("There are a total of ", length(unique(atac.seq.vista.ids.all)), " Vista elements overlapping ATAC-seq peaks across all cell types."))
  print(paste0("There are a total of ", length(unique(atac.seq.vista.ids.all[atac.seq.vista.ids.all %in% cns.ids])), " CNS-related Vista elements overlapping ATAC-seq peaks across all cell types."))
  pdf(paste0(output.prefix, "/atac.seq.vista.elements.venn.pdf"), width=6, height=6) 
  vp <- venn.diagram(atac.seq.vista.ids.by.cell.type, alpha = 0.33, filename = NULL, fill=c(colors[1:length(cell.types)]), 
                     main.fontfamily="sans", cat.fontfamily="sans", fontfamily="sans")
  grid.draw(vp)  
  dev.off()
  
  # Write table of results to file.
  vista.results.compiled <- cbind(vista.data, vista.results)
  write.table(vista.results.compiled, file=paste0(output.prefix, "/vista.results.all.txt"), sep="\t", row.names=F, col.names=F, quote=F)
  vista.results.positive.compiled <- vista.results.compiled[vista.results.compiled$vista_ID %in% interacting.vista.ids.all, ]
  write.table(vista.results.positive.compiled, file=paste0(output.prefix, "/vista.results.positive.txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
}

# Plots # of cases where the Vista element interacts with its nearest gene, etc.
# Notes:
#
printVistaConcordance <- function(vista.results, colors, output.prefix) {
  
  # Get a list of cell types in the results table.
  cell.types <- gsub("_interactions", "", grep("_interactions", colnames(vista.results), value=TRUE))
  
  # Get the positive elements.
  num.elements <- length(vista.data[, 1])
  interacting.vista.ids.by.cell.type <- list()
  interacting.vista.ids.all <- c()
  for (cell.type in cell.types)
    interacting.vista.ids.all <- c(interacting.vista.ids.all, rownames(vista.results)[which(vista.results[, paste0(cell.type, "_interactions")] > 0)])
  vista.results.interacting <- vista.results[rownames(vista.results) %in% interacting.vista.ids.all, ]
  num.interacting <- length(vista.results.interacting[, 1])
  
  # Get the interacting positive elements without same fragment ambiguity.
  print(table(vista.results.interacting$same_fragment_ambiguity))

  # Annotate cases where there is same fragment ambiguity, nearest gene concordance, or no nearest gene concordance, as well as the # of other target genes in each case.
  same.fragment.ambiguity <- vista.results.interacting$same_fragment_ambiguity
  nearest.gene.concordance <- rep(FALSE, num.interacting)
  num.other.targets <- rep(0, num.interacting)
  for (i in 1:num.interacting) {
    
    # Get unique target genes across all cell types, and determine whether the nearest gene is one of these target genes.
    target.genes.all <- unique(unlist(strsplit(paste0(vista.results.interacting[i, paste0(cell.types, "_target_genes")], collapse=","), split=",")))
    target.genes.all <- target.genes.all[target.genes.all != ""]
    if (vista.results.interacting$nearest_gene[i] %in% target.genes.all)
      nearest.gene.concordance[i] <- TRUE

    # Count the # of unique target genes that is not the nearest gene
    other.target.genes <- target.genes.all[!(target.genes.all %in% vista.results.interacting$nearest_gene[i])]
    num.other.targets[i] <- length(other.target.genes)
    
  }
  
  # Print cases and their corresponding # of other targets.
  print(paste0("same fragment ambiguity: ", length(which(same.fragment.ambiguity))))
  print(paste0("same fragment ambiguity # other targets: ", sum(num.other.targets[which(same.fragment.ambiguity)])))
  print(paste0("nearest gene only: ", length(which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets == 0)))))
  print(paste0("nearest gene only # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets == 0))])))
  print(paste0("nearest gene plus other: ", length(which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets > 0)))))
  print(paste0("nearest gene plus other # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets > 0))])))
  print(paste0("no nearest gene: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance))))
  print(paste0("no nearest gene # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & !nearest.gene.concordance)])))
  
  # Print more detailed breakdowns of certain cases.
  print(paste0("same fragment ambiguity + no interaction: ", length(which(same.fragment.ambiguity & (num.other.targets == 0)))))
  print(paste0("same fragment ambiguity + interaction: ", length(which(same.fragment.ambiguity & (num.other.targets > 0)))))
  print(paste0("no nearest gene + no interaction: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance & (num.other.targets == 0)))))
  print(paste0("no nearest gene + interaction: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance & (num.other.targets > 0)))))

}

# Print out basic results for Vista element analysis.
# Notes:
#
plotSNPResults <- function(disease.data, all.snp.results, tag.snp.results, colors, disease, output.prefix) {
  
  # Get a list of cell types in the results table.
  cell.types <- gsub("_interactions", "", grep("_interactions", colnames(all.snp.results), value=TRUE))
  num.all.snps <- length(all.snp.results[, 1])
  num.tag.snps <- length(tag.snp.results[, 1])
  
  # Print a table of how many cell types have an imputed interaction for each tag SNP.
  print("# of interacting cell types across all tag SNPs:")
  print(table(rowSums(tag.snp.results[, paste0(cell.types, "_interactions_imputed")])))
  
  # Count the # of interacting tag SNPs per cell type.
  interacting.snp.ids.by.cell.type <- list()
  interacting.snp.ids.all <- c()
  for (cell.type in cell.types) {
    
    interacting.snp.ids.by.cell.type[[cell.type]] <- rownames(tag.snp.results)[which(tag.snp.results[, paste0(cell.type, "_interactions_imputed")] == TRUE)]
    interacting.snp.ids.all <- c(interacting.snp.ids.all, interacting.snp.ids.by.cell.type[[cell.type]])
    
  }
  
  # Print the # of interacting tag SNPs across all cell types and plot a Venn diagram showing the intersections between sets.
  num.interacting.snps <- length(unique(interacting.snp.ids.all))
  print(paste0("There are a total of ", num.interacting.snps, "/", num.tag.snps, " interacting tag SNPs across all cell types."))
  pdf(paste0(output.prefix, "/interacting.snps.venn.", disease, ".pdf"), width=6, height=6) 
  vp <- venn.diagram(interacting.snp.ids.by.cell.type, alpha = 0.33, filename = NULL, fill=c(colors[1:length(cell.types)]), 
                     main.fontfamily="sans", cat.fontfamily="sans", fontfamily="sans", 
                     main=paste0(num.interacting.snps, "/", num.tag.snps, " interacting tag SNPs for disease: ", disease))
  grid.draw(vp)  
  dev.off()
  
  # Write table of results to file.
  if (!all(disease.data$rsid == rownames(all.snp.results)))
    print("Error")
  all.snp.results.compiled <- cbind(disease.data[, 1:9], all.snp.results)
  write.table(all.snp.results.compiled, file=paste0(output.prefix, "/all.snp.results.", disease, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  tag.snp.data <- disease.data[disease.data$is_query_snp == TRUE, ]
  if (!all(tag.snp.data$rsid == rownames(tag.snp.results)))
    print("Error")
  tag.snp.results.compiled <- cbind(tag.snp.data[, 1:9], tag.snp.results)
  write.table(tag.snp.results.compiled, file=paste0(output.prefix, "/tag.snp.results.", disease, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
  # Write results with ATAC-seq peaks only.
  selector <- apply(all.snp.results.compiled[, c(paste0(cell.types, "_promoter_atac_seq_peak"), paste0(cell.types, "_distal_atac_seq_peak"))], 1, any)  
  print(paste0("# of SNPs overlapping promoter or distal ATAC-seq peaks: ", length(which(selector)), "/", num.all.snps))
  all.snp.results.compiled.filtered <- all.snp.results.compiled[selector, ]
  write.table(all.snp.results.compiled.filtered, file=paste0(output.prefix, "/all.snp.results.filtered.", disease, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
}

# Get lists of target genes across all cell types for each disease.
# Notes:
#
getTargetGenes <- function(all.snp.results, promoters, expression.data, expression.threshold, output.prefix) {
  
  # Write lists of gene targets for each disease to file.
  cell.types <- gsub("_interactions", "", grep("_interactions", colnames(all.snp.results), value=TRUE))
  target.genes.all <- unique(unlist(strsplit(paste0(as.vector(as.matrix(all.snp.results[, paste0(cell.types, "_target_genes")])), collapse=","), split=",")))
  target.genes.all <- target.genes.all[target.genes.all != ""]
  print(paste0("Total # of unique gene targets: ", length(target.genes.all)))
  target.genes.filtered <- c()
  for (i in 1:length(target.genes.all)) {
    
    if (target.genes.all[i] %in% promoters$gene_name)
      if (any(expression.data[unique(promoters$gene_id[promoters$gene_name == target.genes.all[i]]), ] > expression.threshold))
        target.genes.filtered <- c(target.genes.filtered, target.genes.all[i])
    
  }
  print(paste0("Total # of unique gene targets that are expressed: ", length(target.genes.filtered)))
  write.table(target.genes.filtered, file=paste0(output.prefix, "/target.genes.all.", disease, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
}

# Plots # of cases where a tag SNP interacts with its nearest gene, etc. (may be through a linked SNP).
# Notes:
#
printSNPConcordance <- function(tag.snp.results) {
  
  # Get a list of cell types in the results table.
  cell.types <- gsub("_interactions_imputed", "", grep("_interactions_imputed", colnames(tag.snp.results), value=TRUE))
  
  # Get the interacting tag SNPs.
  num.tag.snps <- length(tag.snp.results[, 1])
  interacting.snp.ids.by.cell.type <- list()
  interacting.snp.ids.all <- c()
  for (cell.type in cell.types)
    interacting.snp.ids.all <- c(interacting.snp.ids.all, rownames(tag.snp.results)[which(tag.snp.results[, paste0(cell.type, "_interactions_imputed")] == TRUE)])
  tag.snp.results.interacting <- tag.snp.results[rownames(tag.snp.results) %in% interacting.snp.ids.all, ]
  num.interacting <- length(tag.snp.results.interacting[, 1])
  
  # Annotate cases where there is same fragment ambiguity, nearest gene concordance, or no nearest gene concordance, as well as the # of other target genes in each case.
  same.fragment.ambiguity <- tag.snp.results.interacting$same_fragment_ambiguity
  nearest.gene.concordance <- rep(FALSE, num.interacting)
  num.other.targets <- rep(0, num.interacting)
  for (i in 1:num.interacting) {
    
    # Get unique target genes across all cell types, and determine whether the nearest gene is one of these target genes.
    target.genes.all <- unique(unlist(strsplit(paste0(tag.snp.results.interacting[i, paste0(cell.types, "_target_genes_imputed")], collapse=","), split=",")))
    target.genes.all <- target.genes.all[target.genes.all != ""]
    if (tag.snp.results.interacting$nearest_gene[i] %in% target.genes.all)
      nearest.gene.concordance[i] <- TRUE
    
    # Count the # of unique target genes that is not the nearest gene
    other.target.genes <- target.genes.all[!(target.genes.all %in% tag.snp.results.interacting$nearest_gene[i])]
    num.other.targets[i] <- length(other.target.genes)
    
  }
  
  # Print cases and their corresponding # of other targets.
  print(paste0("same fragment ambiguity: ", length(which(same.fragment.ambiguity))))
  print(paste0("same fragment ambiguity # other targets: ", sum(num.other.targets[which(same.fragment.ambiguity)])))
  print(paste0("nearest gene only: ", length(which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets == 0)))))
  print(paste0("nearest gene only # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets == 0))])))
  print(paste0("nearest gene plus other: ", length(which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets > 0)))))
  print(paste0("nearest gene plus other # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & nearest.gene.concordance & (num.other.targets > 0))])))
  print(paste0("no nearest gene: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance))))
  print(paste0("no nearest gene # other targets: ", sum(num.other.targets[which(!same.fragment.ambiguity & !nearest.gene.concordance)])))
  
  # Print more detailed breakdowns of certain cases.
  print(paste0("same fragment ambiguity + no interaction: ", length(which(same.fragment.ambiguity & (num.other.targets == 0)))))
  print(paste0("same fragment ambiguity + interaction: ", length(which(same.fragment.ambiguity & (num.other.targets > 0)))))
  print(paste0("no nearest gene + no interaction: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance & (num.other.targets == 0)))))
  print(paste0("no nearest gene + interaction: ", length(which(!same.fragment.ambiguity & !nearest.gene.concordance & (num.other.targets > 0)))))
  
}

# Prints out statistics about the #s of unique genes interacting in each cell type.
# Notes:
#
printInteractingGenes <- function(hub.results, promoters, colors, output.prefix) {
  
  # Get lists of protein coding and noncoding RNA gene IDs.
  protein.coding.ids <- unique(promoters$gene_id[promoters$gene_type == "protein_coding"])
  noncoding.rna.ids <- unique(promoters$gene_id[promoters$gene_type != "protein_coding"])
  
  # Print # of interacting genes per cell type.
  interacting.gene.ids.by.cell.type <- list()
  interacting.gene.ids.all <- c()
  for (cell.type in cell.types) {
    
    interacting.gene.ids <- rownames(hub.results[[cell.type]])[hub.results[[cell.type]]$num_interactions_all > 0]
    print(paste0(length(interacting.gene.ids), " genes interacting in cell type: ", cell.type))
    print(paste0(length(interacting.gene.ids[interacting.gene.ids %in% protein.coding.ids]), " protein coding genes interacting in cell type: ", cell.type))
    print(paste0(length(interacting.gene.ids[interacting.gene.ids %in% noncoding.rna.ids]), " noncoding RNA genes interacting in cell type: ", cell.type))
    interacting.gene.ids.by.cell.type[[cell.type]] <- interacting.gene.ids
    interacting.gene.ids.all <- c(interacting.gene.ids.all, interacting.gene.ids)
    
  }
  
  # Print # of interacting genes across all cell types.
  print(paste0(length(unique(interacting.gene.ids.all)), " genes interacting across all cell types."))
  print(paste0(length(unique(interacting.gene.ids.all[interacting.gene.ids.all %in% protein.coding.ids])), " protein coding genes interacting across all cell types"))
  print(paste0(length(unique(interacting.gene.ids.all[interacting.gene.ids.all %in% noncoding.rna.ids])), " noncoding RNA genes interacting across all cell types"))
  
  # Print Venn diagrams of interacting genes.
  pdf(paste0(output.prefix, "/interacting.genes.venn.pdf"), width=6, height=6) 
  vp <- venn.diagram(interacting.gene.ids.by.cell.type, alpha = 0.33, filename = NULL, fill=c(colors[1:length(cell.types)]), 
                     main.fontfamily="sans", cat.fontfamily="sans", fontfamily="sans")
  grid.draw(vp)  
  dev.off()
  
}

# Write interactions and interaction annotations to file.
# Notes:
#
writeInteractions <- function(interactions.sig, interactions.ann, features, promoters.all, specificity.entries, omit.list, output.prefix) {
  
  # Prepare a dictionary of feature types/names for translating feature IDs.
  cell.types <- names(interactions.sig)
  feature.list <- c()
  for (cell.type in cell.types)
    for (feature.type in names(features[[cell.type]]))
      for (feature.name in names(features[[cell.type]][[feature.type]]))
        feature.list <- rbind(feature.list, c(cell.type, feature.type, feature.name))
      
  # Define order of columns.
  annotations <- c("promoter_lhs", "promoter_lhs_ids", "promoter_other_lhs", "promoter_other_lhs_ids", "promoter_ATAC-seq_lhs",
                   "promoter_rhs", "promoter_rhs_ids", "promoter_other_rhs", "promoter_other_rhs_ids", "promoter_ATAC-seq_rhs",
                    "distal_ATAC-seq_lhs", "distal_ATAC-seq_rhs",
                   paste0(names(features[[cell.types[1]]]$vista), "_lhs"), paste0(names(features[[cell.types[1]]]$vista), "_lhs_ids"),
                   paste0(names(features[[cell.types[1]]]$vista), "_rhs"), paste0(names(features[[cell.types[1]]]$vista), "_rhs_ids"),
                   c(rbind(paste0(names(features[[cell.types[1]]]$snp), "_lhs"), paste0(names(features[[cell.types[1]]]$snp), "_lhs_ids"),
                           paste0(names(features[[cell.types[1]]]$snp), "_rhs"), paste0(names(features[[cell.types[1]]]$snp), "_rhs_ids"))),
                   paste0(names(features[[cell.types[1]]]$hars), "_lhs"), paste0(names(features[[cell.types[1]]]$hars), "_lhs_ids"),
                   paste0(names(features[[cell.types[1]]]$hars), "_rhs"), paste0(names(features[[cell.types[1]]]$hars), "_rhs_ids"))
  for (to.omit in omit.list)
    annotations <- annotations[annotations != to.omit]
  
  # Get a list of included diseases.
  included.diseases <- names(features[[cell.types[1]]]$snp)[names(features[[cell.types[1]]]$snp) %in% gsub("_lhs_ids", "", annotations[grep("_lhs_ids", annotations)])]
  
  # Make a mapping of gene IDs to gene names.
  promoters.reference <- unique(promoters.all[, c("gene_id", "gene_name")])
  
  # Write data for each cell type.
  merged.data.all <- list()
  for (cell.type in cell.types) {
    
    # Assemble interaction and annotation data.
    merged.data <- interactions.sig[[cell.type]][, c(1:10)]
    for (annotation in annotations) {
      
      merged.data <- cbind(merged.data, interactions.ann[[cell.type]][[annotation]])
      colnames(merged.data)[length(merged.data[1, ])] <- annotation
      if (grepl("_ids", annotation))
        merged.data[, annotation] <- as.character(merged.data[, annotation])
      
    }
    
    # Translate feature and promoter IDs.
    num.interactions <- length(merged.data[, 1])
    print(paste0(num.interactions, " interactions to translate IDs for..."))
    for (i in 1:num.interactions) {
      
      # Track progress.
      if (i %% 10000 == 0)
        print(paste0(i, " interactions processed"))
      
      # Translate promoter IDs.
      if (merged.data[i, "promoter_lhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "promoter_lhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, promoters.reference$gene_name[promoters.reference$gene_id == id])
        merged.data[i, "promoter_lhs_ids"] <- paste0(sort(unique(replacement.ids)), collapse=",")
        
      }
      if (merged.data[i, "promoter_rhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "promoter_rhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, promoters.reference$gene_name[promoters.reference$gene_id == id])
        merged.data[i, "promoter_rhs_ids"] <- paste0(sort(unique(replacement.ids)), collapse=",")
        
      }
      if (merged.data[i, "promoter_other_lhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "promoter_other_lhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, promoters.reference$gene_name[promoters.reference$gene_id == id])
        merged.data[i, "promoter_other_lhs_ids"] <- paste0(sort(unique(replacement.ids)), collapse=",")
        
      }
      if (merged.data[i, "promoter_other_rhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "promoter_other_rhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, promoters.reference$gene_name[promoters.reference$gene_id == id])
        merged.data[i, "promoter_other_rhs_ids"] <- paste0(sort(unique(replacement.ids)), collapse=",")
        
      }
      
      # Translate Vista IDs.
      if (merged.data[i, "vista_lhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "vista_lhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, features[[cell.type]][["vista"]][["vista"]]$vista_ID[features[[cell.type]][["vista"]][["vista"]]$feature_ID == id])
        merged.data[i, "vista_lhs_ids"] <- paste0(sort(replacement.ids), collapse=",")

      }
      if (merged.data[i, "vista_rhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "vista_rhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, features[[cell.type]][["vista"]][["vista"]]$vista_ID[features[[cell.type]][["vista"]][["vista"]]$feature_ID == id])
        merged.data[i, "vista_rhs_ids"] <- paste0(sort(replacement.ids), collapse=",")
        
      }
      
      # Translate HAR IDs.
      if (merged.data[i, "all_lhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "all_lhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, features[[cell.type]][["hars"]][["all"]]$V4[features[[cell.type]][["hars"]][["all"]]$feature_ID == id])
        merged.data[i, "all_lhs_ids"] <- paste0(sort(replacement.ids), collapse=",")
        
      }
      if (merged.data[i, "all_rhs_ids"] != "") {
        
        exploded.ids <- unlist(strsplit(merged.data[i, "all_rhs_ids"], split=","))
        replacement.ids <- c()
        for (id in exploded.ids)
          replacement.ids <- c(replacement.ids, features[[cell.type]][["hars"]][["all"]]$V4[features[[cell.type]][["hars"]][["all"]]$feature_ID == id])
        merged.data[i, "all_rhs_ids"] <- paste0(sort(replacement.ids), collapse=",")
        
      }
      
      # Translate rsids.
      for (disease in included.diseases) {
        
        # There may be multiple rsids per feature ID so we have to split and recombine twice.
        if (merged.data[i, paste0(disease, "_lhs_ids")] != "") {

          exploded.ids <- unlist(strsplit(merged.data[i, paste0(disease, "_lhs_ids")], split=","))
          replacement.ids <- c()
          for (id in exploded.ids)
            replacement.ids <- c(replacement.ids, features[[cell.type]][["snp"]][[disease]]$rsid[features[[cell.type]][["snp"]][[disease]]$feature_ID == id])
          merged.data[i, paste0(disease, "_lhs_ids")] <- paste0(sort(replacement.ids), collapse=",")
          merged.data[i, paste0(disease, "_lhs_ids")] <- paste0(sort(unique(unlist(strsplit(merged.data[i, paste0(disease, "_lhs_ids")], split=",")))), collapse=",")

        }
        if (merged.data[i, paste0(disease, "_rhs_ids")] != "") {

          exploded.ids <- unlist(strsplit(merged.data[i, paste0(disease, "_rhs_ids")], split=","))
          replacement.ids <- c()
          for (id in exploded.ids)
            replacement.ids <- c(replacement.ids, features[[cell.type]][["snp"]][[disease]]$rsid[features[[cell.type]][["snp"]][[disease]]$feature_ID == id])
          merged.data[i, paste0(disease, "_rhs_ids")] <- paste0(sort(replacement.ids), collapse=",")
          merged.data[i, paste0(disease, "_rhs_ids")] <- paste0(sort(unique(unlist(strsplit(merged.data[i, paste0(disease, "_rhs_ids")], split=",")))), collapse=",")

        }
        
      }
      
    }
    
    # Add category column.
    ann <- interactions.ann[[cell.type]]
    Pr_Pr.selector <- ((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & 
      ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0))
    Pr_Np.selector <- (((ann[["promoter_lhs"]] > 0) | (ann[["promoter_ATAC-seq_lhs"]] > 0)) & ((ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0)))
    Np_Pr.selector <- (((ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0)) & ((ann[["promoter_rhs"]] > 0) | (ann[["promoter_ATAC-seq_rhs"]] > 0)))
    Np_Np.selector <- (ann[["promoter_lhs"]] == 0) & (ann[["promoter_ATAC-seq_lhs"]] == 0) & 
      (ann[["promoter_rhs"]] == 0) & (ann[["promoter_ATAC-seq_rhs"]] == 0)
    merged.data <- cbind(merged.data[, 1:10], rep("", num.interactions), merged.data[, 11:length(merged.data)])
    merged.data[, 11] <- as.character(merged.data[, 11])
    colnames(merged.data)[11] <- "category"
    merged.data[Pr_Pr.selector, 11] <- "Pr_Pr"
    merged.data[Pr_Np.selector, 11] <- "Pr_Np"
    merged.data[Np_Pr.selector, 11] <- "Np_Pr"
    merged.data[Np_Np.selector, 11] <- "other"
    
    # Make specificity column.
    merged.data <- cbind(merged.data[, 1:11], rep(".", num.interactions), merged.data[, 12:length(merged.data)])
    merged.data[, 12] <- as.character(merged.data[, 12])
    colnames(merged.data)[12] <- "specificity"
    if (!is.null(specificity.entries)) {
      
      specificity.entries.cell.type <- specificity.entries[specificity.entries[, 1] %in% interactions.sig[[cell.type]]$interaction_ID, ]
      for (i in 1:num.interactions)
        merged.data[i, 12] <- specificity.entries.cell.type[specificity.entries.cell.type[, 1] == interactions.sig[[cell.type]]$interaction_ID[i], 2]
      
    }
    
    # Write table to file.
    write.table(sortInteractions(merged.data), file=paste0(output.prefix, "/significant.interactions.", cell.type, ".txt"), sep="\t", row.names=F, col.names=T, quote=F)
    merged.data.all[[cell.type]] <- merged.data
    
  }
  
  # Return table with merged interactions and annotations.
  return(merged.data.all)

}


# End ---------------------------------------------------------------------

