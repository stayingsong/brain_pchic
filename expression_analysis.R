#
# Inputs:
# Outputs:
getAllDifferentiallyExpressedGenes <- function(expression.counts, threshold, output.prefix) {
  
  # Compile expression counts for each cell type into a single data fraame.
  all.results <- rep(NA, length(expression.counts[[cell.types[1]]]))
  groups <- c()
  for (i in 1:length(cell.types)) {
    
    all.results <- cbind(all.results, expression.counts[[cell.types[i]]])
    groups <- c(groups, rep(cell.types[i], length(expression.counts[[cell.types[i]]])))
    
  }
  all.results <- all.results[, -1]
  rownames(all.results) <- rownames(expression.counts[[cell.types[1]]])
  
  # Read data into edgeR format and perform TMM normalization.
  groups <- factor(groups, ordered=T, levels=cell.types)
  y <- DGEList(counts=all.results, group=groups)
  counts.per.mi <- cpm(y)
  cpm.filter <- counts.per.mi > 1
  cpm.keep <- which(rowSums(cpm.filter) >= 2)
  y <- y[cpm.keep, ]
  y <- calcNormFactors(y, method="TMM")
  
  # Prepare MDS plot.
  colors <- c(cortical="#BFEFFF", hippocampal="#FFD39B", astrocyte="#BDFCC9")
  points <- c(cortical=20, hippocampal=20, astrocyte=20)
  pdf(paste0(output.prefix, "/gene_expression_MDS_plot.pdf"), width=6, height=6)
  plotMDS(y, col=colors[groups], main="MDS Plot", xlim=c(-8, 8), ylim=c(-8, 8), pch=points[groups])
  legend("topleft", legend=levels(groups), pch=points, col=colors, ncol=1)
  dev.off()
  
  # Prepare MD plot.
  pdf(paste0(output.prefix, "/gene_expression_MD_plot.pdf"), width=6, height=6)
  plotMD(cpm(y, log=TRUE), column=1)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
  
  #
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  
  #
  y <- estimateDisp(y, design, robust=TRUE)
  y$common.dispersion
  pdf(paste0(output.prefix, "/gene_expression_BCV_plot.pdf"), width=6, height=6)
  plotBCV(y)
  dev.off()
  
  #
  fit <- glmQLFit(y, design, robust=TRUE)
  head(fit$coefficients)
  pdf(paste0(output.prefix, "/gene_expression_QLDisp_plot.pdf"), width=6, height=6)
  plotQLDisp(fit)
  dev.off()
  
  #
  con <- makeContrasts(cortical.hippocampal=cortical-hippocampal, 
                       hippocampal.astrocyte=hippocampal-astrocyte,
                       cortical.astrocyte=cortical-astrocyte, 
                       levels=design)
  anov <- glmQLFTest(fit, contrast=con)
  d <- decideTests(anov, p.value=threshold)
  num.sig <- table(d)[["1"]]
  sig.results <- data.frame(topTags(anov, num.sig))
  
  #
  return(as.character(rownames(sig.results)))
  
}

#
# Inputs:
# Outputs:
getSpecificDifferentiallyExpressedGenes <- function(expression.counts, threshold, output.prefix) {
  
  #
  all.results <- rep(NA, length(expression.counts[[cell.types[1]]]))
  groups <- c()
  for (i in 1:length(cell.types)) {
    
    all.results <- cbind(all.results, expression.counts[[cell.types[i]]])
    groups <- c(groups, rep(cell.types[i], length(expression.counts[[cell.types[i]]])))
    
  }
  all.results <- all.results[, -1]
  rownames(all.results) <- rownames(expression.counts[[cell.types[1]]])
  
  # Read data into edgeR format and perform TMM normalization.
  groups <- factor(groups, ordered=T, levels=cell.types)
  y <- DGEList(counts=all.results, group=groups)
  counts.per.mi <- cpm(y)
  cpm.filter <- counts.per.mi > 1
  cpm.keep <- which(rowSums(cpm.filter) >= 2)
  y <- y[cpm.keep, ]
  y <- calcNormFactors(y, method="TMM")
  y <- estimateDisp(y)
  
  #
  de.genes <- list()
  de.genes.results <- list()
  for (cell.type in cell.types) {
    
    #
    other.types <- cell.types[!(cell.types %in% cell.type)]
    e1 <- exactTest(y, pair=c(cell.type, other.types[1]))
    e2 <- exactTest(y, pair=c(cell.type, other.types[2]))
    d1 <- data.frame(decideTests(e1, p.value=threshold))
    d2 <- data.frame(decideTests(e2, p.value=threshold))
    s1 <- c(rownames(d1)[which(d1 == 1)], rownames(d1)[which(d1 == -1)])
    s2 <- c(rownames(d2)[which(d2 == 1)], rownames(d2)[which(d2 == -1)])
    de.genes[[cell.type]] <- as.character(s1[s1 %in% s2])
    print(paste0("Discovered ", length(de.genes[[cell.type]]), " differentially expressed genes for cell type: ", cell.type))
    
    #
    de.genes.results[[cell.type]] <- data.frame(matrix(NA, length(de.genes[[cell.type]]), 7), stringsAsFactors=F)
    colnames(de.genes.results[[cell.type]]) <- c("gene_id", "logFC1", "logCPM1", "p-value1", "logFC2", "logCPM2", "p-value2")
    de.genes.results[[cell.type]][, 1] <- as.character(de.genes.results[[cell.type]][, 1])
    de.genes.results[[cell.type]][, 2] <- as.numeric(de.genes.results[[cell.type]][, 2])
    de.genes.results[[cell.type]][, 3] <- as.numeric(de.genes.results[[cell.type]][, 3])
    de.genes.results[[cell.type]][, 4] <- as.numeric(de.genes.results[[cell.type]][, 4])
    de.genes.results[[cell.type]][, 5] <- as.numeric(de.genes.results[[cell.type]][, 5])
    de.genes.results[[cell.type]][, 6] <- as.numeric(de.genes.results[[cell.type]][, 6])
    de.genes.results[[cell.type]][, 7] <- as.numeric(de.genes.results[[cell.type]][, 7])
    for (i in 1:length(de.genes[[cell.type]])) {
      
      de.genes.results[[cell.type]][i, 1] <- de.genes[[cell.type]][i]
      de.genes.results[[cell.type]][i, 2] <- e1$table[de.genes[[cell.type]][i], 1]
      de.genes.results[[cell.type]][i, 3] <- e1$table[de.genes[[cell.type]][i], 2]
      de.genes.results[[cell.type]][i, 4] <- e1$table[de.genes[[cell.type]][i], 3]
      de.genes.results[[cell.type]][i, 5] <- e2$table[de.genes[[cell.type]][i], 1]
      de.genes.results[[cell.type]][i, 6] <- e2$table[de.genes[[cell.type]][i], 2]
      de.genes.results[[cell.type]][i, 7] <- e2$table[de.genes[[cell.type]][i], 3]
      
    }
    
    #
    print(table((de.genes.results[[cell.type]][, 2] > 0) & (de.genes.results[[cell.type]][, 5] > 0)))
    print(table((de.genes.results[[cell.type]][, 2] > 0) & (de.genes.results[[cell.type]][, 5] < 0)))
    print(table((de.genes.results[[cell.type]][, 2] < 0) & (de.genes.results[[cell.type]][, 5] < 0)))
    
  }
  
  #
  return(de.genes)
  
}

