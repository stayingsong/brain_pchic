#
# allelic_bias_analysis.R
# Author: Michael Song
# Last modified: 2019-02-05
# This script analyzes contact matrices from HiC-Pro/Fit-Hi-C for allelic bias across significantly interacting loci.
#


# Load dependencies.
require(dplyr)


# Set up environment and load settings ------------------------------------


# Clear workspace before running script.
rm(list=ls())

# Turn off scientific notation for writing output.
options(scipen=999)

# Convert warnings into errors.
options(warn=2)

# Process CLI arguments (or override if running interactively).
cli.args = commandArgs(trailingOnly=TRUE)
print(cli.args[1]) # Set the home directory for the pipeline.
if ((exists("cli.args") == FALSE) | (length(cli.args) == 0)) {
  
  cli.args <- c()
  cli.args[1] <- "/Users/Michael/Box Sync/MS_general/analysis/brain_pchic_revision"
  
}

# Set home directory.
home.dir <- cli.args[1]


# Set up analysis environment ---------------------------------------------


# Source files should be in ~/scripts.
scripts.dir <- paste0(home.dir, "/scripts")
source(paste0(scripts.dir, "/utilities.R"))
source(paste0(scripts.dir, "/preprocess_data.R"))

# Define replicates and cell types.
replicates <- c("MS001", "MS002", "MS003", "MS137", "MS142", "MS136", "MS140", "MS141")
identities <- c("cortical", "cortical", "cortical", "hippocampal", "hippocampal", "motor", "motor", "motor")
cell.types <- unique(identities)
selector <- 1

# Set resolution and score cutoff.
resolution <- 10000
overlap.threshold <- 10^-9
score.cutoff <- 3

# Flags to determine which part of the script to run.
merge.replicates <- FALSE
filter.counts <- FALSE

# Set run ID.
run.id <- "brain_revision"

# Set folder containing processed input files from HiC-Pro/Fit-Hi-C.
processed.dir <- "/other/allelic_bias/unfiltered/processed"

# Set run folders.
output.dir <- paste0(home.dir, "/figures")
dir.list <- c("", 
              "/allelic_bias", "/allelic_bias/filtered", "/allelic_bias/results", "/allelic_bias/merged")
for (dir in dir.list)
  if (!file.exists(paste0(output.dir, "/", run.id, dir)))
    dir.create(paste0(output.dir, "/", run.id, dir))

# Read in VCF file.
vcf.data <- read.table(paste0(home.dir, "/resources/vcf/wtc_PASS_hg19.phased_v1.filtered.vcf"), sep="\t", header=F, stringsAsFactors=F)
vcf.bed <- cbind(vcf.data[, 1], vcf.data[, 2]-1, vcf.data[, 2])


# Merge replicates --------------------------------------------------------


# Merge data from different replicates. Data is in the format (chr1 pos1 chr2 pos2 N_reads) and was prepared using HiC-Pro and Fit-Hi-C.
if (merge.replicates) {
  
  # Merge replicates by cell type.
  for (cell.type in cell.types) {
    
    # Get replicates corresponding to the current cell type.
    print(paste0("Merging data from replicates for cell type: ", cell.type))
    current.replicates <- replicates[identities == cell.type]
    
    # Combine the count data across all replicates for the current cell type.
    G1.all <- c()
    G2.all <- c()
    for (current.replicate in current.replicates) {
      
      # Read in files of the format (chr1 pos1 chr2 pos2 N_reads) for the current replicate.
      print(paste0("Reading in data for: ", current.replicate))
      G1 <- read.table(paste0(home.dir, "/", processed.dir, "/", current.replicate, "_", resolution, "_G1.txt"), sep="\t", header=F, stringsAsFactors=F,
                       colClasses=c("character", "integer", "character", "integer", "integer"))
      G2 <- read.table(paste0(home.dir, "/", processed.dir, "/", current.replicate, "_", resolution, "_G2.txt"), sep="\t", header=F, stringsAsFactors=F,
                       colClasses=c("character", "integer", "character", "integer", "integer"))
      
      # Retain only cis interactions and concatenate the results.
      G1 <- G1[G1$V1 == G1$V3, ]
      G2 <- G2[G2$V1 == G2$V3, ]
      G1.all <- rbind(G1.all, G1)
      G2.all <- rbind(G2.all, G2)
      print(paste0(length(G1[, 1]), " cis interactions read in for G1"))
      print(paste0(length(G2[, 1]), " cis interactions read in for G2"))
      
    }
    
    # Combine duplicate rows, summing their counts.
    G1.merged <- G1.all %>% group_by_(.dots=c("V1", "V2", "V3", "V4")) %>% summarize(V5=sum(V5))
    G2.merged <- G2.all %>% group_by_(.dots=c("V1", "V2", "V3", "V4")) %>% summarize(V5=sum(V5))
    G1.merged <- data.frame(G1.merged)
    G2.merged <- data.frame(G2.merged)
    print(paste0(length(G1.merged[, 1]), " unique interactions across all replicates for G1"))
    print(paste0(length(G2.merged[, 1]), " unique interactions across all replicates for G2"))
    
    # Save results.
    save(G1.merged, file=paste0(output.dir, "/", run.id, "/allelic_bias/merged/G1.merged.", cell.type, ".", resolution, ".Rdata"))
    save(G2.merged, file=paste0(output.dir, "/", run.id, "/allelic_bias/merged/G2.merged.", cell.type, ".", resolution, ".Rdata"))
    
  }
  
}


# Filter the counts data by significantly interacting bins ----------------


# Retain only the count data for bins that are significantly interacting for each cell type.
if (filter.counts) {
  
  # Load and expand significant interactions to the specified resolution.
  load(file=paste0(home.dir, "/figures/", run.id, "/saved/interactions.raw.processed.Rdata"))
  interactions.sig.res <- list()
  for (cell.type in cell.types)
    interactions.sig.res[[cell.type]] <- expandInteractions(filterInteractionsBySignificance(interactions.raw.processed[[cell.type]], score.cutoff), 
                                                            resolution, neighboring.fragments = NULL)
  rm(interactions.raw.processed)
  gc()
  
  # Select cell type to be processed.
  cell.type <- cell.types[selector]
  print(paste0("Filtering count data for cell type: ", cell.type))
  current.interactions <- interactions.sig.res[[cell.type]]
  print(paste0("There are ", length(current.interactions[, 1]), " interactions for cell type ", cell.type, " with score cutoff ", score.cutoff, " and resolution ", resolution, "..."))
  
  # Load count data.
  load(file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/merged/G1.merged.", cell.type, ".", resolution, ".Rdata"))
  load(file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/merged/G2.merged.", cell.type, ".", resolution, ".Rdata"))
  
  # Convert count data to BED format.
  G1.bed <- data.frame(cbind(G1.merged[, 1], G1.merged[, 2] - resolution/2, G1.merged[, 2] + resolution/2, 
                             G1.merged[, 3], G1.merged[, 4] - resolution/2, G1.merged[, 4] + resolution/2, 
                             G1.merged[, 5]), stringsAsFactors=F)
  G2.bed <- data.frame(cbind(G2.merged[, 1], G2.merged[, 2] - resolution/2, G2.merged[, 2] + resolution/2, 
                             G2.merged[, 3], G2.merged[, 4] - resolution/2, G2.merged[, 4] + resolution/2, 
                             G2.merged[, 5]), stringsAsFactors=F)
  colnames(G1.bed) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")
  colnames(G2.bed) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")
  G1.bed[, 2] <- as.integer(G1.bed[, 2])
  G1.bed[, 3] <- as.integer(G1.bed[, 3])
  G1.bed[, 5] <- as.integer(G1.bed[, 5])
  G1.bed[, 6] <- as.integer(G1.bed[, 6])
  G1.bed[, 7] <- as.integer(G1.bed[, 7])
  G2.bed[, 2] <- as.integer(G2.bed[, 2])
  G2.bed[, 3] <- as.integer(G2.bed[, 3])
  G2.bed[, 5] <- as.integer(G2.bed[, 5])
  G2.bed[, 6] <- as.integer(G2.bed[, 6])
  G2.bed[, 7] <- as.integer(G2.bed[, 7])
  
  # Filter out interactions between the same bins.
  G1.bed <- G1.bed[!((G1.bed$start1 == G1.bed$start2) & (G1.bed$end1 == G1.bed$end2)), ]
  G2.bed <- G2.bed[!((G2.bed$start1 == G2.bed$start2) & (G2.bed$end1 == G2.bed$end2)), ]
  
  # Check that the count data is ordered. Should be TRUE.
  print(table(G1.bed$end1 <= G1.bed$start2))
  print(table(G2.bed$end1 <= G2.bed$start2))
  
  # Get all unique interacting bins from the count data and assign them IDs.
  unique.intervals <- unique(rbind(G1.bed[, 1:6], G2.bed[, 1:6]))
  unique.intervals <- cbind(unique.intervals[, 1:6], paste0("id", 1:length(unique.intervals[, 1])))
  colnames(unique.intervals)[7] <- "ID"
  
  # Grab the lhs and rhs bins for the unique interacting bins and the significant interactions.
  ids.1 <- as.character(current.interactions[, 11])
  ids.2 <- as.character(unique.intervals[, 7])
  interactions.lhs.1 <- current.interactions[, c(1:3, 11)]
  interactions.lhs.2 <- unique.intervals[, c(1:3, 7)]
  interactions.rhs.1 <- current.interactions[, c(5:7, 11)]
  interactions.rhs.2 <- unique.intervals[, c(4:6, 7)]
  
  # Filter out all unique interacting bins that don't overlap their respective ends of a significant interaction.
  filter.lhs <- bedTools.2in(bed1=interactions.lhs.2, bed2=interactions.lhs.1, opt.string="-c")
  filter.rhs <- bedTools.2in(bed1=interactions.rhs.2, bed2=interactions.rhs.1, opt.string="-c")
  interactions.lhs.2 <- interactions.lhs.2[filter.lhs[, 5] > 0, ]
  interactions.rhs.2 <- interactions.rhs.2[filter.rhs[, 5] > 0, ]
  
  # Intersect the unique interacting bins and the significant interactions and get the IDs of the former that do overlap significant interactions.
  overlap.results.lhs <- bedTools.2in(bed1=interactions.lhs.1, bed2=interactions.lhs.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
  overlap.results.rhs <- bedTools.2in(bed1=interactions.rhs.1, bed2=interactions.rhs.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
  interacting.bin.ids <- c()
  num.matches <- 0
  print(paste0(length(ids.1), " significant interactions to be processed..."))
  for (i in 1:length(ids.1)) {

    # Track progress.
    if (i %% 10000 == 0)
      print(paste0(i, " interactions processed"))

    # For the current significant interaction, get the IDs of unique interacting bins overlapping the lhs and rhs bins separately, then take the intersection.
    intersections.lhs <- as.character(overlap.results.lhs[overlap.results.lhs[, 4] == ids.1[i], 8])
    intersections.rhs <- as.character(overlap.results.rhs[overlap.results.rhs[, 4] == ids.1[i], 8])
    intersections.both <- intersections.lhs[intersections.lhs %in% intersections.rhs]
    if (length(intersections.both) > 0) {

      interacting.bin.ids <- c(interacting.bin.ids, intersections.both)
      num.matches <- num.matches + 1

    }

  }
  
  # Print out how many significant interactions were supported by count data at the current resolution (should be high unless very short interactions).
  print(paste0(num.matches, " of ", length(ids.1), " significant interactions had supporting count data"))
  
  # Print and save results.
  interacting.bin.ids <- unique(interacting.bin.ids)
  print(paste0(length(interacting.bin.ids), " unique interacting bins found to overlap significant interactions"))
  save(interacting.bin.ids, file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".interacting.bin.ids.Rdata"))
  load(file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".interacting.bin.ids.Rdata"))
  
  # Retrieve counts for filtered unique interacting bins.
  unique.intervals <- unique.intervals[unique.intervals$ID %in% interacting.bin.ids, 1:6]
  num.intervals <- length(unique.intervals[, 1])
  G1.counts <- cbind(unique.intervals, rep(NA, num.intervals))
  colnames(G1.counts)[7] <- "count"
  G2.counts <- cbind(unique.intervals, rep(NA, num.intervals))
  colnames(G2.counts)[7] <- "count"
  print(paste0(num.intervals, " unique intervals to retrieve count data for..."))
  for (i in 1:num.intervals) {
    
    # Track progress.
    if (i %% 500 == 0)
      print(paste0(i, " unique intervals processed"))
    
    # Grab the count data for each allele corresponding to the current unique interval.
    G1.value <- G1.bed[(G1.bed[, 1] == unique.intervals[i, 1]) & (G1.bed[, 2] == unique.intervals[i, 2]) & (G1.bed[, 3] == unique.intervals[i, 3]) &
                         (G1.bed[, 4] == unique.intervals[i, 4]) & (G1.bed[, 5] == unique.intervals[i, 5]) & (G1.bed[, 6] == unique.intervals[i, 6]), 7]
    if (length(G1.value) > 0)
      G1.counts[i, 7] <- G1.value
    else
      G1.counts[i, 7] <- 0
    G2.value <- G2.bed[(G2.bed[, 1] == unique.intervals[i, 1]) & (G2.bed[, 2] == unique.intervals[i, 2]) & (G2.bed[, 3] == unique.intervals[i, 3]) &
                         (G2.bed[, 4] == unique.intervals[i, 4]) & (G2.bed[, 5] == unique.intervals[i, 5]) & (G2.bed[, 6] == unique.intervals[i, 6]), 7]
    if (length(G2.value) > 0)
      G2.counts[i, 7] <- G2.value
    else
      G2.counts[i, 7] <- 0
    
    # Check to make sure that the count data for at least one allele is nonzero.
    if ((length(G1.value) == 0) && (length(G2.value) == 0))
      print("Error: # of supporting reads should be > 0 for at least one allele")
    
  }
  
  # Save results.
  save(G1.counts, file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".G1.counts.Rdata"))
  save(G2.counts, file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".G2.counts.Rdata"))
  
}


# Statistical analysis ----------------------------------------------------


# Load significantly interacting counts data for each allele.
cell.type <- cell.types[selector]
name <- paste0(cell.type, ".", resolution, ".", score.cutoff)
print(paste0("Analyzing data for cell type and resolution: ", cell.type, ", ", resolution))
load(file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".G1.counts.Rdata"))
load(file=paste0(home.dir, "/figures/", run.id, "/allelic_bias/filtered/", cell.type, ".", resolution, ".", score.cutoff, ".G2.counts.Rdata"))

# Plot number of counts across both alleles for each bin pair.
print(paste0("Mean counts per bin pair: ", mean(G1.counts$count + G2.counts$count)))
pdf(file=paste0(output.dir, "/", run.id, "/allelic_bias/results/", name,".counts.per.bin.pair.pdf"), height=6, width=6)
hist(G1.counts$count + G2.counts$count, xlim=c(0, 50), breaks=seq(0, 5000, 2))
abline(v=mean(G1.counts$count + G2.counts$count), col="red")
dev.off()

# Filter bin pairs by a minimum # of counts across both alleles.
min.counts <- 10
keep <- ((G1.counts$count + G2.counts$count) >= min.counts)
print(paste0("Number of bin pairs with at least ", min.counts, " reads: ", table(keep)[["TRUE"]], " (", 
             signif(table(keep)[["TRUE"]]*100/(length(which(keep)) + length(which(!keep))), 3), "%)"))
print(paste0("Number of bin pairs filtered out with fewer than ", min.counts, " reads: ", length(which(!keep))))
G1.filtered <- G1.counts[keep, ]
G2.filtered <- G2.counts[keep, ]
bins.filtered <- G1.filtered[, 1:6]
G1.counts.filtered <- G1.filtered$count
G2.counts.filtered <- G2.filtered$count

# Prepare values for qq plot.
pvalues <- 0
for (j in 1:length(G1.counts.filtered))
  pvalues[j] <- binom.test(G1.counts.filtered[j], G1.counts.filtered[j] + G2.counts.filtered[j], 0.5, alternative="two.sided")$p.value
observed <- -log10(sort(pvalues, decreasing=F))
expected <- -log10(ppoints(length(G1.counts.filtered)))

# Plot FDR adjusted plot.
pvalues.fdr <- p.adjust(pvalues, "BH")
significant.fdr <- (pvalues.fdr <= 0.05)
observed.fdr <- -log10(sort(pvalues.fdr, decreasing=F))
colors <- rep("black", length(significant.fdr))
colors[1:length(significant.fdr[significant.fdr == TRUE])] <- "purple"
pdf(file=paste0(output.dir, "/", run.id, "/allelic_bias/results/", name,".fdr.plot.pdf"), height=6, width=6)
plot(expected, observed.fdr, xlim=c(0, 5), ylim=c(0, 15), main=paste0(name, " allelic bias fdr plot"),
     sub=paste0(signif(length(which(significant.fdr))*100/length(significant.fdr), 5), "%, ", table(significant.fdr)[["TRUE"]], "/", length(significant.fdr), " bin pairs with >= ", min.counts, " counts"), col=colors)
abline(h=-log10(0.05), lty=2, col="red")
dev.off()

# Prepare qq plot, highlighting the points which pass the FDR cutoff above.
colors <- rep("black", length(significant.fdr))
colors[1:length(significant.fdr[significant.fdr])] <- "red"
jpeg(file=paste0(output.dir, "/", run.id, "/allelic_bias/results/", name,".qq.plot.with.fdr.jpg"), height=6, width=6, units="in", res=300)
plot(expected, observed, xlim=c(0, 5), ylim=c(0, 15), main=paste0(name, " allelic bias qq plot"),
     sub=paste0(signif(length(which(significant.fdr))*100/length(significant.fdr), 5), "%, ", table(significant.fdr)[["TRUE"]], "/", length(significant.fdr), " bin pairs with >= ", min.counts, " counts"), col=colors)
abline(c(0, 0), c(1, 1))
dev.off()

# Plot counts.
to.reverse <- (G1.counts.filtered > G2.counts.filtered)
higher <- G2.counts.filtered
higher[to.reverse] <- G1.counts.filtered[to.reverse]
lower <- G1.counts.filtered
lower[to.reverse] <- G2.counts.filtered[to.reverse]
colors <- rep("blue", length(to.reverse))
colors[to.reverse] <- "purple"
types <- rep(1, length(to.reverse))
types[to.reverse] <- 3
epsilon <- 10^-9
pdf(file=paste0(output.dir, "/", run.id, "/allelic_bias/results/", name,".allelic.counts.pdf"), height=6, width=6)
plot(lower + epsilon, higher + epsilon, xlim=c(0, 125), ylim=c(0, 125), col=colors, pch=types,
     xlab="lower frequency allele (# reads)", ylab="higher frequency allele (# reads)", main=paste0(cell.type, " allelic bias counts plot"))
points(lower[significant.fdr] + epsilon, higher[significant.fdr] + epsilon, col="red", pch=1, cex=2)
dev.off()

# Write significant bin interactions to file (use less stringent p-value cutoff).
significant <- pvalues <= 0.001
table(significant)
distances <- bins.filtered$start2[significant] - bins.filtered$end1[significant]
snps1 <- bedTools.2in(bed1=bins.filtered[significant, 1:3], bed2=vcf.bed, opt.string="-c")
snps2 <- bedTools.2in(bed1=bins.filtered[significant, 4:6], bed2=vcf.bed, opt.string="-c")
sig.bins <- cbind(bins.filtered[significant, ], distances, snps1[, 4], snps2[, 4],
                  G1.counts.filtered[significant], G2.counts.filtered[significant], pvalues[significant], pvalues.fdr[significant])
colnames(sig.bins) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "distance", "snps1", "snps2", "count1", "count2", "p-value", "q-value")
sig.bins <- sig.bins[order(sig.bins$'q-value'), ]
write.table(sig.bins, file=paste0(output.dir, "/", run.id, "/allelic_bias/results/", name, ".significant.interactions.pval.tsv"), sep="\t", row.names=F, col.names=F, quote=F)

# Save in interaction format.
allelic.interactions.sig <- cbind(sig.bins[, 1:3], rep(".", length(sig.bins[, 1])), sig.bins[, 4:6], rep(".", length(sig.bins[, 1])),
                                  apply(sig.bins[, 10:11], 1, paste0, collapse="|"), -log10(sig.bins$`p-value`), paste0(cell.type, "_", 1:length(sig.bins[, 1])))
colnames(allelic.interactions.sig) <- c("bait_chr", "bait_start", "bait_end", "bait_name", "otherEnd_chr", "otherEnd_start", "otherEnd_end", "otherEnd_name", "N_reads", "score", "interaction_ID")
save(allelic.interactions.sig, file=paste0(output.dir, "/", run.id, "/saved/allelic.interactions.sig.pval.", cell.type, ".Rdata"))


# End ---------------------------------------------------------------------

