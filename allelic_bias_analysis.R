#
# allelic_bias_analysis.R
# Author: Michael Song
# Date created: 2018-01-28
# Last modified: 2018-05-10
# This script analyzes the output from HiC-Pro and evaluates allelic bias in the set of significant interactions.
#


# Miscellaneous settings --------------------------------------------------


# Clear workspace.
rm(list=ls())

# Turn off scientific notation for output/printing.
options(scipen=999)


# Load packages -----------------------------------------------------------


# Install and load packages.
#source("https://bioconductor.org/biocLite.R")
#biocLite("HiTC")
require(HiTC)


# Functions ---------------------------------------------------------------


# Load auxiliary functions.
scripts.dir <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/scripts/pipeline/"
source(paste0(scripts.dir, "/utilities.R"))
source(paste0(scripts.dir, "/read_data.R"))
source(paste0(scripts.dir, "/preprocess_data.R"))
source(paste0(scripts.dir, "/annotation.R"))
source(paste0(scripts.dir, "/interaction_analysis.R"))
source(paste0(scripts.dir, "/expression_analysis.R"))
source(paste0(scripts.dir, "/prepare_plots.R"))


# Preprocessing -----------------------------------------------------------


# Set parameters.
resolution <- 10000
downsampling.depth <- 1000000
overlap.threshold <- 10^-9
score.cutoff <- 5
type <- "raw"
epsilon <- 1e-9
cell.type <- "cortical"
run.id <- "final"
name <- paste0(run.id, ".", resolution, ".", score.cutoff)

# Read in WTC SNPs file and process to BED format for intersections with interactions.
setwd("/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files")
snps <- read.table("wtc_PASS_hg19.haploseq.phased.vcf", sep="\t", skip=59, stringsAsFactors=F, header=T)
colnames(snps) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",	"FILTER",	"INFO",	"FORMAT",	"WTC")
genotype <- snps$WTC
genotype_split <- rep(NA, length(genotype))
for (i in 1:length(genotype)) {
  genotype_split[i] <- unlist(strsplit(genotype[i], split="\\:"))[1]
}
table(genotype_split)
snps.bed <- cbind(snps$CHROM, snps$POS, snps$POS+1, snps$ID, snps$QUAL, snps$FILTER, snps$FORMAT, snps$WTC)
head(snps.bed)

# Read in interactions.
interaction.file <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/pipeline_files/chicago_interactions/minNPerBaits_50k/hicup_0/cortical_all_350.ibed"
raw.interactions <- read.table(file=interaction.file, sep="\t", stringsAsFactors=F, header=T,
                               colClasses=c("character", "integer", "integer", "character", "character", "integer", "integer", "character", "integer", "numeric"))
filtered.interactions <- raw.interactions[raw.interactions[, 10] >= score.cutoff, ]
print(dim(filtered.interactions))
head(filtered.interactions)
expanded.interactions <- expandInteractions(filtered.interactions, resolution)
downsampled.interactions <- downsampleInteractions(expanded.interactions, downsampling.depth)
downsampled.interactions <- cbind(downsampled.interactions[, 1:10], paste0("cortical", 1:length(downsampled.interactions[, 1])))
colnames(downsampled.interactions)[11] <- "ID"

#
setwd("/Users/michael/Box Sync/MS_analysis/brain_pchic/figures/common/allelic_bias")

#
G1.path <- paste0("/Users/michael/Box Sync/MS_analysis/brain_pchic/other/allelic_bias_files/final/processed/", 
                  cell.type, "_", resolution, "_G1.txt")
G2.path <- paste0("/Users/michael/Box Sync/MS_analysis/brain_pchic/other/allelic_bias_files/final/processed/", 
                  cell.type, "_", resolution, "_G2.txt")
G1 <- read.table(G1.path, sep="\t", stringsAsFactors=F, header=F, colClasses=c("character", "integer", "character", "integer", "integer"))
G2 <- read.table(G2.path, sep="\t", stringsAsFactors=F, header=F, colClasses=c("character", "integer", "character", "integer", "integer"))
G1.bed <- data.frame(cbind(G1[, 1], G1[, 2]-resolution/2, G1[, 2]+resolution/2, G1[, 3], G1[, 4]-resolution/2, G1[, 4]+resolution/2, G1[, 5]), stringsAsFactors=F)
G2.bed <- data.frame(cbind(G2[, 1], G2[, 2]-resolution/2, G2[, 2]+resolution/2, G2[, 3], G2[, 4]-resolution/2, G2[, 4]+resolution/2, G2[, 5]), stringsAsFactors=F)
colnames(G1.bed) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")
colnames(G2.bed) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "count")
G1.bed[, 1] <- as.character(G1.bed[, 1])
G1.bed[, 2] <- as.integer(G1.bed[, 2])
G1.bed[, 3] <- as.integer(G1.bed[, 3])
G1.bed[, 4] <- as.character(G1.bed[, 4])
G1.bed[, 5] <- as.integer(G1.bed[, 5])
G1.bed[, 6] <- as.integer(G1.bed[, 6])
G1.bed[, 7] <- as.integer(G1.bed[, 7])
G2.bed[, 1] <- as.character(G2.bed[, 1])
G2.bed[, 2] <- as.integer(G2.bed[, 2])
G2.bed[, 3] <- as.integer(G2.bed[, 3])
G2.bed[, 4] <- as.character(G2.bed[, 4])
G2.bed[, 5] <- as.integer(G2.bed[, 5])
G2.bed[, 6] <- as.integer(G2.bed[, 6])
G2.bed[, 7] <- as.integer(G2.bed[, 7])

# Filter out trans interactions.
table(G1.bed$chr1 != G1.bed$chr2)
table(G2.bed$chr1 != G2.bed$chr2)
G1.bed <- G1.bed[G1.bed$chr1 == G1.bed$chr2, ]
G2.bed <- G2.bed[G2.bed$chr1 == G2.bed$chr2, ]

# Filter out interactions between the same bins.
G1.bed <- G1.bed[!((G1.bed$start1 == G1.bed$start2) & (G1.bed$end1 == G1.bed$end2)), ]
G2.bed <- G2.bed[!((G2.bed$start1 == G2.bed$start2) & (G2.bed$end1 == G2.bed$end2)), ]

# Check that all lhs intervals are smaller than rhs intervals.
table(G1.bed$end1 <= G1.bed$start2)
table(G2.bed$end1 <= G2.bed$start2)

# Get all unique interacting intervals and intersect them with interactions of interest.
all.intervals <- rbind(G1.bed[, 1:6], G2.bed[, 1:6])
unique.intervals <- unique(all.intervals)
unique.intervals <- cbind(unique.intervals[, 1:6], paste0("id", 1:length(unique.intervals[, 1])))
colnames(unique.intervals)[7] <- "id"
ids.2 <- as.character(unique.intervals[, 7])
ids.1 <- as.character(downsampled.interactions[, 11])
lhs.interactions.2 <- unique.intervals[, c(1:3, 7)]
lhs.interactions.1 <- downsampled.interactions[, c(1:3, 11)]
rhs.interactions.2 <- unique.intervals[, c(4:6, 7)]
rhs.interactions.1 <- downsampled.interactions[, c(5:7, 11)]
lhs.overlap.results <- bedTools.2in(bed1=lhs.interactions.1, bed2=lhs.interactions.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
rhs.overlap.results <- bedTools.2in(bed1=rhs.interactions.1, bed2=rhs.interactions.2, opt.string=paste0('-wa -wb -f ', overlap.threshold))
interacting.bin.ids <- c()
print(length(ids.1))
for (i in 1:length(ids.1)) {
  
  #
  if (i %% 10000 == 0)
    print(i)
  
  # Find interactions in other cell type that intersect both ends of current interaction.
  lhs.intersections <- as.character(lhs.overlap.results[lhs.overlap.results[, 4] == ids.1[i], 8])
  rhs.intersections <- as.character(rhs.overlap.results[rhs.overlap.results[, 4] == ids.1[i], 8])
  both.intersections <- lhs.intersections[lhs.intersections %in% rhs.intersections]
  if (length(both.intersections) > 0) {
    interacting.bin.ids <- c(interacting.bin.ids, both.intersections)
  }
  
}
interacting.bin.ids <- unique(interacting.bin.ids)
length(interacting.bin.ids)

# Populate counts for unique interacting intervals.
unique.intervals <- unique.intervals[unique.intervals$id %in% interacting.bin.ids, 1:6]
G1.contacts <- cbind(unique.intervals, rep(0, length(unique.intervals[, 1])))
colnames(G1.contacts)[7] <- "count"
G2.contacts <- cbind(unique.intervals, rep(0, length(unique.intervals[, 1])))
colnames(G2.contacts)[7] <- "count"
print(dim(unique.intervals))
for (i in 1:length(unique.intervals[, 1])) {
  
  #
  if (i %% 10000 == 0)
    print(i)
  
  #
  G1.value <- G1.bed[(G1.bed[, 1] == unique.intervals[i, 1]) & (G1.bed[, 2] == unique.intervals[i, 2]) & (G1.bed[, 3] == unique.intervals[i, 3]) &
                       (G1.bed[, 4] == unique.intervals[i, 4]) & (G1.bed[, 5] == unique.intervals[i, 5]) & (G1.bed[, 6] == unique.intervals[i, 6]), 7]
  if (length(G1.value) > 0)
    G1.contacts[i, 7] <- G1.value 
  G2.value <- G2.bed[(G2.bed[, 1] == unique.intervals[i, 1]) & (G2.bed[, 2] == unique.intervals[i, 2]) & (G2.bed[, 3] == unique.intervals[i, 3]) &
                       (G2.bed[, 4] == unique.intervals[i, 4]) & (G2.bed[, 5] == unique.intervals[i, 5]) & (G2.bed[, 6] == unique.intervals[i, 6]), 7]
  if (length(G2.value) > 0)
    G2.contacts[i, 7] <- G2.value
  
}


# Prepare plots -----------------------------------------------------------


# Quick overview of data.
table((G1.contacts$count > 0) | (G2.contacts$count > 0))
mean(G1.contacts$count + G2.contacts$count)
pdf(file=paste0("allelic.bias.total.counts.distribution.", name,".pdf"), height=8, width=8)
hist(G1.contacts$count + G2.contacts$count, xlim=c(0, 50), breaks=seq(0, 500, 5))
abline(v=mean(G1.contacts$count + G2.contacts$count), col="red")
dev.off()

# Filter interacting bins by minimum counts.
min.counts <- 15
keep <- ((G1.contacts$count + G2.contacts$count) >= min.counts)
G1.contacts.filtered <- G1.contacts[keep, ]
G2.contacts.filtered <- G2.contacts[keep, ]
plot(G1.contacts.filtered$count, G2.contacts.filtered$count)

#
filtered.bins <- G1.contacts.filtered[, 1:6]
filtered.G1 <- G1.contacts.filtered$count
filtered.G2 <- G2.contacts.filtered$count

# Prepare qq plot.
pvalues <- 0
for (j in 1:length(filtered.G1)) {
  
  pvalues[j] <- binom.test(filtered.G1[j], filtered.G1[j]+filtered.G2[j], 0.5, alternative="two.sided")$p.value
  
}
observed <- -log10(sort(pvalues, decreasing=F))
expected <- -log10(ppoints(length(filtered.G1)))
significant <- (pvalues < 0.01)
significant.percentage <- (length(which(significant))/length(significant))*100
table(significant)

# Plot qq plot.
sig.colors <- rep("black", length(expected))
sig.colors[1:length(significant.fdr[significant.fdr == TRUE])] <- "red"

pdf(file=paste0("allelic.bias.qq.fdr.plot.", name,".pdf"), height=6, width=6)
plot(expected, observed, xlim=c(0, 4), ylim=c(0, 12), main=paste0(cell.type, " interactions allelic bias qq plot"), pch=19,
     sub=paste0(signif(100*table(significant.fdr)[["TRUE"]]/length(significant.fdr), 5), "%, ", table(significant.fdr)[["TRUE"]], "/", length(significant.fdr), " interactions with >= ", min.counts, " counts"), col=sig.colors)
abline(c(0, 0), c(1, 1))
dev.off()

# Plot FDR-controlled qq plot.
pvalues.fdr <- p.adjust(pvalues, "fdr")
significant.fdr <- (pvalues.fdr < 0.05)
significant.fdr.percentage <- (length(which(significant.fdr))/length(significant.fdr))*100
observed <- -log10(sort(pvalues.fdr, decreasing=F))
sig.colors <- rep("black", length(expected))
sig.colors[1:length(significant.fdr[significant.fdr == TRUE])] <- "purple"
pdf(file=paste0("allelic.bias.fdr.plot.", name,".pdf"), height=8, width=8)
plot(expected, observed, xlim=c(0, 5), ylim=c(0, 12), main=paste0(cell.type, " interactions allelic bias fdr plot"), 
     sub=paste0(signif(significant.fdr.percentage, 5), "%, ", table(significant.fdr)[["TRUE"]], "/", length(significant), " interactions with >= ", min.counts, " counts"), col=sig.colors)
abline(h=-log10(0.05), lty=2, col="red")
abline(h=3, lty=2, col="purple")
#abline(c(0, 0), c(1, 1))
dev.off()

# Plot counts.
reverse <- (filtered.G1 > filtered.G2)
abs_higher <- filtered.G2
abs_higher[reverse] <- filtered.G1[reverse]
abs_lower <- filtered.G1
abs_lower[reverse] <- filtered.G2[reverse]
shades <- rep("blue", length(reverse))
shades[reverse] <- "purple"
types <- rep(1, length(reverse))
types[reverse] <- 3
pdf(file=paste0("allelic.bias.allelic.counts.", name,".pdf"), height=8, width=8)
plot(abs_lower+epsilon, abs_higher+epsilon, xlim=c(0, 100), ylim=c(0, 100), col=shades, main=paste0(cell.type, " interactions allelic bias allelic counts"), pch=types,
     xlab="lower frequency allele (# reads)", ylab="higher frequency allele (# reads)")
points(abs_lower[significant]+epsilon, abs_higher[significant]+epsilon, col="red", pch=1, cex=2)
dev.off()

# Write significant bin interactions to file.
significant <- significant
distances <- filtered.bins$start2[significant] - filtered.bins$end1[significant]
sig.bins <- cbind(filtered.bins[significant, ], paste0("id", 1:length(distances)), 
                  filtered.G1[significant], filtered.G2[significant], distances, pvalues[significant], pvalues.fdr[significant], 
                  rep(NA, length(which(significant))), rep(NA, length(which(significant))),
                  rep(NA, length(which(significant))), rep(NA, length(which(significant))))
colnames(sig.bins) <- c("chr1", "start1", "end1", "chr2", "start2", "end2", "id", "count1", "count2", "distance", "p-value", "q-value", 
                        "lhs.snps.count", "rhs.snps.count", "lhs.snps.info", "rhs.snps.info")

#
lhs.snp.overlaps <- bedTools.2in(bed1=sig.bins[, c(1:3, 7)], bed2=snps.bed, opt.string="-wb")
rhs.snp.overlaps <- bedTools.2in(bed1=sig.bins[, c(4:6, 7)], bed2=snps.bed, opt.string="-wb")
for (i in 1:length(sig.bins[, 1])) {
  sig.bins$lhs.snps.count[i] <- length(lhs.snp.overlaps[as.character(lhs.snp.overlaps[, 4]) == sig.bins$id[i], 4])
  sig.bins$rhs.snps.count[i] <- length(rhs.snp.overlaps[as.character(rhs.snp.overlaps[, 4]) == sig.bins$id[i], 4])
}

#
write.table(sig.bins, file=paste0("allelic.bias.significant.interactions.", name, ".tsv"), sep="\t", quote=F, row.names=F, col.names=F)

