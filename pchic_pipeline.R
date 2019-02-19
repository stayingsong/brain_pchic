#
# pchic_pipeline.R
# Author: Michael Song
# Last modified: 2019-01-15
# This script takes sets of PCHi-C interactions, ATAC-seq peaks, and RNA-seq data, performs integrative analysis, and generates tables/figures.
#


# Load dependencies.
require(ggplot2)
require(dplyr)
require(edgeR)
require(plotrix)
require(gtools)
require(gplots)
require(eulerr)
require(vioplot)
require(VennDiagram)


# Set up environment and load settings ------------------------------------


# Clear workspace before running script.
rm(list=ls())

# Turn off scientific notation for writing output.
options(scipen=999)

# Convert warnings into errors.
options(warn=1)

# Process CLI arguments (or override if running interactively).
cli.args = commandArgs(trailingOnly=TRUE)
print(cli.args[1]) # Set the home directory for the pipeline.
if ((exists("cli.args") == FALSE) | (length(cli.args) == 0)) {
  
  cli.args <- c()
  cli.args[1] <- "/Users/Michael/Box Sync/MS_general/analysis/brain_pchic_revision"

}

# Set home directory.
home.dir <- cli.args[1]
start.index <- cli.args[2]
end.index <- cli.args[3]
disease <- cli.args[4]

# Home directory should be structured as such:
# ~/scripts           - contains all scripts used in the analysis
# ~/scripts/config.R  - configuration file
# ~/pipeline          - contains all raw input files used in the analysis
# ~/figures           - default output folder
# ~/figures/saved     - contains saved intermediate files used in the analysis
# ~/resources         - contains supplemmentary/reference files used in the analysis


# Set up analysis environment ---------------------------------------------


# Source files should be in ~/scripts.
scripts.dir <- paste0(home.dir, "/scripts")
source(paste0(scripts.dir, "/interaction_analysis.R"))
source(paste0(scripts.dir, "/prepare_plots.R"))
source(paste0(scripts.dir, "/preprocess_data.R"))
source(paste0(scripts.dir, "/read_data.R"))
source(paste0(scripts.dir, "/utilities.R"))

# Load configuration file.
source(paste0(scripts.dir, "/config.R"))

# Create output directories.
output.dir <- paste0(home.dir, "/figures")
dir.list <- c("", 
              "/atac_seq_results",
              "/feature_results", "/feature_results/snps", "/feature_results/vista",
              "/interaction_results", "/interaction_results/general", "/interaction_results/hubs", "/interaction_results/other", "/interaction_results/specificity",
              "/rna_seq_results",
              "/saved", 
              "/visualization")
for (dir in dir.list)
  if (!file.exists(paste0(output.dir, "/", run.id, dir)))
    dir.create(paste0(output.dir, "/", run.id, dir))


# Read input data ---------------------------------------------------------


# Read in interactions.
interactions.raw <- read.interactions(interactions.list, paste0(home.dir, "/pipeline"))
save(interactions.raw, file=paste0(output.dir, "/", run.id, "/saved/interactions.raw.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.raw.Rdata"))

# Read in ATAC-seq peaks.
atac.seq.peaks <- read.atac.seq.data(atac.seq.peak.list, paste0(home.dir, "/pipeline"), allowed.chrs, max.atac.seq.peaks)
save(atac.seq.peaks, file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.Rdata"))

# Read in features.
features <- read.features(features.list, paste0(home.dir, "/pipeline"))
save(features, file=paste0(output.dir, "/", run.id, "/saved/features.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/features.Rdata"))

# Read in RNA-seq transcript quantification results.
transcript.lengths.individual <- read.rna.seq.data(rna.seq.output.list.replicates, paste0(home.dir, "/pipeline"), "length")
transcript.counts.individual <- read.rna.seq.data(rna.seq.output.list.replicates, paste0(home.dir, "/pipeline"), "expected_count")
transcript.tpm.individual <- read.rna.seq.data(rna.seq.output.list.replicates, paste0(home.dir, "/pipeline"), "TPM")
save(transcript.lengths.individual, file=paste0(output.dir, "/", run.id, "/saved/transcript.lengths.individual.Rdata"))
save(transcript.counts.individual, file=paste0(output.dir, "/", run.id, "/saved/transcript.counts.individual.Rdata"))
save(transcript.tpm.individual, file=paste0(output.dir, "/", run.id, "/saved/transcript.tpm.individual.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/transcript.lengths.individual.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/transcript.counts.individual.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/transcript.tpm.individual.Rdata"))

# Get cell types.
if (exists("interactions.raw")) {
  cell.types <- names(interactions.raw)
} else {
  cell.types <- names(atac.seq.peaks)
}

# Get lists of promoter intervals.
allowed.types <- c("protein_coding", "miRNA", "miscRNA", "rRNA", "snoRNA", "snRNA", "lincRNA")
promoters.all <- read.promoters(paste0(home.dir, "/resources/", annotation.file), allowed.chrs, tss.upstream, tss.downstream, NULL)
promoters.strict <- read.promoters(paste0(home.dir, "/resources/", annotation.file), allowed.chrs, tss.upstream, tss.downstream, allowed.types)
other.types <- unique(promoters.all$gene_type)[!(unique(promoters.all$gene_type) %in% allowed.types)]
promoters.other <- read.promoters(paste0(home.dir, "/resources/", annotation.file), allowed.chrs, tss.upstream, tss.downstream, other.types)

# Read in exonic and intronic intervals.
exons <- read.table(paste0(home.dir, "/resources/", exons.file), sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
colnames(exons) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")                      
introns <- read.table(paste0(home.dir, "/resources/", introns.file), sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
colnames(introns) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")     

# Read in RMAP and BAITMAP CHiCAGO design files.
rmap <- read.table(paste0(home.dir, "/resources/", rmap.file), sep="\t", stringsAsFactors=F, header=F)
colnames(rmap) <- c("chrom", "chromStart", "chromEnd", "rmap_ID")
if (all(grepl("chr", rmap$chrom) == FALSE)) {
  rmap$chrom <- paste0("chr", rmap$chrom)
}
baitmap <- read.table(paste0(home.dir, "/resources/", baitmap.file), sep="\t", stringsAsFactors=F, header=F)
colnames(baitmap) <- c("chrom", "chromStart", "chromEnd", "rmap_ID", "baitmap_ID")
if (all(grepl("chr", baitmap$chrom) == FALSE)) {
  baitmap$chrom <- paste0("chr", baitmap$chrom)
}

# Read in chromosome sizes.
chrom.sizes <- read.table(paste0(home.dir, "/resources/", chrom.sizes.file), sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer"))

# Read in TAD data.
tads <- read.table(file=paste0(home.dir, "/resources/", tads.file), sep="\t", header=F, stringsAsFactors=F)

# Read in positive vista elements.
vista.data <- read.table(file=paste0(home.dir, "/resources/", vista.file), sep="\t", header=F, stringsAsFactors=F)
colnames(vista.data) <- c("chrom", "chromStart", "chromEnd", "vista_ID", "annotation")
  

# Preprocess expression data ----------------------------------------------


# Calculate mean TPM or TMM-normalized FPKM expression values for each gene and cell type.
tpm.data.individual <- summarizeExpressionResults(transcript.tpm.individual, type="TPM")
rpkm.data.individual <- summarizeExpressionResults(transcript.counts.individual, transcript.lengths.individual, type="RPKM")

# Write results to file.
write.table(tpm.data.individual, file=paste0(output.dir, "/", run.id, "/rna_seq_results/tpm.data.individual.txt"), sep="\t", row.names=T, col.names=T, quote=F)
write.table(rpkm.data.individual, file=paste0(output.dir, "/", run.id, "/rna_seq_results/rpkm.data.individual.txt"), sep="\t", row.names=T, col.names=T, quote=F)


# Preprocess interactions -------------------------------------------------


# Preproces interactions.
interactions.raw.processed <- list()
for (cell.type in cell.types) {

  # First filter for cis interactions on the main chromosomes, then remove interactions with large bins and reorder interactions based on the positions of their lhs and rhs bins.
  num.all <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- removeTransInteractions(interactions.raw[[cell.type]])
  num.cis <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- filterInteractionsByChromosome(interactions.raw[[cell.type]], allowed.chrs)
  num.fil <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- removeLowResInteractions(interactions.raw[[cell.type]], max.bin.size)
  num.res <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- reorderInteractionBins(interactions.raw[[cell.type]])

  # Deduplicate interactions.
  interactions.raw.processed[[cell.type]] <- deduplicateInteractions(interactions.raw[[cell.type]], cell.type)
  num.dup <- length(interactions.raw.processed[[cell.type]][, 1])

  # Print out preprocessing summary.
  print(paste0("Preprocessing interactions for cell type: ", cell.type))
  print(paste0("Retaining ", num.cis, " cis interactions from ", num.all, " total interactions"))
  print(paste0("Retaining ", num.fil, " interactions on main chromosomes"))
  print(paste0("Retaining ", num.res, " interactions with bins not larger than: ", max.bin.size, " bp"))
  print(paste0("Retaining ", num.dup, " interactions with deduplicated coordinates"))

}

# Save results.
save(interactions.raw.processed, file=paste0(output.dir, "/", run.id, "/saved/interactions.raw.processed.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.raw.processed.Rdata"))

# Plot #s of interactions as a function of score cutoff.
plotNumInteractionsBySignificance(interactions.raw.processed, seq(0, 10, 0.5), xlim=c(1, 10), ylim=c(0, 1000000), NULL, paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotNumInteractionsBySignificance(interactions.raw.processed, seq(0, 10, 0.5), xlim=c(2, 10), ylim=c(0, 500000), score.cutoff, paste0(output.dir, "/", run.id, "/interaction_results/general"))

# Filter interactions by a score cutoff.
interactions.sig <- list()
if (!is.null(score.cutoff)) {

  print(paste0("Applying score cutoff of ", score.cutoff, ":"))
  for (cell.type in cell.types) {

    interactions.sig[[cell.type]] <- filterInteractionsBySignificance(interactions.raw.processed[[cell.type]], score.cutoff)
    num.raw <- length(interactions.raw.processed[[cell.type]][, 1])
    num.sig <- length(interactions.sig[[cell.type]][, 1])
    print(paste0("Retaining ", num.sig, " significant interactions from ", num.raw, " interactions for cell type: ", cell.type))

  }

}

# Save results.
save(interactions.sig, file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.Rdata"))

# Write interactions for display on the WashU Epigenome Browser.
for (cell.type in cell.types) {

  washu.old <- convertToWashUFormat(interactions.sig[[cell.type]], center=T)
  write.table(washu.old, file=paste0(output.dir, "/", run.id, "/visualization/", cell.type, ".cutoff.", score.cutoff, ".washU.txt"),
              sep="\t", row.names=F, col.names=F, quote=F)
  washu.new <- updateWashUFormat(washu.old)
  write.table(sortBed(washu.new), file=paste0(output.dir, "/", run.id, "/visualization/", cell.type, ".cutoff.", score.cutoff, ".longrange.bed"),
              sep="\t", row.names=F, col.names=F, quote=F)

}

# Plot CDFs of interaction distances for each cell type.
plotDistances(interactions.sig, type="hist", 10, 250, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotDistances(interactions.sig, type="cdf", 10, 500, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotDistances(interactions.sig, type="both", 10, 500, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/general"))


# Analyze, downsample and expand interactions -----------------------------


# Expand interactions to the minimum resolution.
interactions.sig.res <- list()
interactions.all.res <- list()
for (cell.type in cell.types) {

  interactions.sig.res[[cell.type]] <- expandInteractions(interactions.sig[[cell.type]], interaction.resolution, neighboring.fragments)
  interactions.all.res[[cell.type]] <- expandInteractions(interactions.raw.processed[[cell.type]], interaction.resolution, neighboring.fragments)

}

# Save results.
save(interactions.sig.res, file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.res.Rdata"))
save(interactions.low.res, file=paste0(output.dir, "/", run.id, "/saved/interactions.low.res.Rdata"))
save(interactions.all.res, file=paste0(output.dir, "/", run.id, "/saved/interactions.all.res.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.res.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.low.res.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.all.res.Rdata"))


# Expand and annotate ATAC-seq peaks --------------------------------------


# Expand ATAC-seq peaks to minimum peak width and annotate with promoters/interactions.
atac.seq.peaks.res <- list()
atac.seq.peaks.ann <- list()
for (cell.type in cell.types) {

  print(paste0("Expanding and annotating ATAC-seq peaks for cell type: ", cell.type))
  atac.seq.peaks.res[[cell.type]] <- expandFeatures(atac.seq.peaks[[cell.type]], atac.seq.peak.width)
  atac.seq.peaks.ann[[cell.type]] <- annotateATACSeqPeaks(atac.seq.peaks.res[[cell.type]], atac.seq.peaks[[cell.type]], interactions.sig.res[[cell.type]], promoters.strict, promoters.other, exons, introns)

}
  
# Save results.
save(atac.seq.peaks.res, file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.res.Rdata"))
save(atac.seq.peaks.ann, file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.ann.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.res.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.ann.Rdata"))

# Plot/print ATAC-seq peak annotation results.
# plotATACSeqPeakAnnotations(atac.seq.peaks.ann, paste0(output.dir, "/", run.id, "/atac_seq_results"))


# Annotate interactions with ATAC-seq peaks and features ------------------


# Annotate significant and nonsignificant interactions.
interactions.sig.ann <- list()
interactions.low.ann <- list()
for (cell.type in cell.types) {

  print(paste0("Annotating interactions for cell type: ", cell.type))
  interactions.sig.ann[[cell.type]] <- annotateInteractions(interactions.sig.res[[cell.type]], atac.seq.peaks.res[[cell.type]], atac.seq.peaks.ann[[cell.type]], features[[cell.type]], promoters.strict, promoters.other)
  interactions.low.ann[[cell.type]] <- annotateInteractions(interactions.low.res[[cell.type]], atac.seq.peaks.res[[cell.type]], atac.seq.peaks.ann[[cell.type]], features[[cell.type]], promoters.strict, promoters.other)

}

# Save results.
save(interactions.sig.ann, file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.ann.Rdata"))
save(interactions.low.ann, file=paste0(output.dir, "/", run.id, "/saved/interactions.low.ann.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.sig.ann.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.low.ann.Rdata"))

# Plot interaction annotation results.
plotInteractionAnnotations(interactions.sig.ann, paste0(output.dir, "/", run.id, "/interaction_results/general"))

# Downsample nonsignificant interactions to same level as significant interactions.
interactions.low.res.sampled <- list()
interactions.low.ann.sampled <- list()
for (cell.type in cell.types) {

  # Count # of significant interactions and sample random indices.
  num.significant <- length(interactions.sig.res[[cell.type]][, 1])
  num.nonsignificant <- length(interactions.low.res[[cell.type]][, 1])
  sampled.indices <- sample(1:num.nonsignificant, num.significant, replace=F)

  # Sample nonsignificant interactions and their corresponding annotations.
  interactions.low.res.sampled[[cell.type]] <- interactions.low.res[[cell.type]][sampled.indices, ]
  interactions.low.ann.sampled[[cell.type]] <- interactions.low.ann[[cell.type]]
  for (feature in names(interactions.low.ann[[cell.type]]))
    interactions.low.ann.sampled[[cell.type]][[feature]] <- interactions.low.ann.sampled[[cell.type]][[feature]][sampled.indices]

}

# Save results.
save(interactions.low.res.sampled, file=paste0(output.dir, "/", run.id, "/saved/interactions.low.res.sampled.Rdata"))
save(interactions.low.ann.sampled, file=paste0(output.dir, "/", run.id, "/saved/interactions.low.ann.sampled.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.low.res.sampled.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.low.ann.sampled.Rdata"))

# Perform gene promoter hub analysis.
hub.results.sig <- analyzeGenePromoterHubs(interactions.sig.res, interactions.sig.ann, promoters.strict, rpkm.data.individual)
hub.results.low <- analyzeGenePromoterHubs(interactions.low.res.sampled, interactions.low.ann.sampled, promoters.strict, rpkm.data.individual)

# Save results.
save(hub.results.sig, file=paste0(output.dir, "/", run.id, "/saved/hub.results.sig.Rdata"))
save(hub.results.low, file=paste0(output.dir, "/", run.id, "/saved/hub.results.low.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/hub.results.sig.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/hub.results.low.Rdata"))

# Filter hub results for genes of interest. Plot gene promoter hubs histograms of # of interactions per gene.
hub.results.sig.filtered <- filterHubResults(hub.results.sig, unique(promoters.strict$gene_id), -1, 1)
plotGenePromoterHubsHistogram(hub.results.sig.filtered, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/hubs"), "strict")
hub.results.sig.filtered <- filterHubResults(hub.results.sig, unique(promoters.strict$gene_id[promoters.strict$gene_type == "protein_coding"]), -1, 1)
plotGenePromoterHubsHistogram(hub.results.sig.filtered, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/hubs"), "protein_coding")

# Print # of interacting genes by cell type.
printInteractingGenes(hub.results.sig, promoters.strict, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/general"))


# Specificity analysis ----------------------------------------------------


# Fill in heatmap of interaction scores at all significant interactions across all cell types.
specificity.heatmap <- makeLociSpecificityHeatmap(interactions.sig.res, interactions.all.res, overlap.threshold)

# Stitch together components of specificity matrix (generated using parallelized version of makeLociSpecificityHeatmap() function).
load(file=paste0(output.dir, "/", run.id, "/saved/specificity.heatmap.(", cell.types[1], ").Rdata"))
specificity.matrix <- specificity.heatmap
for (cell.type in cell.types) {

  load(file=paste0(output.dir, "/", run.id, "/saved/specificity.heatmap.(", cell.type, ").Rdata"))
  specificity.matrix[, cell.type] <- specificity.heatmap[, cell.type]

}
head(specificity.matrix)
specificity.table <- specificity.matrix >= score.cutoff

# Make mapping of individual interaction IDs to cell type specificity.
specificity.entries <- c()
for (i in 1:length(specificity.table[, 1])) {

  interaction.ids <- unlist(strsplit(rownames(specificity.table)[i], split=","))
  for (interaction.id in interaction.ids)
    specificity.entries <- rbind(specificity.entries, c(interaction.id, paste0(cell.types[specificity.table[i, ]], collapse=",")))

}
specificity.entries[, 1] <- as.character(specificity.entries[, 1])
specificity.entries[, 2] <- as.character(specificity.entries[, 2])

# Save results.
save(specificity.matrix, file=paste0(output.dir, "/", run.id, "/saved/specificity.matrix.Rdata"))
save(specificity.table, file=paste0(output.dir, "/", run.id, "/saved/specificity.table.Rdata"))
save(specificity.entries, file=paste0(output.dir, "/", run.id, "/saved/specificity.entries.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/specificity.matrix.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/specificity.table.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/specificity.entries.Rdata"))

# Plot specificity Venn diagrams.
plotSpecificityVenn(specificity.table, cell.type.colors, paste0(output.dir, "/", run.id, "/interaction_results/specificity"))


# GO analysis -------------------------------------------------------------


# Set output prefix for writing input files for GO analysis.
output.prefix <- paste0(output.dir, "/", run.id, "/interaction_results/specificity")

# Find genes participating in interactions shared across all cell types.
shared.terms <- getSpecificInteractingGenes(interactions.sig.res, interactions.sig.ann, atac.seq.peaks.ann, specificity.entries, c("cortical,hippocampal,motor,astrocyte"), promoters.strict, rpkm.data.individual, 0.5, mode="both_atac")
write.table(shared.terms, file=paste0(output.prefix, "/go.shared.txt"), sep="\t", row.names=F, col.names=F, quote=F)

# Find genes participating in specific interactions for each cell type.
for (cell.type in cell.types) {
  
  specific.terms <- getSpecificInteractingGenes(interactions.sig.res, interactions.sig.ann, atac.seq.peaks.ann, specificity.entries, cell.type, promoters.strict, rpkm.data.individual, 0.5, mode="both_atac")
  specific.terms <- specific.terms[!(specific.terms %in% shared.terms)]
  write.table(specific.terms, file=paste0(output.prefix, "/go.", cell.type, ".txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
}


# Write interactions ------------------------------------------------------


# Write interactions and their annotations to file.
omit.list <- c("expanded_lhs", "expanded_lhs_ids", "expanded_rhs", "expanded_rhs_ids")
for (control.disease in c("AMD", "BH", "BMI", "RA", "T1D", "T2D")) {

  omit.list <- c(omit.list, paste0(control.disease, "_lhs"))
  omit.list <- c(omit.list, paste0(control.disease, "_lhs_ids"))
  omit.list <- c(omit.list, paste0(control.disease, "_rhs"))
  omit.list <- c(omit.list, paste0(control.disease, "_rhs_ids"))

}
merged.data.all <- writeInteractions(interactions.sig, interactions.sig.ann, features, promoters.all, specificity.entries, omit.list, paste0(output.dir, "/", run.id, "/interaction_results/general"))

# Save results.
save(merged.data.all, file=paste0(output.dir, "/", run.id, "/saved/merged.data.all.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/merged.data.all.Rdata"))


# Vista element analysis --------------------------------------------------


# Get list of unique annotations for Vista elements.
all.annotations <- paste0(vista.data$annotation, collapse=",")
all.annotations <- gsub("ganglion, cranial", "ganglion|cranial", all.annotations)
all.annotations <- unlist(strsplit(all.annotations, split=","))
for (i in 1:length(all.annotations))
  all.annotations[i] <- unlist(strsplit(all.annotations[i], split="\\["))[1]
unique.annotations <- unique(all.annotations)
cns.annotations <- c("neural tube", "hindbrain (rhombencephalon)", "cranial nerve", "midbrain (mesencephalon)", "forebrain", 
                     "mesenchyme derived from neural crest", "dorsal root ganglion", "trigeminal V (ganglion|cranial)")

# Analyze overlap between Vista elements, interactions, and gene promoters.
vista.results <- analyzeVistaOverlap(vista.data, interactions.sig.res, interactions.sig.ann, atac.seq.peaks.res, promoters.strict, promoters.all, rmap, interaction.resolution)

# Save results.
save(vista.results, file=paste0(output.dir, "/", run.id, "/saved/vista.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/vista.results.Rdata"))

# Print out results for Vista element analysis.
plotVistaResults(vista.data, vista.results, cns.annotations, cell.type.colors, paste0(output.dir, "/", run.id, "/feature_results/vista"))
printVistaConcordance(vista.results, cell.type.colors, paste0(output.dir, "/", run.id, "/feature_results/vista"))


# GWAS SNP analysis -------------------------------------------------------


# Perform overlap analysis for all SNPs.
all.snp.results <- list()
tag.snp.results <- list()
for (disease in diseases) {

  print(paste0("Analyzing SNPs for disease: ", disease))
  disease.data <- features[[cell.types[1]]][["snp"]][[disease]]
  all.snp.results[[disease]] <- analyzeAllSNPOverlap(disease.data, interactions.sig.res, interactions.sig.ann, atac.seq.peaks.res, atac.seq.peaks.ann,
                                                     promoters.strict, promoters.all, rmap, 1000)
  tag.snp.results[[disease]] <- analyzeTagSNPOverlap(disease.data, all.snp.results[[disease]], promoters.strict, promoters.all, rmap, 1000)

}

# Save results.
save(all.snp.results, file=paste0(output.dir, "/", run.id, "/saved/all.snp.results.Rdata"))
save(tag.snp.results, file=paste0(output.dir, "/", run.id, "/saved/tag.snp.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/all.snp.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/tag.snp.results.Rdata"))

# Print out results for SNP analysis.
for (disease in diseases) {

  print(paste0("Processing results for disease: ", disease))
  disease.data <- features[[cell.types[1]]][["snp"]][[disease]]
  plotSNPResults(disease.data, all.snp.results[[disease]], tag.snp.results[[disease]], cell.type.colors, disease, paste0(output.dir, "/", run.id, "/feature_results/snps"))
  printSNPConcordance(tag.snp.results[[disease]])
  getTargetGenes(all.snp.results[[disease]], promoters.all, rpkm.data.individual, 0, paste0(output.dir, "/", run.id, "/feature_results/snps"))

}

# Print out results combined across all diseases.
merged.results <- combineTagSNPResults(tag.snp.results, features, diseases)
printSNPConcordance(merged.results)


# Allelic bias annotation -------------------------------------------------


# Annotate and write interactions with allelic bias for each cell type.
allelic.interactions.sig.ann <- list()
for (cell.type in c("cortical", "motor")) {

  load(file=paste0(output.dir, "/", run.id, "/saved/allelic.interactions.sig.pval.", cell.type, ".Rdata"))
  print(length(allelic.interactions.sig[, 1]))
  temp <- allelic.interactions.sig
  allelic.interactions.sig <- list()
  allelic.interactions.sig[[cell.type]] <- temp
  allelic.interactions.sig.ann[[cell.type]] <- annotateInteractions(allelic.interactions.sig[[cell.type]], atac.seq.peaks.res[[cell.type]], atac.seq.peaks.ann[[cell.type]], features[[cell.type]], promoters.strict, promoters.other)
  writeInteractions(allelic.interactions.sig, allelic.interactions.sig.ann, features, promoters.all, NULL, omit.list, paste0(output.dir, "/", run.id, "/allelic_bias/results"))

}


# End ---------------------------------------------------------------------

