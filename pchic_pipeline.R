#
# pchic_pipeline.R
# Author: Michael Song
# Date created: 2018-01-29
# Last modified: 2018-05-10
# This script performs integrative analysis of PCHi-C interactions, ATAC-seq peaks, RNA-seq data, and other (epi)genomic features.
#


# Miscellaneous settings --------------------------------------------------


# Clear workspace.
rm(list=ls())

# Turn off scientific notation for output/printing.
options(scipen=999)


# Load packages -----------------------------------------------------------


# Install and load packages.
#source("https://bioconductor.org/biocLite.R")
#biocLite("plotrix")
require(plotrix)
#biocLite("gtools")
require(gtools)
#biocLite("eulerr")
require(eulerr)
#biocLite("ggplot2")
require(ggplot2)
#biocLite("DiffBind")
require(DiffBind)
#biocLite("eulerr")
require(eulerr)
#biocLite("edgeR")
require(edgeR)
#biocLite("preprocessCore")
require(preprocessCore)
#biocLite("statmod")
require(statmod)
#biocLite("vioplot")
require(vioplot)
#install.packages("UpSetR")
require(UpSetR)
#biocLite("topGO")
require(topGO)
#biocLite("VennDiagram")
require(VennDiagram)
#biocLite("gplots")
require(gplots)
#biocLite("RColorBrewer")
require(RColorBrewer)
#biocLite("GMD")
require(GMD)


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


# Script options ----------------------------------------------------------


tss.upstream <- 500            # Defines how far a promoter extends upstream from TSS for intersection with ATAC-seq peaks/interactions.
tss.downstream <- 500          # Defines how far a promoter extends downstream from TSS for intersection with ATAC-seq peaks/interactions.
atac.peak.width <- 500         # Defines minimum width of ATAC-seq peak.
score.cutoff <- 5              # Keep interactions at/above this score cutoff.
permissive.score.cutoff <- 3   # More permissive version of the score cutoff for specific analyses.
min.bin.size <- 50000          # All interacting bins/fragments exceeding this width will be filtered out.
downsampling.depth <- NULL     # Defines how many interactions to downsample each CHiCAGO interaction set to before further analysis.  Or set to NULL
interaction.resolution <- 5000 # Defines resolution of interactions (all restriction fragments within this window are considered as interacting).
neighboring.fragments <- FALSE # Defines whether or not to consider neighboring fragments as interacting, regardless of interaction resolution.
clustering.resolution <- 5000  # Defines resolution of interactions for clustering (not implemented yet).
overlap.threshold <- 10^-9     # Defines minimum overlap threshold for two interactions to cluster together.
tpm.cutoff <- 1                # Cutoff for determining whether or not a gene is significantly expressed.


# Define inputs -----------------------------------------------------------


# # Read in and print out CLI arguments.
# args = commandArgs(trailingOnly=TRUE)
# print(args[1]) # File containing significance-filtered interaction filepaths (CHiCAGO IBED format).
# print(args[2]) # File containing ATAC-seq peak filepaths (narrowPeak format).
# print(args[3]) # File containing genome-wide TSS list (chrom \t chromStart \t chromEnd \t strand \t gene ID \t gene type \t gene name).
# print(args[4]) # File containing RNA-seq output filepaths for each replicate (RSEM output format).
# print(args[5]) # Output filepath.
# print(args[6]) # File containing restriction fragment coordinates (chrom \t chromStart \t chromEnd \t rmap ID).
# print(args[7]) # File containing supplementary feature data (filepath \t cell type \t category).
# print(args[8]) # File containing unfiltered interaction filepaths (CHiCAGO IBED format).
# print(args[9]) # Name of subfolder for organizing different runs.
# print(args[10]) # Chromosome sizes file for reference genome used.
# print(args[11]) # BED file of baited restriction fragments.
# print(args[12]) # BED file of all introns.
# print(args[13]) # BED file of all exons.
# print(args[14]) # File containing RNA-seq output filepaths for each sample (RSEM output format).
# print(args[15]) # BED file of TADs.

# Override CLI arguments for running the pipeline interactively.
home.folder <- "/Users/michael/Box Sync/MS_analysis/brain_pchic"
args <- c()
args[1] <- paste0(home.folder, "/pipeline_files/interaction_input_list_50k_hicup_5.txt")
args[2] <- paste0(home.folder, "/pipeline_files/atac-seq_input_list_pooled_macs2_parameters_0.05.txt")
args[3] <- paste0(home.folder, "/supplementary_files/gencode.v19.annotation.bed")
args[4] <- paste0(home.folder, "/pipeline_files/rna-seq_input_list_individual.txt")
args[5] <- paste0(home.folder, "/figures")
args[6] <- paste0(home.folder, "/supplementary_files/hg19.rmap")
args[7] <- paste0(home.folder, "/pipeline_files/features_list.txt")
args[8] <- paste0(home.folder, "/pipeline_files/interaction_input_list_50k_hicup_0.txt")
args[9] <- "final"
args[10] <- paste0(home.folder, "/supplementary_files/hg19.chrom.sizes")
args[11] <- paste0(home.folder, "/supplementary_files/hg19.baitmap")
args[12] <- paste0(home.folder, "/supplementary_files/gencode_v19_introns.bed")
args[13] <- paste0(home.folder, "/supplementary_files/gencode_v19_exons.bed")
args[14] <- paste0(home.folder, "/pipeline_files/rna-seq_input_list_pooled.txt")
args[15] <- paste0(home.folder, "/supplementary_files/hESC.tad")

# Convert to more readable variable names.
interaction.list.sig <- args[1]
atac.seq.list <- args[2]
tss.annotation.file <- args[3]
rna.seq.list.individual <- args[4]
output.dir <- args[5]
rmap.file <- args[6]
features.list <- args[7]
interaction.list.all <- args[8]
run.id <- args[9]
chrom.sizes.file <- args[10]
baitmap.file <- args[11]
introns.file <- args[12]
exons.file <- args[13]
rna.seq.list.pooled <- args[14]
tads.file <- args[15]


# Read input data ---------------------------------------------------------


# Read in interactions, ATAC-seq data, RNA-seq data, and features.
interactions.raw <- read.interactions(interaction.list.all)
atac.seq.peaks <- read.atac.seq.data(atac.seq.list)
expression.counts.individual <- read.expression.counts.data(rna.seq.list.individual)
gene.lengths.individual <- read.expression.gene.lengths(rna.seq.list.individual)
expression.tpm.individual <- read.expression.tpm.data(rna.seq.list.individual)
expression.counts.pooled <- read.expression.counts.data(rna.seq.list.pooled)
gene.lengths.pooled <- read.expression.gene.lengths(rna.seq.list.pooled)
expression.tpm.pooled <- read.expression.tpm.data(rna.seq.list.pooled)
features <- read.features(features.list)

# Get list of cell types.
cell.types <- names(interactions.raw)

# Read in TSS coordinates.
tss.annotation <- read.table(tss.annotation.file, sep="\t", header=F, stringsAsFactors=F, 
                             colClasses=c("character", "integer", "integer", "character", "character", "character", "character"))
colnames(tss.annotation) <- c("chrom", "chromStart", "chromEnd", "strand", "gene_id", "gene_type", "gene_name")

# Expand TSS coordinates to get set of promoters.
tss.plus <- tss.annotation[tss.annotation[, 4] == "+", ]
tss.plus[, 2] <- tss.plus[, 2] - tss.upstream
tss.plus[tss.plus[, 2] < 1, 2] <- 1
tss.plus[, 3] <- tss.plus[, 3] + tss.downstream
tss.minus <- tss.annotation[tss.annotation[, 4] == "-", ]
tss.minus[, 2] <- tss.minus[, 2] - tss.downstream
tss.minus[tss.minus[, 2] < 1, 2] <- 1
tss.minus[, 3] <- tss.minus[, 3] + tss.upstream
promoters <- rbind(tss.plus, tss.minus)

# Filter only protein coding and noncoding RNA promoters.
promoters.all <- promoters
allowed.types <- c("protein_coding", "miRNA", "miscRNA", "rRNA", "snoRNA", "snRNA", "lincRNA")
promoters <- promoters[promoters$gene_type %in% allowed.types, ]
#dim(promoters)
#length(unique(promoters$gene_id))
coding.promoters <- promoters[promoters$gene_type == "protein_coding", ]
noncoding.promoters <- promoters[promoters$gene_type != "protein_coding", ]

# Read in BED files of introns and exons.
introns <- read.table(introns.file, sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
exons <- read.table(exons.file, sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
colnames(introns) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")     
colnames(exons) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")                      

# Read in file containing restriction fragment coordinates.
rmap <- read.table(rmap.file, sep="\t", stringsAsFactors=F, header=F)
rmap[, 1] <- paste0("chr", rmap[, 1])
colnames(rmap) <- c("chrom", "chromStart", "chromEnd", "rmap_ID")
baitmap <- read.table(baitmap.file, sep="\t", stringsAsFactors=F, header=F)
baitmap[, 1] <- paste0("chr", baitmap[, 1])
colnames(baitmap) <- c("chrom", "chromStart", "chromEnd", "rmap_ID", "baitmap_ID")

# Read in chromosome sizes file.
chrom.sizes <- read.table(chrom.sizes.file, sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer"))

# Set and create output directories.
setwd(output.dir)
dir.create(paste0(output.dir, "/", run.id))
dir.create(paste0(output.dir, "/", run.id, "/saved"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results/features"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results/general"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results/heatmap"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results/hubs"))
dir.create(paste0(output.dir, "/", run.id, "/interaction_results/specificity"))
dir.create(paste0(output.dir, "/", run.id, "/visualization"))
dir.create(paste0(output.dir, "/", run.id, "/atac-seq_results"))
dir.create(paste0(output.dir, "/", run.id, "/rna-seq_results"))
dir.create(paste0(output.dir, "/", run.id, "/SNPs"))
dir.create(paste0(output.dir, "/", run.id, "/vista"))

# Read in TADs information.
tads <- read.table(file=tads.file, sep="\t", header=F, stringsAsFactors=F)


# Preprocess expression data ----------------------------------------------


# Summarize TPM and TMM-normalized RPKMs for each cell type in data frames.
tpm.data.individual <- summarizeExpressionResults(expression.tpm.individual, type="TPM")
tpm.data.pooled <- summarizeExpressionResults(expression.tpm.pooled, type="TPM")
rpkm.data.individual <- summarizeExpressionResults(expression.counts.individual, gene.lengths=gene.lengths.individual, type="RPKM")
rpkm.data.pooled <- summarizeExpressionResults(expression.counts.pooled, gene.lengths=gene.lengths.pooled, type="RPKM")

# Analyze correlation between replicates for each expression metric.
for (i in 1:length(cell.types)) {
  
  print(cell.types[i])
  results <- cbind(tpm.data.individual[, i], tpm.data.pooled[, i], rpkm.data.individual[, i], rpkm.data.pooled[, i])
  #print(head(results))
  colnames(results) <- c("tpm.individual", "tpm.pooled", "rpkm.individual", "rpkm.pooled")
  print(cor(results, method="pearson", use="complete.obs"))
  
}

# Get differentially expressed genes for each cell type.
de.genes.all <- getAllDifferentiallyExpressedGenes(expression.counts.individual, 0.01, paste0(output.dir, "/", run.id, "/rna-seq_results"))
de.genes.specific <- getSpecificDifferentiallyExpressedGenes(expression.counts.individual, 0.01, paste0(output.dir, "/", run.id, "/rna-seq_results"))

# Checkpoint.
#save(de.genes.all, file=paste0(output.dir, "/", run.id, "/saved/de.genes.all.Rdata"))
#save(de.genes.specific, file=paste0(output.dir, "/", run.id, "/saved/de.genes.specific.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/de.genes.all.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/de.genes.specific.Rdata"))


# Preprocess interactions -------------------------------------------------


# Retain only cis interactions.
for (cell.type in cell.types) {
  
  num.all <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- removeTransInteractions(interactions.raw[[cell.type]])
  num.cis <- length(interactions.raw[[cell.type]][, 1])
  print(paste0("Retaining ", num.cis, " cis interactions from ", num.all, " (", signif((num.cis/num.all)*100, 3), "%) total interactions for cell type: ", cell.type))
  
}

# Filter out interactions for which either end exceeds the specified bin size.
for (cell.type in cell.types) {
  
  num.all <- length(interactions.raw[[cell.type]][, 1])
  interactions.raw[[cell.type]] <- removeLowResInteractions(interactions.raw[[cell.type]], min.bin.size)
  num.res <- length(interactions.raw[[cell.type]][, 1])
  print(paste0("Retaining ", num.res, " resolution-filtered interactions from ", num.all, " (", signif((num.res/num.all)*100, 3), "%) interactions for cell type: ", cell.type))
  
}

# Reverse order of bins if rhs coordinates are smaller than lhs coordinates.
for (cell.type in cell.types) {
  
  incorrect.start <- interactions.raw[[cell.type]]$bait_start > interactions.raw[[cell.type]]$otherEnd_start
  incorrect.end <- interactions.raw[[cell.type]]$bait_end > interactions.raw[[cell.type]]$otherEnd_end
  temp.start <- interactions.raw[[cell.type]][incorrect.start, 1:4]
  interactions.raw[[cell.type]][incorrect.start, 1:4] <- interactions.raw[[cell.type]][incorrect.start, 5:8]
  interactions.raw[[cell.type]][incorrect.start, 5:8] <- temp.start
  
}

# Plot # of interactions as a function of score.
plotNumInteractionsBySignificance(interactions.raw, seq(0, 100, 0.2), xlim=c(1, 10), ylim=c(0, 500000), NULL, paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotNumInteractionsBySignificance(interactions.raw, seq(0, 100, 0.2), xlim=c(2, 10), ylim=c(0, 200000), score.cutoff, paste0(output.dir, "/", run.id, "/interaction_results/general"))

# Plot bin width distributions.
plotBinWidths(interactions.raw, min.bin.size, paste0(output.dir, "/", run.id, "/interaction_results/general"), "raw")

# Filter interactions by score cutoff.
interactions.sig <- list()
interactions.low <- list()
if (!is.null(score.cutoff)) {
  
  print(paste0("Applying score cutoff of ", score.cutoff))
  for (cell.type in cell.types) {
    
    num.raw <- length(interactions.raw[[cell.type]][, 1])
    interactions.sig[[cell.type]] <- filterInteractionsBySignificance(interactions.raw[[cell.type]], score.cutoff)
    interactions.low[[cell.type]] <- getNonSignificantInteractions(interactions.raw[[cell.type]], 0, 1, 100000)
    num.sig <- length(interactions.sig[[cell.type]][, 1])
    print(paste0("Retaining ", num.sig, " significant interactions from ", num.raw, " (", signif((num.sig/num.raw)*100, 3), "%) interactions for cell type: ", cell.type))
    
  }
  
}

# Write WashU browser pairwise interaction format.
interactions.con <- interactions.sig
for (cell.type in cell.types) {
  
  print(length(interactions.sig[[cell.type]]$score))
  washu.sig <- convertToWashUFormat(interactions.sig, cell.type)
  write.table(washu.sig, file=paste0(output.dir, "/", run.id, "/visualization/", cell.type, ".full.5.washU.txt"), sep="\t", row.names=F, col.names=F, quote=F)
  
}

#
plotDistances(interactions.sig, type="hist", paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotDistances(interactions.sig, type="cdf", paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotDistances(interactions.sig, type="hybrid", paste0(output.dir, "/", run.id, "/interaction_results/general"))


# Analyze, downsample and expand interactions -----------------------------


# Downsample and expand interactions.
interactions.sam <- list()
interactions.res <- list()
interactions.res.low <- list()
for (cell.type in cell.types) {
  
  print(paste0("Downsampling and expanding interactions for cell type: ", cell.type))
  if (!is.null(downsampling.depth)) {
    interactions.sam[[cell.type]] <- downsampleInteractions(interactions.sig[[cell.type]], downsampling.depth)
  } else {
    interactions.sam[[cell.type]] <- interactions.sig[[cell.type]]
  }
  interactions.res[[cell.type]] <- expandInteractions(interactions.sam[[cell.type]], interaction.resolution, neighboring.fragments)
  interactions.res.low[[cell.type]] <- expandInteractions(interactions.low[[cell.type]], interaction.resolution, neighboring.fragments)
  
}

# Expand unfiltered interactions as well for use later on.
interactions.res.all <- list()
for (cell.type in cell.types) {
  
  interactions.res.all[[cell.type]] <- expandInteractions(interactions.raw[[cell.type]], interaction.resolution, neighboring.fragments=F)
  
}

# Plot width of interaction bins after expanding to minimum resolution.
plotBinWidths(interactions.res, min.bin.size, paste0(output.dir, "/", run.id, "/interaction_results/general"), "expanded")


# Expand and annotate ATAC-seq peaks --------------------------------------


# Expand ATAC-seq peaks to minimum peak width and annotate with promoters/interactions.
atac.seq.peaks.res <- list()
atac.seq.peaks.ann <- list()
for (cell.type in cell.types) {
  
  print(paste0("Expanding and annotating ATAC-seq peaks for cell type: ", cell.type))
  atac.seq.peaks.res[[cell.type]] <- expandATACPeaks(atac.seq.peaks[[cell.type]], atac.peak.width)
  atac.seq.peaks.ann[[cell.type]] <- annotateATACPeaks(atac.seq.peaks.res[[cell.type]], atac.seq.peaks[[cell.type]], interactions.res[[cell.type]], promoters, introns, exons)
  
}

# Checkpoint.
#save(atac.seq.peaks.ann, file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.ann.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/atac.seq.peaks.ann.Rdata"))

# Plot/print ATAC-seq peak annotation results.
plotATACPeakAnnotations(atac.seq.peaks.ann, paste0(output.dir, "/", run.id, "/atac-seq_results"))
distanceATACPeaksInteractionLoci(atac.seq.peaks.res, atac.seq.peaks.ann, interactions.raw, paste0(output.dir, "/", run.id, "/atac-seq_results"))
printATACPeakAnnotations(atac.seq.peaks.ann, paste0(output.dir, "/", run.id, "/atac-seq_results"))


# Annotate interactions with ATAC-seq peaks and features ------------------


# Annotate interactions with ATAC-seq peaks + all features.
interactions.ann <- list()
interactions.ann.low <- list()
for (cell.type in cell.types) {
  
  print(paste0("Annotating interactions for cell type: ", cell.type))
  interactions.ann[[cell.type]] <- annotateInteractions(interactions.res[[cell.type]], promoters, atac.seq.peaks.res[[cell.type]], atac.seq.peaks.ann[[cell.type]], features[[cell.type]])
  interactions.ann.low[[cell.type]] <- annotateInteractions(interactions.res.low[[cell.type]], promoters, atac.seq.peaks.res[[cell.type]], atac.seq.peaks.ann[[cell.type]], features[[cell.type]])
  
}

# Checkpoint.
#save(interactions.ann, file=paste0(output.dir, "/", run.id, "/saved/interactions.ann.Rdata"))
#save(interactions.ann.low, file=paste0(output.dir, "/", run.id, "/saved/interactions.ann.low.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.ann.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/interactions.ann.low.Rdata"))

# Plot interaction annotation results.
plotAllInteractionClasses(interactions.ann, paste0(output.dir, "/", run.id, "/interaction_results/general"))
plotPromoterOtherInteractionTranscriptTypes(interactions.ann, promoters, paste0(output.dir, "/", run.id, "/interaction_results/general"))

# Analyze # of interactions at each promoter.
hub.results <- analyzePromoterHubs(interactions.res, interactions.ann, promoters, rpkm.data.pooled)
hub.results.low <- analyzePromoterHubs(interactions.res.low, interactions.ann.low, promoters, rpkm.data.pooled)

# Checkpoint.
#save(hub.results, file=paste0(output.dir, "/", run.id, "/saved/hub.results.Rdata"))
#save(hub.results.low, file=paste0(output.dir, "/", run.id, "/saved/hub.results.low.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/hub.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/hub.results.low.Rdata"))

# Analyze general features of interactions at each promoter.
plotPromoterHubsGeneral(hub.results, unique(promoters$gene_id), 0.5, paste0(output.dir, "/", run.id, "/interaction_results/hubs"), "all.genes")

# Analyze gene expression grouped into # of interactions with enhancer states at each promoter for each cell type.
plotPromoterHubsSpecific(hub.results=hub.results, subset.gene.ids=unique(promoters$gene_id), min.expression=1, bin.size=2, min.points.per.bin=5, y.max=c(125, 125, 125), max.num.bins=7, 
                         category='num_interactions_PE', output.prefix=paste0(output.dir, "/", run.id, "/interaction_results/hubs"), output.suffix="all.genes.num_interactions_PE")

# Analyze gene expression grouped into # of interactions with enhancer states at each promoter for all cell types.
hub.results.merged <- list()
hub.results.merged[["combined"]] <- rbind(hub.results[["cortical"]], hub.results[["hippocampal"]], hub.results[["astrocyte"]])
plotPromoterHubsSpecific(hub.results=hub.results.merged, subset.gene.ids=unique(promoters$gene_id), min.expression=1, bin.size=2, min.points.per.bin=10, y.max=c(100, 100, 100), max.num.bins=8, 
                          category='num_interactions_PE', output.prefix=paste0(output.dir, "/", run.id, "/interaction_results/hubs"), output.suffix="all.genes.num_interactions_PE.combined")

# Analyzes mean expression of promoters interacting with enhancer versus repressive states.
plotPromoterHubsComparisons(hub.results=hub.results, hub.results.control=hub.results.low, subset.gene.ids=unique(promoters$gene_id), min.expression=1, filter.percent=0.05, y.max=c(160, 110, 160),
                            output.prefix=paste0(output.dir, "/", run.id, "/interaction_results/hubs"))

# Analyzes one-to-one versus one-to-many interactions.
analyzeInteractionHubs(interactions.res, paste0(output.dir, "/", run.id, "/interaction_results/hubs"))


# Write results to file ---------------------------------------------------


# Write all interaction overlap results to file.
writeFeatureOverlapResults(interactions.sig, interactions.ann, output.prefix=paste0(output.dir, "/", run.id, "/interaction_results/general"))


# Cell type specific overlap analysis -------------------------------------


# For each interaction, determine which cell type(s) it appears in.
specificity.matrix <- getInteractionSpecificity(interactions.res, overlap.threshold)
specificity.categories <- getSpecificityCategories(specificity.matrix)

# Checkpoint.
#save(specificity.categories, file=paste0(output.dir, "/", run.id, "/saved/specificity.categories.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/specificity.categories.Rdata"))

# Plot specificity Venn diagram.
plotSpecificityVenn(specificity.categories, interactions.res, paste0(output.dir, "/", run.id, "/interaction_results/specificity"))

# Write interactions for preparing the specificity heatmap.
writeHeatmapFiles(interactions.res, interactions.res.all, chunk.size=1000, output.folder=paste0(output.dir, "/", run.id, "/interaction_results/heatmap"))

# Draw the heatmap.
plotHeatmapResults(interactions.res, paste0(output.dir, "/", run.id, "/interaciton_results/heatmap/results"), paste0(output.dir, "/", run.id, "/interaciton_results/heatmap"))
setwd(output.dir)


# Motif analysis ----------------------------------------------------------


# Get list of cell type specific distal ATAC-seq peaks involved in interactions.
distal.ATAC.ids <- list()
genes.distal.ATAC.ids <- list()
for (cell.type in cell.types) {
  
  # Select interaction IDs with distal ATAC-seq peaks.
  rhs.DA.interactions.selector <- ((((interactions.ann[[cell.type]]$lhs_promoter > 0) | (interactions.ann[[cell.type]]$'lhs_promoter_ATAC-seq' > 0))
                                    & ((interactions.ann[[cell.type]]$rhs_promoter == 0) & (interactions.ann[[cell.type]]$'rhs_promoter_ATAC-seq' == 0)))
                                   & ((interactions.ann[[cell.type]]$'rhs_distal_ATAC-seq' > 0)))
  lhs.DA.interactions.selector <- ((((interactions.ann[[cell.type]]$rhs_promoter > 0) | (interactions.ann[[cell.type]]$'rhs_promoter_ATAC-seq' > 0))
                                    & ((interactions.ann[[cell.type]]$lhs_promoter == 0) & (interactions.ann[[cell.type]]$'lhs_promoter_ATAC-seq' == 0)))
                                   & ((interactions.ann[[cell.type]]$'lhs_distal_ATAC-seq' > 0)))
  rhs.ids <- as.character(interactions.res[[cell.type]]$ID[rhs.DA.interactions.selector])
  lhs.ids <- as.character(interactions.res[[cell.type]]$ID[lhs.DA.interactions.selector])
  DA.interactions.ids <- unique(c(rhs.ids, lhs.ids))
  
  #
  lhs.promoter.ids <- as.character(interactions.ann[[cell.type]]$'lhs_promoter_ATAC-seq_ids'[interactions.res[[cell.type]]$ID %in% DA.interactions.ids])
  rhs.promoter.ids <- as.character(interactions.ann[[cell.type]]$'rhs_promoter_ATAC-seq_ids'[interactions.res[[cell.type]]$ID %in% DA.interactions.ids])
  all.promoter.ids <- unique(c(lhs.promoter.ids, rhs.promoter.ids))
  split.promoter.ids <- c()
  for (i in 1:length(all.promoter.ids)) {
    split.promoter.ids <- c(split.promoter.ids, unlist(strsplit(all.promoter.ids[i], split=",")))
  }
  split.promoter.ids <- unique(split.promoter.ids)
  genes <- atac.seq.peaks.ann[[cell.type]][split.promoter.ids, "promoter_genes"]
  genes.split <- c()
  for (i in 1:length(genes)) {
    genes.split <- c(genes.split, unlist(strsplit(genes[i], split=",")))
  }
  split.promoter.ids <- genes.split
  
  #
  length(split.promoter.ids)
  expression.selector <- (rpkm.data.pooled[split.promoter.ids, cell.type] > 1)
  table(expression.selector)
  genes.distal.ATAC.ids[[cell.type]] <- split.promoter.ids[expression.selector]
  
  # Get ATAC-seq peak IDs from these interactions.
  rhs.ATAC.ids <- as.character(interactions.ann[[cell.type]]$'rhs_distal_ATAC-seq_ids'[rhs.DA.interactions.selector])
  lhs.ATAC.ids <- as.character(interactions.ann[[cell.type]]$'lhs_distal_ATAC-seq_ids'[lhs.DA.interactions.selector])
  distal.ATAC.ids[[cell.type]] <- unique(c(rhs.ATAC.ids, lhs.ATAC.ids))
  distal.ATAC.ids.split <- c()
  for (i in length(distal.ATAC.ids[[cell.type]])) {
    distal.ATAC.ids.split <- c(distal.ATAC.ids.split, unlist(strsplit(distal.ATAC.ids[[cell.type]], split=",")))
  }
  distal.ATAC.ids[[cell.type]] <- unique(distal.ATAC.ids.split)
  print(paste0(length(DA.interactions.ids), " interactions with distal ATAC-seq peaks for cell type: ", cell.type))
  
  #
  specific.ids <- getSpecificityCategoryInteractionIDs(specificity.categories, cell.type)
  print(paste0(length(DA.interactions.ids[DA.interactions.ids %in% specific.ids]), " of these interactions are specific to cell type: ", cell.type))
  rhs.specific.ids <- rhs.ids[rhs.ids %in% specific.ids]
  lhs.specific.ids <- lhs.ids[lhs.ids %in% specific.ids]
  distal.ATAC.seq.ids <- c(interactions.ann[[cell.type]]$'rhs_distal_ATAC-seq_ids'[interactions.res[[cell.type]]$ID %in% rhs.specific.ids],
                           interactions.ann[[cell.type]]$'lhs_distal_ATAC-seq_ids'[interactions.res[[cell.type]]$ID %in% lhs.specific.ids])
  distal.ATAC.seq.ids <- unique(unlist(strsplit(paste0(distal.ATAC.seq.ids, collapse=","), split=",")))
  print(paste0(length(distal.ATAC.seq.ids), " distal ATAC-seq peaks participating in cell type specific interactions for: ", cell.type))
  distal.ATAC.seq.peaks <- atac.seq.peaks[[cell.type]][atac.seq.peaks[[cell.type]]$ID %in% distal.ATAC.seq.ids, ]
  plus.strand <- distal.ATAC.seq.peaks[, c(1:6)] # For GREAT analysis
  negative.strand <- cbind(distal.ATAC.seq.peaks[, c(1:3, 11)], rep("-", length(distal.ATAC.seq.peaks[, 1])))
  write.table(plus.strand, file=paste0(run.id, "/interaction_results/specificity/distal_ATAC-seq_peaks.plus.", cell.type, ".bed"), quote=F, row.names=F, col.names=F, sep="\t")

}

# Draw the motif enrichment plot.
plotMotifs(paste0(output.dir, "/", run.id, "/interaction_results/specificity/motifs_data.csv"), paste0(output.dir, "/", run.id, "/interaction_results/specificity"))


# eQTL analysis -----------------------------------------------------------


#
setwd(paste0(output.dir, "/", run.id, "/eQTLs"))
eqtl.cell.type <- "cortical"
interaction.cell.type <- "cortical"
eQTL.files <- list()
eQTL.files[["blood"]] <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/Whole_Blood.v7.egenes.txt"
eQTL.files[["liver"]] <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/Liver.v7.egenes.txt"
eQTL.files[["cortical"]] <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/Brain_Cortex.v7.egenes.txt"
eQTL.files[["hippocampal"]] <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/Brain_Hippocampus.v7.egenes.txt"
eQTL.files[["lung"]] <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/Lung.v7.egenes.txt"
eQTL.data <- read.table(eQTL.files[[eqtl.cell.type]], sep="\t", stringsAsFactors=F, header=T)
eQTL.bed <- cbind(eQTL.data$chr, eQTL.data$pos-1, eQTL.data$pos, eQTL.data$rs_id_dbSNP147_GRCh37p13, eQTL.data$qval, eQTL.data$gene_id, eQTL.data$gene_name, eQTL.data$tss_distance)
print(paste0("Total eQTL entries: ", length(eQTL.bed[, 1])))
eQTL.bed <- eQTL.bed[eQTL.bed[, 5] < 0.05, ]
print(paste0("Significant eQTL entries: ", length(eQTL.bed[, 1])))
eQTL.bed[, 1] <- paste0("chr", eQTL.bed[, 1])
eQTL.bed[, 6] <- unlist(strsplit(eQTL.bed[, 6], split="\\."))[(1:length(eQTL.bed[, 6]))*2-1]
eQTL.bed <- data.frame(eQTL.bed, stringsAsFactors=F)
eQTL.bed[, 1] <- as.character(eQTL.bed[, 1])
eQTL.bed[, 2] <- as.numeric(eQTL.bed[, 2])
eQTL.bed[, 3] <- as.numeric(eQTL.bed[, 3])
eQTL.bed[, 4] <- as.character(eQTL.bed[, 4])
eQTL.bed[, 5] <- as.numeric(eQTL.bed[, 5])
eQTL.bed[, 6] <- as.character(eQTL.bed[, 6])
eQTL.bed[, 7] <- as.character(eQTL.bed[, 7])
eQTL.bed[, 8] <- as.numeric(eQTL.bed[, 8])
eQTL.bed <- cbind(eQTL.bed, paste0("eqtl", 1:length(eQTL.bed[, 1])))
colnames(eQTL.bed) <- c("chrom", "start", "end", "rsid", "qval", "gene_id", "gene_name", "tss_distance", "ID")
write.table(eQTL.bed, file=paste0("/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/", cell.type, ".eQTL.processed.bed"),
            quote=F, row.names=F, col.names=F, sep="\t")


# SNP analysis ------------------------------------------------------------


#
load(file=paste0(output.dir, "/", run.id, "/saved/reported.gene.table.all.Rdata"))

#
setwd(paste0(output.dir, "/", run.id, "/SNPs"))
diseases <- c("AD", "ADHD", "ASD", "BD", "EP", "MDD", "PD", "SCZ")
cell.types <- names(interactions.ann)
snp.resolution <- 5000
all.tag.snp.results <- list()
all.imputed.snp.results <- list()
for (d in 1:length(diseases)) {
  
  #
  current.disease <- diseases[d]
  print(current.disease)
  
  # Filter SNPs for current disease by overlapping with exons.
  snps.bed <- features[[cell.type[1]]][["SNP"]][[diseases[d]]]
  num.all.tag.snps <- length(unique(snps.bed[snps.bed$is_query_snp, 4]))
  num.all.imputed.snps <- length(unique(snps.bed$rsID))
  exon.overlap.results <- bedTools.2in(bed1=snps.bed, bed2=exons, opt.string="-c")
  snps.bed <- snps.bed[exon.overlap.results[, 10] == 0, ]
  num.noncoding.tag.snps <- length(unique(snps.bed[snps.bed$is_query_snp, 4]))
  num.noncoding.imputed.snps <- length(unique(snps.bed$rsID))
  print(paste0("Retained ", num.noncoding.tag.snps, " (", signif(num.noncoding.tag.snps/num.all.tag.snps*100, 4), 
               "%) noncoding tag SNPs from ", num.all.tag.snps, " total tag SNPs."))
  print(paste0("Retained ", num.noncoding.imputed.snps, " (", signif(num.noncoding.imputed.snps/num.all.imputed.snps*100, 4), 
               "%) noncoding imputed SNPs from ", num.all.imputed.snps, " total imputed SNPs."))
  
  # Create translation table between imputed SNPs and tag SNPs.
  translation <- c()
  imputed.snps <- as.character(unique(snps.bed$rsID))
  for (j in 1:length(imputed.snps)) {
    
    #
    query.snps <- features[[cell.types[1]]][["SNP"]][[diseases[d]]]$query_snp_rsid[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID == imputed.snps[j]]
    if (length(query.snps) != 1)
      print("Error")
    
    #
    exploded.query.splits <- unlist(strsplit(query.snps, split=","))
    for (k in 1:length(exploded.query.splits)) {
      translation <- rbind(translation, c(imputed.snps[j], exploded.query.splits[k]))
    }
    
  }
  colnames(translation) <- c("rsID", "query_snp_rsid")
  
  # Create results tables for tag and imputed SNPs.
  tag.snp.results <- data.frame(matrix(NA, num.noncoding.tag.snps, 7 + length(cell.types)*4), stringsAsFactors=F)
  rownames(tag.snp.results) <- as.character(unique(snps.bed[snps.bed$is_query_snp, 4]))
  colnames(tag.snp.results) <- c("reported_genes", "nearest_gene", "nearest_gene_dist", "nearest_expressed_gene", "nearest_expressed_gene_dist", "min_resolution_genes", 
                                 cell.types, paste0(cell.types, "_genes"), paste0(cell.types, "_imputed"), paste0(cell.types, "_imputed_genes"), "min_resolution_ambiguity")
  tag.snp.results[, 1] <- as.character(tag.snp.results[, 1])
  tag.snp.results[, 2] <- as.character(tag.snp.results[, 2])
  tag.snp.results[, 3] <- as.numeric(tag.snp.results[, 3])
  tag.snp.results[, 4] <- as.character(tag.snp.results[, 4])
  tag.snp.results[, 5] <- as.numeric(tag.snp.results[, 5])
  tag.snp.results[, 6] <- as.character(tag.snp.results[, 6])
  tag.snp.results[, 19] <- as.logical(tag.snp.results[, 19])
  for (i in 1:length(cell.types)) {
    tag.snp.results[, 6+i] <- as.logical(tag.snp.results[, 6+i])
    tag.snp.results[, 6+i+length(cell.types)*1] <- as.character(tag.snp.results[, 6+i]+length(cell.types)*1)
    tag.snp.results[, 6+i+length(cell.types)*2] <- as.logical(tag.snp.results[, 6+i]+length(cell.types)*2)
    tag.snp.results[, 6+i+length(cell.types)*3] <- as.character(tag.snp.results[, 6+i]+length(cell.types)*3)
  }
  imputed.snp.results <- data.frame(matrix(NA, num.noncoding.imputed.snps, 6 + length(cell.types)*6), stringsAsFactors=F)
  rownames(imputed.snp.results) <- as.character(unique(snps.bed$rsID))
  colnames(imputed.snp.results) <- c("nearest_gene", "nearest_gene_dist", "nearest_expressed_gene", "nearest_expressed_gene_dist", "min_resolution_genes", 
                                     cell.types, paste0(cell.types, "_genes"), paste0(cell.types, "_atac-seq_distal"), paste0(cell.types, "_atac-seq_distal_distance"),
                                     paste0(cell.types, "_atac-seq_promoter"), paste0(cell.types, "_atac-seq_promoter_distance"), "is_query_snp")
  imputed.snp.results[, 1] <- as.character(imputed.snp.results[, 1])
  imputed.snp.results[, 2] <- as.numeric(imputed.snp.results[, 2])
  imputed.snp.results[, 3] <- as.character(imputed.snp.results[, 3])
  imputed.snp.results[, 4] <- as.numeric(imputed.snp.results[, 4])
  imputed.snp.results[, 5] <- as.character(imputed.snp.results[, 5])
  for (i in 1:length(cell.types)) {
    imputed.snp.results[, 5+i] <- as.logical(imputed.snp.results[, 5+i])
    imputed.snp.results[, 5+i+length(cell.types)*1] <- as.character(imputed.snp.results[, 5+i]+length(cell.types)*1)
    imputed.snp.results[, 5+i+length(cell.types)*2] <- as.logical(imputed.snp.results[, 5+i+length(cell.types)*2])
    imputed.snp.results[, 5+i+length(cell.types)*3] <- as.numeric(imputed.snp.results[, 5+i+length(cell.types)*3])
    imputed.snp.results[, 5+i+length(cell.types)*4] <- as.logical(imputed.snp.results[, 5+i+length(cell.types)*4])
    imputed.snp.results[, 5+i+length(cell.types)*5] <- as.numeric(imputed.snp.results[, 5+i+length(cell.types)*5])
  }
  
  # Calculate nearest and same fragment genes for each SNP.
  snp.positions <- unique(snps.bed[, 1:4])
  snp.positions.sorted <- sortBed(snp.positions)
  all.promoters.sorted <- sortBed(promoters)
  max.expression.per.gene <- apply(rpkm.data.pooled[all.promoters.sorted$gene_id, ], 1, max)
  expressed.promoters.sorted <- all.promoters.sorted[max.expression.per.gene > 0, ]
  query.chr <- unique(snp.positions[, 1])
  all.promoters.sorted.filtered <- all.promoters.sorted[all.promoters.sorted[, 1] %in% query.chr, ]
  expressed.promoters.sorted.filtered <- expressed.promoters.sorted[expressed.promoters.sorted[, 1] %in% query.chr, ]
  closest.results.all <- bedTools.2in(functionstring="bedtools closest", bed1=snp.positions.sorted, bed2=all.promoters.sorted.filtered, opt.string="-t first -d")
  closest.results.expressed <- bedTools.2in(functionstring="bedtools closest", bed1=snp.positions.sorted, bed2=expressed.promoters.sorted.filtered, opt.string="-t first -d")
  
  # Calculate same fragment gene intersections.
  snp.positions.expanded <- snp.positions
  widths <- abs(snp.positions.expanded[, 2] - snp.positions.expanded[, 3])
  centers <- round((snp.positions.expanded[, 2] + snp.positions.expanded[, 3])/2)
  snp.positions.expanded[, 2] <- centers[widths < snp.resolution] - snp.resolution/2
  snp.positions.expanded[, 3] <- centers[widths < snp.resolution] + snp.resolution/2
  snp.positions.expanded[snp.positions.expanded[, 2] < 1, 2] <- 1
  promoters.with.id <- cbind(promoters, paste0("promoter", 1:length(promoters[, 1])))
  snp.min.resolution.intersections <- bedTools.2in(bed1=snp.positions.expanded, bed2=promoters.with.id, opt.string="-wb")
  
  # Get overlapped SNP rsIDs for each cell type for reference later.
  overlapping.snp.ids <- list()
  for (cell.type in cell.types) {
    
    #
    lhs.snp.ids <- unique(interactions.ann[[cell.type]][[paste0("lhs_", diseases[d], "_ids")]])
    expanded.lhs.snp.ids <- c()
    for (j in 1:length(lhs.snp.ids)) {
      expanded.lhs.snp.ids <- c(expanded.lhs.snp.ids, unlist(strsplit(lhs.snp.ids[j], split=",")))
    }
    lhs.rsids <- unique(as.character(features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$ID %in% expanded.lhs.snp.ids]))
    lhs.rsids <- lhs.rsids[lhs.rsids != ""]
    
    #
    rhs.snp.ids <- unique(interactions.ann[[cell.type]][[paste0("rhs_", diseases[d], "_ids")]])
    expanded.rhs.snp.ids <- c()
    for (j in 1:length(rhs.snp.ids)) {
      expanded.rhs.snp.ids <- c(expanded.rhs.snp.ids, unlist(strsplit(rhs.snp.ids[j], split=",")))
    }
    rhs.rsids <- unique(as.character(features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$ID %in% expanded.rhs.snp.ids]))
    rhs.rsids <- rhs.rsids[rhs.rsids != ""]
    
    #
    overlapping.snp.ids[[cell.type]] <- c(lhs.rsids, rhs.rsids)
    
  }
  
  # Also create a way to figure out which interactions intersect which rsIDs. Column 1 is the rsID and column 2 is the interaction ID.
  lhs.translation <- list()
  rhs.translation <- list()
  for (cell.type in cell.types) {
    
    #
    lhs.translation[[cell.type]] <- c(-1, -1)
    rhs.translation[[cell.type]] <- c(-1, -1)
    
    #
    lhs.snp.ann.ids <- interactions.ann[[cell.type]][[paste0("lhs_", diseases[d], "_ids")]]
    for (i in 1:length(lhs.snp.ann.ids)) {
      
      #
      current.interaction.id <- as.character(interactions.res[[cell.type]]$ID[i])
      
      #
      if (lhs.snp.ann.ids[i] != "") {
        
        #
        exploded.ids <- unlist(strsplit(lhs.snp.ann.ids[i], split=","))
        rsids <- unique(as.character(features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$ID %in% exploded.ids]))
        if (length(exploded.ids) != length(rsids))
          print("No")
        for (j in 1:length(rsids)) {
          lhs.translation[[cell.type]] <- rbind(lhs.translation[[cell.type]], c(rsids[j], current.interaction.id))
        }
        
      }
      
    }
    
    #
    rhs.snp.ann.ids <- interactions.ann[[cell.type]][[paste0("rhs_", diseases[d], "_ids")]]
    for (i in 1:length(rhs.snp.ann.ids)) {
      
      #
      current.interaction.id <- as.character(interactions.res[[cell.type]]$ID[i])
      
      #
      if (rhs.snp.ann.ids[i] != "") {
        
        #
        exploded.ids <- unlist(strsplit(rhs.snp.ann.ids[i], split=","))
        rsids <- unique(as.character(features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$ID %in% exploded.ids]))
        if (length(exploded.ids) != length(rsids))
          print("No")
        for (j in 1:length(rsids)) {
          rhs.translation[[cell.type]] <- rbind(rhs.translation[[cell.type]], c(rsids[j], current.interaction.id))
        }
        
      }
      
    }
    
    #
    lhs.translation[[cell.type]] <- lhs.translation[[cell.type]][-1, ]
    rhs.translation[[cell.type]] <- rhs.translation[[cell.type]][-1, ]
    
  }
  
  # Intersect SNPs with distal ATAC-seq peaks and get distance to nearest peaks.
  atac.intersections <- list()
  atac.intersections.distal <- list()
  atac.distances <- list()
  atac.distances.distal <- list()
  for (cell.type in cell.types) {
    
    #
    atac.distal.selector <- atac.seq.peaks.ann[[cell.type]]$annotation == "distal"
    distal.atac.peaks <- atac.seq.peaks[[cell.type]][atac.distal.selector, ]
    distal.atac.peaks.res <- atac.seq.peaks.res[[cell.type]][atac.distal.selector, ]
    
    #
    snp.positions.sorted <- sortBed(snp.positions)
    
    #
    atac.intersections[[cell.type]] <- bedTools.2in(bed1=snp.positions.sorted, bed2=distal.atac.peaks.res, opt.string="-c")
    table(atac.intersections[[cell.type]][, 5] > 0)
    
    #
    query.chr <- unique(snp.positions.sorted[, 1])
    atac.seq.peaks.filtered <- distal.atac.peaks[distal.atac.peaks[, 1] %in% query.chr, ]
    atac.seq.peaks.filtered <- sortBed(atac.seq.peaks.filtered)
    atac.distances[[cell.type]] <- bedTools.2in(functionstring="bedtools closest", bed1=snp.positions.sorted, bed2=atac.seq.peaks.filtered, opt.string="-t all -d")
    
    #
    atac.distances.distal <- atac.distances
    atac.intersections.distal <- atac.intersections
    
  }
  
  # Intersect SNPs with promoter ATAC-seq peaks and get distance to nearest peaks.
  atac.intersections <- list()
  atac.intersections.promoter <- list()
  atac.distances <- list()
  atac.distances.promoter <- list()
  for (cell.type in cell.types) {
    
    #
    atac.distal.selector <- atac.seq.peaks.ann[[cell.type]]$annotation == "promoter"
    distal.atac.peaks <- atac.seq.peaks[[cell.type]][atac.distal.selector, ]
    distal.atac.peaks.res <- atac.seq.peaks.res[[cell.type]][atac.distal.selector, ]
    
    #
    snp.positions.sorted <- sortBed(snp.positions)
    
    #
    atac.intersections[[cell.type]] <- bedTools.2in(bed1=snp.positions.sorted, bed2=distal.atac.peaks.res, opt.string="-c")
    table(atac.intersections[[cell.type]][, 5] > 0)
    
    #
    query.chr <- unique(snp.positions.sorted[, 1])
    atac.seq.peaks.filtered <- distal.atac.peaks[distal.atac.peaks[, 1] %in% query.chr, ]
    atac.seq.peaks.filtered <- sortBed(atac.seq.peaks.filtered)
    atac.distances[[cell.type]] <- bedTools.2in(functionstring="bedtools closest", bed1=snp.positions.sorted, bed2=atac.seq.peaks.filtered, opt.string="-t all -d")
    
    #
    atac.distances.promoter <- atac.distances
    atac.intersections.promoter <- atac.intersections
    
  }
  
  # Fill out result tables (imputed SNPs first).
  print(length(imputed.snp.results[, 1]))
  for (i in 1:length(imputed.snp.results[, 1])) {
    
    #
    if (i %% 1000 == 0)
      print(i)
    
    #
    current.rsid <- rownames(imputed.snp.results)[i]
    
    #
    imputed.snp.results[i, 1] <- paste0(unique(as.character(closest.results.all[as.character(closest.results.all[, 4]) == current.rsid, 11])), collapse=",")
    imputed.snp.results[i, 2] <- unique(as.numeric(closest.results.all[as.character(closest.results.all[, 4]) == current.rsid, 12]))
    imputed.snp.results[i, 3] <- paste0(unique(as.character(closest.results.expressed[as.character(closest.results.expressed[, 4]) == current.rsid, 11])), collapse=",")
    imputed.snp.results[i, 4] <- unique(as.numeric(closest.results.expressed[as.character(closest.results.expressed[, 4]) == current.rsid, 12]))
    imputed.snp.results[i, 5] <- paste0(unique(as.character(snp.min.resolution.intersections[as.character(snp.min.resolution.intersections[, 4]) == current.rsid, 11])), collapse=",")
    
    #
    if (snps.bed[snps.bed[, 4] == current.rsid, 8] == TRUE) {
      imputed.snp.results[i, "is_query_snp"] <- TRUE
    } else {
      imputed.snp.results[i, "is_query_snp"] <- FALSE
    }
    
    #
    for (cell.type in cell.types) {
      
      # If the SNP overlaps interactions in a cell type, retrieve any potentially interacting genes.
      if (current.rsid %in% overlapping.snp.ids[[cell.type]]) {
        
        #
        imputed.snp.results[i, cell.type] <- TRUE
        interacting.genes <- c()
        
        # If SNP is in lhs of interaction, then fetch all genes associated with rhs of interaction.
        if (current.rsid %in% lhs.translation[[cell.type]][, 1]) {
          
          interaction.ids <- lhs.translation[[cell.type]][lhs.translation[[cell.type]][, 1] == current.rsid, 2]
          selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
          rhs.genes <- as.character(interactions.ann[[cell.type]]$rhs_promoter_ids[selector])
          rhs.genes <- paste0(rhs.genes, collapse=",")
          rhs.genes <- unique(unlist(strsplit(rhs.genes, split=",")))
          rhs.coords <- interactions.res[[cell.type]][selector, 5:7]
          interacting.genes <- c(interacting.genes, rhs.genes)
          
        }
        
        # If SNP is in rhs of interaction, then fetch all genes associated with lhs of interaction.
        if (current.rsid %in% rhs.translation[[cell.type]][, 1]) {
          
          #
          interaction.ids <- rhs.translation[[cell.type]][rhs.translation[[cell.type]][, 1] == current.rsid, 2]
          selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
          lhs.genes <- as.character(interactions.ann[[cell.type]]$lhs_promoter_ids[selector])
          lhs.genes <- paste0(lhs.genes, collapse=",")
          lhs.genes <- unique(unlist(strsplit(lhs.genes, split=",")))
          lhs.coords <- interactions.res[[cell.type]][selector, 1:3]
          interacting.genes <- c(interacting.genes, lhs.genes)
          
        }
        
        #
        interacting.genes <- promoters$gene_name[promoters$gene_id %in% interacting.genes]
        interacting.genes <- unique(interacting.genes)
        interacting.genes <- interacting.genes[interacting.genes != ""]
        imputed.snp.results[i, paste0(cell.type, "_genes")] <- paste0(interacting.genes, collapse=",")
        
      } else {
        
        #
        imputed.snp.results[i, cell.type] <- FALSE
        imputed.snp.results[i, paste0(cell.type, "_genes")] <- ""
        
      }
      
      # Check out distal ATAC-seq peak overlaps.
      if (atac.intersections.distal[[cell.type]][atac.intersections.distal[[cell.type]][, 4] == current.rsid, 5] > 0) {
        
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_distal")] <- TRUE
        
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_distal_distance")] <- atac.distances.distal[[cell.type]][atac.distances.distal[[cell.type]][, 4] == current.rsid, 16]
        
      } else {
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_distal")] <- FALSE
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_distal_distance")] <- NA
        
      }
      
      # Check out promoter ATAC-seq peak overlaps.
      if (atac.intersections.promoter[[cell.type]][atac.intersections.promoter[[cell.type]][, 4] == current.rsid, 5] > 0) {
        
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_promoter")] <- TRUE
        
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_promoter_distance")] <- atac.distances.promoter[[cell.type]][atac.distances.promoter[[cell.type]][, 4] == current.rsid, 16]
        
      } else {
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_promoter")] <- FALSE
        imputed.snp.results[i, paste0(cell.type, "_atac-seq_promoter_distance")] <- NA
        
      }
      
    }
    
  }
  
  # Fill out result tables (now tag SNPs).
  print(length(tag.snp.results[, 1]))
  for (i in 1:length(tag.snp.results[, 1])) {
    
    #
    if (i %% 50 == 0)
      print(i)
    
    #
    current.rsid <- rownames(tag.snp.results)[i]
    
    #
    reported.genes <- reported.gene.table.all[[current.disease]][reported.gene.table.all[[current.disease]][, 1] == current.rsid, 2]
    reported.genes <- unique(reported.genes)
    reported.genes <- reported.genes[reported.genes != "NR"]
    if (length(reported.genes) == 0)
      reported.genes <- "N/A"
    tag.snp.results[i, 1] <- paste0(reported.genes, collapse=",")
    tag.snp.results[i, 2] <- paste0(unique(as.character(closest.results.all[as.character(closest.results.all[, 4]) == current.rsid, 11])), collapse=",")
    tag.snp.results[i, 3] <- unique(as.numeric(closest.results.all[as.character(closest.results.all[, 4]) == current.rsid, 12]))
    tag.snp.results[i, 4] <- paste0(unique(as.character(closest.results.expressed[as.character(closest.results.expressed[, 4]) == current.rsid, 11])), collapse=",")
    tag.snp.results[i, 5] <- unique(as.numeric(closest.results.expressed[as.character(closest.results.expressed[, 4]) == current.rsid, 12]))
    tag.snp.results[i, 6] <- paste0(unique(as.character(snp.min.resolution.intersections[as.character(snp.min.resolution.intersections[, 4]) == current.rsid, 11])), collapse=",")
    
    #
    related.imputed.snps <- translation[translation[, 2] == current.rsid, 1]
    related.results <- imputed.snp.results[related.imputed.snps, ]
    if ((current.rsid %in% related.imputed.snps) == FALSE)
      print("Error")
    
    # 
    nearest.gene.coordinates <- promoters[promoters$gene_name == tag.snp.results[i, 2], 1:3]
    gene.rmap.expanded <- bedTools.2in(bed1=nearest.gene.coordinates, bed2=rmap, opt.string="-wb")
    imputed.snp.coordinates <- snp.positions[snp.positions[, 4] %in% related.imputed.snps, ]
    snp.rmap.expanded <- bedTools.2in.compact(bed1=imputed.snp.coordinates, bed2=rmap, opt.string="-wb")
    gene.rmap.ids <- unique(gene.rmap.expanded[, 7])
    snp.rmap.ids <- unique(snp.rmap.expanded[, 8])
    if ((length(gene.rmap.ids) == 1) & (length(snp.rmap.ids) == 1) & (gene.rmap.ids[1] == snp.rmap.ids[1])) {
      tag.snp.results[i, 19] <- TRUE
    } else {
      tag.snp.results[i, 19] <- FALSE
    }
    
    #
    for (cell.type in cell.types) {
      
      # If the SNP overlaps interactions in a cell type, retrieve any potentially interacting genes.
      if (current.rsid %in% overlapping.snp.ids[[cell.type]]) {
        
        #
        tag.snp.results[i, cell.type] <- TRUE
        interacting.genes <- c()
        
        # If SNP is in lhs of interaction, then fetch all genes associated with rhs of interaction.
        if (current.rsid %in% lhs.translation[[cell.type]][, 1]) {
          
          interaction.ids <- lhs.translation[[cell.type]][lhs.translation[[cell.type]][, 1] == current.rsid, 2]
          selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
          rhs.genes <- as.character(interactions.ann[[cell.type]]$rhs_promoter_ids[selector])
          rhs.genes <- paste0(rhs.genes, collapse=",")
          rhs.genes <- unique(unlist(strsplit(rhs.genes, split=",")))
          interacting.genes <- c(interacting.genes, rhs.genes)
          
        }
        
        # If SNP is in rhs of interaction, then fetch all genes associated with lhs of interaction.
        if (current.rsid %in% rhs.translation[[cell.type]][, 1]) {
          
          #
          interaction.ids <- rhs.translation[[cell.type]][rhs.translation[[cell.type]][, 1] == current.rsid, 2]
          selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
          lhs.genes <- as.character(interactions.ann[[cell.type]]$lhs_promoter_ids[selector])
          lhs.genes <- paste0(lhs.genes, collapse=",")
          lhs.genes <- unique(unlist(strsplit(lhs.genes, split=",")))
          interacting.genes <- c(interacting.genes, lhs.genes)
          
        }
        
        #
        interacting.genes <- promoters$gene_name[promoters$gene_id %in% interacting.genes]
        interacting.genes <- unique(interacting.genes)
        interacting.genes <- interacting.genes[interacting.genes != ""]
        tag.snp.results[i, paste0(cell.type, "_genes")] <- paste0(interacting.genes, collapse=",")
        
      } else {
        
        #
        tag.snp.results[i, cell.type] <- FALSE
        tag.snp.results[i, paste0(cell.type, "_genes")] <- ""
        
      }
      
      # Process results for imputed SNPs related to this tag SNP.
      tag.snp.results[i, paste0(cell.type, "_imputed")] <- any(related.results[, cell.type])
      related.interacting.genes <- unique(unlist(strsplit(paste0(related.results[, paste0(cell.type, "_genes")], collapse=","), split=",")))
      related.interacting.genes <- related.interacting.genes[related.interacting.genes != ""]
      tag.snp.results[i, paste0(cell.type, "_imputed_genes")] <- paste0(related.interacting.genes, collapse=",")
      
    }
    
  }
  
  #
  all.tag.snp.results[[current.disease]] <- tag.snp.results
  all.imputed.snp.results[[current.disease]] <- imputed.snp.results
  
}

# Checkpoint.
save(all.tag.snp.results, file=paste0(output.dir, "/", run.id, "/saved/all.tag.snp.results.Rdata"))
save(all.imputed.snp.results, file=paste0(output.dir, "/", run.id, "/saved/all.imputed.snp.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/all.tag.snp.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/all.imputed.snp.results.Rdata"))

# Write results to file.
for (d in 1:length(diseases)) {
  
  #
  current.disease <- diseases[d]
  print(current.disease)
  
  # Filter SNPs for current disease by overlapping with exons.
  snps.bed <- features[[cell.type[1]]][["SNP"]][[diseases[d]]]
  num.all.tag.snps <- length(unique(snps.bed[snps.bed$is_query_snp, 4]))
  num.all.imputed.snps <- length(unique(snps.bed$rsID))
  exon.overlap.results <- bedTools.2in(bed1=snps.bed, bed2=exons, opt.string="-c")
  snps.bed <- snps.bed[exon.overlap.results[, 10] == 0, ]
  
  # Create translation table between imputed SNPs and tag SNPs.
  translation <- c()
  imputed.snps <- as.character(unique(snps.bed$rsID))
  for (j in 1:length(imputed.snps)) {
    
    query.snps <- features[[cell.types[1]]][["SNP"]][[diseases[d]]]$query_snp_rsid[features[[cell.types[1]]][["SNP"]][[diseases[d]]]$rsID == imputed.snps[j]]
    if (length(query.snps) != 1)
      print("What")
    
    exploded.query.splits <- unlist(strsplit(query.snps, split=","))
    for (k in 1:length(exploded.query.splits)) {
      translation <- rbind(translation, c(imputed.snps[j], exploded.query.splits[k]))
    }
    
  }
  colnames(translation) <- c("rsID", "query_snp_rsid")
  
  #
  print(table(rownames(all.imputed.snp.results[[current.disease]]) == snps.bed$rsID))
  all.imputed.snp.results[[current.disease]] <- cbind(snps.bed[, 1:6], all.imputed.snp.results[[current.disease]])
  
  #
  snps.bed.filtered <- snps.bed[snps.bed$is_query_snp == TRUE, ]
  print(table(rownames(all.tag.snp.results[[current.disease]]) == snps.bed.filtered$rsID))
  all.tag.snp.results[[current.disease]] <- cbind(snps.bed.filtered[, 1:6], all.tag.snp.results[[current.disease]])
  all.tag.snp.results[[current.disease]] <- cbind(all.tag.snp.results[[current.disease]], rep("", length(all.tag.snp.results[[current.disease]][, 1])))
  colnames(all.tag.snp.results[[current.disease]])[26] <- "imputed_SNPs"
  all.tag.snp.results[[current.disease]][, 26] <- as.character(all.tag.snp.results[[current.disease]][, 26])
  for (i in 1:length(all.tag.snp.results[[current.disease]][, 1])) {
    
    all.tag.snp.results[[current.disease]][i, 26] <- as.character(paste0(as.character(translation[as.character(translation[, 2]) == as.character(all.tag.snp.results[[current.disease]]$rsID[i]), 1]), collapse=","))
    
  }
  
  #
  atac.supported.entries <- c()
  master.selector <- rep(FALSE, length(all.imputed.snp.results[[current.disease]][, 1]))
  for (cell.type in cell.types) {
    selector <- (all.imputed.snp.results[[current.disease]][, cell.type] == TRUE) & 
      (all.imputed.snp.results[[current.disease]][, paste0(cell.type, "_genes")] != "") &
      ((all.imputed.snp.results[[current.disease]][, paste0(cell.type, "_atac-seq_promoter")] == TRUE) | (all.imputed.snp.results[[current.disease]][, paste0(cell.type, "_atac-seq_distal")] == TRUE))
    master.selector <- master.selector | selector
  }
  table(master.selector)
  atac.supported.entries <- all.imputed.snp.results[[current.disease]][master.selector, ]
  atac.supported.entries <- atac.supported.entries[, c(1:8, 30, 12, 15, 24, 18, 13, 16, 25, 19, 14, 17, 26, 20)]
  write.table(atac.supported.entries, file=paste0("tables/", current.disease, ".atac-seq.supported.interactions.txt"), sep="\t", quote=F, col.names=T, row.names=F)
  
  all.tag.snp.results[[current.disease]] <- all.tag.snp.results[[current.disease]][, c(1:9, 25:26, 13, 16, 19, 22, 14, 17, 20, 23, 15, 18, 21, 24)]
  all.imputed.snp.results[[current.disease]] <- all.imputed.snp.results[[current.disease]][, c(1:8, 30, 12, 15, 13, 16, 14, 17)] 
  
  #
  write.table(all.tag.snp.results[[current.disease]], file=paste0("tables/", current.disease, ".all.tag.snp.results.txt"), sep="\t", quote=F, col.names=T, row.names=F)
  write.table(all.imputed.snp.results[[current.disease]], file=paste0("tables/", current.disease, ".all.imputed.snp.results.txt"), sep="\t", quote=F, col.names=T, row.names=F)
  
  
}

# Print numbers of  interacting SNPs.
for (d in 1:length(diseases)) {
  
  #
  print(diseases[d])
  print(length(all.tag.snp.results[[diseases[d]]][, 1]))
  print(length(all.imputed.snp.results[[diseases[d]]][, 1]))
  
  #
  count.single <- 0
  count.multiple <- 0
  count.cortical <- 0
  count.hippocampal <- 0
  count.astrocyte <- 0
  for (i in 1:length(all.tag.snp.results[[diseases[d]]][, 1])) {
    num.cell.types.overlap <- length(which(as.logical(all.tag.snp.results[[diseases[d]]][i, 13:15])))
    if (num.cell.types.overlap == 1) {
      count.single <- count.single + 1
      if (all.tag.snp.results[[diseases[d]]][i, 13] == TRUE)
        count.cortical <- count.cortical + 1
      if (all.tag.snp.results[[diseases[d]]][i, 14] == TRUE)
        count.hippocampal <- count.hippocampal + 1
      if (all.tag.snp.results[[diseases[d]]][i, 15] == TRUE)
        count.astrocyte <- count.astrocyte + 1
    } else if (num.cell.types.overlap > 1) {
      count.multiple <- count.multiple + 1
    }
  }
  
  #
  print(length(all.tag.snp.results[[diseases[d]]][, 1]) - count.single - count.multiple)
  print(count.single)
  print(count.multiple)
  print(count.cortical)
  print(count.hippocampal)
  print(count.astrocyte)
  
}

# Plot correspondence plots.
union.total <- 0
gwas.nearest.total <- 0
new.targets.total <- 0
for (d in 1:length(diseases)) {
  
  #
  print(diseases[d])
  
  #
  interacting <- rep(FALSE, length(all.tag.snp.results[[diseases[d]]][, 1]))
  nearest.agreement <- rep(FALSE, length(all.tag.snp.results[[diseases[d]]][, 1]))
  nearest.expressed.agreement <- rep(FALSE, length(all.tag.snp.results[[diseases[d]]][, 1]))
  reported.agreement <- rep(FALSE, length(all.tag.snp.results[[diseases[d]]][, 1]))
  num.targets <- rep(NA, length(all.tag.snp.results[[diseases[d]]][, 1]))
  targeted.genes <- c()
  min.resolution.ambiguity <- rep(FALSE, length(all.tag.snp.results[[diseases[d]]][, 1]))
  for (i in 1:length(all.tag.snp.results[[diseases[d]]][, 1])) {
    
    #
    tag.interacting.genes <- unique(unlist(strsplit(paste0(all.tag.snp.results[[diseases[d]]][i, paste0(cell.types, "_genes")], collapse=","), split=",")))
    imputed.interacting.genes <- unique(unlist(strsplit(paste0(all.tag.snp.results[[diseases[d]]][i, paste0(cell.types, "_imputed_genes")], collapse=","), split=",")))
    interacting.genes <- imputed.interacting.genes
    interacting.genes <- interacting.genes[interacting.genes != ""]
    num.targets[i] <- length(interacting.genes)
    targeted.genes <- c(targeted.genes, interacting.genes)
    
    #
    min.resolution.ambiguity[i] <- all.tag.snp.results[[diseases[d]]][i, 19]
    
    #
    if (all.tag.snp.results[[diseases[d]]][i, 2] %in% interacting.genes)
      nearest.agreement[i] <- TRUE
    if (all.tag.snp.results[[diseases[d]]][i, 4] %in% interacting.genes)
      nearest.expressed.agreement[i] <- TRUE
    if (length(interacting.genes) > 0)
      interacting[i] <- TRUE
    if (all.tag.snp.results[[diseases[d]]][i, 1] %in% interacting.genes)
      reported.agreement[i] <- TRUE
    
  }
  
  #
  union <- length(which((interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == TRUE)))
  gwas.nearest <- length(which((interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == FALSE)))
  new.targets <- sum(num.targets[(interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == TRUE)] - 1) +
    sum(num.targets[(interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == FALSE)])
  
  # Draw Venn diagram for current disease.
  pdf(paste0(output.dir, "/", run.id, "/SNPs/agreement_venn.", diseases[d], ".pdf"), width=5, height=5)
  values <- c(gwas.nearest, union, new.targets)
  names(values) <- c("Nearest", "Nearest&New interaction targets", "New interaction targets")
  fit <- euler(values)
  colors <- c("#FFE4C4", "#00FF7F")
  plot(fit, fill_opacity = 0.01, shape="circle", title=paste0("Agreement for ", diseases[d]),
       fill=colors, border="transparent", fontsize=12, quantities=list(fontsize=12))
  dev.off()
  
  #
  union.total <- union.total + union
  gwas.nearest.total <- gwas.nearest.total + gwas.nearest
  new.targets.total <- new.targets.total + new.targets
  
}

# Draw Venn diagram for all diseases.
pdf(paste0(output.dir, "/", run.id, "/SNPs/agreement_venn.all.pdf"), width=5, height=5)
values <- c(gwas.nearest.total, union.total, new.targets.total)
names(values) <- c("Nearest", "Nearest&New interaction targets", "New interaction targets")
fit <- euler(values)
colors <- c("#FFE4C4", "#00FF7F")
plot(fit, fill_opacity = 0.01, shape="circle", title=paste0("Agreement for all diseases"),
     fill=colors, border="transparent", fontsize=12, quantities=list(fontsize=12))
dev.off()

#
plotDiseaseEnrichment(paste0(output.dir, "/", run.id, "/SNPs/disease_cell_type_enrichment.csv"), paste0(output.dir, "/", run.id, "/SNPs"))


# Vista analysis ----------------------------------------------------------


#
vista.all.data <- read.table("/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/vista_processing_new/vista_all.csv", 
                             sep=",", stringsAsFactors=F, header=F, colClasses=rep("character", 5))
length(grep("mm", vista.all.data[, 1]))
length(grep("hs", vista.all.data[, 1]))

#
table(vista.all.data[grep("hs", vista.all.data[, 1]), 2] == "")
table(vista.all.data[grep("mm", vista.all.data[, 1]), 2] == "")
human.mapped.ids <- vista.all.data[vista.all.data[, 2] != "", 1]
length(human.mapped.ids)

#
human.file <- "/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/vista_processing/human_positive_vista.txt"
human.con <- file(human.file, open="r")
lines <- readLines(human.con)
num.lines <- length(lines)
headers <- lines[substr(lines, 0, 1) == ">"]
headers <- headers[grep(">Human\\|", headers)]
length(headers)
hs.positive.ids <- c()
human.regions <- c()
for (i in 1:length(headers)) {
  
  current <- headers[i]
  loc <- unlist(strsplit(current, split=">Human\\|"))[2]
  loc <- unlist(strsplit(loc, split=" "))[1]
  chr <- unlist(strsplit(loc, split=":"))[1]
  coords <- unlist(strsplit(loc, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  end <- unlist(strsplit(coords, split="-"))[2]
  other <- trimws(unlist(strsplit(current, split="\\|")))
  id <- other[3]
  id <- gsub("element ", "hs", id)
  associations <- paste0(other[5:length(other)], collapse=",")
  human.regions <- rbind(human.regions, c(chr, start, end, id, associations))
  
}
close(human.con)
dim(human.regions)
table(human.regions[, 4] %in% human.mapped.ids)

#
both.file <- "/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/vista_processing/both_positive_vista.txt"
both.con <- file(both.file, open="r")
lines <- readLines(both.con)
num.lines <- length(lines)
headers <- lines[substr(lines, 0, 1) == ">"]
headers <- headers[grep(">Mouse\\|", headers)]
length(headers)
mm.positive.ids <- c()
mouse.regions <- c()
for (i in 1:length(headers)) {
  
  current <- headers[i]
  loc <- unlist(strsplit(current, split=">Mouse\\|"))[2]
  loc <- unlist(strsplit(loc, split=" "))[1]
  chr <- unlist(strsplit(loc, split=":"))[1]
  coords <- unlist(strsplit(loc, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  end <- unlist(strsplit(coords, split="-"))[2]
  other <- trimws(unlist(strsplit(current, split="\\|")))
  id <- other[3]
  id <- gsub("element ", "mm", id)
  associations <- paste0(other[5:length(other)], collapse=",")
  mouse.regions <- rbind(mouse.regions, c(chr, start, end, id, associations))
  
}
close(both.con)
dim(mouse.regions)
table(mouse.regions[, 4] %in% human.mapped.ids)

#
positive.element.ids <- c(human.regions[, 4], mouse.regions[mouse.regions[, 4] %in% human.mapped.ids, 4])
length(unique(positive.element.ids))

#
all.regions <- rbind(human.regions, mouse.regions)
positive.element.bed <- c()
for (i in 1:length(positive.element.ids)) {
  
  #
  current.element.id <- positive.element.ids[i]
  position <- vista.all.data[vista.all.data[, 1] == current.element.id, 2]
  if (position == "")
    print("Error")
  
  #
  chr <- unlist(strsplit(position, split=":"))[1]
  coords <- unlist(strsplit(position, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  end <- unlist(strsplit(coords, split="-"))[2]
  
  #
  annotation <- all.regions[all.regions[, 4] == current.element.id, 5]
  
  #
  positive.element.bed <- rbind(positive.element.bed, c(chr, start, end, current.element.id, annotation))
  
}

#
positive.element.bed <- data.frame(positive.element.bed, stringsAsFactors=F)
positive.element.bed[, 1] <- as.character(positive.element.bed[, 1])
positive.element.bed[, 2] <- as.numeric(positive.element.bed[, 2])
positive.element.bed[, 3] <- as.numeric(positive.element.bed[, 3])
positive.element.bed[, 4] <- as.character(positive.element.bed[, 4])
positive.element.bed[, 5] <- as.character(positive.element.bed[, 5])
write.table(positive.element.bed, file="/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/vista_processing_new/vista.positive.human.elements.bed",
            sep="\t", quote=F, row.names=F, col.names=F)

# For each element, figure out its nearest gene, whether it's interacting in each cell type, and what genes it is interacting with.
vista.results <- matrix(NA, length(positive.element.bed[, 1]), 4+length(cell.types)*2)
rownames(vista.results) <- positive.element.bed[, 4]
colnames(vista.results) <- c("nearest_gene", "nearest_gene_dist", "min_resolution_genes", cell.types, paste0(cell.types, "_genes"), "ambiguity")
vista.results <- data.frame(vista.results, stringsAsFactors=F)
vista.results[, 1] <- as.character(vista.results[, 1])
vista.results[, 2] <- as.character(vista.results[, 2])
vista.results[, 3] <- as.character(vista.results[, 3])
vista.results[, 4] <- as.logical(vista.results[, 4])
vista.results[, 5] <- as.logical(vista.results[, 5])
vista.results[, 6] <- as.logical(vista.results[, 6])
vista.results[, 7] <- as.character(vista.results[, 7])
vista.results[, 8] <- as.character(vista.results[, 8])
vista.results[, 9] <- as.character(vista.results[, 9])
vista.results[, 10] <- as.logical(vista.results[, 10])

# Calculate nearest and same fragment genes for each Vista region.
vista.positions <- unique(positive.element.bed[, 1:4])
vista.positions.sorted <- sortBed(vista.positions)
all.promoters.sorted <- sortBed(promoters)
query.chr <- unique(vista.positions[, 1])
all.promoters.sorted.filtered <- all.promoters.sorted[all.promoters.sorted[, 1] %in% query.chr, ]
closest.results.all <- bedTools.2in(functionstring="bedtools closest", bed1=vista.positions.sorted, bed2=all.promoters.sorted.filtered, opt.string="-t first -d")

# Calculate same fragment gene intersections.
vista.sizes <- vista.positions[, 3] - vista.positions[, 2]
vista.resolution <- 5000
vista.positions.expanded <- vista.positions
widths <- abs(vista.positions.expanded[, 2] - vista.positions.expanded[, 3])
centers <- round((vista.positions.expanded[, 2] + vista.positions.expanded[, 3])/2)
vista.positions.expanded[widths < vista.resolution, 2] <- centers[widths < vista.resolution] - vista.resolution/2
vista.positions.expanded[widths < vista.resolution, 3] <- centers[widths < vista.resolution] + vista.resolution/2
vista.positions.expanded[vista.positions.expanded[, 2] < 1, 2] <- 1
table((vista.positions.expanded[, 3] - vista.positions.expanded[, 2]) > vista.resolution)
promoters.with.id <- cbind(promoters, paste0("promoter", 1:length(promoters[, 1])))
vista.min.resolution.intersections <- bedTools.2in(bed1=vista.positions.expanded, bed2=promoters.with.id, opt.string="-wb")

#
print(length(vista.results[, 1]))
for (i in 1:length(vista.results[, 1])) {
  
  #
  if (i %% 100 == 0)
    print(i)
  
  #
  current.id <- rownames(vista.results)[i]
  
  #
  vista.results[i, 1] <- paste0(unique(as.character(closest.results.all[as.character(closest.results.all[, 4]) == current.id, 11])), collapse=",")
  vista.results[i, 2] <- unique(as.numeric(closest.results.all[as.character(closest.results.all[, 4]) == current.id, 12]))
  vista.results[i, 3] <- paste0(unique(as.character(vista.min.resolution.intersections[as.character(vista.min.resolution.intersections[, 4]) == current.id, 11])), collapse=",")
  
  #
  for (cell.type in cell.types) {
    
    interacting.genes <- c()
    interacting.promoter.ids <- c()
    
    lhs.intersection <- bedTools.2in.compact(bed1=interactions.res[[cell.type]][, c(1:3, 11)], bed2=positive.element.bed[i, 1:4], opt.string="-wb")
    if (length(lhs.intersection[, 1]) > 0) {
      
      vista.results[i, cell.type] <- TRUE
      interaction.ids <- as.character(lhs.intersection[, 4])
      selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
      rhs.genes <- paste0(as.character(interactions.ann[[cell.type]]$rhs_promoter_ids[selector]), collapse=",")
      interacting.genes <- c(interacting.genes, rhs.genes)
      
      rhs.coords <- interactions.res[[cell.type]][selector, 5:7]
      promoter.results <- bedTools.2in(bed1=promoters.with.id, bed2=rhs.coords, opt.string="-c")
      interacting.promoter.ids <- c(interacting.promoter.ids, as.character(promoter.results[promoter.results[, 9] > 0, 8]))
      
    }
    
    rhs.intersection <- bedTools.2in.compact(bed1=interactions.res[[cell.type]][, c(5:7, 11)], bed2=positive.element.bed[i, 1:4], opt.string="-wb")
    if (length(rhs.intersection[, 1]) > 0) {
      
      vista.results[i, cell.type] <- TRUE
      interaction.ids <- as.character(rhs.intersection[, 4])
      selector <- interactions.res[[cell.type]]$ID %in% interaction.ids
      lhs.genes <- paste0(as.character(interactions.ann[[cell.type]]$lhs_promoter_ids[selector]), collapse=",")
      interacting.genes <- c(interacting.genes, lhs.genes)
      
      lhs.coords <- interactions.res[[cell.type]][selector, 1:3]
      promoter.results <- bedTools.2in(bed1=promoters.with.id, bed2=lhs.coords, opt.string="-c")
      interacting.promoter.ids <- c(interacting.promoter.ids, as.character(promoter.results[promoter.results[, 9] > 0, 8]))
      
    }
    
    if ((length(lhs.intersection[, 1]) == 0) & (length(rhs.intersection[, 1]) == 0))
      vista.results[i, cell.type] <- FALSE
    
    interacting.genes <- unlist(strsplit(paste0(interacting.genes, collapse=","), split=","))
    interacting.genes <- promoters$gene_name[promoters$gene_id %in% interacting.genes]
    interacting.genes <- unique(interacting.genes)
    interacting.genes <- interacting.genes[interacting.genes != ""]
    vista.results[i, paste0(cell.type, "_genes")] <- paste0(interacting.genes, collapse=",")
    
  }
  
  # Evalute minimum resolution ambiguity.
  vista.coordinates <- positive.element.bed[i, 1:4]
  nearest.gene.coordinates <- promoters[promoters$gene_name == vista.results[i, 1], 1:3]
  gene.rmap.expanded <- bedTools.2in(bed1=nearest.gene.coordinates, bed2=rmap, opt.string="-wb")
  vista.rmap.expanded <- bedTools.2in.compact(bed1=vista.coordinates, bed2=rmap, opt.string="-wb")
  gene.rmap.ids <- unique(gene.rmap.expanded[, 7])
  vista.rmap.ids <- unique(vista.rmap.expanded[, 8])
  #print(length(vista.rmap.ids))
  if ((length(gene.rmap.ids) == 1) & (length(vista.rmap.ids) == 1) & (gene.rmap.ids[1] == vista.rmap.ids[1])) {
    vista.results[i, 10] <- TRUE
  } else {
    vista.results[i, 10] <- FALSE
  }
  
}

# Checkpoint.
save(vista.results, file=paste0(output.dir, "/", run.id, "/saved/vista.results.Rdata"))
load(file=paste0(output.dir, "/", run.id, "/saved/vista.results.Rdata"))

# Count how many elements appear in interactions for each cell type.
vista.ids.per.cell.type <- list()
for (cell.type in cell.types) {
  
  vista.ids.per.cell.type[[cell.type]] <- as.character(rownames(vista.results)[which(vista.results[, cell.type])])
  vista.results[is.na(vista.results[, cell.type]), cell.type] <- FALSE
  
}

#
pdf(file=paste0(output.dir, "/", run.id, "/vista/interacting_vista_regions_by_cell_type.pdf"), width=5, height=5)
venn.plot <- venn.diagram(vista.ids.per.cell.type , NULL, fill=c("#BFEFFF", "#FFD39B", "#98FB98"), alpha=c(0.5, 0.5, 0.5), 
                          cex=2, category.names=c("cortical", "hippocampal", "astrocyte"), 
                          main="Overlap of Vista regions with significant interactions")
grid.draw(venn.plot)
dev.off()

#
all.vista.ids <- unique(c(vista.ids.per.cell.type[["cortical"]], vista.ids.per.cell.type[["hippocampal"]], vista.ids.per.cell.type[["astrocyte"]]))
length(all.vista.ids)
compiled.results <- cbind(positive.element.bed, vista.results)
write.table(compiled.results, file=paste0(output.dir, "/", run.id, "/vista/vista_results_all.txt"), quote=F, row.names=F, col.names=F, sep="\t")
positive.results <- compiled.results[apply(compiled.results[, cell.types], 1, any), ]
write.table(positive.results, file=paste0(output.dir, "/", run.id, "/vista/vista_results_positive.txt"), quote=F, row.names=F, col.names=F, sep="\t")

#
interacting <- rep(FALSE, length(vista.results[, 1]))
nearest.agreement <- rep(FALSE, length(vista.results[, 1]))
num.targets <- rep(NA, length(vista.results[, 1]))
min.resolution.ambiguity <- rep(FALSE, length(vista.results[, 1]))
targeted.genes <- c()
for (i in 1:length(vista.results[, 1])) {
  
  #
  interacting.genes <- unique(unlist(strsplit(paste0(vista.results[i, paste0(cell.types, "_genes")], collapse=","), split=",")))
  interacting.genes <- interacting.genes[interacting.genes != ""]
  num.targets[i] <- length(interacting.genes)
  targeted.genes <- c(targeted.genes, interacting.genes)
  
  #
  min.resolution.ambiguity[i] <- vista.results[i, 10]
  
  #
  if (vista.results[i, 1] %in% interacting.genes)
    nearest.agreement[i] <- TRUE
  if (length(interacting.genes) > 0)
    interacting[i] <- TRUE
  
}

#
table((interacting == TRUE) & (min.resolution.ambiguity == TRUE))[["TRUE"]]
union <- length(which((interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == TRUE)))
gwas.nearest <- length(which((interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == FALSE)))
new.targets <- sum(num.targets[(interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == TRUE)] - 1) +
  sum(num.targets[(interacting == TRUE) & (min.resolution.ambiguity == FALSE) & (nearest.agreement == FALSE)])

# Draw Venn diagram.
pdf(paste0(output.dir, "/", run.id, "/vista/vista_agreement_venn.pdf"), width=5, height=5)
values <- c(gwas.nearest, union, new.targets)
names(values) <- c("Nearest", "Nearest&New interaction targets", "New interaction targets")
fit <- euler(values)
colors <- c("#FFE4C4", "#00FF7F")
plot(fit, fill_opacity = 0.01, shape="circle", title=paste0("Agreement for Vista regions"),
     fill=colors, border="transparent", fontsize=12, quantities=list(fontsize=12))
dev.off()


# Allelic bias annotation -------------------------------------------------


#
filename <- paste0("allelic.bias.significant.interactions.", run.id, ".", interaction.resolution, ".", score.cutoff, ".tsv")
custom.allelic <- read.table(paste0(output.dir, "/common/allelic_bias/", filename), sep="\t", header=F, stringsAsFactors=F)
allelic.interactions <- data.frame(matrix(NA, length(custom.allelic[, 1]), 11), stringsAsFactors=F)
allelic.interactions[, 1] <- as.character(allelic.interactions[, 1])
allelic.interactions[, 2] <- as.numeric(allelic.interactions[, 2])
allelic.interactions[, 3] <- as.numeric(allelic.interactions[, 3])
allelic.interactions[, 4] <- as.character(allelic.interactions[, 4])
allelic.interactions[, 5] <- as.character(allelic.interactions[, 5])
allelic.interactions[, 6] <- as.numeric(allelic.interactions[, 6])
allelic.interactions[, 7] <- as.numeric(allelic.interactions[, 7])
allelic.interactions[, 8] <- as.character(allelic.interactions[, 8])
allelic.interactions[, 9] <- as.numeric(allelic.interactions[, 9])
allelic.interactions[, 10] <- as.numeric(allelic.interactions[, 10])
allelic.interactions[, 11] <- as.character(allelic.interactions[, 11])
allelic.interactions[, 1] <- as.character(custom.allelic[, 1])
allelic.interactions[, 2] <- as.numeric(custom.allelic[, 2])
allelic.interactions[, 3] <- as.numeric(custom.allelic[, 3])
allelic.interactions[, 4] <- rep(".", length(allelic.interactions[, 1]))
allelic.interactions[, 5] <- as.character(custom.allelic[, 4])
allelic.interactions[, 6] <- as.numeric(custom.allelic[, 5])
allelic.interactions[, 7] <- as.numeric(custom.allelic[, 6])
allelic.interactions[, 8] <- rep(".", length(allelic.interactions[, 1]))
allelic.interactions[, 9] <- as.numeric(custom.allelic[, 8]) + as.numeric(custom.allelic[, 9])
allelic.interactions[, 10] <- as.numeric(custom.allelic[, 11])
allelic.interactions[, 11] <- as.character(custom.allelic[, 7])
colnames(allelic.interactions) <- colnames(interactions.res[["cortical"]])

#
allelic.annotation <- annotateInteractions(allelic.interactions, promoters, atac.seq.peaks.res[["cortical"]], atac.seq.peaks.ann[["cortical"]], 
                                           features[["cortical"]])
allelic.res <- list()
allelic.res[["cortical"]] <- allelic.interactions
allelic.ann <- list()
allelic.ann[["cortical"]] <- allelic.annotation
writeFeatureOverlapResults(interactions.res=allelic.res, interactions.ann=allelic.ann, output.prefix=paste0(output.dir, "/", run.id, "/interaction_results/"))


# Miscellaneous -----------------------------------------------------------


#
plotAstrocyteMarkers(paste0(output.dir, "/", run.id, "/RNA-seq_results/astrocyte_markers.csv"), paste0(output.dir, "/", run.id, "/RNA-seq_results"))

