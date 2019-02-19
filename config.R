#
# config.R
# Author: Michael Song
# Last modified: 2018-12-17
# This file contains settings for the main analysis pipeline.
#


# Environmental settings --------------------------------------------------


#
run.id <- "brain_revision"
interactions.list <- "interactions_fastp_250_0.txt"
atac.seq.peak.list <- "atac_seq_peaks_optimal.txt"
rna.seq.output.list.replicates <- "rna_seq_output_replicates.txt"
features.list <- "features_list.txt"
annotation.file <- "gencode/gencode_v19_annotation.bed"
exons.file <- "gencode/gencode_v19_exons.bed"
introns.file <- "gencode/gencode_v19_introns.bed"
chrom.sizes.file <- "hg19.chrom.sizes"
rmap.file <- "hg19.rmap"
baitmap.file <- "hg19.baitmap"
tads.file <- "tads/tads.CP.bed"
cell.type.colors <- c("#1667A8", "#D4542D", "#278A2A", "#F0AD0E")
vista.file <- "vista/all.elements.bed"


# Analysis settings -------------------------------------------------------


# Chromosomes included in the analysis.
allowed.chrs <- paste0("chr", c(1:22, "X", "Y"))

# Score cutoff for filtering significant interactions.
score.cutoff <- 5

# Defines how far a promoter extends upstream/downstream from each TSS for intersection with ATAC-seq peaks/interactions.
tss.upstream <- 500            
tss.downstream <- 500

# Defines minimum width of an ATAC-seq peak.
atac.seq.peak.width <- 500

# Defines minimum resolution of interacting bins (all bins smaller than this are expanded to the minimum resolution for intersection with features).
interaction.resolution <- 5000 

# Defines whether or not to consider neighboring bins as interacting, regardless of 'interaction.resolution'.
neighboring.bins <- FALSE

# Defines maximum bin size for annotating interactions between promoters and PIRs.
max.bin.size <- 100000

# Defines minimum overlap threshold in bedtools for intersections between bins or between bins and features.
overlap.threshold <- 10^-9

# Defines cutoffs for determining whether or not a gene is significantly expressed.
tpm.cutoff <- 0.5
rpkm.cutoff <- 0.5

# Defines how many interactions to downsample each interaction set to before further analysis.
downsampling.depth <- NULL

# Defines how many ATAC-seq peaks to downsample each ATAC-seq peak set to before further analysis.
max.atac.seq.peaks <- NULL

# List of diseases to analyze.
diseases <- c("AD", "ADHD", "ALS", "ASD", "BD", "EP", "FTD", "MP", "PD", "SCZ", "UD")


# End ---------------------------------------------------------------------

