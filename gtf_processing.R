#
# gtf_processing.R
# Author: Michael Song
# Last modified: 2019-01-15
# This script processes a GTF file for use with the main analysis script.
#


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


# Main --------------------------------------------------------------------


# Define GTF file to be processed.
gtf.file <- paste0(home.dir, "/resources/gencode/gencode.v29lift37.annotation.gtf")
gtf.version <- "V29"

# Read in GTF file downloaded from GENCODE.
gtf.raw <- read.table(paste0(gtf.file), stringsAsFactors=F, sep="\t", header=F)
gencode.transcripts <- gtf.raw[gtf.raw[, 3] == "transcript", ]

# Parse transcripts in GTF file to get metadata for each entry.
num.transcripts <- length(gencode.transcripts[, 1])
key.value.pairs <- c("gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type", "transcript_status", "transcript_name", "exon_number", "exon_id", "level",
                     "tag", "ccdsid", "havana_gene", "havana_transcript", "protein_id", "ont", "transcript_support_level",
                     "remap_status", "remap_original_id", "remap_original_location", "remap_num_mappings", "remap_target_status", "remap_substituted_missing_target")
processed.transcripts <- cbind(gencode.transcripts[, 1:8], matrix("", num.transcripts, 24))
colnames(processed.transcripts) <- c("chrom", "source", "feature_type", "start", "end", "score", "strand", "phase", key.value.pairs)
for (key.value.pair in key.value.pairs)
  processed.transcripts[, key.value.pair] <- as.character(processed.transcripts[, key.value.pair])

# Process key value pairs for each transcript one by one.
print(paste0(num.transcripts, " entries to process in the GTF file..."))
for (i in 1:num.transcripts) {
  
  # Track progress.
  if (i %% 10000 == 0)
    print(paste0(i, " transcripts processed"))
  
  # Explode 9th column with key-value pairs.
  if (substr(gencode.transcripts[i, 9], nchar(gencode.transcripts[i, 9]), nchar(gencode.transcripts[i, 9])) == ";")
    gencode.transcripts[i, 9] <- substr(gencode.transcripts[i, 9], 1, nchar(gencode.transcripts[i, 9]) - 1)
  fields <- trimws(unlist(strsplit(gencode.transcripts[i, 9], split=";")))
  expanded.keys <- unlist(strsplit(fields, split=" "))[(1:length(fields))*2 - 1]
  expanded.values <- unlist(strsplit(fields, split=" "))[(1:length(fields))*2]
  
  # Fill in the table with the expanded fields.
  for (j in 1:length(expanded.keys)) {
    
    # Append if there is already an existing entry.
    if (processed.transcripts[i, expanded.keys[j]] == "") {
      processed.transcripts[i, expanded.keys[j]] <- expanded.values[j]
    } else {
      processed.transcripts[i, expanded.keys[j]] <- paste0(processed.transcripts[i, expanded.keys[j]], ";", expanded.values[j])
    }
    
    # Identify unlisted key-value pairs.
    # if ((expanded.keys[j] %in% colnames(processed.transcripts) == FALSE) && (!is.na(expanded.keys[j])))
    #   print(expanded.keys[j])
    
  }
  
}

# Save results.
save(processed.transcripts, file=paste0(home.dir, "/resources/gencode/processed.transcripts.", gtf.version, ".Rdata"))
load(file=paste0(home.dir, "/resources/gencode/processed.transcripts.", gtf.version, ".Rdata"))

# Stitch together BED file based on the strand information.
tss.bed.plus <- processed.transcripts[processed.transcripts$strand == "+", c("chrom", "start", "start", "strand", "gene_id", "gene_type", "gene_name")]
tss.bed.minus <- processed.transcripts[processed.transcripts$strand == "-", c("chrom", "end", "end", "strand", "gene_id", "gene_type", "gene_name")]
tss.bed <- unique(rbind(tss.bed.plus, setNames(tss.bed.minus, names(tss.bed.plus))))
colnames(tss.bed) <- c("chrom", "start", "end", "strand", "gene_id", "gene_type", "gene_name")

# Filter BED entries.
allowed.chrs <- paste0("chr", c(as.character(c(1:22)), "X", "Y"))
tss.bed <- tss.bed[tss.bed[, 1] %in% allowed.chrs, ]

# Remove version number for the Ensembl gene ID.
tss.bed$gene_id <- unlist(strsplit(tss.bed$gene_id, split="\\."))[(1:length(tss.bed$gene_id))*2 - 1]

# Write BED file.
write.table(tss.bed, file=paste0(gtf.file, ".bed"), sep="\t", row.names=F, col.names=F, quote=F)


# End ---------------------------------------------------------------------

