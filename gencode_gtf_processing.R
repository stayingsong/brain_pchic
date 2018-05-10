#
# gencode_gtf_processing.R
# Author: Michael Song
# Date created: 2018-01-05
# Last modified: 2018-01-05
# This script processes a GENCODE GTF file to a more easily readable BED format for downstream processing.
#


# Miscellaneous settings --------------------------------------------------


# Clear workspace.
rm(list=ls())

# Turn off scientific notation for output/printing.
options(scipen=999)


# Processing --------------------------------------------------------------


# Process GENCODE GTF file for translating gene names and Ensembl gene IDs.
setwd("/Users/Michael/Box Sync/MS_analysis/brain_pchic/supplementary_files")
gtf_file <- "gencode.v19.annotation"
gencode_gtf <- read.table(paste0(gtf_file, ".gtf"), stringsAsFactors=F, sep="\t", header=F)
gencode_genes <- gencode_gtf[gencode_gtf[, 3] == "gene", ]
gencode_transcripts <- gencode_gtf[gencode_gtf[, 3] == "transcript", ]

# Parse GENCODE transcripts to get TSS, Ensembl gene ID, Ensembl transcript IDs, and HGNC symbol for each entry.
num_transcripts <- length(gencode_transcripts[, 1])
gtf_transcripts <- cbind(gencode_transcripts[, 1:8], matrix("", num_transcripts, 18))
gtf_transcripts[, 9] <- as.character(gtf_transcripts[, 9])
gtf_transcripts[, 10] <- as.character(gtf_transcripts[, 10])
gtf_transcripts[, 11] <- as.character(gtf_transcripts[, 11])
gtf_transcripts[, 12] <- as.character(gtf_transcripts[, 12])
gtf_transcripts[, 13] <- as.character(gtf_transcripts[, 13])
gtf_transcripts[, 14] <- as.character(gtf_transcripts[, 14])
gtf_transcripts[, 15] <- as.character(gtf_transcripts[, 15])
gtf_transcripts[, 16] <- as.character(gtf_transcripts[, 16])
gtf_transcripts[, 17] <- as.character(gtf_transcripts[, 17])
gtf_transcripts[, 18] <- as.character(gtf_transcripts[, 18])
gtf_transcripts[, 19] <- as.character(gtf_transcripts[, 19])
gtf_transcripts[, 20] <- as.character(gtf_transcripts[, 20])
gtf_transcripts[, 21] <- as.character(gtf_transcripts[, 21])
gtf_transcripts[, 22] <- as.character(gtf_transcripts[, 22])
gtf_transcripts[, 23] <- as.character(gtf_transcripts[, 23])
gtf_transcripts[, 24] <- as.character(gtf_transcripts[, 24])
gtf_transcripts[, 25] <- as.character(gtf_transcripts[, 25])
gtf_transcripts[, 26] <- as.character(gtf_transcripts[, 26])
colnames(gtf_transcripts) <- c("chrom", "source", "feature_type", "start", "end", "score", "strand", "phase",
                               "gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", "transcript_type", 
                               "transcript_status", "transcript_name", "exon_number", "exon_id", "level",
                               "tag", "ccdsid", "havana_gene", "havana_transcript", "protein_id", "ont", "transcript_support_level")
print(num_transcripts)
for (i in 1:num_transcripts) {
  
  #
  if (i %% 10000 == 0)
    print(i)
  
  #
  annotations <- trimws(unlist(strsplit(gencode_transcripts[i, 9], split=";")))
  expanded_names <- unlist(strsplit(annotations, split=" "))[(1:10)*2-1]
  expanded_values <- unlist(strsplit(annotations, split=" "))[(1:10)*2]
  
  #
  for (j in 1:length(expanded_names)) {
    gtf_transcripts[i, which(colnames(gtf_transcripts) == expanded_names[j])] <- expanded_values[j]
    if ((expanded_names[j] %in% colnames(gtf_transcripts) == FALSE) && (!is.na(expanded_names[j]))) {
      print(expanded_names[j])
    }
    
  }
  
}

#
save(gtf_transcripts, file=paste0(gtf_file, ".parsed.Rdata"))

# Write BED file.
setwd("/Users/Michael/Box Sync/MS_analysis/brain_pchic/supplementary_files")
gtf_file <- "gencode.v19.annotation"
load(file=paste0(gtf_file, ".parsed.Rdata"))
TSS_bed_plus <- gtf_transcripts[gtf_transcripts$strand == "+", c(1, 4, 4, 7, 9, 11, 13)]
colnames(TSS_bed_plus) <- c("chrom", "start", "end", "strand", "gene_id", "gene_type", "gene_name")
TSS_bed_minus <- gtf_transcripts[gtf_transcripts$strand == "-", c(1, 5, 5, 7, 9, 11, 13)]
colnames(TSS_bed_minus) <- c("chrom", "start", "end", "strand", "gene_id", "gene_type", "gene_name")
TSS_bed <- rbind(TSS_bed_plus, TSS_bed_minus)
allowed.chrs <- c(as.character(c(1:22)), "X", "Y")
allowed.chrs <- paste0("chr", allowed.chrs)
TSS_bed <- TSS_bed[!((TSS_bed[, 2] == 0) | (TSS_bed[, 3] == 0)), ]
TSS_bed <- TSS_bed[TSS_bed[, 1] %in% allowed.chrs, ]
TSS_bed$gene_id <- unlist(strsplit(TSS_bed$gene_id, split="\\."))[(1:length(TSS_bed$gene_id))*2-1]
write.table(TSS_bed, file=paste0(gtf_file, ".bed"), quote=F, sep="\t", row.names=F, col.names=F)
