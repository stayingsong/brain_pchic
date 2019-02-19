#
# process_vista_elements.R
# Author: Michael Song
# Last modified: 2019-01-14
# This script processes data (all positive elements with human homologues) from the Vista Enhancer Browser.
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


# Read in data (all entries from Vista Enhancer Browser). CSV was prepared by copying table in "Experimental Data" page.
vista.data <- read.table(paste0(home.dir, "/resources/vista/vista_all.csv"),
                         sep=",", stringsAsFactors=F, header=F, colClasses=rep("character", 5))
print(paste0("Read in ", length(grep("hs", vista.data[, 1])), " human elements."))
print(paste0("Read in ", length(grep("mm", vista.data[, 1])), " mouse elements."))
print(paste0(length(which(vista.data[grep("hs", vista.data[, 1]), 2] != "")), " of the human elements have human homologues."))
print(paste0(length(which(vista.data[grep("mm", vista.data[, 1]), 2] != "")), " of the mouse elements have human homologues."))
human.mapped.ids <- vista.data[vista.data[, 2] != "", 1]
print(paste0("There are a total of ", length(human.mapped.ids), " vista elements with coordinates mapped in humans."))

# Read in and process headers for all positive elements with human homologues.
both.elements.file <- paste0(home.dir, "/resources/vista/vista_both_positive.txt")
both.elements.con <- file(both.elements.file, open="r")
both.elements.lines <- readLines(both.elements.con)
close(both.elements.con)
both.elements.headers <- both.elements.lines[substr(both.elements.lines, 0, 1) == ">"]
human.elements.headers <- both.elements.headers[grep(">Human\\|", both.elements.headers)]
mouse.elements.headers <- both.elements.headers[grep(">Mouse\\|", both.elements.headers)]
print(paste0("Processing details for ", length(human.elements.headers), " human positive elements."))
print(paste0("Processing details for ", length(mouse.elements.headers), " mouse positive elements."))

# For each human positive element, process its header and grab the BED coordinates, ID, and annotation.
human.elements.bed <- c()
for (i in 1:length(human.elements.headers)) {
  
  # Get BED coordinates.
  current.header <- human.elements.headers[i]
  pos <- unlist(strsplit(unlist(strsplit(current.header, split=">Human\\|"))[2], split=" "))[1]
  chr <- unlist(strsplit(pos, split=":"))[1]
  coords <- unlist(strsplit(pos, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  end <- unlist(strsplit(coords, split="-"))[2]
  
  # Get ID and annotation.
  metadata <- trimws(unlist(strsplit(current.header, split="\\|")))
  id <- metadata[3]
  id <- gsub("element ", "hs", id)
  annotation <- paste0(metadata[5:length(metadata)], collapse=",")
  human.elements.bed <- rbind(human.elements.bed, c(chr, start, end, id, annotation))
  
}

# For each mouse positive element with human homologues, process its header and grab the BED coordinates, ID, and annotation.
mouse.elements.bed <- c()
for (i in 1:length(mouse.elements.headers)) {

  # Get BED coordinates.
  current.header <- mouse.elements.headers[i]
  pos <- unlist(strsplit(unlist(strsplit(current.header, split=">Mouse\\|"))[2], split=" "))[1]
  chr <- unlist(strsplit(pos, split=":"))[1]
  coords <- unlist(strsplit(pos, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  end <- unlist(strsplit(coords, split="-"))[2]
  
  # Get ID and annotation.
  metadata <- trimws(unlist(strsplit(current.header, split="\\|")))
  id <- metadata[3]
  id <- gsub("element ", "mm", id)
  annotation <- paste0(metadata[5:length(metadata)], collapse=",")
  mouse.elements.bed <- rbind(mouse.elements.bed , c(chr, start, end, id, annotation))

}

# Print results.
print(paste0(length(which(human.elements.bed[, 4] %in% human.mapped.ids)), " human positive elements in list of human mapped IDs."))
print(paste0(length(which(mouse.elements.bed[, 4] %in% human.mapped.ids)), " mouse positive elements in list of human mapped IDs."))

# Combine BED files.
unmapped.elements.bed <- mouse.elements.bed[!(mouse.elements.bed[, 4] %in% human.mapped.ids), ]
mouse.elements.bed <- mouse.elements.bed[mouse.elements.bed[, 4] %in% human.mapped.ids, ]
all.elements.bed <- data.frame(rbind(human.elements.bed, mouse.elements.bed), stringsAsFactors=F)
num.vista.elements <- length(all.elements.bed[, 1])
print(paste0(num.vista.elements, " total positive elements with coordinates in hg19."))
colnames(all.elements.bed) <- c("chrom", "chromStart", "chromEnd", "vista_ID", "annotation")
all.elements.bed$chrom <- as.character(all.elements.bed$chrom)
all.elements.bed$chromStart <- as.integer(all.elements.bed$chromStart)
all.elements.bed$chromEnd <- as.integer(all.elements.bed$chromEnd)
all.elements.bed$vista_ID <- as.character(all.elements.bed$vista_ID)
all.elements.bed$annotation <- as.character(all.elements.bed$annotation)

# Check hg19 coordinates for all elements (also grab hg19 coordinates for mouse elements with mm9 coordinates).
for (i in 1:num.vista.elements) {
  
  # Get current element.
  current.element <- all.elements.bed[i, ]
  
  # Grab reference hg19 coordinates.
  reference.entry <- vista.data[vista.data[, 1] %in% current.element[4], 2]
  chr <- unlist(strsplit(reference.entry, split=":"))[1]
  coords <- unlist(strsplit(reference.entry, split=":"))[2]
  start <- unlist(strsplit(coords, split="-"))[1]
  start <- as.integer(gsub(",", "", start))
  end <- unlist(strsplit(coords, split="-"))[2]
  end <- as.integer(gsub(",", "", end))
  
  # Make sure the hg19 coordinates agree, fixing if needed.
  if ((current.element[1] != chr) | (current.element[2] != start) | (current.element[3] != end)) {
    
    all.elements.bed[i, 1] <- chr
    all.elements.bed[i, 2] <- start
    all.elements.bed[i, 3] <- end
    
  }
  
}

# Write results to file.
write.table(all.elements.bed, file=paste0(home.dir, "/resources/vista/all.elements.bed"), sep="\t", row.names=F, col.names=F, quote=F)
save(all.elements.bed, file=paste0(home.dir, "/resources/vista/all.elements.bed.Rdata"))

# Write.results for display.
vista.display.bed <- cbind(all.elements.bed[, 1:4], rep(".", length(all.elements.bed[, 1])), rep(".", length(all.elements.bed[, 1])))
write.table(vista.display.bed, file=paste0(home.dir, "/resources/vista/vista.display.bed"), sep="\t", row.names=F, col.names=F, quote=F)


# End ---------------------------------------------------------------------

