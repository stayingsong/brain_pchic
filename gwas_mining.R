#
# gwas_mining.R
# Author: Michael Song
# Date created: 2018-03-20
# Last modified: 2018-05-10
# This script mines SNPs from the GWAS Catalog, then imputes them and formats them for downstream analysis.
#


# Miscellaneous settings --------------------------------------------------


# Clear workspace.
rm(list=ls())

# Turn off scientific notation for output/printing.
options(scipen=999)


# Load packages -----------------------------------------------------------


# Install and load packages.
#source("https://bioconductor.org/biocLite.R")
#biocLite(haploR)
require(haploR)
#biocLite("GenomicRanges")
require(GenomicRanges)

#
#source("http://bioconductor.org/workflows.R")
#workflowInstall("liftOver")
require(rtracklayer)


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


# SNP imputation ----------------------------------------------------------


#
translate <- read.table("/Users/Michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/InfiniumPsychArray-24v1-2_A1_b144_rsids.txt", 
                        sep="\t", header=T, quote = "", stringsAsFactors=F)

# Read in pivot tables and SNP tables.
sig.cutoff <- 10^-5
diseases <- c("AD", "ADHD", "ASD", "BD", "EP", "MDD", "PD", "SCZ", "T1D")
pivot.tables <- list()
snp.tables <- list()
pop.sorted.tables <- list()
imputation.results <- list()
pops <- c("AFR", "AMR", "ASN", "EUR")
reported.gene.table.all <- list()
for (i in 1:length(diseases)) {
  
  #
  current.disease <- diseases[i]
  reported.gene.table <- c()
  
  # Read in data.
  setwd("/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/gwas_catalog/raw")
  print(paste0("Reading in data for: ", current.disease))
  pivot.tables[[current.disease]] <- read.table(paste0(current.disease, ".pivot.csv"), sep=",", header=T, quote = "", stringsAsFactors=F, 
                                            colClasses=c("character", "character", "integer", "character", "integer"))
  snp.tables[[current.disease]] <- read.table(paste0("2018-03-20-", current.disease, ".tsv"), sep="\t", header=T, stringsAsFactors=F, quote="", fill=T)

  #
  print(dim(snp.tables[[current.disease]]))
  
  #
  for (i in 1:length(snp.tables[[current.disease]]$REPORTED.GENE.S.)) {
    rsids <- trimws(unlist(strsplit(snp.tables[[current.disease]]$SNPS[i], split=";")))
    rsids <- rsids[grep("rs", rsids)]
    reported.genes <- trimws(unlist(strsplit(snp.tables[[current.disease]]$REPORTED.GENE.S.[i], split=",")))
    if ((length(rsids) > 0) & (length(reported.genes) > 0)) {
      for (j in 1:length(rsids)) {
        for (k in 1:length(reported.genes)) {
          reported.genes[k] <- toupper(reported.genes[k])
          if (!is.na(rsids[j])) {
            if ((reported.genes[k] == "Intergenic") | (reported.genes[k] == "intergenic") | (reported.genes[k] == "INTERGENIC"))
              gene.name.correspondance <- "Intergenic"
            else if (reported.genes[k] == "NR")
              gene.name.correspondance <- "NR"
            else if (reported.genes[k] %in% toupper(promoters$gene_name))
              gene.name.correspondance <- "TRUE"
            else if (!(reported.genes[k] %in% toupper(promoters$gene_name)))
              gene.name.correspondance <- "FALSE"
            else
              print("Error")
            reported.gene.table <- rbind(reported.gene.table, c(rsids[j], reported.genes[k], gene.name.correspondance))
          }
        }
      }
    }
  }
  reported.gene.table[reported.gene.table[, 3] == "FALSE", ]
  reported.gene.table.all[[current.disease]] <- reported.gene.table

  #
  print(paste0("Number of SNPs with study information available: ", 
               table(as.character(snp.tables[[current.disease]]$PUBMEDID) %in% pivot.tables[[current.disease]]$PMID)["TRUE"]))
  
  # Filter out studies that do not relate to current disorder.
  pivot.tables[[current.disease]] <- pivot.tables[[current.disease]][pivot.tables[[current.disease]][, 5] == 1, ]
  print(length(unique(pivot.tables[[current.disease]][, 1])))
  
  # Filter SNPs according to p-value cutoff.
  filtered.snps <- snp.tables[[current.disease]][as.numeric(snp.tables[[current.disease]]$P.VALUE) < sig.cutoff, ]
  print(paste0("Number of SNPs passing significance cutoff of ", sig.cutoff, ": ", length(filtered.snps[, 1])))
  
  # Try to translate SNPs that don't have rsIDs.
  curated.snp.ids <- filtered.snps$SNP_ID_CURRENT
  snps <- filtered.snps$SNPS
  na.positions <- which(is.na(filtered.snps$SNP_ID_CURRENT))
  snps[na.positions]
  curated.snp.ids[na.positions]
  if (length(na.positions) > 0) {
    
    for (j in 1:length(na.positions)) {
      
      test <- grep("PsychChip", snps[na.positions[j]])
      query <- NULL
      if ((length(test) > 0) && (test == TRUE)) {
        query <- unlist(strsplit(snps[na.positions[j]], split="\\|"))[2]
      } else {
        query <- snps[na.positions[j]]
      }
      
      multi.split <- unlist(strsplit(query, split="; "))
      if (query %in% translate$Name) {
        result <- translate$RsID[translate$Name == query]
        curated.snp.ids[na.positions[j]] <- result
      } else if (length(multi.split) > 1) {
        for (k in 1:length(multi.split)) {
          multi.split[k] <- gsub("rs", "", multi.split[k])
        }
        curated.snp.ids[na.positions[j]] <- paste0(multi.split, collapse=",")
      } else {
        print(paste0("Could not translate: ", query, " (", na.positions[j], ")"))
      }
      
    }
    curated.snp.ids[na.positions] <- gsub("rs", "", curated.snp.ids[na.positions])
  
  }

  # Get PMIDs per population, then grab the SNPs for that study.
  pop.sorted.tables[[current.disease]] <- list()
  for (j in 1:length(pops)) {
    
    # For each study, find if it has current population, and if so, store the SNPs.
    pop.sorted.tables[[current.disease]][[pops[j]]] <- c()
    for (k in 1:length(pivot.tables[[current.disease]][, 1])) {
      current.pops <- unlist(strsplit(pivot.tables[[current.disease]][k, 2], split="_"))
      if (pops[j] %in% current.pops) {
        current.snps <- as.character(curated.snp.ids[filtered.snps$PUBMEDID == pivot.tables[[current.disease]][k, 1]])
        for (snp in current.snps) {
          exploded <- unlist(strsplit(snp, split=","))
          if (length(exploded) == 1) {
            pop.sorted.tables[[current.disease]][[pops[j]]] <- c(pop.sorted.tables[[current.disease]][[pops[j]]], exploded)
          } else if (length(exploded) > 1) {
            pop.sorted.tables[[current.disease]][[pops[j]]] <- c(pop.sorted.tables[[current.disease]][[pops[j]]], exploded)
            #print(paste0("Adding: ", paste0(exploded, collapse=" ")))
          } else {
            print("What")
          }
        }
      }
    }
    
  }
  
  #
  setwd("/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/gwas_catalog/population_sorted")
  
  #
  imputation.results[[current.disease]] <- c()
  if (length(pop.sorted.tables[[current.disease]][["AFR"]]) > 0) {
    
    AFR <- pop.sorted.tables[[current.disease]][["AFR"]][!is.na(pop.sorted.tables[[current.disease]][["AFR"]])]
    AFR <- paste0("rs", AFR)
    length(AFR)
    print(length(unique(AFR)))
    write.table(unique(AFR), file=paste0(current.disease, ".AFR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AFR.results <- queryHaploreg(query = unique(AFR), file = NULL, study = NULL, ldThresh = 0.8,
                                 ldPop = "AFR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                 url = "http://archive.broadinstitute.org/mammals/haploreg/haploreg.php",
                                 timeout = 100000, encoding = "UTF-8", verbose = FALSE)
    AFR.results <- data.frame(AFR.results)
    print(paste0("AFR imputation resulted in # rows: ", length(AFR.results[, 1])))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AFR.results)
    
  }
  if (length(pop.sorted.tables[[current.disease]][["AMR"]]) > 0) {
    
    AMR <- pop.sorted.tables[[current.disease]][["AMR"]][!is.na(pop.sorted.tables[[current.disease]][["AMR"]])]
    AMR <- paste0("rs", AMR)
    length(AMR)
    print(length(unique(AMR)))
    write.table(unique(AMR), file=paste0(current.disease, ".AMR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AMR.results <- queryHaploreg(query = unique(AMR), file = NULL, study = NULL, ldThresh = 0.8,
                                 ldPop = "AMR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                 url = "http://archive.broadinstitute.org/mammals/haploreg/haploreg.php",
                                 timeout = 100000, encoding = "UTF-8", verbose = FALSE)
    AMR.results <- data.frame(AMR.results)
    print(paste0("AMR imputation resulted in # rows: ", length(AMR.results[, 1])))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AMR.results)
    
  }
  if (length(pop.sorted.tables[[current.disease]][["ASN"]]) > 0) {
    
    ASN <- pop.sorted.tables[[current.disease]][["ASN"]][!is.na(pop.sorted.tables[[current.disease]][["ASN"]])]
    ASN <- paste0("rs", ASN)
    length(ASN)
    print(length(unique(ASN)))
    write.table(unique(ASN), file=paste0(current.disease, ".ASN.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    ASN.results <- queryHaploreg(query = unique(ASN), file = NULL, study = NULL, ldThresh = 0.8,
                                 ldPop = "ASN", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                 url = "http://archive.broadinstitute.org/mammals/haploreg/haploreg.php",
                                 timeout = 100000, encoding = "UTF-8", verbose = FALSE)
    ASN.results <- data.frame(ASN.results)
    print(paste0("ASN imputation resulted in # rows: ", length(ASN.results[, 1])))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], ASN.results)
    
  }
  if (length(pop.sorted.tables[[current.disease]][["EUR"]]) > 0) {
    
    EUR <- pop.sorted.tables[[current.disease]][["EUR"]][!is.na(pop.sorted.tables[[current.disease]][["EUR"]])]
    EUR <- paste0("rs", EUR)
    length(EUR)
    print(length(unique(EUR)))
    write.table(unique(EUR), file=paste0(current.disease, ".EUR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    EUR.results <- queryHaploreg(query = unique(EUR), file = NULL, study = NULL, ldThresh = 0.8,
                                 ldPop = "EUR", epi = "vanilla", cons = "siphy", genetypes = "gencode",
                                 url = "http://archive.broadinstitute.org/mammals/haploreg/haploreg.php",
                                 timeout = 100000, encoding = "UTF-8", verbose = FALSE)
    EUR.results <- data.frame(EUR.results)
    print(paste0("EUR imputation resulted in # rows: ", length(EUR.results[, 1])))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], EUR.results)
    
  }
  
  #
  print(paste0("Total imputed rows: ", length(imputation.results[[current.disease]][, 1])))
  
  #
  save(imputation.results, file="/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/gwas_catalog/imputation.results.Rdata")

}


# Process HaploReg results ------------------------------------------------


#
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

# Remove exonic SNPs.
exons.file <- "/Users/michael/Box Sync/MS_analysis/brain_pchic/supplementary_files/gencode_v19_exons.bed"
exons <- read.table(exons.file, sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
colnames(exons) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")    

#
load(file="/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/gwas_catalog/imputation.results.Rdata")
diseases <- c("AD", "ADHD", "ASD", "BD", "EP", "MDD", "PD", "SCZ", "T1D")
setwd("/Users/Michael/Box Sync/MS_analysis/brain_pchic/other/gwas_catalog/processed")
merged.translated <- c()
for (i in 1:length(diseases)) {
  
  #
  current.disease <- diseases[i]
  print(current.disease)
  print(length(imputation.results[[current.disease]][, 1]))
  positions <- as.integer(as.character(imputation.results[[current.disease]][, 2]))
  no.coord <- imputation.results[[current.disease]][is.na(positions), ]
  print(length(no.coord[, 1]))
  
  #
  imputation.results[[current.disease]] <- imputation.results[[current.disease]][!(is.na(positions)), ]
  positions <- as.integer(as.character(imputation.results[[current.disease]][, 2]))
  print(length(imputation.results[[current.disease]][, 1]))

  # Lift over coordinates.
  gr <- GRanges(seqnames=as.character(imputation.results[[current.disease]][, 1]),
          ranges=IRanges(positions-1, end=positions,
                         names=as.character(imputation.results[[current.disease]][, 6])),
          strand=rep("*", length(imputation.results[[current.disease]][, 1])),
          rsID=as.character(imputation.results[[current.disease]][, 6]),
          ref=as.character(imputation.results[[current.disease]][, 7]),
          alt=as.character(imputation.results[[current.disease]][, 8]),
          query_snp=as.character(imputation.results[[current.disease]][, 33]),
          is_query_snp=as.character(imputation.results[[current.disease]][, 5]))
  
  #
  seqlevelsStyle(gr) = "UCSC"  # necessary
  hg19.coords <- unlist(liftOver(gr, ch))
  genome(hg19.coords) = "hg19"
  hg19.coords <- data.frame(hg19.coords)

  #
  print(length(hg19.coords[, 1]))
  
  #
  translated <- hg19.coords[, c(1, 2, 3, 6, 7, 8, 9, 10)]
  translated <- unique(translated[, c(1:6)])
  print(length(translated[, 1]))
  print(length(unique(translated$rsID)))
  
  #
  if (current.disease == "SCZ") {
    
    translated$rsID[duplicated(translated$rsID)]
    translated[translated$rsID == "rs11811290", ]
    translated <- translated[!(translated$rsID == "rs11811290"), ]
    translated <- rbind(translated, c("chr1", 147114284, 147114286, "rs11811290", "A", "T"))
    
  }
  
  #
  translated <- cbind(translated[, c(1:6)], query_snp_rsid=rep("", length(translated[, 1])), is_query_snp=rep(FALSE, length(translated[, 1])))
  translated[, 7] <- as.character(translated[, 7])
  translated[, 8] <- as.logical(translated[, 8])
  query_snps <- unique(hg19.coords$rsID[hg19.coords$is_query_snp == 1])
  for (j in 1:length(translated[, 1])) {
    
    if (translated[j, 4] %in% query_snps)
      translated[j, 8] <- TRUE
    
    query.snps <- paste0(unique(as.character(imputation.results[[current.disease]]$query_snp_rsid[imputation.results[[current.disease]]$rsID == translated[j, 4]])), collapse=",")
    translated[j, 7] <- query.snps

  }
  
  # Fix zero width intervals.
  translated[, 2] <- as.numeric(translated[, 2])
  translated[, 3] <- as.numeric(translated[, 3])
  zero <- translated$start == translated$end
  print("Zero width intervals:")
  translated[zero, 2] <- translated[zero, 2] - 1
  zero_check <- translated$start == translated$end

  #
  print(paste0("unique tag SNPs:", length(unique(translated[translated[, 8] == TRUE, 4]))))
  print(paste0("unique imputed SNPs:", length(unique(translated[, 4]))))
  
  #
  intersect <- bedTools.2in(bed1=translated, bed2=exons, opt.string="-c")
  print(table(intersect[, 9] > 0))
  print(1-table(intersect[, 9] > 0)[["TRUE"]]/table(intersect[, 9] > 0)[["FALSE"]])
  translated.noncoding <- translated[intersect[, 9] == 0, ]
  dim(translated.noncoding)
  dim(unique(translated.noncoding[, 1:4]))
  
  #
  print(paste0("unique tag SNPs:", length(unique(translated.noncoding[translated.noncoding[, 8] == TRUE, 4]))))
  print(paste0("unique imputed SNPs:", length(unique(translated.noncoding[, 4]))))
  
  #
  write.table(translated.noncoding, file=paste0(current.disease, ".feature.noncoding.snps.bed"), quote=F, sep="\t", row.names=F, col.names=F)
  write.table(sortBed(unique(translated.noncoding[, 1:4])), quote=F, sep="\t", row.names=F, col.names=F, file=paste0(current.disease, ".display.noncoding.snps.bed"))
  merged.translated <- rbind(merged.translated, translated.noncoding[, 1:4])

}

#
write.table(sortBed(unique(merged.translated)), quote=F, sep="\t", row.names=F, col.names=F, file=paste0("all.display.noncoding.snps.bed"))

