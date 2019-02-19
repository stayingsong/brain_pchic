#
# gwas_mining.R
# Author: Michael Song
# Last modified: 2019-01-21
# This script mines SNPs from the GWAS Catalog, then imputes them and formats them for downstream analysis.
#


# Load dependencies.
require(haploR)
require(GenomicRanges)
require(rtracklayer)
require(rsnps)


# Set up environment and load settings ------------------------------------


# Clear workspace before running script.
rm(list=ls())

# Turn off scientific notation for writing output.
options(scipen=999)

# Convert warnings into errors.
# options(warn=2)

# Process CLI arguments (or override if running interactively).
cli.args = commandArgs(trailingOnly=TRUE)
print(cli.args[1]) # Set the home directory for the pipeline.
if ((exists("cli.args") == FALSE) | (length(cli.args) == 0)) {
  
  cli.args <- c()
  cli.args[1] <- "/Users/Michael/Box Sync/MS_general/analysis/brain_pchic_revision"
  
}

# Set home directory.
home.dir <- cli.args[1]

# Load auxiliary functions.
scripts.dir <- paste0(home.dir, "/scripts")
source(paste0(scripts.dir, "/utilities.R"))

# Define the significance cutoff for filtering SNPs from the GWAS Catalog.
sig.cutoff <- 10^-6

# Define the LD threshold for performing imputation.
ld.threshold <- 0.8

# Define the populations to be used for imputation.
populations <- c("AFR", "AMR", "ASN", "EUR")

# List of diseases to mine SNPs for, as well ass the dates where the associations and study metadata were downloaded from the GWAS Catalog.
controls <- c("AMD", "BH", "BMI", "RA", "T1D", "T2D")
diseases <- c("AD", "ADHD", "ALS", "ASD", "BD", "EP", "FTD", "MP", "PD", "SCZ", "UD", controls)
dates <- c(rep("2018-10-30", 2), "2018-12-02", rep("2018-10-30", 4), "2018-12-10", rep("2018-10-30", 4), rep("2018-12-09", 2), rep("2018-10-30", 3))


# Load resources ----------------------------------------------------------


# Read in tables for translating rsids.
translation.psycharray <- read.table(paste0(home.dir, "/resources/InfiniumPsychArray-24v1-2_A1_b144_rsids.txt"), 
                                     sep="\t", header=T, stringsAsFactors=F, colClasses=c("character", "character"))
translation.custom <- read.table(paste0(home.dir, "/resources/custom_rsids.txt"), 
                                 sep="\t", header=T, stringsAsFactors=F, colClasses=c("character", "character"))
translation.all <- rbind(translation.psycharray, translation.custom)

# Read in exonic intervals.
exons.file <- "gencode/gencode_v19_exons.bed"
exons <- read.table(paste0(home.dir, "/resources/", exons.file), sep="\t", header=F, stringsAsFactors=F, colClasses=c("character", "integer", "integer", "character", "numeric", "character"))
colnames(exons) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")     

# Prepare chain files for lifting over to hg19.
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)

# Read in the RsMergeArch table for updating rsids.
rs.merge.file <- paste0(home.dir, "/resources/RsMergeArch.bcp")
rs.merge <- read.table(rs.merge.file, sep="\t", header=F, stringsAsFactors=F)
rs.merge[, 1] <- paste0("rs", rs.merge[, 1])
rs.merge[, 2] <- paste0("rs", rs.merge[, 2])


# SNP imputation ----------------------------------------------------------


# Read in study and association data for each disease and perform imputation by population(s).
studies <- list()
associations <- list()
snps.by.population <- list()
imputation.results <- list()
for (i in 1:length(diseases)) {
  
  # Set the current disease.
  current.disease <- diseases[i]
  
  # Read in study and association data for each disease.
  print(paste0("Imputing SNPs for: ", current.disease))
  studies[[current.disease]] <- read.table(paste0(home.dir, "/other/gwas_mining/raw/", current.disease, "_", dates[i], "_studies.txt"), 
                                           sep="\t", header=T, stringsAsFactors=F, quote="",
                                           colClasses=c(rep("character", 9), "integer", "character", "integer"))
  associations[[current.disease]] <- read.table(paste0(home.dir, "/other/gwas_mining/raw/", current.disease, "_", dates[i], "_withChildTraits.tsv"), 
                                                sep="\t", header=T, stringsAsFactors=F, fill=T, quote="")
  
  # Print #s of studies and associations.
  print(paste0("Number of studies: ", length(studies[[current.disease]][, 1])))
  print(paste0("Number of associations: ", sum(studies[[current.disease]]$Association.count)))
  print("Studies with associations that also have metadata available:")
  print(table(unique(associations[[current.disease]]$STUDY.ACCESSION) %in% studies[[current.disease]]$Study.accession))
  # print("Studies with associations without metadata available:")
  # print(unique(associations[[current.disease]]$STUDY.ACCESSION[!(associations[[current.disease]]$STUDY.ACCESSION %in% studies[[current.disease]]$Study.accession)]))
  
  # Filter out studies unrelated to the current disease.
  # studies[[current.disease]] <- studies[[current.disease]][studies[[current.disease]]$Include == 1, ]
  # print(paste0("Number of studies used for imputation: ", length(unique(studies[[current.disease]]$Study.accession))))
  
  # Filter SNPs according to the significance cutoff.
  snps.filtered <- associations[[current.disease]][as.numeric(associations[[current.disease]]$P.VALUE) < sig.cutoff, ]
  print(paste0("Number of SNPs passing significance cutoff of ", sig.cutoff, ": ", length(snps.filtered[, 1])))
  
  # Count how many SNPs do not have a current rsid ('SNP_ID_CURRENT' field is NA).
  snps.current.ids <- snps.filtered$SNP_ID_CURRENT
  snps.na.ids <- which(is.na(snps.filtered$SNP_ID_CURRENT))
  print(paste0("Number of SNPs without a current rsid ('SNP_ID_CURRENT' field is NA): ", length(snps.na.ids)))
  
  # Detect instances where the rsid in the 'SNPS' field doesn't agree with the 'SNP_ID_CURRENT' field. Print out all such instances.
  # 'SNP_ID_CURRENT' takes precedence over 'SNPS', so use the former as long as it is available.
  snps.reported <- snps.filtered$SNPS
  num.disagreements <- length(which(paste0("rs", snps.current.ids) != snps.reported))
  print(paste0("Number of instances where the rsid in the 'SNPS' field doesn't agree with the 'SNP_ID_CURRENT' field: ", num.disagreements))
  if (num.disagreements > 0) {
    
    which.disagreements <- !(paste0("rs", snps.current.ids) == snps.reported)
    mismatches <- cbind(snps.current.ids[which.disagreements], snps.reported[which.disagreements], is.na(snps.filtered$SNP_ID_CURRENT)[which.disagreements])
    colnames(mismatches) <- c("SNP_ID_CURRENT", "SNPS", "IS.NA")
    # print(mismatches)
    
  }

  # Try to use the information in the 'SNPS' field to resolve instances where the 'SNP_CURRENT_ID' is NA. Keep track of which instances were resolved.
  resolved <- rep(FALSE, length(snps.na.ids))
  if (length(snps.na.ids) > 0) {
    
    # Try to resolve each NA case one by one.
    for (j in 1:length(snps.na.ids)) {
      
      # Set current query to try to resolve.
      query <- snps.reported[snps.na.ids[j]]
      
      # First check if it is a delimited lists of multiple rsids by splitting via comma, then by semicolon.
      query.split <- unlist(strsplit(query, split=", "))
      if (length(query.split) == 1) {
        query.split <- unlist(strsplit(query, split="; "))
      }
      
      # If splitting resulted in multiple entries, process them--otherwise, turn to the translation table.
      if (length(query.split) > 1) {
        
        # Remove "rs" prefix for conformity, then update the 'SNP_ID_CURRENT' field (the "#" character is used as a delimiter in this script).
        query.split <- gsub("rs", "", query.split)
        snps.current.ids[snps.na.ids[j]] <- paste0(query.split, collapse="#")
        resolved[j] <- TRUE
        
      } else {
        
        # Try to translate through one of the translation tables.
        result <- unique(translation.all$RsID[translation.all$Name == query])
        if (length(result) == 1) {
          
          # Remove "rs" prefix for conformity, then update the 'SNP_ID_CURRENT' field.
          result <- gsub("rs", "", result)
          # results <- unlist(strsplit(result, split=";"))
          snps.current.ids[snps.na.ids[j]] <- result
          resolved[j] <- TRUE
          
        }
        
      }
      
    }

  }
  
  # Print unresolved and resolved cases.
  print(paste0("Unresolved cases (", length(which(!resolved)), "):"))
  print(data.frame(mismatches[mismatches[, 3] == "TRUE", ], stringsAsFactors=F)[!resolved, ])
  # print(paste0("Resolved cases (", length(which(resolved)), "):"))
  # print(cbind(as.data.frame(mismatches[mismatches[, 3] == "TRUE", ], stringsAsFactors=F)[resolved, ], NEW_SNP_ID_CURRENT=snps.current.ids[snps.na.ids][resolved]))

  # For each population, get a list of SNPs associated with that population.
  snps.by.population[[current.disease]] <- list()
  for (j in 1:length(populations)) {
    
    # For each study, look up what populations are associated with it. If that includes the current population, add all the SNPs from that study to the master list.
    snps.by.population[[current.disease]][[populations[j]]] <- c()
    for (k in 1:length(studies[[current.disease]][, 1])) {
      
      # See what populations are associated with the current study.
      current.populations <- unlist(strsplit(studies[[current.disease]]$Population[k], split="\\|"))
      
      # If those populationss include the current population, add all the SNPs from that study to the master list.
      if (populations[j] %in% current.populations) {
        
        # Retrieve the rsids for the current study.
        snps.study <- as.character(snps.current.ids[snps.filtered$STUDY.ACCESSION == studies[[current.disease]]$Study.accession[k]])
        for (snp in snps.study) {
          
          # Explode each SNP by the "#" delimiter and add them to the master list.
          exploded <- unlist(strsplit(snp, split="#"))
          if (length(exploded) >= 1)
            snps.by.population[[current.disease]][[populations[j]]] <- c(snps.by.population[[current.disease]][[populations[j]]], exploded)
          
        }
        
      }
      
    }
    
    # Remove NA values, add "rs" prefix to rsids, and retain only unique entries.
    if (!is.null(snps.by.population[[current.disease]][[populations[j]]])) {
      
      snps.by.population[[current.disease]][[populations[j]]] <- snps.by.population[[current.disease]][[populations[j]]][!is.na(snps.by.population[[current.disease]][[populations[j]]])]
      snps.by.population[[current.disease]][[populations[j]]] <- unique(snps.by.population[[current.disease]][[populations[j]]])
      snps.by.population[[current.disease]][[populations[j]]] <- paste0("rs", snps.by.population[[current.disease]][[populations[j]]])
      
    }
    
  }
  
  # Perform imputation for each population using HaploReg.
  imputation.results[[current.disease]] <- c()
  chunk.size <- 500
  if (length(snps.by.population[[current.disease]][["AFR"]]) > 0) {
    
    AFR <- snps.by.population[[current.disease]][["AFR"]]
    print(paste0("AFR unique SNPs: ", length(AFR)))
    write.table(AFR, file=paste0(home.dir, "/other/gwas_mining/sorted/", current.disease, ".AFR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AFR.results <- queryHaploreg(query=AFR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="AFR", epi="vanilla", cons="siphy", genetypes="gencode",
                                 url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
    AFR.results <- data.frame(AFR.results, stringsAsFactors=F)
    print(paste0("AFR imputation resulted in ", length(AFR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AFR.results)
    
  }
  if (length(snps.by.population[[current.disease]][["AMR"]]) > 0) {
    
    AMR <- snps.by.population[[current.disease]][["AMR"]]
    if (current.disease == "UD") {
      AMR <- AMR[-which(AMR == "rs748441912")]
      AMR <- AMR[-which(AMR == "rs782472239")]
    }
    print(paste0("AMR unique SNPs: ", length(AMR)))
    write.table(AMR, file=paste0(home.dir, "/other/gwas_mining/sorted/", current.disease, ".AMR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    AMR.results <- queryHaploreg(query=AMR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="AMR", epi="vanilla", cons="siphy", genetypes="gencode",
                                 url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
    AMR.results <- data.frame(AMR.results, stringsAsFactors=F)
    print(paste0("AMR imputation resulted in ", length(AMR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], AMR.results)
    
  }
  if (length(snps.by.population[[current.disease]][["ASN"]]) > 0) {
    
    ASN <- snps.by.population[[current.disease]][["ASN"]]
    print(paste0("ASN unique SNPs: ", length(ASN)))
    write.table(ASN, file=paste0(home.dir, "/other/gwas_mining/sorted/", current.disease, ".ASN.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    if (length(ASN) > chunk.size) {
      ASN.results <- c()
      for (i in 1:floor(length(ASN)/chunk.size)) {
        ASN.results.temp <- queryHaploreg(query=ASN[(1+(i-1)*chunk.size):(i*chunk.size)], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                          url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
        ASN.results <- rbind(ASN.results, data.frame(ASN.results.temp))
      }
      ASN.results.temp <- queryHaploreg(query=ASN[(1+floor(length(ASN)/chunk.size)*chunk.size):(length(ASN))], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                        url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      ASN.results <- rbind(ASN.results, data.frame(ASN.results.temp))
    } else {
      ASN.results <- queryHaploreg(query=ASN, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="ASN", epi="vanilla", cons="siphy", genetypes="gencode",
                                   url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      ASN.results <- data.frame(ASN.results, stringsAsFactors=F)
    }
    print(paste0("ASN imputation resulted in ", length(ASN.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], ASN.results)
    
  }
  if (length(snps.by.population[[current.disease]][["EUR"]]) > 0) {
    
    EUR <- snps.by.population[[current.disease]][["EUR"]]
    print(paste0("EUR unique SNPs: ", length(EUR)))
    write.table(EUR, file=paste0(home.dir, "/other/gwas_mining/sorted/", current.disease, ".EUR.txt"), sep="\t", row.names=F, col.names=F, quote=F)
    if (length(EUR) > chunk.size) {
      EUR.results <- c()
      for (i in 1:floor(length(EUR)/chunk.size)) {
        EUR.results.temp <- queryHaploreg(query=EUR[(1+(i-1)*chunk.size):(i*chunk.size)], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                     url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
        EUR.results <- rbind(EUR.results, data.frame(EUR.results.temp))
      }
      EUR.results.temp <- queryHaploreg(query=EUR[(1+floor(length(EUR)/chunk.size)*chunk.size):(length(EUR))], file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                        url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      EUR.results <- rbind(EUR.results, data.frame(EUR.results.temp))
    } else {
      EUR.results <- queryHaploreg(query=EUR, file=NULL, study=NULL, ldThresh=ld.threshold, ldPop="EUR", epi="vanilla", cons="siphy", genetypes="gencode",
                                   url="http://pubs.broadinstitute.org/mammals/haploreg/haploreg.php", timeout=500000, encoding="UTF-8", verbose=FALSE)
      EUR.results <- data.frame(EUR.results, stringsAsFactors=F)
    }
    print(paste0("EUR imputation resulted in ", length(EUR.results[, 1]), " entries."))
    imputation.results[[current.disease]] <- rbind(imputation.results[[current.disease]], EUR.results)
    
  }
  
  # Print total # of imputed entries.
  print(paste0("Total number of imputed entries: ", length(imputation.results[[current.disease]][, 1])))

}

# Save results.
save(studies, file=paste0(home.dir, "/other/gwas_mining/saved/studies.Rdata"))
save(associations, file=paste0(home.dir, "/other/gwas_mining/saved/associations.Rdata"))
save(snps.by.population, file=paste0(home.dir, "/other/gwas_mining/saved/snps.by.population.Rdata"))
save(imputation.results, file=paste0(home.dir, "/other/gwas_mining/saved/imputation.results.Rdata"))


# Process HaploReg results ------------------------------------------------


# Load results.
load(file=paste0(home.dir, "/other/gwas_mining/saved/imputation.results.Rdata"))

# Process imputed data.
processed.snps.noncoding <- list()
processed.snps.merged <- c()
processed.snps.all <- list()
chunk.size <- 200
for (i in 1:length(diseases)) {
  
  # Set the current disease and summarize imputation results.
  current.disease <- diseases[i]
  print(paste0("Processing SNPs for disease: ", current.disease))
  print(paste0("# of entries from imputation: ", length(imputation.results[[current.disease]][, 1])))
  print(paste0("# of unique rsids: ", length(unique(imputation.results[[current.disease]]$rsID))))
  
  # Count the # of unique SNPs w/ and w/o positional info.
  # print(table(which(imputation.results[[current.disease]][, 2] == "") %in% which(is.na(imputation.results[[current.disease]][, 1])))) # Should be TRUE.
  has.pos <- imputation.results[[current.disease]][!is.na(imputation.results[[current.disease]][, 1]), ]
  missing.pos <- imputation.results[[current.disease]][is.na(imputation.results[[current.disease]][, 1]), ]
  print(paste0("# of unique rsids with positional info: ", length(unique(has.pos$rsID))))
  print(paste0("# of unique rsids without positional info: ", length(unique(missing.pos$rsID))))
  # print(table(unique(missing.pos$rsID) %in% unique(has.pos$rsID))) # Should be FALSE.
  
  # Update rsids w/o positional info. First update any merged rsids then query NCBI's dbSNP for information.
  num.updated <- 0
  num.recovered <- 0
  original.rsids <- missing.pos$rsID
  for (j in 1:length(missing.pos$rsID)) {
    
    missing.pos$rsID[j] <- update.rsid(missing.pos$rsID[j], rs.merge)
    if (missing.pos$rsID[j] != original.rsids[j])
      num.updated <- num.updated + 1
    
  }
  if (length(missing.pos$rsID) > chunk.size) {
    
    results <- c()
    for (j in 1:floor(length(missing.pos$rsID)/chunk.size)) {
      results <- rbind(results, ncbi_snp_query(SNPs=missing.pos$rsID[(1 + (j - 1)*chunk.size):(j*chunk.size)]))
    }
    results <- rbind(results, ncbi_snp_query(SNPs=missing.pos$rsID[(1 + (ceiling(length(missing.pos$rsID)/chunk.size) - 1)*chunk.size):length(missing.pos$rsID)]))

  } else {
    
    results <- ncbi_snp_query(SNPs=missing.pos$rsID)
    
  }
  
  # Update the information for SNPs w/o positional info using the query results.
  unresolved.rsids <- c()
  for (j in 1:length(missing.pos[, 1])) {
    
    # There should be one unique match per entry.
    match <- unique(results[results$Marker == missing.pos$rsID[j], ])
    if (length(match[, 1]) == 1) {
    
      if (!is.na(match$Chromosome) & !is.na(match$BP)) {
        
        missing.pos$chr[j] <- match$Chromosome
        missing.pos$pos_hg38[j] <- match$BP
        num.recovered <- num.recovered + 1
        
      } else {
        
        # print(paste0(missing.pos$rsID[j], " not recovered"))
        unresolved.rsids <- c(unresolved.rsids, missing.pos$rsID[j])
        
      }
        
    } else if (length(match[, 1]) == 0) {
      
      print(paste0(missing.pos$rsID[j], " did not return any results"))
      unresolved.rsids <- c(unresolved.rsids, missing.pos$rsID[j])
      
    }
    
  }
  
  # Print how many rsids were updated and how many had their positional info recovered.
  print(paste0(num.updated, " of ", length(missing.pos[, 1]), " rsids updated."))
  print(paste0(num.recovered, " of ", length(missing.pos[, 1]), " rsids had their positional info recovered."))
  
  # Use the old rsids for now (we can add the new rsids back later anytime).
  mapping <- unique(cbind(original.rsids, missing.pos$rsID))
  if (!(length(original.rsids) == dim(missing.pos)[1]))
    print("Inconsistency between original and updated rsids")
  missing.pos$rsID <- original.rsids
  
  # Recombine the recovered rsids with the original rsids with positional info. Also print any rsids that were not recovered.
  combined.pos <- rbind(has.pos, missing.pos[!is.na(missing.pos$chr), ])
  print(paste0("# of unique rsids recovered: ", length(unique(missing.pos[!is.na(missing.pos$chr), ]))))
  print(paste0("# of unique rsids which could not be recovered: ", length(unique(missing.pos$rsID[is.na(missing.pos$chr)]))))
  # print(table(is.na(combined.pos[, 1]))) # Should be FALSE.
  print(paste0("Check: ", length(unique(unresolved.rsids))))
  print(unique(unresolved.rsids))
  
  # Show table of converted rsids in case one of the rsids to be recovered was updated.
  unresolved.mapping <- mapping[(mapping[, 1] %in% unresolved.rsids) | (mapping[, 2] %in% unresolved.rsids), ]
  print(table(unresolved.mapping[, 1] == unresolved.mapping[, 2]))
  print(unresolved.mapping[unresolved.mapping[, 1] != unresolved.mapping[, 2], ])

  # Create GenomicRanges object with rsids to be lifted over.
  gr <- GRanges(seqnames=as.character(combined.pos$chr),
          ranges=IRanges(as.numeric(combined.pos$pos_hg38) - 1, end=as.numeric(combined.pos$pos_hg38)),
          strand=rep("*", length(combined.pos[, 1])),
          rsid=as.character(combined.pos$rsID),
          ref=as.character(combined.pos$ref),
          alt=as.character(combined.pos$alt),
          query_snp=as.character(combined.pos$query_snp_rsid),
          is_query_snp=as.character(combined.pos$is_query_snp))
  
  # Lift over coordinates from hg38 to hg19.
  seqlevelsStyle(gr) = "UCSC"  # necessary
  hg19.coords <- unlist(liftOver(gr, ch))
  genome(hg19.coords) = "hg19"
  hg19.coords <- data.frame(hg19.coords, stringsAsFactors=F)
  hg19.coords <- hg19.coords[, c(1, 2, 3, 6, 7, 8, 9, 10)]
  print(paste0("# of unique positions to lift over: ", length(unique(combined.pos[, 1:2])[, 1])))
  print(paste0("# of unique positions successfully lifted over: ", length(unique(hg19.coords[, 1:3])[, 1])))
  print(table(hg19.coords$query_snp %in% hg19.coords$rsid))

  # Determine which entries were not lifted over correctly.
  print(table(combined.pos$rsID %in% hg19.coords$rsid))
  failed.rsids <- combined.pos$rsID[!(combined.pos$rsID %in% hg19.coords$rsid)]
  failed.pos <- combined.pos[combined.pos$rsID %in% failed.rsids, c(1:10, 33)]
  print(failed.pos)
  missing.queries <- unique(hg19.coords$query_snp[!(hg19.coords$query_snp %in% hg19.coords$rsid)])
  print(missing.queries)

  # Manually recover rsids that lack positional info in hg38 but are reported in hg19. Also fix cases where the same rsid has multiple positions.
  {
    
    if (current.disease == "AD") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr19", "45219014", "45219015", "rs3201345", "C", "T", "rs79638902", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr6", "107307825", "107307826", "rs9486479", "G", "A", "rs140633572", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr6", "107307830", "107307831", "rs9486480", "C", "T", "rs140633572", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr19", "45221172", "45221173", "rs55710026", "G", "T", "rs79638902", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr19", "45217788", "45217789", "rs78358160", "G", "C", "rs79638902", "0"))
    } else if (current.disease == "ADHD") {
    } else if (current.disease == "ALS") {
    } else if (current.disease == "ASD") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr22", "42537126", "42537127", "rs3020736", "A", "G", "rs3020736", "1")) # Query SNP.
      # dbSNP:
      hg19.coords <- rbind(hg19.coords, c("chr1", "163646390", "163646391", "rs5778340", "TA", "T", "rs7521492", "0"))
    } else if (current.disease == "BD") {
    } else if (current.disease == "EP") {
    } else if (current.disease == "FTD") {
    } else if (current.disease == "MP") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr17", "34869154", "34869155", "rs8882", "G", "A", "rs3736166", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "34869154", "34869155", "rs8882", "G", "A", "rs8070260", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "35000571", "35000572", "rs211985", "A", "G", "rs211988", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42528857", "42528858", "rs28439297", "C", "T", "rs133383", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42528857", "42528858", "rs28439297", "C", "T", "rs2142694", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42528857", "42528858", "rs28439297", "C", "T", "rs2743461", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42528857", "42528858", "rs28439297", "C", "T", "rs5758605", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42528857", "42528858", "rs28439297", "C", "T", "rs2743465", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr6", "119204628", "119204629", "rs36166681", "A", "G", "rs9401101", "0"))
      # dbSNP:
      hg19.coords <- hg19.coords[-which(hg19.coords$rsid == "rs5884966"), ]
      hg19.coords <- rbind(hg19.coords, c("chr7", "75145265", "75145266", "rs5884966", "TA", "T", "rs1167827", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr2", "172833138", "172833139", "rs13004237", "A", "G", "rs2165934", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150158000", "150158001", "rs10708380", "TG", "T", "rs55802315", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150092238", "150092239", "rs11310019", "TA", "T", "rs7554367", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150092238", "150092239", "rs11310019", "TA", "T", "rs55802315", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44126472", "44126473", "rs11318581", "AG", "A,AC", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44126472", "44126473", "rs11318581", "AG", "A,AC", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44104409", "44104410", "rs35654420", "TCTC", "T,TC,TGAG", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44104409", "44104410", "rs35654420", "TCTC", "T,TC,TGAG", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr6", "119204624", "119204625", "rs36139342", "A", "AC", "rs9374757", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs4327090", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs17426174", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs17563986", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs62057107", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44201357", "44201358", "rs56359268", "TG", "T,TC", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44201357", "44201358", "rs56359268", "TG", "T,TC", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs4327090", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44058698", "44058699", "rs66635171", "TAGG", "T,TCCT", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44058698", "44058699", "rs66635171", "TAGG", "T,TCCT", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43849571", "43849572", "rs66940642", "AT", "A,AA", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43849571", "43849572", "rs66940642", "AT", "A,AA", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs4327090", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44094132", "44094133", "rs67084886", "CA", "C,CT", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44094132", "44094133", "rs67084886", "CA", "C,CT", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43994507", "43994508", "rs67180752", "GAGA", "G,GTCT", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43994507", "43994508", "rs67180752", "GAGA", "G,GTCT", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs4327090", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44028468", "44028469", "rs67600861", "CA", "C,CT", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44028468", "44028469", "rs67600861", "CA", "C,CT", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44177303", "44177304", "rs67668514", "TA", "T,TT", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44177303", "44177304", "rs67668514", "TA", "T,TT", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44125276", "44125277", "rs67892405", "CT", "C,CA", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44125276", "44125277", "rs67892405", "CT", "C,CA", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43733706", "43733707", "rs68054598", "GC", "G,GG", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44061607", "44061608", "rs71910877", "GT", "G,GA", "rs242559", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44061607", "44061608", "rs71910877", "GT", "G,GA", "rs8080583", ""))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44161090", "44161091", "rs71927936", "GCTCCCCTGGTAAGTCCTAAA", "G,GTTTAGGACTTACCAGGGGAG", "rs242559", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44161090", "44161091", "rs71927936", "GCTCCCCTGGTAAGTCCTAAA", "G,GTTTAGGACTTACCAGGGGAG", "rs8080583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr3", "85588062", "85588063", "rs72250906", "A", "AAA,AAC", "rs17518584", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs4327090", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr3", "107835258", "107835259", "rs78660433", "T", "TA", "rs2117762", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs4327090", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs17563986", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43778314", "43778315", "rs112067021", "CA", "C,CT", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs17426174", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs62057107", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "145598991", "145598992", "rs138163379", "TTTTAA", "T", "rs61816194", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "35766562", "35766563", "rs140611363", "CA", "C", "rs1016678", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr6", "163852838", "163852839", "rs141016793", "A", "AA,AC", "rs1737329", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42531063", "42531064", "rs141606817", "AAAC", "A", "rs133383", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42531063", "42531064", "rs141606817", "AAAC", "A", "rs2142694", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42531063", "42531064", "rs141606817", "AAAC", "A", "rs2743461", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42531063", "42531064", "rs141606817", "AAAC", "A", "rs2743465", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42531063", "42531064", "rs141606817", "AAAC", "A", "rs5758605", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150159379", "150159380", "rs142725249", "CCT", "C", "rs55802315", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr2", "172833619", "172833620", "rs147918589", "AAAAC", "A", "rs2165934", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150157999", "150158000", "rs201783277", "TTG", "T", "rs55802315", "0"))
    } else if (current.disease == "PD") {
      # dbSNP:
      hg19.coords <- rbind(hg19.coords, c("chr17", "44126472", "44126473", "rs11318581", "AG", "A,AC", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44126472", "44126473", "rs11318581", "AG", "A,AC", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43796395", "43796396", "rs35343144", "TG", "T,TC", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44104409", "44104410", "rs35654420", "TCTC", "T,TC,TGAG", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44104409", "44104410", "rs35654420", "TCTC", "T,TC,TGAG", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43898361", "43898362", "rs45502095", "AAGG", "A,AAAG,ACCT", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43707837", "43707838", "rs56131117", "TG", "T,TC,TGTT", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43707837", "43707838", "rs56131117", "TG", "T,TC,TGTT", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43707837", "43707838", "rs56131117", "TG", "T,TC,TGTT", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44201357", "44201358", "rs56359268", "TG", "T,TC", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44201357", "44201358", "rs56359268", "TG", "T,TC", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43908475", "43908476", "rs56396707", "GT", "G,GA", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44058698", "44058699", "rs66635171", "TAGG", "T,TCCT", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44058698", "44058699", "rs66635171", "TAGG", "T,TCCT", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43819596", "43819597", "rs66857688", "AC", "A,AG", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43849571", "43849572", "rs66940642", "AT", "A,AA", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43954981", "43954982", "rs67040965", "TTCTCTTCCTCA", "T,TTGAGGAAGAGA", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44094132", "44094133", "rs67084886", "CA", "C,CT", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44094132", "44094133", "rs67084886", "CA", "C,CT", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43994507", "43994508", "rs67180752", "GAGA", "G,GTCT", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43994507", "43994508", "rs67180752", "GAGA", "G,GTCT", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43897025", "43897026", "rs67584582", "TGGAG", "T,TCTCC", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44028468", "44028469", "rs67600861", "CA", "C,CT", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44028468", "44028469", "rs67600861", "CA", "C,CT", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44177303", "44177304", "rs67668514", "TA", "T,TT", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44177303", "44177304", "rs67668514", "TA", "T,TT", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44125276", "44125277", "rs67892405", "CT", "C,CA", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44125276", "44125277", "rs67892405", "CT", "C,CA", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43733706", "43733707", "rs68054598", "GC", "G,GG", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43733706", "43733707", "rs68054598", "GC", "G,GG", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43733706", "43733707", "rs68054598", "GC", "G,GG", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44061607", "44061608", "rs71910877", "GT", "G,GA", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44061607", "44061608", "rs71910877", "GT", "G,GA", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44161090", "44161091", "rs71927936", "GCTCCCCTGGTAAGTCCTAAA", "G,GTTTAGGACTTACCAGGGGAG", "rs8070723", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44161090", "44161091", "rs71927936", "GCTCCCCTGGTAAGTCCTAAA", "G,GTTTAGGACTTACCAGGGGAG", "rs17577094", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43913359", "43913360", "rs72498449", "TGAGGGCAGGAGG", "T,TCCTCCTGCCCTC", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43814369", "43814370", "rs111581254", "AC", "A,AG", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43971205", "43971206", "rs111907488", "CTAATT", "C,CAATTA", "rs17649553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43778314", "43778315", "rs112067021", "CA", "C,CT", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43778314", "43778315", "rs112067021", "CA", "C,CT", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43777834", "43777835", "rs112694589", "TATCCTGCTTCCC", "T,TGGGAAGCAGGAT", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs365825", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs393152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs2942168", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs12185268", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "43827155", "43827156", "rs112702678", "GA", "G,GT", "rs17649553", "0"))
    } else if (current.disease == "SCZ") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr1", "147089346", "147089347", "rs2644555", "A", "T", "rs583583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "147089345", "147089346", "rs2644556", "G", "A,T", "rs583583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr22", "42537126", "42537127", "rs3020736", "A", "G", "rs3020736", "1"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143736136", "143736137", "rs35985277", "C", "T", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143737689", "143737690", "rs4567014", "G", "C", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143741380", "143741381", "rs4736362", "G", "A", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143728386", "143728387", "rs6990126", "T", "C", "rs11785400", "0"))
      # dbSNP:
      hg19.coords <- hg19.coords[-which(hg19.coords$rsid == "rs11811290"), ]
      hg19.coords <- hg19.coords[-which(hg19.coords$rsid == "rs35328046"), ]
      hg19.coords <- rbind(hg19.coords, c("chr8", "143743880", "143743881", "rs3839865", "TG", "T", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "163646390", "163646391", "rs5778340", "TA", "T", "rs7521492", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "150092238", "150092239", "rs11310019", "TA", "T", "rs140505938", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "147114284", "147114285", "rs11811290", "A", "T", "rs583583", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143724549", "143724550", "rs202112237", "TTG", "T", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143736999", "143737000", "rs35328046", "ACT", "A", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr2", "68445860", "68445861", "rs71907403", "T", "TGT,TTG,TTTT", "rs12052801", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143736301", "143736302", "rs74994924", "GAGA", "G", "rs11785400", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr8", "143726067", "143726068", "rs143712250", "CGAG", "C", "rs11785400", "0"))
    } else if (current.disease == "UD") {
    } else if (current.disease == "AMD") {
    } else if (current.disease == "BH") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr20", "34105974", "34105975", "rs224390", "G", "A", "rs2236164", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr20", "34109718", "34109718", "rs1810742", "A", "G", "rs2236164", "0"))
    } else if (current.disease == "BMI") {
      # dbSNP:
      hg19.coords <- hg19.coords[-which(hg19.coords$rsid == "rs5884966"), ]
      hg19.coords <- rbind(hg19.coords, c("chr7", "75145265", "75145266", "rs5884966", "TA", "T", "rs1167827", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "34932497", "34932498", "rs12952492", "C", "G", "rs12150665", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr1", "27967192", "27967193", "rs201326079", "AAC", "A", "rs2076463", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr12", "31438352", "31438353", "rs201532982", "CAT", "C", "rs80234489", "0"))
    } else if (current.disease == "RA") {
      # dbSNP:
      hg19.coords <- rbind(hg19.coords, c("chr1", "147205875", "147205876", "rs6593804", "C", "T", "rs6593803", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "123643423", "123643424", "rs11412415", "T", "TT,TTTAA", "rs881375", "0"))
    } else if (current.disease == "T1D") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr9", "136142999", "136143000", "rs543040", "A", "T", "rs657152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136143120", "136143121", "rs543968", "T", "C", "rs657152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136142216", "136142217", "rs644234", "T", "G", "rs657152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136149097", "136149098", "rs116552240", "T", "A", "rs657152", "0"))
      # dbSNP:
      hg19.coords <- rbind(hg19.coords, c("chr9", "136149094", "136149095", "rs8176646", "CA", "C", "rs657152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44126472", "44126473", "rs11318581", "AG", "A,AC", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr4", "123343276", "123343277", "rs34753082", "T", "TC,TTT", "rs2069762", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44104409", "44104410", "rs35654420", "TCTC", "T,TC,TGAG", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44201357", "44201358", "rs56359268", "TG", "T,TC", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136146067", "136146068", "rs59602812", "A", "AT", "rs657152", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44058698", "44058699", "rs66635171", "TAGG", "T,TCCT", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44094132", "44094133", "rs67084886", "CA", "C,CT", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44028468", "44028469", "rs67600861", "CA", "C,CT", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44177303", "44177304", "rs67668514", "TA", "T,TT", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44125276", "44125277", "rs67892405", "CT", "C,CA", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44061607", "44061608", "rs71910877", "GT", "G,GA", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr17", "44161090", "44161091", "rs71927936", "GCTCCCCTGGTAAGTCCTAAA", "G,GTTTAGGACTTACCAGGGGAG", "rs1052553", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136146447", "136146448", "rs139840563", "TAAGAC", "T", "rs657152", "0"))
    } else if (current.disease == "T2D") {
      # Liftover:
      hg19.coords <- rbind(hg19.coords, c("chr9", "136142999", "136143000", "rs543040", "A", "T", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136143120", "136143121", "rs543968", "T", "C", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136142216", "136142217", "rs644234", "T", "G", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136149097", "136149098", "rs116552240", "T", "A", "rs505922", "0"))
      # dbSNP:
      hg19.coords <- rbind(hg19.coords, c("chr9", "136149094", "136149095", "rs8176646", "CA", "C", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136146067", "136146068", "rs59602812", "A", "AT", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136146447", "136146448", "rs139840563", "TAAGAC", "T", "rs505922", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr9", "136149708", "136149709", "rs146518861", "AC", "A", "rs635634", "0"))
      hg19.coords <- rbind(hg19.coords, c("chr12", "31438352", "31438353", "rs201532982", "CAT", "C", "rs147538848", "0"))
    } 
    
  }
  
  # Check liftover success.
  table(combined.pos$rsID %in% hg19.coords$rsid)
  
  # Determine the difference between the # of unique rsids and the # of unique positions.
  print(paste0("Difference between the # of unique rsids and the # of unique positions: ", dim(unique(hg19.coords[, 1:6]))[1] - dim(unique(hg19.coords[, 1:3]))[1]))
  
  # Collapse entries by rsid, ref, alt, and query_snp fields so that each entry has a unique position.
  unique.coords <- hg19.coords %>% group_by_(.dots=c("seqnames", "start", "end")) %>% 
    summarize(rsid=paste0(unique(rsid), collapse=","), ref=paste0(ref, collapse=","), alt=paste0(alt, collapse=","), query_snp=paste0(unique(query_snp), collapse=","))
  unique.coords <- data.frame(unique.coords, stringsAsFactors=F)
  num.pos <- length(unique.coords[, 1])
  for (j in 1:num.pos) {
    
    unique.coords$ref[j] <- paste0(unique(unlist(strsplit(unique.coords$ref[j], split=","))), collapse=",")
    unique.coords$alt[j] <- paste0(unique(unlist(strsplit(unique.coords$alt[j], split=","))), collapse=",")
    
  }
  
  # Print entries where two rsids have the same coordinates, and the reconciled entry under the new approach.
  unique.coords.old <- unique(hg19.coords[, 1:6])
  duplicate.coords <- unique.coords.old[duplicated(unique.coords.old[, 1:3]), 1:3]
  if (length(duplicate.coords[, 1]) > 0) {

    for (j in 1:length(duplicate.coords[, 1])) {

      duplicate.indices <- which((unique.coords.old[, 1] == duplicate.coords[j, 1]) & (unique.coords.old[, 2] == duplicate.coords[j, 2]) & (unique.coords.old[, 3] == duplicate.coords[j, 3]))
      # print(unique.coords.old[duplicate.indices, ])

    }
    
  }
  unique.coords[grep(",", unique.coords$rsid), ]
  
  # Confirm that the #s of unique rsids and unique positions agree. See above comment.
  print(length(unique(unique.coords[, 1:3])[, 1]) == length(unique(unique.coords$rsid))) # Should be TRUE.
  
  ### OLD APPROACH ###
  
  # Print out the positions with duplicate rsids. We used to resolve cases where multiple rsids had the same position, but now we are retaining these cases (just merge the last 3 fields).
  # unique.coords <- unique(hg19.coords[, 1:6])
  # duplicate.coords <- unique.coords[duplicated(unique.coords[, 1:3]), 1:3]
  # if (length(duplicate.coords[, 1]) > 0) {
  #   
  #   for (j in 1:length(duplicate.coords[, 1])) {
  #     
  #     duplicate.indices <- which((unique.coords[, 1] == duplicate.coords[j, 1]) & (unique.coords[, 2] == duplicate.coords[j, 2]) & (unique.coords[, 3] == duplicate.coords[j, 3]))
  #     print(unique.coords[duplicate.indices, ])
  #     duplicate.saved <- unique.coords[duplicate.indices, ]
  #     duplicate.merged <- duplicate.saved %>% group_by_(.dots=c("seqnames", "start", "end")) %>% 
  #       summarize(rsid=paste0(rsid, collapse=";"), ref=paste0(ref, collapse=";"), alt=paste0(alt, collapse=";"))
  #     duplicate.merged <- data.frame(duplicate.merged, stringsAsFactors=F)
  #     unique.coords <- unique.coords[-duplicate.indices, ]
  #     unique.coords <- rbind(unique.coords, duplicate.merged)
  #     print(duplicate.merged)
  #     
  #   }
  #   
  # }  
  
  # Confirm that the #s of unique rsids and unique positions agree. See above comment.
  # print(length(unique(unique.coords[, 1:3])[, 1]) == length(unique(unique.coords$rsid))) # Should be TRUE.
  # num.pos <- length(unique.coords[, 1])
  # print(paste0("# of processed unique rsids: ", num.pos))
  
  # Add back query SNP information to each entry.
  # processed.snps <- cbind(unique.coords, query_snp_rsid=rep("", num.pos), is_query_snp=rep(FALSE, num.pos))
  # processed.snps$seqnames <- as.character(processed.snps$seqnames)
  # processed.snps$start <- as.integer(processed.snps$start)
  # processed.snps$end <- as.integer(processed.snps$end)
  # processed.snps$rsid <- as.character(processed.snps$rsid)
  # processed.snps$ref <- as.character(processed.snps$ref)
  # processed.snps$alt <- as.character(processed.snps$alt)
  # processed.snps$query_snp_rsid <- as.character(processed.snps$query_snp_rsid)
  # processed.snps$is_query_snp <- as.logical(processed.snps$is_query_snp)
  # for (j in 1:num.pos) {
  #   
  #   # For each SNP, look up its query SNP info and summarize it in the processed results table.
  #   selector <- hg19.coords$rsid == processed.snps$rsid[j]
  #   processed.snps$query_snp_rsid[j] <- paste0(unique(hg19.coords$query_snp[selector]), collapse=",")
  #   if (processed.snps$rsid[j] %in% unique(hg19.coords$query_snp[selector]))
  #     processed.snps$is_query_snp[j] <- TRUE
  #   
  # }
  
  ### OLD APPROACH ###
  
  # Add logical field determining whether entry is a query SNP.
  processed.snps <- cbind(unique.coords, is_query_snp=rep(FALSE, num.pos))
  colnames(processed.snps)[7] <- "query_snp_rsid"
  processed.snps$seqnames <- as.character(processed.snps$seqnames)
  processed.snps$seqnames <- as.character(processed.snps$seqnames)
  processed.snps$start <- as.integer(processed.snps$start)
  processed.snps$end <- as.integer(processed.snps$end)
  processed.snps$rsid <- as.character(processed.snps$rsid)
  processed.snps$ref <- as.character(processed.snps$ref)
  processed.snps$alt <- as.character(processed.snps$alt)
  processed.snps$query_snp_rsid <- as.character(processed.snps$query_snp_rsid)
  processed.snps$is_query_snp <- as.logical(processed.snps$is_query_snp)
  for (j in 1:num.pos) {
    
    rsids <- unlist(strsplit(processed.snps$rsid[j], split=","))
    query.rsids <- unlist(strsplit(processed.snps$query_snp_rsid[j], split=","))
    if (length(rsids) == 1) {
      
      if (rsids %in% query.rsids)
        processed.snps$is_query_snp[j] <- TRUE
      
    } else {
      
      print(processed.snps[j, ])
      if (any(rsids %in% query.rsids)) {
        
        print("Query SNP detected.")
        processed.snps$is_query_snp[j] <- TRUE
        
      }
      
    }
    
  }
  
  # Set each SNP to have a width of at least 1.
  interval.width <- processed.snps$end - processed.snps$start
  processed.snps$start[interval.width == 0] <- processed.snps$start[interval.width == 0] - 1

  # Determine which SNPs intersect exonic intervals and which are noncoding.
  intersect <- bedTools.2in(bed1=processed.snps, bed2=exons, opt.string="-c")
  processed.snps.noncoding[[current.disease]] <- processed.snps[intersect[, 9] == 0, ]

  # Store all processed SNPs for reference later.
  processed.snps.all[[current.disease]] <- processed.snps

  # Print final #s of processed SNPs.
  print(paste0("# of unique tag SNPs: ", length(unique(processed.snps[processed.snps[, 8] == TRUE, 4]))))
  print(paste0("# of unique imputed SNPs: ", length(unique(processed.snps[, 4]))))
  print(paste0("# of unique noncoding imputed SNPs: ", length(processed.snps.noncoding[[current.disease]][, 1])))

  # Add "exon_overlap" column as well as any redundant SNPs (rsids at the same position).
  intersect <- bedTools.2in(bed1=processed.snps.all[[current.disease]][, 1:8], bed2=exons, opt.string="-c")
  processed.snps.all[[current.disease]] <- cbind(processed.snps.all[[current.disease]][, 1:8], intersect[, 9] > 0)
  if (!(current.disease %in% controls))
    processed.snps.merged <- rbind(processed.snps.merged, processed.snps.all[[current.disease]][, 1:4])

  # Write processed SNPs to files for downstream analysis.
  write.table(processed.snps.all[[current.disease]], file=paste0(home.dir, "/other/gwas_mining/features/", current.disease, ".complete.snps.features.bed"),
              sep="\t", row.names=F, col.names=F, quote=F)
  # write.table(processed.snps.noncoding[[current.disease]], file=paste0(home.dir, "/other/gwas_mining/features/", current.disease, ".noncoding.snps.features.bed"),
  #             sep="\t", row.names=F, col.names=F, quote=F)
  write.table(sortBed(unique(processed.snps.all[[current.disease]][, 1:4])), file=paste0(home.dir, "/other/gwas_mining/display/", current.disease, ".complete.snps.display.bed"),
              sep="\t", row.names=F, col.names=F, quote=F)
  # write.table(sortBed(unique(processed.snps.noncoding[[current.disease]][, 1:4])), file=paste0(home.dir, "/other/gwas_mining/display/", current.disease, ".noncoding.snps.display.bed"),
  #             sep="\t", row.names=F, col.names=F, quote=F)
  
  # Final check to make sure all query SNPs are also entries.
  print(table(hg19.coords$query_snp %in% hg19.coords$rsid))
  
}

# Write processed SNPs from all non-control diseases for display at the same time.
write.table(sortBed(unique(processed.snps.merged)), file=paste0(home.dir, "/other/gwas_mining/display/all.complete.snps.display.bed"), sep="\t", row.names=F, col.names=F, quote=F)

# Save results.
save(processed.snps.all, file=paste0(home.dir, "/other/gwas_mining/saved/processed.snps.all.Rdata"))


# End ---------------------------------------------------------------------

