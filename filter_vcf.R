#
# filter_vcf.R
# Author: Michael Song
# Last modified: 2018-12-23
# This script processes a VCF file for use with HiC-Pro.
#


# Filtering ---------------------------------------------------------------


# Read in raw VCF file.
vcf.file <- "/Users/michael/Box Sync/MS_general/analysis/brain_pchic_revision/resources/vcf/wtc_PASS_hg19.phased_v1.vcf"
vcf.data <- read.table(vcf.file, sep="\t", header=F, skip=60, stringsAsFactors=F)
colnames(vcf.data) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "WTC")

# Get the GT and GQ fields for each entry. Filter out entries with only a GT entry.
vcf.data <- vcf.data[vcf.data$FORMAT != "GT", ]
num.entries <- length(vcf.data[, 1])
format <- vcf.data$FORMAT
wtc <- vcf.data$WTC
GT <- rep(NA, num.entries)
GQ <- rep(NA, num.entries)
print(num.entries)
for (i in 1:num.entries) {
  
  # Determine which parts of the string contain the GT and GQ info.
  fields <- unlist(strsplit(format[i], split="\\:"))
  GT[i] <- unlist(strsplit(wtc[i], split="\\:"))[which(fields == "GT")]
  GQ[i] <- unlist(strsplit(wtc[i], split="\\:"))[which(fields == "GQ")]
  
}

# Define criteria for filtering VCF entries by GT and GQ info.
informative.GT <- c("0|1", "1|0")
min.GQ <- 20

# Print filtering statistics.
print(paste0("Total # of VCF entries: ", num.entries))
print("Filter flag results:")
print(table(vcf.data$FILTER))
print("Formats distribution:")
print(table(vcf.data$FORMAT))
print("GTs distribution:")
print(table(GT))
print("# GQs passing filter:")
print(table(GQ > min.GQ))
print("# informative GTs:")
print(table(GT %in% informative.GT))
print("Total usable entries:")
print(table((vcf.data$FILTER == "PASS") & (GT %in% informative.GT) & (GQ >= min.GQ)))

# Filter VCF file according to various metrics.
usable.entries <- (vcf.data$FILTER == "PASS") & (GT %in% informative.GT) & (GQ >= min.GQ)
vcf.data.filtered <- vcf.data[usable.entries, ]
GT.filtered <- GT[usable.entries]
num.filtered <- length(vcf.data.filtered[, 1])
print(num.filtered)

# Switch around REF and ALT bases such that all entries are 0|1 (0 is REF, 1 is ALT). HiC-Pro also requires the GT field to be 0/1.
reversed <- GT.filtered == "1|0"
ref.temp <- vcf.data.filtered$REF[reversed]
alt.temp <- vcf.data.filtered$ALT[reversed]
vcf.data.filtered$REF[reversed] <- alt.temp
vcf.data.filtered$ALT[reversed] <- ref.temp
vcf.data.filtered$FORMAT <- rep("GT", num.filtered)
vcf.data.filtered$WTC <- rep("0/1", num.filtered)

# Write filtered table to file.
write.table(vcf.data.filtered, "/Users/michael/Box Sync/MS_general/analysis/brain_pchic_revision/resources/vcf/wtc_PASS_hg19.phased_v1.filtered.vcf", row.names=F, col.names=F, quote=F, sep="\t")


# End ---------------------------------------------------------------------

