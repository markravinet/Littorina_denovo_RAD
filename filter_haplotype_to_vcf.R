# Filter Stacks haplotype to match the vcf 
rm(list = ls())

source("~/source/vcf_functions.R")

# Specify command line arguments
args <- commandArgs(TRUE)
infile_hap <- read.delim(args[1])
infile_vcf <- readVCF(args[2])
outfile <- args[3]

# read haplotype file
stacks_hap <- infile_hap

# read vcf file
stacks_vcf <- infile_vcf

# extract vcf loci
vcf_loci <- unique(stacks_vcf$data$ID)
# extract hap file loci
hap_loci <- unique(stacks_hap$Catalog.ID)

# then filter the haplotype data for the loci in the vcf - generate an output
filtered_stacks_hap <- stacks_hap[which(hap_loci %in% vcf_loci), ]

# write out the filtered output
write.table(filtered_stacks_hap, outfile, quote = FALSE, sep = "\t")

