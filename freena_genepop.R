### Rewrite genepop file to account for null alleles ####
# New version to split the files into chunks of 5000 loci
# 14th May 2014 - Mark Ravinet, National Institute of Genetics

rm(list = ls())
library(diveRsity)
source("~/source/custom_diversity_09052014.R")

args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]
header <- args[3]
use.n <- args[4]
n <- as.numeric(args[5])

myData <- rgp(infile)
# if not filtering by n, set n to length of loci
if(use.n == FALSE) n <- length(myData$locs)

# code below provides the indices to write multiple genepop files
# split the loci into blocks of 5000
loc_bloc <- split(myData$locs, ceiling(seq_along(myData$locs)/n))
# get the first and last loci of each
loc_bloc <- lapply(loc_bloc, function(z){
  c(z[1], z[length(z)])
})
# find the indexes of these loci
loc_bloc <- lapply(loc_bloc, function(z){
  which(myData$locs %in% z)
})

# run custom statCalc
inputs <- mr_statCalc(rsDat = myData$genos, idx = NULL,
                      al = myData$af, fst = TRUE, bs = FALSE)

# convert for FREENA
inputs <- null_p(inputs)
# get total ind per pop
nind <- sapply(myData$indnms, length)
# get number of expected null homozygotes
inputs <- null_n(inputs, nind)
# get number of null/missing genotypes
inputs$missing <- missing_geno(myData$genos)
# retype genotypes as either nulls or technical errors
new_genos <- retype_missing(myData$genos, inputs)


loci_names <- lapply(loc_bloc, function(z){
  # prep to write out a genepop file
  myData$locs[z[1]:z[2]]
  #pop_data <- lapply(1:length(new_genos), function(x){
   # cbind(myData$indnms[[x]], new_genos[[x]][, z[1]:z[2]])
  })

pop_data <- lapply(loc_bloc, function(z){
  pop_data <- lapply(1:length(new_genos), function(x){
  cbind(myData$indnms[[x]], new_genos[[x]][, z[1]:z[2]])
})
})

sapply(1:length(loc_bloc), function(z){
  # Now write the output 
  newfile <- file(paste0(outfile, "_", z, ".genepop"), "w")
  cat(header, "\n", file = newfile, append = TRUE)
  cat(loci_names[[z]], sep = "\n", file = newfile, append = TRUE)
  for(i in 1:length(pop_data[[z]])){
    cat("pop\n", file = newfile, append = TRUE)
    write.table(pop_data[[z]][[i]], file = newfile, append = TRUE, col.names = FALSE,
                row.names = FALSE, sep = " ", quote = FALSE)
  }
  close(newfile)
})
