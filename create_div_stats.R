### Generate locus statistics - including FST, null allele corrections and hom het counts
# 5th June 2014 - Mark Ravinet, National Institute of Genetics

rm(list = ls())
library(diveRsity)
library(Rcpp)
source("~/source/custom_diversity.R")
sourceCpp("~/source/glbWC.cpp")

args <- commandArgs(TRUE)
infile <- args[1]
sex <- read.delim(args[2], header = F)
outfile <- args[3]

names(sex) <- c("id", "sex")
levels(sex$sex)[3] <- NA

# read in data
myData <- rgp(infile)

pop_He <- sapply(1:length(myData$locs), function(y){
  1 - sum(apply(myData$af[[y]], 1, function(z) z[1]*z[2]))
})

# pop_Ho <- sapply(1:length(myData$locs), function(z){
#   hom <- c(myData$genos[[1]][, z, 1] == myData$genos[[1]][, z ,2],
#            myData$genos[[2]][, z, 1] == myData$genos[[2]][, z ,2])
#   (1 - sum(hom, na.rm = TRUE)/length(na.omit(hom)))
# })

# run custom statCalc
inputs <- mr_statCalc(rsDat = myData$genos, idx = NULL,
                      al = myData$af, fst = TRUE, bs = FALSE)

# derive number of alleles
n_all <- unlist(lapply(myData$af, nrow))

# calculate major allele frequency
maf <- major_allele_freq(myData)

# calculate nind 
nind <- unlist(lapply(myData$ps, sum))

# derive max_r
n_miss <- missing_geno(myData$genos)
h_sum <- nhoms_calc(myData$genos)
n <- 24
r_max <- (h_sum + (2*n_miss))/(2*n)

# correct r_hat to match r_max
inputs <- r_hat_correct(inputs, r_max)

# derive mean_r_hat
inputs$mean_r_hat <- lapply(inputs$r_hat, function(x) c(mean(x), mean(x)))

# format hsum output
hsum <- lapply(inputs$hsum, tabMerge)

### calculate uncorrected W+C statistics ###
allVarComp <- mapply(FUN = "glbWCcpp", hsum = hsum, af = inputs$alOut,
                     indtyp = inputs$indtyp, SIMPLIFY = FALSE)
# fst
fst <- sapply(allVarComp, function(x){
  return(fstCalc(a = x$a, b = x$b, cdat = x$c))
})
# fit
fit <- sapply(allVarComp, function(x){
  return(fitCalc(a = x$a, b = x$b, cdat = x$c))
})
# fis
fis <- sapply(allVarComp, function(x){
  return(fisCalc(b = x$b, cdat = x$c))
})
# he
wc_he <- sapply(allVarComp, function(x){
  return(x$b + x$c)
})
# new he
new_he <- sapply(allVarComp, function(x){
  return(x$c/(x$b + x$c))
})
# create a dataframe output
WC <- data.frame(fst, fit, fis, wc_he)

### calculated corrected W+C statistics with r_hat per popn
null <- suppressWarnings(null_corrections(inputs, mean_r_hat = FALSE))
# loop over all loci
allVarComp <- mapply(FUN = "glbWCcpp", hsum = null$hsums, af = null$alOut,
                     indtyp = null$indtyp, SIMPLIFY = FALSE)
# fst
corr_fst <- sapply(allVarComp, function(x){
  return(fstCalc(a = x$a, b = x$b, cdat = x$c))
})
# fit
corr_fit <- sapply(allVarComp, function(x){
  return(fitCalc(a = x$a, b = x$b, cdat = x$c))
})
# fis
corr_fis <- sapply(allVarComp, function(x){
  return(fisCalc(b = x$b, cdat = x$c))
})
# he
corr_wc_he <- sapply(allVarComp, function(x){
  return(x$b + x$c)
})
# create a dataframe output
corr_WC <- data.frame(corr_fst, corr_fit, corr_fis, corr_wc_he)

### calculated corrected W+C statistics with mean_r_hat
null <- suppressWarnings(null_corrections(inputs, mean_r_hat = TRUE))
# loop over all loci
allVarComp <- mapply(FUN = "glbWCcpp", hsum = null$hsums, af = null$alOut,
                     indtyp = null$indtyp, SIMPLIFY = FALSE)
# fst
corr_fst_mr <- sapply(allVarComp, function(x){
  return(fstCalc(a = x$a, b = x$b, cdat = x$c))
})
# fit
corr_fit_mr <- sapply(allVarComp, function(x){
  return(fitCalc(a = x$a, b = x$b, cdat = x$c))
})
# fis
corr_fis_mr <- sapply(allVarComp, function(x){
  return(fisCalc(b = x$b, cdat = x$c))
})
# he
corr_wc_he_mr <- sapply(allVarComp, function(x){
  return(x$b + x$c)
})
# create a dataframe output
corr_WC_mr <- data.frame(corr_fst_mr, corr_fit_mr, corr_fis_mr, corr_wc_he_mr)

# calculate numbers of het/hom males and females
het_hom <- het_hom_mat(myData$genos, sex)

# perform g test for autosomal and sex-linkage
sex_link <- sex.g.test(myData, sex)

### Final output ###
he <- do.call(rbind, inputs$he)
colnames(he) <- paste0("he", 1:ncol(he))
ho <- do.call(rbind, inputs$ho)
colnames(ho) <- paste0("ho", 1:ncol(ho))
r_hat <- do.call(rbind, inputs$r_hat)
colnames(r_hat) <- paste0("r_hat", 1:ncol(ho))
mean_r_hat <- do.call(rbind, inputs$mean_r_hat)[, 1]
colnames(r_max) <- paste0("r_max", 1:ncol(ho))

out <- data.frame(he, ho, n_all, r_hat, mean_r_hat, r_max, WC, corr_WC, corr_WC_mr, 
                  het_hom, sex_link, maf, nind)
rownames(out) <- myData$locs
write.table(out, outfile, quote = FALSE, sep = "\t")

