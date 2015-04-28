#### Custom Diversity functions - 9th May 2014
## Now updated to take advantage of Rcpp functions

# updated statCalc to produce r_hat and heterozygosity estimates
mr_statCalc <- function (rsDat, idx = NULL, al, fst, bs = TRUE) {

#    rsDat = myData$genos
#     idx = NULL
#     al = myData$af
#     fst = TRUE
#     bs = FALSE
  
  if (bs) {
    rsFun <- function(x, y) {
      return(x[y, , ])
    }
    rsDat <- mapply(rsFun, x = rsDat, y = idx, SIMPLIFY = FALSE)
  }
  alf <- lapply(rsDat, function(x) {
    apply(x, 2, function(y) {
      if (all(is.na(y))) {
        return(NA)
      }
      else {
        y <- as.vector(na.omit(y))
        nms <- unique(y)[order(unique(y))]
        ot <- myTab(y)
        names(ot) <- nms
        return(ot)
      }
    })
  })
  nloci <- length(al)
  alf <- lapply(1:length(al), function(i) {
    lapply(alf, "[[", i)
  })
  alSort <- function(x, y) {
    idx <- lapply(x, function(z) {
      match(names(z), rownames(y))
    })
    for (i in 1:length(idx)) {
      y[idx[[i]], i] <- x[[i]]
    }
    return(y)
  }
  alOut <- mapply(alSort, x = alf, y = al, SIMPLIFY = FALSE)
  popSizes <- lapply(rsDat, function(x) {
    lgths <- apply(x, 2, function(y) {
      nrow(na.omit(y))
    })
    return(lgths)
  })
  ps <- do.call("cbind", popSizes)
  indtyp <- lapply(1:nloci, function(i) {
    vapply(rsDat, function(x) {
      op <- length(x[!is.na(x[, i, 1]), i, 1])
      if (op == 0L) {
        return(NA)
      }
      else {
        return(op)
      }
    }, FUN.VALUE = numeric(1))
  })
  if (fst) {
    hsums <- lapply(rsDat, function(x) {
      hts <- lapply(1:dim(x)[2], function(i) {
        out <- x[, i, 1] == x[, i, 2]
        return(out)
      })
      gts <- lapply(1:dim(x)[2], function(i) {
        x[, i, ]
      })
      alls <- lapply(gts, function(y) {
        unique(y[!is.na(y)])
      })
      htCount <- function(gts, hts, alls) {
        if (all(is.na(gts))) {
          out <- 0
          names(out) <- "NA"
          return(out)
        }
        else {
          gts <- gts[!hts, ]
          ht <- sapply(alls, function(al) {
            return(sum((gts == al), na.rm = TRUE))
          })
          names(ht) <- alls
          return(ht)
        }
      }
      htcounts <- mapply(htCount, gts = gts, hts = hts, 
                         alls = alls, SIMPLIFY = FALSE)
      return(htcounts)
    })
    hsums <- lapply(1:nloci, function(i) {
      lapply(hsums, "[[", i)
    })
    
    # calculate observed heterozygosity 
    ho <- lapply(rsDat, function(x){
      apply(x, 2, function(y){
        1 - (sum(na.omit(y[ ,1] == y[ ,2])) / length(na.omit(y[,1])))
      })
    })
    ho <- do.call(cbind, ho)
    
    # calculate expected heterozygosity 
    he <- t(sapply(al, function(x){
      apply(x, 2, function(y){
        return(1 - sum(y^2))
      })
    }))
    
    # calculate Chakraborty's (1992) r_hat for null alleles
    r_hat <- sapply(1:ncol(he), function(x){
      (he[, x] - ho[, x])/(he[, x] + ho[, x])
    })
    
    # correct r_hat
    # set NAN values to NA
    r_hat[is.na(r_hat)] <- 0.0001
    r_hat[r_hat == 1] <- 0.9999
    # set to negatives to zero
    r_hat[sign(r_hat) == -1] <- 0.0001
    
    # prep for output
    # ho
    ho <- split(ho, row(ho))
    names(ho) <- NULL
    # he
    he <- split(he, row(he))
    names(he) <- NULL
    # r_hat
    r_hat <- split(r_hat, row(r_hat))
    names(r_hat) <- NULL
    
    list(alOut = alOut, ps = ps, hsums = hsums, indtyp = indtyp, 
         ho = ho, he = he, r_hat = r_hat)
  }
  else {
    list(alOut = alOut, ps = ps, indtyp = indtyp)
  }
}

null_corrections <- function(inputs, mean_r_hat = FALSE){
  nloc <- length(inputs$indtyp)
  npop <- length(inputs$indtyp[[1]])
  n <- inputs$indtyp
  # determine whether to use mean r_hat or standard
  if(mean_r_hat) r_hat <- inputs$mean_r_hat else
    r_hat <- inputs$r_hat
  pj <- inputs$alOut
  hj <- inputs$hsums
  
  # derive corrected n
  ncorr <- lapply(1:nloc, function(x){
    n[[x]]/(1 - r_hat[[x]]^2)
  })
  
  # derive corrected pj
  pjcorr <- lapply(1:nloc, function(x){
    sapply(1:npop, function(y){
      pj[[x]][, y] - (pj[[x]][, y] * r_hat[[x]][y])
    })
  })
  
  # derive corrected hj
  # warnings can be ignored
  hjcorr <- lapply(1:nloc, function(x){
    sapply(1:npop, function(y){
      ((hj[[x]][[y]]*n[[x]][y]) + ((2*pjcorr[[x]][, y])*r_hat[[x]][y]*ncorr[[x]][y]))/ncorr[[x]][y]
    })
  })
  
  # alter structure of hjcorr for input into WC_rcpp function
  hjcorr<- lapply(hjcorr, function(x){
    z <- apply(x, 1, sum)
    return(z)
  })
  
  # alter input for output
  inputs$indtyp <- ncorr
  inputs$alOut <- pjcorr
  inputs$hsums <- hjcorr
  return(inputs)
}

# Fst
fstCalc <- function(a, b, cdat){
  return(a/(a+b+cdat))
}
# Fit
fitCalc <- function(a, b, cdat){
  return(1 - (cdat/(a+b+cdat)))
}
# Fis
fisCalc <- function(b, cdat){
  return(1 - (cdat/(b+cdat)))
}

# tabMerge function to process hsums
tabMerge <- function(...){
  ip <- unlist(list(...))
  idx <- names(ip) != "NA"
  ip <- ip[idx]
  out <- sapply(split(ip, names(ip)), sum)
  if(length(out) == 0L){
    ot <- NA
    names(ot) <- "NA"
    return(ot)
  } else {
    return(out)
  }
}

# calculate the number of missing genos per locus
missing_geno <- function(genos){
  npop <- length(genos)
  # how many genotypes per locus are missing?
  missing <- lapply(1:npop, function(x){
    eval <- is.na(genos[[x]][, , 1]) & is.na(genos[[x]][, , 2])
    colSums(eval)
  })
  return(do.call(cbind, missing))
}

# function to estimate expected proportion of null 
# homozygotes from r_hat
null_p <- function(inputs){
  r_hat <- inputs$r_hat
  r_hat <- do.call(rbind, r_hat)
  p_null <- r_hat^2 
  inputs$p_null <- p_null
  return(inputs)
}

# function to calculate number of null homozygotes 
# based on sample size

null_n <- function(inputs, nind){
  p_null <- inputs$p_null
  n_null <- sapply(1:length(nind), function(x){
    round(p_null[, x]*nind[x])
  })
  inputs$n_null <- n_null
  return(inputs)
}

# large function to retype missing genotypes as either nulls or technical
# errors
retype_missing <- function(genotypes, inputs){
  n_null <- inputs$n_null
  missing <- inputs$missing
  genos <- genotypes
  npop <- length(genos)
  nind <- lapply(genos, nrow)
  # is number of expected nulls greater than missing?
  decision <- n_null > missing
  
  # recreate and correct genotypes
  genotypes <- lapply(1:npop, function(z){
    out <- paste0(genos[[z]][, , 1], genos[[z]][, , 2])
    out <- matrix(out, nrow = nind[[z]])
    out[out == "NANA"] <- NA
    return(out)
  })
  
  new_genotypes <- lapply(1:length(genotypes), function(y){
    ev <- decision[, y]
    geno_df <- genotypes[[y]]
    nn <- n_null[, y]
    sapply(1:length(ev), function(z){
      if(ev[z] == FALSE){
        convert_to_null <- sample(which(is.na(geno_df[, z])), nn[z])
        geno_df[convert_to_null, z] <- "999999"
        geno_df[, z][is.na(geno_df[, z])] <- "000000"
        geno_df[, z]
      } else {
        geno_df[, z][is.na(geno_df[, z])] <- "999999"
        geno_df[, z]
      }
    })
  })
  return(new_genotypes) 
}

# calculate the number of homozygotes in a sample
nhoms_calc <- function(genos){
  npop <- length(genos)
  # how many genotypes per locus are missing?
  nhoms <- lapply(1:npop, function(x){
    eval <- genos[[x]][, , 1] == genos[[x]][, , 2]
    colSums(eval, na.rm = TRUE)
  })
  return(do.call(cbind, nhoms))
}

# create a het_hom matrix
het_hom_mat <- function(genos, sex){
  npop <- length(genos)
  # how many genotypes per locus are homozygotes?
  homs <- lapply(1:npop, function(x){
    eval <- genos[[x]][, , 1] == genos[[x]][, , 2]
  })
  homs <- do.call(rbind, homs)
  # convert to an easier to understand matrix
  homs <- ifelse(homs == TRUE, "hom", "het")
  # add in sex data
  homs <- cbind(sex, homs)
  
  # get counts of heterozygotes and homozygotes of both sexes
  het_homs <- apply(homs[, 3:ncol(homs)], 2, function(x){
    x <- factor(x, levels = c("het", "hom"))
    z <- table(x, homs$sex)
    z
  })
  rownames(het_homs) <- c("het_f", "hom_f", "het_m", "hom_m")
  return(t(het_homs))
}

# a function to extract a set of loci from a genepop file
extract_loci_gp <- function(myData, loci_to_extract){
  # identify indexes of loci to extract
  extract <- which(myData$locs %in% loci_to_extract)
  # get loci names
  new_locs <- myData$locs[myData$locs %in% loci_to_extract]
  # identify variables
  genos <- myData$genos
  npop <- length(genos)
  nind <- lapply(genos, nrow)
  # subset genotypes array
  genos <- lapply(1:npop, function(z){
    genos[[z]][, extract, ]
  })
  # recreate genotypes
  genotypes <- lapply(1:npop, function(z){
    out <- paste0(genos[[z]][, , 1], genos[[z]][, , 2])
    out <- matrix(out, nrow = nind[[z]])
    out[out == "NANA"] <- "000000"
    return(out)
  })
  output <- list(locs = new_locs, genos = genotypes, indnms = myData$indnms)
  return(output)
}

# calculate all possible genotype levels 
genotype_levels <- function(alleles){
  hets <- apply(combn(alleles, 2), 2, function(z){
    paste0(z[1], z[2])
  })
  homs <- paste0(unique(alleles), unique(alleles))
  loc_geno <- c(hets, homs)
  return(loc_geno[order(loc_geno)])
}

# a long and quite complicated function to perform g_tests for 
# sex linkage
sex.g.test <- function(myData, sex){
  require(Deducer)
  # identify variables
  genos <- myData$genos
  npop <- length(genos)
  nind <- lapply(genos, nrow)
  npop <- length(genos)
  # recreate genotypes
  genotypes <- lapply(1:npop, function(z){
    out <- paste0(genos[[z]][, , 1], genos[[z]][, , 2])
    out <- matrix(out, nrow = nind[[z]])
    out[out == "NANA"] <- "000000"
    return(out)
  })
  # combine to a single pop
  genotypes <- do.call(rbind, genotypes)
  # add in locus names
  colnames(genotypes) <- myData$locs
  # create counts for loci under autosomal exp
  autosomal <- sapply(1:ncol(genotypes), function(z){
    gen_lev <- genotype_levels(row.names(myData$af[[z]]))
    table(factor(genotypes[, z], levels = gen_lev), sex$sex)
  })
  
  # now generate counts expected under sex linkage
  sex_pop <- split(sex, factor(substr(sex$id, 1, 2)))
  # create a dataframe with sex included
  sex_alleles <- sapply(1:length(myData$genos), function(z){
    a <- rbind(myData$genos[[z]][, , 1], myData$genos[[z]][, , 2])
    b <- rbind(sex_pop[[z]], sex_pop[[z]])
    cbind(b, a)
  }, simplify = FALSE)
  sex_alleles <- do.call(rbind, sex_alleles)
  # count the occurrence of alleles in females
  female_counts <- apply(sex_alleles[, 3:ncol(sex_alleles)], 2, function(z){
    table(z, sex_alleles$sex)[, 'f']
  })
  # change the names to match genotype counts
  female_geno_counts <- lapply(female_counts, function(z){
    names(z) <- paste0(names(z), names(z))
    return(z)
  })
  # swap in to female counts
  sex_chr <- sapply(1:length(autosomal), function(z){
    autosomal[[z]][, 1] <- 0
    index <- which(names(autosomal[[z]][, 1]) %in% names(female_geno_counts[[z]]))
    autosomal[[z]][index] <- female_geno_counts[[z]]
    autosomal[[z]]
  })
  
  # perform g.tests
  sex_chr_con <- t(sapply(1:length(sex_chr), function(z){
    a <- likelihood.test(sex_chr[[z]], conservative = FALSE)
    c(as.numeric(a$parameter), as.numeric(a$statistic))
  }))
  autosomal_con <- t(sapply(1:length(autosomal), function(z){
    a <- likelihood.test(autosomal[[z]], conservative = FALSE)
    c(as.numeric(a$parameter), as.numeric(a$statistic))
  }))
  
  gstats <- data.frame(autosomal_con, sex_chr_con)
  names(gstats) <- c("auto_df", "auto_g", "sex_df", "sex_g")
  rownames(gstats) <- myData$locs
  return(gstats)
}

# function to correct r_hat by setting it to r_max if greater than
# r_max
r_hat_correct <- function(inputs, r_max){
  r_hat <- inputs$r_hat
  corr_r_hat <- sapply(1:length(r_hat), function(z){
    sapply(1:ncol(r_max), function(y){
      if(r_hat[[z]][y] > r_max[z, y]) r_max[z, y] else
        r_hat[[z]][y]
    })
  }, simplify = FALSE)
  inputs$r_hat <- corr_r_hat
  return(inputs)  
}

# calculate minor allele frequency
minor_allele_freq <- function(myData){
  nloc <- length(myData$locs)
  genos <- myData$genos
  maf <- sapply(1:nloc, function(z){
    x <- table(c(genos[[1]][, z, 1], genos[[1]][, z, 2],
                 genos[[2]][, z, 1], genos[[2]][, z, 2]))
    min(x/sum(x))
  })
  return(maf)
}

# calculate major allele frequency
major_allele_freq <- function(myData){
  nloc <- length(myData$locs)
  genos <- myData$genos
  maf <- sapply(1:nloc, function(z){
    x <- table(c(genos[[1]][, z, 1], genos[[1]][, z, 2],
                 genos[[2]][, z, 1], genos[[2]][, z, 2]))
    max(x/sum(x))
  })
  return(maf)
}

# quick function to caclulate the expected heterozygosity from
# freena allele frequency estimates
calculate_freena_het <- function(freena_het){
  z <- freena_het$est_fr
  names(z) <-  freena_het$pop 
  z <- split(z, as.factor(freena_het$locus))
  z <- lapply(z, function(y){
    a <- 1 - sum(y[names(y) == "1"]^2)
    b <- 1 - sum(y[names(y) == "2"]^2)
    c(a, b)
  })
  het <- do.call(rbind, z)
  mean_het <- apply(het, 1, mean)
  het <- data.frame(het, mean_het)
  names(het) <- c("freena_h1", "freena_h2", "freena_mean_h")
  return(het)
}