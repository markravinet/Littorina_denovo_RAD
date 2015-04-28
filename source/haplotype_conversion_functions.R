### Stacks haplotype conversion - functions ###
# Mark Ravinet - October 2013

# wrapper function to produce filtered haplotype output
# added 11th December 2014
haplo_processor <- function(infile, depth, vcf_coord){
  
  loc_depth <- depth[!duplicated(depth[, 1]), ]
  rad_depth <- apply(loc_depth[, -c(1:3)], 1, sum)
  names(rad_depth) <- loc_depth$vcf_coord
  
  # Assign RAD_loci names
  loci <- infile$Catalog.ID
  
  # Create a haplotype matrix
  haplo <- file_convert(infile)
  
  # remove inds with > 2 alleles
  # identify individual with > 2 aleles
  allele_check <- apply(haplo, 2, function(z){
    str_count(z, "/") > 1
  })
  haplo[allele_check] <- NA
  # remove all blank calls
  haplo[haplo == "-"] <- NA
  
  # calc expected He
  he <- apply(haplo, 2, function(z){
    count <- table(unlist(strsplit(z, "/")))
    freq <- count/sum(count)
    1-sum(freq^2)
  })
  
  # filter where He is < 0.1
  haplo <- haplo[, !he < 0.1]
  
  # work out the number of individuals typed
  ind_typ <- apply(haplo, 2, function(z) length(na.omit(z)))
  
  # placeholder for now
  rad_depth <- head(rad_depth)
  depth_ind <- rad_depth/ind_typ
  # remove loci with depth outside this quantile
  depth_lim <- quantile(depth_ind, c(0.025, 0.975))
  haplo <- haplo[, depth_ind >= depth_lim[1] & depth_ind <= depth_lim[2]]
  
  # now create a presence absence matrix
  haplo_pa <- !is.na(haplo)
  # filter based on too large or too small p/a scores
  pa_ind <- apply(haplo_pa, 1, sum) 
  pa_lim <- quantile(pa_ind, c(0.025, 0.975))
  haplo_pa <- haplo_pa[pa_ind >= pa_lim[1] & pa_ind <= pa_lim[2], ]
  
  return(haplo_pa)
  
}

# file_convert - converts from original TSV to a haplotype matrix
file_convert <- function(infile){
  
  haplo <- t(infile[, -1:-2])
  colnames(haplo) <- infile$Catalog.ID

  return(haplo)
}

# Calculate nucleotide diversity using Nei & Li's (1979) method
nuc.div <- function(x, seq_length) {
  
  haplotype_freq <- table(x)/length(x)
  
  if(length(unique(x)) <= 1){
    nd <- 0
  } else{  
    
    haplotype_comb <- combn(unique(x), 2)
    # Calculate pairwise nucleotide mismatches
    pairwise_mm <- apply(haplotype_comb, 2, function(z) {
      mapply(function(x, y) sum(x!=y), strsplit(z[1], ""), strsplit(z[2], ""))
    })
    # Calculate nucleotide diversity from frequency and mismatch data
    multi_freq <- apply(haplotype_comb, 2, function(x) haplotype_freq[x[1]] * haplotype_freq[x[2]])
    nd <- 2*sum(multi_freq*pairwise_mm)/seq_length
  }
  
  return(nd)
}

# Basic method
# # Based on formula in Hartl and Clark - poss computationally slower
# L <- 92
# pair_comp <- combn(seq, 2)
# pair_mm <- apply(pair_comp, 2, function(z) {
#   mapply(function(x, y) sum(x!=y), strsplit(z[1], ""), strsplit(z[2], ""))
# })
# 
# nuc_mm <- sum(pair_mm)/2346
# 
# ((length(seq)^2)-length(seq))
# 
# ncol(pair_comp)
# nuc_div <- nuc_mm/L 
