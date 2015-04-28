### Catalogue diagnosis functions ###
# 13th October 2014
# Mark Ravinet

# make sex linked plot
sex_linkage_plot <- function(pop){
  sex_plot <- ggplot(pop, aes(het_m ,het_f)) + geom_point(aes(colour = chr), alpha = .6, position = "jitter") + 
    xlab("No. of male heterozygotes") + ylab("No. of female heterozygotes")
  return(sex_plot)
}

# a function to create the output requited to combine freena data with diversity stats
create_freena_output <- function(fna_fst, fna_rhat, fna_het){
  
  # read in freena input files
  freena_fst <- read.delim(fna_fst, skip = 16, header = F, sep = "")
  freena_rhat <-  read.delim(fna_rhat, skip = 11, header = F, sep = "")
  freena_het <- read.delim(fna_het, skip = 10, header = F, sep = "")
  
  # organise freena data.frame
  freena <- cbind(freena_rhat[freena_rhat[, 2] == 1, 3], freena_rhat[freena_rhat[, 2] == 2, 3],
                  freena_fst[, 2:3])
  names(freena) <- c("fna_rhat1", "fna_rhat2", "raw_fst", "ena_fst")
  
  # calculate freena heterozygosity
  names(freena_het) <- c("locus", "pop", "allele", "allele_size",
                         "obs_fr", "est_fr", "act_hom", "est_hom")
  freena_het <- calculate_freena_het(freena_het)
  
  # combine freena data
  freena_out <- cbind(freena, freena_het)
  return(freena_out)
  
}

#extract legend
#https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# a function to create a hom-het plot
hh_depth_plot <- function(pop){
  outplot <- ggplot(pop, aes(exp_hom_n, nhom, colour = loc_depth)) + geom_point(position = "jitter") +
    geom_abline() + xlab("Expected homozygotes") + ylab("Observed homozygotes") + 
    scale_colour_gradient(name = "Mean\n depth", high = "red", low = "blue") +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          legend.text = element_text(size = 11),
          legend.title = element_text(size = 11),
          legend.background = element_blank(),
          legend.key = element_blank())
  return(outplot)
}

# a function to create an FIS distribution plot
fis_dist_plot <- function(pop){
  outplot <- ggplot(pop, aes(n_fis)) + geom_histogram(colour = "black", fill = "white", binwidth = 0.02) +
    xlab(expression(paste("F"["IS"]))) +
    theme(axis.text.x = element_text(size = 11),
          axis.text.y = element_text(size = 11),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12))
  return(outplot)
}

# a function to create haplotype output for diagnostic plots
create_hap_output <- function(prefix, div_stats_suffix, depth_mat_suffix,
                              vcf_coords_suffix, haplotypes_suffix, all = TRUE){
  
  # read in custom diveRsity data
  pop <- read.delim(paste0(prefix, div_stats_suffix))
  # read in depth matrix - nb this for SNPs
  depth_mat <- read.delim(paste0(prefix, depth_mat_suffix))
  #depth_mat <- depth_mat[1:(nrow(depth_mat)-2), ]
  # read in vcf coords to get the RAD loci identities
  vcf_coord <- read.delim(paste0(prefix, vcf_coords_suffix))[, 3]
  # combine the two
  depth_mat <- data.frame(vcf_coord, depth_mat)
  # remove duplicates to get depth per RAD loci
  RAD_depth <- depth_mat[!duplicated(depth_mat[, 1]), ]
  # calculate mean depth per loci
  loc_depth <- apply(RAD_depth[, 4:ncol(RAD_depth)], 1, function(z){
    mean(z[z > 0])
  })
  names(loc_depth) <- RAD_depth$vcf_coord
  
  # now read in haplotype matrix and convert depth
  haplo <- read.delim(paste0(prefix, haplotypes_suffix))
  new_haplo <- t(haplo[, 3:ncol(haplo)])
  new_haplo[new_haplo == "-"] <- NA
  het_hom <- apply(new_haplo, 2, function(z){
    z[grepl("/", z)] <- "het"
    z[z != "het"] <- "hom"  
    return(z)
  })
  colnames(het_hom) <- haplo$Catalog.ID
  # transform RAD depth data
  haplo_depth <- t(RAD_depth[, 4:ncol(RAD_depth)])
  colnames(haplo_depth) <- RAD_depth$vcf_coord
  
  # get mean het hom depth
  hh_depth <- t(sapply(1:ncol(het_hom), function(z){
    hom_depth <- mean(haplo_depth[which(het_hom[, z] == "hom"), z])
    het_depth <- mean(haplo_depth[which(het_hom[, z] == "het"), z])
    nhom <- sum(het_hom[, z] == "hom", na.rm = TRUE)
    nhet <- sum(het_hom[, z] == "het", na.rm = TRUE)
    return(c(hom_depth, het_depth, nhom, nhet))
  }))
  hh_depth <- data.frame(hh_depth)
  rownames(hh_depth) <- RAD_depth$vcf_coord
  colnames(hh_depth) <- c("hom_depth", "het_depth", "nhom", "nhet")
  # deal with NANs
  hh_depth$hom_depth[is.nan(hh_depth$hom_depth)] <- 0
  hh_depth$het_depth[is.nan(hh_depth$het_depth)] <- 0
  # calc hom_het_ratio
  hh_depth$hom_het_ratio <- hh_depth$hom_depth/hh_depth$het_depth
  
  # now filter and add to population data.frame
  # first combine loc depth 
  hh_depth <- data.frame(hh_depth, loc_depth)
  # add to pop
  pop <- data.frame(pop, hh_depth[which(rownames(pop) %in% rownames(hh_depth)), ])
  
  # calc additional FIS values for pop
  if(all == TRUE){
    pop$he <- apply(cbind(pop$he1, pop$he2, pop$he3), 1, mean)
    pop$ho <- apply(cbind(pop$ho1, pop$ho2, pop$ho3), 1, mean)
  } else{
    pop$he <- apply(cbind(pop$he1, pop$he2), 1, mean)
    pop$ho <- apply(cbind(pop$ho1, pop$ho2), 1, mean)
  }
  # add to dataframe
  pop$n_fis <- (pop$he - pop$ho)/pop$he
  # calculate expected number of homozygotes
  pop$nind_a <- pop$nhom + pop$nhet
  pop$exp_hom_n<- pop$nind_a*(1-pop$he)
  # calculate hom het obs ratio
  pop$obs_exp_hom <- pop$nhom/pop$exp_hom_n
  return(pop)
}