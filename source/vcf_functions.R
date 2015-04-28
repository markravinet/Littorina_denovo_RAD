# VCF functions

# basic VCF reading function, returns a list with meta, header and data
readVCF <- function(infile){
  # use readLines to read into R
  vcf <- readLines(infile)
  # locate header format
  header_loc <- grep("#CHROM", vcf)
  # assign meta data
  meta <- vcf[1:header_loc-1]
  # assign header
  header <- unlist(strsplit(vcf[header_loc], "\t"))
  # assign data and remove vcf object
  data_lines <- vcf[-c(1:header_loc)]
  remove(vcf)
  # read vcf data as dataframe - maybe, it might be easier to parse line by line
  data <- read.delim(infile, header = TRUE, skip = header_loc-1)
  # output a list with meta, header and data as elements
  out = list(meta = meta, header = header, data = data)
  return(out)
}

# Simple function to parse VCF info section
# NOTE - SPECIFIC TO STACKS ONLY AT THIS STAG
parse_vcf_info <- function(data){
  info <- data$INFO
  # use sub to remove the INFO string and add sep
  info <- sub("NS=", "", info)
  info <- sub("AF=", "", info)
  info <- sub(",", ";", info)
  # return a 3 col matrix of seq depth, af1 and af2
  info <- strsplit(info, ";", )
  info <- do.call(rbind, info)
class(info) <- "numeric"
  return(info)
}

# simple function to parse genotype and depth data per ind
# for now this function does not output likelihood data
# does not format data
parse_vcf_geno <- function(data, bases = FALSE){
  ind_names <- names(data[, 10:ncol(data)])
  x <- as.matrix(data[, 10:ncol(data)])
  nind <- ncol(x)
  nloc <- nrow(x)
  split_x <- strsplit(x, ":")
  split_x <- do.call(rbind, split_x)
  geno <- matrix(split_x[,1], nrow = nloc, ncol = nind)
  if(bases == TRUE){
    bases <- data[, c('REF', 'ALT')]
    bases_geno <- sapply(1:nrow(geno), function(z){
      # first replace the ref allele
      y <- gsub("0", as.character(bases[z, 1]), geno[z, ])
      # then replace alt allele
      y <- gsub("1", as.character(bases[z, 2]), y)
      # remove all unecessary chars
      y[y == "./."] <- NA
      y
    })
    geno <- t(bases_geno)
  } else{
    geno <- data.frame(geno)
    names(geno) <- ind_names
    depth <- matrix(split_x[,2], nrow = nloc, ncol = nind)
    depth <- data.frame(depth)
    names(depth) <- ind_names
    out <- list(geno = geno, depth = depth)
  }
  return(out)
}

# convert VCF genotype data to binary code
geno_to_binary <- function(x){
  x <- as.character(x)
  x[x == "./."] <- NA
  split <- strsplit(x, "/")
  binary <- do.call(rbind, lapply(split, function(x) sum(as.numeric(x))))
  binary[is.na(binary)] <- "?"
  return(binary)
}

# a larger wrapper function to create haplotypes per RAD locus
haplo_by_locus <- function(data, onesnp = FALSE){
  # create a vector of loci Id
  loci <- unique(data$ID)
  # use sapply to get a per locus
  haplo_mat <- sapply(loci, function(x){
    # create a vector of loci IDs
    loci_subset <- data[data$ID == x, ]
    if(onesnp == TRUE)  loci_subset <- loci_subset[1, ] else
      loci_subset <- loci_subset
    # extract the real genotypes
    real_gt <- t(apply(loci_subset, 1, geno_extractor))
    # now generate haplotype
    haplo_gen(real_gt)
  })
  colnames(haplo_mat) <- loci
  return(haplo_mat)
}

# a function wrapper for the haplotype generator apply function
haplo_gen <- function(real_gt){
  apply(real_gt, 2, function(y){
    if(NA %in% y) NA else
      paste0(paste(sapply(y, function(x) substr(x, 1, 1)), collapse = ""), "/",
             paste(sapply(y, function(x) substr(x, 3, 3)), collapse = ""))
  })
}

# a function that will take a line of vcf data
# and extract the real genotype
geno_extractor <- function(x){
  # set up ref and alt base
  ref <- x[4] # if 0, this base
  alt <- x[5] # if 1, this base
  # extract genotype
  geno <- substr(unlist(x[10:length(x)]), 1, 3)
  # convert to base
  real_geno <- sapply(geno, function(x){
    if(x == "./.") NA else
      if(x == "0/0") paste0(ref, "/", ref) else
        if(x == "1/1") paste0(alt, "/", alt) else
          if(x == "0/1") paste0(ref, "/", alt) else
            if(x == "1/0") paste0(alt, "/", ref)    
  }, USE.NAMES = FALSE)
  return(real_geno)
}

# custom functions to shift to main source code
# write_seqs - a wrapper function to write sequences from subset data
# write_seqs <- function(loc_seq, subset_data, subset_sample, outfile){
#   lapply(1:length(loc_seq), function(x){
#     # write out sequence
#     cat(">", sub(".bam", "", names(subset_data)[ncol(subset_data)]),
#         ":", as.character(subset_data$X.CHROM[1]), ":", sep ="", 
#         file = outfile, append = TRUE)
#     # Target coords
#     cat(subset_sample$start[x], "-", subset_sample$stop[x],
#         ":", sep = "", file = outfile, append = TRUE)
#     # True coords
#     cat(subset_data[min(names(loc_seq[[x]])),]$POS, "-", 
#         subset_data[max(names(loc_seq[[x]])),]$POS,
#         ":", sep = "", file = outfile, append = TRUE)
#     # true bp
#     cat(sum(nchar(loc_seq[[x]])), "\n", file = outfile, append = TRUE)
#     # Actual sequence
#     cat(loc_seq[[x]], "\n", sep = "", file = outfile, append = TRUE)
#   })
# }

# write_seqs - a wrapper function to write sequences from subset data
write_seqs <- function(loc_seq, subset_data, subset_sample, outfile){
  lapply(1:length(loc_seq), function(x){
    # write out sequence
    cat(">", sub(".bam", "", names(subset_data)[ncol(subset_data)]),
        ";", sep ="", 
        file = outfile, append = TRUE)
    # write locus target and actual coords
    cat(names(loc_seq[x]), file = outfile, append = TRUE)
    # true bp
    cat(";", sum(nchar(loc_seq[[x]])), "\n", sep = "", file = outfile, append = TRUE)
    # Actual sequence
    cat(loc_seq[[x]], "\n", sep = "", file = outfile, append = TRUE)
  })
}

# parse_loci - a wrapper function to parse loci on a subset chromosome
parse_loci <- function(subset_data){
  # first loop through and find diffs between POS
  POS_diff <- sapply(2:length(subset_data$POS), function(x){
    y <- subset_data$POS[x] - subset_data$POS[x - 1]
  })
  # look for position diffs greater that are greater than the specified interval
  # this gives the start locations of the loci (note first one not included) so add
  # 1 to get all the starting positions
  POS_diff_loc <- which(POS_diff >= 100000)+1
  loc_start <- c(1, POS_diff_loc)
  # then get the locus stop positions
  loc_end <- c((POS_diff_loc - 1), nrow(subset_data))
  # parse sequences for loci on chromosome
  loc_seq <- lapply(1:length(loc_end), function(x){
    y <- c(loc_start[x], loc_end[x])
    apply(subset_data[y[1]:y[2], ], 1, function(x){
      base_call(x, q_cutoff = 30, depth_cutoff = 30,
                mapq_cutoff = 60)
    })
  })
  return(loc_seq)
}

# base_call - wrapper function to call base per line
base_call <- function(vcf_line, q_cutoff,
                      depth_cutoff, mapq_cutoff){
  
  # store ref and alt bases, extract q value
  ref <- as.character(unlist(vcf_line["REF"]))
  alt <- as.character(unlist(vcf_line["ALT"]))
  q <- as.character(unlist(vcf_line["QUAL"]))
  
  # if call quality is not high, reference base is printed
  if(q < q_cutoff) print(ref) else{
    if(alt == ".") print(ref) else{
      # if call quality cutoff is satisfied, rest of data is extracted
      # read per line and process for depth, mapping quality, call quality (FQ)
      info <- as.character(unlist(vcf_line["INFO"]))
      info <- unlist(strsplit(info, ";"))
      # extract read depth and vcf_line threshold
      read_depth <- as.numeric(unlist(strsplit(info[grep("DP=", info)], "="))[2])
      # extract mapping quality
      mapq <- as.numeric(unlist(strsplit(info[grep("MQ=", info)], "="))[2])
      # make sure data passes depth and mapping qualities, print ref if not
      if(read_depth < depth_cutoff | mapq < mapq_cutoff) print(ref) else{
        # if so 
        parse_genotype_data(vcf_line[names(vcf_line)[10]], ref, alt)
      } 
    }
  }  
}

# write a single base (for FASTA)
parse_genotype_data <- function(sample_genotype, ref, alt){
  # first extract info
  geno_info <- as.character(sample_genotype)
  geno_info <- unlist(strsplit(geno_info, ":"))
  genotype <- geno_info[1]
  
  # call base based on genotype - randomly sample heterozygotes
  if(genotype == "0/0") print(ref) else
    if(genotype == "0/1" | genotype == "1/0") print(sample(c(ref, alt), 1)) else
      if(genotype == "1/1") print(alt)
}

# same as above but it writes out both bases rather than one
# and considers cut offs at genotype level too
parse_full_genotype_data <- function(sample_genotype, ref, alt,
                                     site_depth_cutoff, genotype_q_cutoff){
  
  # first extract info
  geno_info <- as.character(sample_genotype)
  geno_info <- unlist(strsplit(geno_info, ":"))
  genotype <- geno_info[1]
  # phred_L <- geno_info[2] redundant info for now
  depth <- as.numeric(geno_info[3])
  geno_q <- as.numeric(geno_info[4])
  
  if (depth < site_depth_cutoff | geno_q < genotype_q_cutoff) out <- NA else
    # call base based on genotype - randomly sample heterozygotes
    if(genotype == "0/0") out <- c(ref, ref) else
      if(genotype == "0/1" | genotype == "1/0") out <- c(ref, alt) else
        if(genotype == "1/1") out <- c(alt, alt)
  
  return(out)
}

# quickly estimate phredq probability
phredq_convert <- function(q){
  return(10^(q/10)) 
}


# extract_locus_name - a function to extract locus name 
extract_locus_name <- function(locus, subset_data, subset_sample){
  # identify and assign target start and stop sites
  target_start <- subset_data[min(names(locus)),]$POS
  target_stop <- subset_data[max(names(locus)),]$POS
  # find which line in the sample the target matches
  y <- subset_sample[target_start >= subset_sample$start & target_stop <= subset_sample$stop, ]
  z <- paste0(y[1], ":", y[2], "-", y[3], ";", target_start, "-", target_stop)
  return(z)
}

# parse vcf allele frequency data and return it as a matrix
vcf_af <- function(myData){
  
  x <- (myData$INFO)
  info <- strsplit(as.character(x), ";")
  af_data <- sapply(info, function(z){
    y <- strsplit(z[2], ",")
  })
  af_data <- lapply(af_data, function(z){
    as.numeric(sub("AF=", "", z))
  })
  af_mat <- do.call(rbind, af_data)
  return(af_mat)
  
}

# convert a genotype matrix to single bp
# using IUPAC codes
convert_to_single_base <- function(gt){
  single_base <- sapply(gt, function(z){
    if(z == "C/C")  "C" else
      if(z == "A/A") "A" else
        if(z == "G/G") "G" else
          if(z == "T/T") "T" else
            if(z == "G/A" | z == "A/G") "R" else
              if(z == "C/T" | z == "T/C") "Y" else
                if(z == "G/C" | z == "C/G") "S" else
                  if(z == "A/T" | z == "T/A") "W" else
                    if(z == "G/T" | z == "T/G") "K" else
                      if(z == "C/A" | z == "A/C") "M" else
                        if(z == "./.") "N" else
                          if(grepl("N", z)) "N" else z
  })
  return(single_base)
}