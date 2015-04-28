### Filtering stacks output prior to outlier analysis ###
# 14/10/14

rm(list = ls())

require(ggplot2)
require(gridExtra)
source("~/source/catalogue_diagnosis_functions.R")
source("~/source/custom_diversity.R")

### HANDLING DATA ###
# read in div stats, add in all depth data, filter down to haplotypes
stacks_dir <- "./example_data/" # set dir here
dirs <- "salto_071014" # set popn here
dirs <- dir(stacks_dir)[c(4, 5, 6)]

# set prefix
stacks_prefix <- paste0(stacks_dir, "/", dirs, "/", dirs)

myData <- lapply(stacks_prefix, function(z){
  
  # preduce pop output 
  pop <- create_hap_output(z, div_stats_suffix = ".div_stats_101014.txt",
                           depth_mat_suffix = "_depth_mat.txt", vcf_coords_suffix = "_vcf_coords.txt",
                           haplotypes_suffix = ".haplotypes_filtered.tsv")
  
  # read in and sort out freena
  free_na <- create_freena_output(fna_fst = paste0(z, ".gFst"),
                                  fna_rhat = paste0(z, ".r"),
                                  fna_het = paste0(z, ".fr"))
  
  # combine pop data and freena data
  pop <- cbind(pop, free_na)
  
  ### ESTIMATING STATS FOR FILTERING ###
  # derive differences in log-likelihoods
  pop$llr <- apply(pop[, c('auto_g', 'sex_g')], 1, function(z){
    (z['auto_g'] - z['sex_g'])/2
  })
  
  # set factor for sex/autsomal
  cut_off <- log(14.23)
  pop$chr <- factor(ifelse(pop$llr >= cut_off, "sex-linked", "autosomal"))
  
  # calculate whether major allele freq exceeds cutoff
  pop$maf_cutoff <- pop$maf*(2*pop$nind) > (2*pop$nind)-2
  
  ### FILTERING THE DATA ###
  nrow(pop)
  # filter loci outside 95% quantiles for mean depth of coverage
  doc_limit <- quantile(pop$loc_depth, c(0.025, 0.975))
  pop_filter <- pop[pop$loc_depth >= doc_limit[1] & pop$loc_depth <= doc_limit[2], ]
  nrow(pop_filter)
  # filter loci above 99% quantile for hom:het read depth ratio
  hhr_limit <- quantile(pop$hom_het_ratio[pop$hom_het_ratio != Inf], 0.99)
  pop_filter <- pop_filter[pop_filter$hom_het_ratio <= hhr_limit, ]
  nrow(pop_filter)
  # filter all loci above quantile for obs_hom:exp_hom ratio
  oheh_limit <- quantile(pop$obs_exp_hom, c(0.01, 0.99))
  pop_filter <- pop_filter[pop_filter$obs_exp_hom >= oheh_limit[1] & pop_filter$obs_exp_hom <= oheh_limit[2], ]
  nrow(pop_filter)
  # remove loci with an allele found in only one individual
  pop_filter <- pop_filter[!pop_filter$maf_cutoff, ]
  nrow(pop_filter)
  # remove loci with mean r_hat > 0.5
  pop_filter <- pop_filter[!pop_filter$mean_r_hat > 0.5, ]
  nrow(pop_filter)
  # remove loci where r_hat is >0.5 in either of the pops
  pop_filter <- pop_filter[pop_filter$r_hat1 <= 0.5, ]
  nrow(pop_filter)
  pop_filter <- pop_filter[pop_filter$r_hat2 <= 0.5, ]
  nrow(pop_filter)
  # remove loci where FIS is less than -0.25
  pop_filter <- pop_filter[!pop_filter$n_fis < -0.25, ]
  nrow(pop_filter)
  # filter out loci identified as sex-linked
  pop_filter <- pop_filter[pop_filter$chr == "autosomal",]
  nrow(pop_filter)
  
  # set up some informative statistics
  nloc <- rbind(c("unfiltered loci", nrow(pop)),
  c("loci outside depth limits", sum(pop$loc_depth < doc_limit[1] | pop$loc_depth > doc_limit[2])),
  c("loci outside hom_het ratio limit", sum(pop$hom_het_ratio > hhr_limit)),
  c("loci outside obs_hom, exp_hom limit", sum(pop$obs_exp_hom < oheh_limit[1] | pop$obs_exp_hom > oheh_limit[2])),
  c("loci with alleles in only a single individual", sum(pop$maf_cutoff)),
  c("loci with mean rhat > 0.5", sum(pop$mean_r_hat > 0.5)),
  c("loci with rhat > 0.5 in pop1", sum(pop$r_hat1 > 0.5)),
  c("loci with rhat > 0.5 in pop2", sum(pop$r_hat2 > 0.5)),
  c("loci with FIS < -0.25", sum(pop$n_fis < -0.25)),
  c("sex-linked loci", sum(pop$chr == "sex-linked")),
  c("total filtered loci", nrow(pop_filter)))
  
  list(pop = pop, pop_filter = pop_filter, nloc = nloc)
})

# set names
names(myData) <- unlist(lapply(strsplit(dirs, "_"), function(x) x[1]))

### DIAGNOSTIC PLOTS ###
## diagnostic plots for unfiltered data
# first extract legend
all_legend <- g_legend(hh_depth_plot(myData$ramso$pop) + theme(legend.position = "bottom", legend.key.width = unit(2, "cm")))
# set up a blank panel
blank <- grid.rect(gp = gpar(col = "white"))
# set plot options
gpar_opt <- gpar(cex = 1, mar = c(0, 0, 0, 0))
text_x_pos <- unit(0.5, "npc")
text_y_pos <- unit(0.5, "npc")
mar1 <- c(0, 0.5, 0, 0)
# arrange plot
grid.arrange(textGrob("Jutholmen", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$jut$pop) + theme(legend.position = "none",
                                                                      plot.margin = unit(mar1, "lines")),                                   
                         fis_dist_plot(myData$jut$pop) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             textGrob("Ramsö", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$ramso$pop) + theme(legend.position = "none",
                                                                        plot.margin = unit(mar1, "lines")),
                         fis_dist_plot(myData$ramso$pop) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             textGrob("Saltö", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$salto$pop) + theme(legend.position = "none",
                                                                        plot.margin = unit(mar1, "lines")),
                         fis_dist_plot(myData$salto$pop) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             arrangeGrob(all_legend, textGrob("Unfiltered data", gp = gpar(font = 2, cex = 1.5)), ncol = 2),
             nrow = 7, heights = c(1.5, 9, 1.5, 9, 1.5, 9, 4))
# Write out
dev.print(device=png, "~/unfiltered_popn_diagnostics.png", height = 800, width = 1200, res = 100)

## diagnostic plots for filtered data
# first extract legend
all_legend <- g_legend(hh_depth_plot(myData$ramso$pop_filter) + theme(legend.position = "bottom", legend.key.width = unit(2, "cm")))
# set up a blank panel
blank <- grid.rect(gp = gpar(col = "white"))
# set plot options
gpar_opt <- gpar(cex = 1, mar = c(0, 0, 0, 0))
text_x_pos <- unit(0.5, "npc")
text_y_pos <- unit(0.5, "npc")
mar1 <- c(0, 0.5, 0, 0)
# arrange plot
grid.arrange(textGrob("Jutholmen", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$jut$pop_filter) + theme(legend.position = "none",
                                                            plot.margin = unit(mar1, "lines")),                                   
                         fis_dist_plot(myData$jut$pop_filter) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             textGrob("Ramsö", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$ramso$pop_filter) + theme(legend.position = "none",
                                                             plot.margin = unit(mar1, "lines")),
                         fis_dist_plot(myData$ramso$pop_filter) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             textGrob("Saltö", gp = gpar_opt, x = text_x_pos, y = text_y_pos),
             arrangeGrob(hh_depth_plot(myData$salto$pop_filter) + theme(legend.position = "none",
                                                            plot.margin = unit(mar1, "lines")),
                         fis_dist_plot(myData$salto$pop_filter) + theme(plot.margin = unit(mar1, "lines")), 
                         ncol = 2),
             arrangeGrob(all_legend, textGrob("Filtered data", gp = gpar(font = 2, cex = 1.5)), ncol = 2),
             nrow = 7, heights = c(1.5, 9, 1.5, 9, 1.5, 9, 4))
# Write out
dev.print(device=png, "~//filtered_popn_diagnostics.png", height = 800, width = 1200, res = 100)

## plot potential sex_linked loci ##
# first extract legend
all_legend <- g_legend(sex_linkage_plot(myData$ramso$pop) + theme(legend.position = "bottom", legend.key.width = unit(2, "cm")))
# set up a blank panel
blank <- grid.rect(gp = gpar(col = "white"))
# set plot options
gpar_opt <- gpar(cex = 1, mar = c(0, 0, 0, 0))
text_x_pos <- unit(0.5, "npc")
text_y_pos <- unit(0.5, "npc")
mar1 <- c(0, 0.5, 0, 0)
# arrange plot
grid.arrange(sex_linkage_plot(myData$jut$pop) + ggtitle("Jutholmen") +
               theme(legend.position = "none", plot.margin = unit(mar1, "lines")),
             sex_linkage_plot(myData$ramso$pop) + ggtitle("Ramsö") +
               theme(legend.position = "none", plot.margin = unit(mar1, "lines")),
             sex_linkage_plot(myData$salto$pop) + ggtitle("Saltö") +
               theme(legend.position = "none", plot.margin = unit(mar1, "lines")),
             all_legend, nrow = 4, heights = c(4, 4, 4, 1))
# Write out
dev.print(device=png, "~/sex_linked_markers.png", height = 800, width = 1200, res = 100)

# set outfiles
outfiles <- paste0(stacks_dir, dirs, "/", names(filtered_data))
# write out data
sapply(1:length(outfiles), function(z){
  write.table(filtered_data[[z]]$pop_filter, paste0(outfiles[z], "_filtered_stats_141014.txt"),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
  write.table(filtered_data[[z]]$pop, paste0(outfiles[z], "_all_stats_141014.txt"),
              quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)
})



