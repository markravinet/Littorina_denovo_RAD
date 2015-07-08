# Littorina_denovo_RAD

Scripts used for analysis of L. saxatilis de novo RAD assembly. The scripts here are used for conversion of stacks outputs into various files that can then be used for downstream analysis. You can use the example data in the example data directory to test them out but the full dataset is also available on <a href="https://datadryad.org/resource/doi:10.5061/dryad.g466r
">Dryad</a>.

## Filter the haplotypes file

First, use the filter_haplotype_to_vcf.R script to filter down the haplotype output produced by Stacks. This will remove all non-polymorphic loci and ensure the haplotype and vcf files match. For example, using R on the command line you would do the following:
```
Rscript filter_haplotype_to_vcf.R in.haplotype in.vcf filtered_hap.out
```

## Convert haplotype to gene pop

Next you want to convert this filtered haplotype file into a genepop file for later analysis. To do so use the haplotype_to_genepop.R script like so:
```
Rscript haplotype_to_genepop.R filtered_hap.out filtered_hap.genepop my_header
```

## Estimate diversity statistics

With your genepop file, you now need to estimate various diversity statistics. The code for this script and it's source are modified versions of some early code written by Kevin Keenan for his <a href="https://github.com/kkeenan02/diveRsity
">diveRsity package</a>. So be sure to check that out for more up-to-date code. To run the custom code here, use the create_div_stats.R script like so:
```
Rscript create_div_stats.R filtered_hap.genepop ind_sex.txt div_stats.txt
```
Here, the second file is a text file with two columns, sample name and sex. The final argument is the output

## Custom filtering

To perform the custom filtering explained in the manuscript, you can use the custom_RAD_loci_filtering.R script. This script is easier to run in the R terminal, altering the arguments within to suit your own needs.

## Prepare for FreeNA analysis

In the manuscript, we used FreeNA to correct for null alleles. This can be achieved using the same genepop file we produced earlier and the freena_genepop.R script like so:
```
Rscript freena_genepop.R filtered_hap.genepop freena.genepop my_header
```
Basically the script will estimate null allele frequencies and then recode missing individuals as 999999 to indicate null alleles. 


