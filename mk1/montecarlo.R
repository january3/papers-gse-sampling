#!/usr/bin/env Rscript

## generate 100 samplings for each number of samples per group
## n must be provided as the number of samples to generate
## ss must be provided as the number of individuals per group
## covar_file is the name of the covariate file
## count_file is the name of the counts file

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 5) {
  warning("Usage: [script] n ss covar_file count_file outputdir")
  stop()
}


n          <- as.numeric(args[1])
ss         <- as.numeric(args[2])
covar_file <- args[3]
count_file <- args[4]
outputdir  <- args[5]

message(sprintf("n=%d, ss=%d, covar_file=%s, count_file=%s, outputdir=%s",
                n, ss, covar_file, count_file, outputdir))

counts <- readRDS(count_file)
covar  <- readRDS(covar_file)

library(DESeq2)
outputfile_fmt <- "deseq_obj_n_%06d_ss_%03d.rds"

groups <- levels(factor(covar$group))
stopifnot(length(groups) == 2)

g1 <- which(covar$group == groups[1])
g2 <- which(covar$group == groups[2])

for(i in 1:n) {
  g1_sel <- sample(g1, ss)
  g2_sel <- sample(g2, ss)
  sel <- c(g1_sel, g2_sel)


  ds <- DESeqDataSetFromMatrix(counts[,sel], colData=covar[sel, ], design=~ group)
  ds <- DESeq(ds)

  save_file <- file.path(outputdir, sprintf(outputfile_fmt, i, ss))

  saveRDS(ds, file=save_file)
}




