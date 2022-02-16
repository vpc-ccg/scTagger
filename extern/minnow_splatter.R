#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(splatter)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]
num_genes <- as.integer(args[2])
num_cells <- as.integer(args[3])
seed <- as.integer(args[4])
rm(args)

set.seed(seed)
sim <- splatSimulate( 
	nGenes=num_genes, 
	batchCells=num_cells, 
	verbose = FALSE
)

dir.create(out_dir, showWarnings=FALSE)
write.table(colnames(sim), file= file.path(out_dir, "quants_mat_cols.txt"), quote=FALSE, col.names=FALSE, row.names=FALSE)
write.table(counts(sim), file= file.path(out_dir, "quants_mat.csv"), quote=FALSE, col.names=FALSE, row.names=FALSE, sep=",")  
