#!/usr/bin/env Rscript

# metadata_duplicate_control.r file.tsv

# this script duplicates the metdata rows and appends "_control" to the duplicate sample IDs within "qiime2_diversity_core.nf" to aid "benchmarking_diversity.nf"
# the first row must be the sample IDs!

# Get params and files from the command line
args <- commandArgs(trailingOnly=TRUE)
FILE <- args[1]

print(paste("Use file",FILE))

df = read.table( FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# make copy and add "_control" to first columnn
con <- df
con[,1] <- paste0( con[,1],"_control" )

# append rows
df <- rbind( df, con )

# Overwrite input
outfile <- paste0(FILE,".final.tsv")
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")
