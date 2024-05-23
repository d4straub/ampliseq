#!/usr/bin/env Rscript

# benchmarking_taxonomy_plot_stats.r file.tsv "F1,recall,precision"

# file.tsv must have headers: match,match_higher,match_higher_nomatch_lower,nomatch,nomatch_higher,not_calculated,method_reftaxdb

#show as horizontal barplot, printing counts on top, colours ranging from blue to red for "Genus & Species expected" (match), "Genus expected, Species not classified" (match_higher), "Genus expected, but Species unexpected" (match_higher_nomatch_lower), "Genus and Species unexpected" (nomatch), "Genus unexpected, Species not classified" (nomatch_higher), but gray for "Not classified" (not_calculated)

# Get params and files from the command line
args <- commandArgs(trailingOnly=TRUE)
FILE <- args[1]
TYPES <- args[2]

# TYPES to list
TYPES <- strsplit(TYPES, ",")[[1]]

# read FILE
df = read.table( FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
#check for expected headers and select only those
exp <- c("tax_level","tax_level_name","type","value","method_reftaxdb")
if ( !all(exp %in% colnames(df)) )  stop( paste("The file",FILE," must have the column headers",paste(exp,collapse=",")) )
df <- subset(df, select = c(exp))

#make sure all TYPES are in FILE
if ( !all(TYPES %in% df$type) )  stop( paste("The file",FILE," must have the column headers",paste(TYPES,collapse=",")) )

#clean data
df <- df[!df$tax_level == "tax_level", ]
df$value <- as.numeric(df$value)

#for each type
for ( type in TYPES ) {
    #select value type
    df_sub <- df[df$type == type, ]

    #plot
    library(ggplot2)
    theme_set(theme_bw())
    p <- ggplot(data = df_sub, aes(x = tax_level, y = value, group = method_reftaxdb))+
        scale_colour_viridis_d( begin = 0, end = 0.9, aesthetics = c("colour", "fill") )+ #end=0.9 prevents yellow, that can be challenging to see
        geom_line( aes(color=method_reftaxdb) )+
        geom_point( aes(color=method_reftaxdb, shape = method_reftaxdb) )+
        xlab("Taxonomic levels") + ylab( type )+
        labs(color="Classifier and database")

    #save plot
    prefix <- sub('\\..*$', '', basename(FILE))

    outfile <- paste0(prefix,".",type,".lineplot.svg")
    print(paste("write",outfile))
    svg(outfile, height = 4, width = 9)
    plot(p)
    invisible(dev.off())

    outfile <- paste0(prefix,".",type,".lineplot.png")
    print(paste("write",outfile))
    png(outfile, height = 400, width = 700)
    plot(p)
    invisible(dev.off())
}
