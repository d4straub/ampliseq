#!/usr/bin/env Rscript

# benchmarking_taxonomy_plot_summary.r file.tsv

# file.tsv must have headers: match,match_higher,match_higher_nomatch_lower,nomatch,nomatch_higher,not_calculated,method_reftaxdb

#show as horizontal barplot, printing counts on top, colours ranging from blue to red for "Genus & Species expected" (match), "Genus expected, Species not classified" (match_higher), "Genus expected, but Species unexpected" (match_higher_nomatch_lower), "Genus and Species unexpected" (nomatch), "Genus unexpected, Species not classified" (nomatch_higher), but gray for "Not classified" (not_calculated)

# Get params and files from the command line
args <- commandArgs(trailingOnly=TRUE)
FILE <- args[1]

df = read.table( FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
#check for expected headers and select only those
exp <- c("match","match_higher","match_higher_nomatch_lower","nomatch","nomatch_higher","not_calculated","method_reftaxdb")
if ( !all(exp %in% colnames(df)) )  stop( paste("The file",FILE," must have the column headers",paste(exp,collapse=",")) )
df <- subset(df, select = c(exp))

df <- df[!df$match == "match", ]
#preserve sequence for classifier & databases
df[, c(1:6)] <- sapply(df[, c(1:6)], as.numeric)
df_ind_factor <- factor(df$method_reftaxdb, levels=df[order(df$match,df$match_higher,df$match_higher_nomatch_lower,df$not_calculated),"method_reftaxdb"])

#wide to long format
df <- cbind(df[7], stack(lapply(df[-c(7)], as.character)))
df$values <- as.numeric(df$values)
#sort categories
df$ind <- factor(df$ind, levels=rev( c("match","match_higher","match_higher_nomatch_lower","not_calculated","nomatch_higher","nomatch") ))
#sort classifier & databases
df$method_reftaxdb <- df_ind_factor

#color
#color_list <- list(match = "navy", match_higher="mediumblue",match_higher_nomatch_lower="magenta", nomatch_higher="orangered", nomatch = "red", not_calculated = "gray")
color_list <- list(match = "#440154FF", match_higher="#39568CFF",match_higher_nomatch_lower="#1F968BFF", nomatch_higher="#B8DE29FF", nomatch = "#FDE725FF", not_calculated = "gray")
#rename
label_list <- list(
	match = "Genus & Species correct",
	match_higher="Genus correct, Species not classified",
	match_higher_nomatch_lower="Genus correct, but Species incorrect",
	nomatch_higher="Genus incorrect, Species not classified",
	nomatch = "Genus and Species incorrect",
	not_calculated = "Not classified")
#plot
library(ggplot2)
theme_set(theme_bw())
p <- ggplot(df, aes(x = method_reftaxdb, y = values))+
    xlab("Classifier and database") + ylab("Sequences")+
    labs(fill = "Classified to expected")+
    geom_col(aes(fill = ind), width = 0.7)+
    coord_flip()+
    scale_fill_manual(labels = label_list, values=color_list)
    #scale_y_continuous(labels = scales::percent) #only if first into proportions!

#save plot
prefix <- sub('\\..*$', '', basename(FILE))

outfile <- paste0(prefix,".barplot.svg")
print(paste("write",outfile))
svg(outfile, height = 4, width = 9)
plot(p)
invisible(dev.off())

outfile <- paste0(prefix,".barplot.png")
print(paste("write",outfile))
png(outfile, height = 400, width = 700)
plot(p)
invisible(dev.off())
