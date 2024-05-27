#!/usr/bin/env Rscript

# benchmarking_diversity.r file.tsv "tag"

# file.tsv must have headers: match,match_higher,match_higher_nomatch_lower,nomatch,nomatch_higher,not_calculated,method_reftaxdb

#show as horizontal barplot, printing counts on top, colours ranging from blue to red for "Genus & Species expected" (match), "Genus expected, Species not classified" (match_higher), "Genus expected, but Species unexpected" (match_higher_nomatch_lower), "Genus and Species unexpected" (nomatch), "Genus unexpected, Species not classified" (nomatch_higher), but gray for "Not classified" (not_calculated)

# Get params and files from the command line
args <- commandArgs(trailingOnly=TRUE)
FILE <- args[1]
TAG <- args[2]

print(paste("Use file",FILE))
print(paste("with tag",TAG))

prefix <- sub('\\..*$', '', basename(FILE))

df = read.table( FILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
colnames(df)[1] = "ID"

#matrix to long format
df <- as.data.frame.table(`row.names<-`(as.matrix(df)[,-1],df$ID))
print(paste("with",nrow(df),"comparisons"))

#remove all lines where Var1 == Var2 (anyway 0)
df <- df[df$Var1 != df$Var2, ]

#remove suffix "_control"
df$Var1 <- sub('_control', '', df$Var1)
df$Var2 <- sub('_control', '', df$Var2)

#choose lines where Var1 == Var2, since this is the observed vs expected comparison
df <- df[df$Var1 == df$Var2, ]
df$Freq <- as.numeric(df$Freq)
df <- unique(df)
print(paste("with",nrow(df),"comparisons to control"))

#clean-up
df <- subset(df, select = c("Var1","Freq"))
df <- df[order(df$Var1),]
colnames(df) <- c("sample","distance_to_control")
print(paste("of",length(unique(df$sample)),"samples:",paste(unique(df$sample),collapse=",")))

#remove duplicate lines
if ( nrow(df) != length(unique(df$sample)) ) {
    print(paste("WARNING:",length(unique(df$sample)), "samples have",nrow(df),"comparisons"))
    print("Remove all duplicated comparisons by producing mean distances!")
    df <- aggregate(.~sample, df, mean)
}

#add TAG
df$tag <- rep(TAG, nrow(df))

# Write detailed output
outfile <- paste0(prefix,".control.tsv")
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Write mean output
MEAN <- mean( as.numeric(df$distance_to_control) )
df_mean <- data.frame(distance_to_control=MEAN, tag=TAG, stringsAsFactors=FALSE)
outfile <- paste0(prefix,".control_mean.tsv")
print(paste("write",outfile))
write.table(df_mean, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

#plot
library(ggplot2)
theme_set(theme_bw())
p <- ggplot(df, aes(x=tag,y=distance_to_control)) +
    geom_violin(trim=TRUE, fill='#A4A4A4') +
    geom_boxplot(width=0.1) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=1)

#save plot
outfile <- paste0(prefix,".violinplot.svg")
print(paste("write",outfile))
svg(outfile, height = 4, width = 5)
plot(p)
invisible(dev.off())

outfile <- paste0(prefix,".violinplot.png")
print(paste("write",outfile))
png(outfile, height = 400, width = 500)
plot(p)
invisible(dev.off())
