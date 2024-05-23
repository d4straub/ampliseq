#!/usr/bin/env Rscript

# benchmark_exact-sequences.r
#
# The script expects the following arguments: expected_seq_file.tsv, calculated_seq_file.tsv

# Get params and files from the command line
args            <- commandArgs(trailingOnly=TRUE)
expFILE         <- args[1] # expected sequences*
resFILE         <- args[2] # calculated sequences*
fbeta           <- args[3] # Fbeta weight, default 2; 1=precision and recall equally weighted (=F1 score), 2=weighs recall higher than precision, 0.5=weighs recall lower than precision.
prefix          <- args[4] # prefix string for output files
id_header       <- args[5] # column header in "resFILE" that contains the first column IDs corresponding to "expFILE"
#*: fasta OR tab separeted file with header: first column with sequences strings, following none or many columns (=samples) with numeric values (=abundance), only presence/absence are used here

# Read params, use defaults if not provided
if ( is.na(fbeta) )         { fbeta <- 2 }
if ( is.na(prefix) )        { prefix <- "" }
if ( is.na(id_header) )     { id_header <- "" }
fbeta <- as.numeric(fbeta)

# Input

# Read stuff
exp = read.table( expFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
res = read.table( resFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)

#In case the two tables are only fasta files, remove lines that start with ">"
if ( length(grep('^>', exp[,1])) > 0 ) {
	exp <- as.data.frame( exp[-grep('^>', exp[,1]),] )
}
if ( length(grep('^>', res[,1])) > 0 ) {
	res <- as.data.frame( res[-grep('^>', res[,1]),] )
}

colnames(exp)[1] <- "seq"
if ( id_header == "" ) {
	colnames(res)[1] <- "seq"
} else {
	colnames(res)[which(colnames(res) == id_header)] <- "seq"
}

# function to produce statistics
get_stats <- function(i_exp,i_res,df,sample) {

	# if none are expected, skip!
	if( length(i_exp)==0 ) {
		print( paste("No expected seq in sample",sample, "- skipping") ); next
	} else {
		print( paste(length(i_exp),"expected seq in sample",sample) )
	}

	# stats
	TP <- intersect(i_exp, i_res)
	FN <- setdiff(i_exp, i_res)
	FP <- setdiff(i_res, i_exp)
	if ( length(TP) > 0 ) {
		Fone <- ( 2 * length(TP)/length(i_res) * length(TP)/length(i_exp) ) / ( length(TP)/length(i_res) + length(TP)/length(i_exp) )
		Fbeta <- ( (1+fbeta^2) * length(TP)/length(i_res) * length(TP)/length(i_exp) ) / ( (fbeta^2) * length(TP)/length(i_res) + length(TP)/length(i_exp) )
	} else {
		Fone <- 0
		Fbeta <- 0
	}

	# save
	types <- c(
		"detected",
		"expected",
		"TP",
		"FN",
		"FP",
		"recall",
		"precision",
		"F1",
		"Fbeta"#,
		#"TP_seq",
		#"FN_seq",
		#"FP_seq"
		)
	values <- c(
		length(i_res),
		length(i_exp),
		length(TP),
		length(FN),
		length(FP),
		length(TP)/length(i_exp),
		length(TP)/length(i_res),
		Fone,
		Fbeta#,
		#paste(TP, collapse=','),
		#paste(FN, collapse=','),
		#paste(FP, collapse=',')
		)
	ids <- rep(sample, length(types))

	df_append <- data.frame(
		sample=ids,
		type= types,
		value= values,
		stringsAsFactors=FALSE)
	df <- rbind( df, df_append )

	return(df)
}

# Output dataframe
df <- data.frame(
	sample=character(),
	type=character(),
	value=character(),
	stringsAsFactors=FALSE)

if( ncol(exp)>=2 && ncol(res)>=2 ) {
	# Find common sample IDs
	exp_samples <- colnames(exp)[2:ncol(exp)]
	res_samples <- colnames(res)[2:ncol(res)]
	samples <- intersect(exp_samples, res_samples)

	if( length(samples)>0 ) {
		print(paste("Found",length(samples),"shared expected and measured samples:"))
		print(paste(samples,collapse=","))
		# for each shared sample
		for (sample in samples) {

			# keep only non-empty columns and sample abundances
			keep_cols <- c("seq",sample)
			s_exp <- subset(exp, select = keep_cols)
			s_res <- subset(res, select = keep_cols)

			# keep only detected (abundance > 0)
			s_exp = s_exp[s_exp[,2] > 0,]
			s_res = s_res[s_res[,2] > 0,]

			df <- get_stats( unique(s_exp$seq), unique(s_res$seq), df, sample)
		}
	} else {
		# if there isnt any shared sample
		print("Found no shared expected and measured samples, compared all sequences.")
		print(paste( "Expected samples:", paste(exp_samples,collapse=",")))
		print(paste( "Measured samples:", paste(res_samples,collapse=",")))
		df <- get_stats( unique(exp$seq), unique(res$seq), df, "unknown")
	}
} else {
	print("Found no abundance information in at least one file, compared all sequences.")
	df <- get_stats( unique(exp$seq), unique(res$seq), df, "unknown")
}

# Write detailed output
outfile <- paste0(prefix,"benchmarking_stats_long.tsv")
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Make & write summary output
df_sum <- data.frame(type=character(), mean=character(), stringsAsFactors=FALSE)
outfile <- paste0(prefix,"benchmarking_stats_mean.tsv")
print("Warnings are expected at that point:")
for (type in unique(df$type)) {
	df_subset <- df[df$type == type, ]
	MEAN <- mean( as.numeric(df_subset$value) )
	df_sum[nrow(df_sum) + 1,] <- list(type, MEAN)
}
print(paste("write",outfile))
write.table(df_sum, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Plot Values
df_subset <- subset(df, type %in% c("recall","precision","F1","Fbeta") )
df_subset$value <- as.numeric(df_subset$value)

outfile <- paste0(prefix,"benchmark_boxplot.svg")
print(paste("write",outfile))
svg(outfile, height = 4, width = 5)
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1))
invisible(dev.off())

outfile <- paste0(prefix,"benchmark_boxplot.png")
print(paste("write",outfile))
png(outfile, height = 400, width = 500)
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1))
invisible(dev.off())

# https://help.displayr.com/hc/en-us/articles/360003262056-How-to-Create-a-Density-Plot-Using-R
df_subset <- subset(df, type == "F1")
df_subset$value <- as.numeric(df_subset$value)
if ( nrow(df_subset) > 1 ) {
	svg(paste0(prefix,"benchmark_density.svg"), height = 4, width = 5)
	plot(density(df_subset$value, from = 0, to = 1, adjust = 0.2), #plot values 0-100 and adjust bandwidth
        main = "F1 score density plot", #this is the title of the plot
        xlab = "F1 score") #this is the title of the x-axis
	invisible(dev.off())
}
