#!/usr/bin/env Rscript

# benchmarking_taxonomy_persequence.r
#
# The script expects the following arguments: expected_taxa_per_sequence.tsv rel-table-ASV_with-*-tax.tsv [from results/qiime2/rel_abundance_tables]

# Get params and files from the command line
args             <- commandArgs(trailingOnly=TRUE)
expFILE          <- args[1] # expected taxa*
resFILE          <- args[2] # calculated taxa*
prefix           <- args[3] # prefix string for output files
id_header        <- args[4] # column header in "resFILE" that contains the first column IDs corresponding to "expFILE"
rm_after_string  <- args[5] # string after that all taxa should be truncated
rm_before_string <- args[6] # string before that the taxonomic string will be deleted
fbeta            <- args[7] # Fbeta weight, default 2; 1=precision and recall equally weighted (=F1 score), 2=weighs recall higher than precision, 0.5=weighs recall lower than precision.
#*: first column with sequences/ids [or $resFILE at $id_header], following none or many columns (=taxonomic levels) with strings (=taxonomy, empty if none)

# Read params, use defaults if not provided
if ( is.na(rm_after_string) )  { rm_after_string <- "" }
if ( is.na(rm_before_string) ) { rm_before_string <- "" }
if ( is.na(fbeta) )          { fbeta <- 2 }
if ( is.na(prefix) )         { prefix <- "" }
if ( is.na(id_header) )      { id_header <- "" }
fbeta <- as.numeric(fbeta)

# Input

# Read stuff
exp = read.table( expFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE)
res = read.table( resFILE, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, strip.white = TRUE, quote = "\"")
res[is.na(res)] <- ""

# make sure that exp has "seq" as first column
if ( colnames(exp)[1] != "seq" ) stop( paste("The file",expFILE," must have 'seq' as first column header.") )

#colnames(exp)[1] <- "seq"
if ( id_header == "" ) {
	colnames(res)[1] <- "seq"
} else {
    if ( "seq" %in% colnames(res) ) stop( paste("The file",resFILE,"contains already the header 'seq', but",id_header,"is supposed to be renamed.") )
	colnames(res)[which(colnames(res) == id_header)] <- "seq"
}

# Keep only columns that intersect between exp and res
print( paste0( "Header in ",expFILE,": ", paste( colnames(exp), collapse=",") ) )
print( paste0( "Header in ",resFILE,": ", paste( colnames(res), collapse=",") ) )
exp_levels <- intersect(colnames(exp), colnames(res))
print( paste("Matched headers:", paste(exp_levels,collapse=",") ) )
if ( length(exp_levels) < 3 ) stop( paste("The file provided by --benchmarking_taxonomy must have at least two column headers matching:",paste(colnames(res),collapse=",")) )
res <- subset(res, select = c(exp_levels))
exp <- subset(exp, select = c(exp_levels))

# retain only everything before string (here: space [" "])! because SILVA uses e.g. "Nodosilinea PCC-7104" instead of "Nodosilinea PCC-7104" or "Phormidesmis ANT.L52.6", "Oxynema BDU 92071", "Roseofilum AO1-A", etc.
if ( rm_after_string != "" ) {
	res <- as.data.frame( sapply(res, function(x) sub( paste0(rm_after_string,".*"), "", x) ) )
}
if ( rm_before_string != "" ) {
	res <- as.data.frame( sapply(res, function(x) sub( paste0(".*",rm_before_string), "", x) ) )
}

# Figure out what sequences are identical & merge
print( paste("Found", length(exp$seq),"expected sequences in",expFILE ) )
print( paste("Found", length(res$seq),"calculated sequences in",resFILE ) )
df <- merge(exp, res, by="seq", all=FALSE)
print( paste("Found", nrow(df),"expected sequences to be calculated" ) )

# Go through levels
df_result <- subset(df, select = "seq")
for ( i in 2:length(exp_levels) ) {
	header_col=exp_levels[i]
	exp_col=i
	res_col=i+length(exp_levels)-1
	print(paste("Compare",colnames(df)[exp_col],"and",colnames(df)[res_col],"as",header_col))
	#print( paste(colnames(df)[2:exp_col]) )
	#print( paste(colnames(df)[(length(exp_levels)+1):res_col]) )
	if (exp_col>2) {
		#expected: df[,2:exp_col]
		exp_list <- do.call(paste, c( df[,2:exp_col], sep=";"))
		#calculated: df[,(length(exp_levels)+1):res_col]
		res_list <- do.call(paste, c( df[,(length(exp_levels)+1):res_col], sep=";"))
	} else {
		exp_list <- df[,2:exp_col]
		res_list <- df[,(length(exp_levels)+1):res_col]
	}
	# iterate through entries
	list_result = c()
	for ( j in 1:length(exp_list) ) {
		e = exp_list[j]
		r = res_list[j]
		x = "no_match" #this is the final result, by default no match
		#if expected is empty
		if ( grepl( ".*;$", e) ) { x = "not_expected"
		} else if ( e == "" ) { x = "not_expected"
		#if calculated is empty
		} else if ( grepl( ".*;$", r) ) { x = "not_calculated"
		} else if ( r == "" ) { x = "not_calculated"
		#if match
		} else if ( e == r ) { x = "match" }
		list_result <- c(list_result,x)
	}
	df_result <- cbind(df_result,list_result )
}
colnames(df_result) <- exp_levels

df <- merge(df_result, df, by="seq", all=TRUE)
# write df to get complete results table
outfile <- paste0(prefix,"_benchmarking_taxonomy_persequence.tsv")
print(paste("write",outfile))
write.table(df, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Make summary
df_summary <- data.frame(
    match = numeric(),
    match_higher = numeric(),
    match_higher_nomatch_lower = numeric(),
    nomatch = numeric(),
    nomatch_higher = numeric(),
    not_calculated = numeric(),
    method_reftaxdb = character(),
    taxa = character(),
    stringsAsFactors=FALSE)
if ( "Genus" %in% colnames(df_result) && "Species" %in% colnames(df_result) ) {
    print("Calculate Genus & Species")
    df_sub = subset(df_result, select = c("Genus","Species") )
    df_table <- as.data.frame(table(df_sub))
    # higher = Genus, lower = Species

    summary <- list(
        match = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'match'),]$Freq, 0),
        match_higher = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        match_higher_nomatch_lower = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch_higher = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        not_calculated = ifelse( nrow(df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),]$Freq, 0),
        method_reftaxdb = paste0(prefix,"_species"),
        taxa = paste(colnames(df_sub),collapse=",")
    )
    df_summary <- rbind(df_summary,summary)
}
if ( "Genus" %in% colnames(df_result) && "Species_exact" %in% colnames(df_result) ) {
    print("Calculate Genus & Species_exact")
    df_sub = subset(df_result, select = c("Genus","Species_exact") )
    df_table <- as.data.frame(table(df_sub))
    # higher = Genus, lower = Species_exact

    summary <- list(
        match = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'match'),]$Freq, 0),
        match_higher = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        match_higher_nomatch_lower = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch_higher = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        not_calculated = ifelse( nrow(df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),]$Freq, 0),
        method_reftaxdb = paste0(prefix,"_speciesexact"),
        taxa = paste(colnames(df_sub),collapse=",")
    )
    df_summary <- rbind(df_summary,summary)
}
if ( nrow(df_summary) == 0 ) {
    print("WARNING: Will extract last two taxonomic levels only!")
    df_sub = df_result[,(ncol(df_result)-1):(ncol(df_result))]
    print(paste("Reduced to",paste(colnames(df_sub),collapse=",")))
    df_table <- as.data.frame(table(df_sub))

    summary <- list(
        match = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'match'),]$Freq, 0),
        match_higher = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        match_higher_nomatch_lower = ifelse( nrow(df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'no_match'),]$Freq, 0),
        nomatch_higher = ifelse( nrow(df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'no_match' & df_table[2] == 'not_calculated'),]$Freq, 0),
        not_calculated = ifelse( nrow(df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),])>0, df_table[(df_table[1] == 'not_calculated' & df_table[2] == 'not_calculated'),]$Freq, 0),
        method_reftaxdb = prefix,
        taxa = paste(colnames(df_sub),collapse=",")
    )
    #add a column with relative numbers
    #df_summary <- rbind(df_summary,df_summary[1,]/rowSums(df_summary)*100)
    #df_summary$type <- c("count","percent")
    df_summary <- rbind(df_summary,summary)
}

# write df_summary to get summary results table
outfile <- paste0(prefix,"_benchmarking_taxonomy_persequence_summary.tsv")
print(paste("write",outfile))
write.table(df_summary, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Produce recall/precision/F-score
df_stats <- data.frame(
	tax_level=character(),
    tax_level_name=character(),
	type=character(),
	value=character(),
	stringsAsFactors=FALSE)
for ( i in 2:ncol(df_result) ) {
    tax_level <- i-1
    tax_level_name <- colnames(df_result)[i]
    columni <- df_result[,i]
	# stats
	TP <- length(which(columni == "match"))
	FN <- length(which(columni == "not_calculated"))
	FP <- length(which(columni == "no_match"))
    nres <- length(which(columni != "not_calculated"))
    nexp <- length(columni)
	if ( TP > 0 ) {
		Fone <- ( 2 * TP/nres * TP/nexp ) / ( TP/nres + TP/nexp )
		Fbeta <- ( (1+fbeta^2) * TP/nres * TP/nexp ) / ( (fbeta^2) * TP/nres + TP/nexp )
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
		"Fbeta"
		)
	values <- c(
		nres,
		nexp,
		TP,
		FN,
		FP,
		TP/nexp,
		TP/nres,
		Fone,
		Fbeta
		)
    tax_level <- rep(tax_level, length(types))
	tax_level_name <- rep(tax_level_name, length(types))

	df_append <- data.frame(
		tax_level=tax_level,
        tax_level_name=tax_level_name,
		type= types,
		value= values,
		stringsAsFactors=FALSE)
	df_stats <- rbind( df_stats, df_append )
}
df_stats$method_reftaxdb <- rep(prefix, nrow(df_stats))

# Write detailed output
outfile <- paste0(prefix,"_benchmarking_taxonomy_persequence_stats_long.tsv")
print(paste("write",outfile))
write.table(df_stats, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Make & write summary output
df_sum <- data.frame(type=character(), mean=character(), stringsAsFactors=FALSE)
outfile <- paste0(prefix,"_benchmarking_taxonomy_persequence_stats_mean.tsv")
for (type in unique(df_stats$type)) {
	df_subset <- df_stats[df_stats$type == type, ]
	MEAN <- mean( as.numeric(df_subset$value) )
	df_sum[nrow(df_sum) + 1,] <- list(type, MEAN)
}
print(paste("write",outfile))
write.table(df_sum, file = outfile, row.names = FALSE, col.names = TRUE, quote = FALSE, na = '', sep="\t")

# Plot Values
df_subset <- subset(df_stats, type %in% c("recall","precision","F1","Fbeta") )
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"_benchmarking_taxonomy_persequence_boxplot.svg"), height = 4, width = 5)
boxplot(value~type, data=df_subset, xlab="Type", ylab="Value", ylim = c(0, 1))
invisible(dev.off())

# Bad looking plot, but potentially helpful
df_subset <- subset(df_stats, type == "F1")
df_subset$value <- as.numeric(df_subset$value)
svg(paste0(prefix,"_benchmarking_taxonomy_persequence_F1.svg"), height = 4, width = 10)
barplot(df_subset$value,
    names.arg = df_subset$tax_level_name,
    ylab = "F1 score",
    xlab = "Taxonomic level")
invisible(dev.off())

