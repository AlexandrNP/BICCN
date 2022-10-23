if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")

if (!require("stringr", quietly = TRUE))
	install.packages("stringr")

if (!require("hash", quietly = TRUE))
	install.packages("hash", repos = "http://cran.us.r-project.org")

if (!require("foreach", quietly = TRUE))
	install.packages("foreach", repos = "http://cran.us.r-project.org")

if (!require("doParallel", quietly = TRUE))
	install.packages("doParallel", repos = "http://cran.us.r-project.org")

if (!require("ensembledb", quietly = TRUE))
	BiocManager::install("ensembldb")

if (!require("EnsDb.Mmusculus.v79", quietly = TRUE))
	BiocManager::install("EnsDb.Mmusculus.v79")


#if (!require("dplyr", quietly = TRUE))
#	install.packages("dplyr")

library("hash")
library("stringr")
library("ensembldb")
library(foreach)
library(doParallel)
library("EnsDb.Mmusculus.v79")
#library("dplyr")

get_evex_data <- function() 
{
	df <- read.table('EVEx_4axes', sep='\t', header=TRUE)
	#print(colnames(df))
	return(df)
}

get_highly_variable <- function()
{
	return(read.table('highlyVariableExons_groupedByVariability', sep='\t', header=TRUE))
}

get_low_variable <- function()
{
	return(read.table('lowVariabilityAltExons_nullSet', sep='\t', header=TRUE))
}

get_developmental_modalities <- function()
{
	return(read.table('exonList_developmentalModalities', sep='\t', header=TRUE))
}

get_cons_set <- function()
{
	return(read.table('consExons_nullSet', sep='\t', header=TRUE))
}


evex_data <- get_evex_data()
highly_variable <- get_highly_variable()
low_variable <- get_low_variable()
developmental_modalities <- get_developmental_modalities()
conserved_set <- get_cons_set()
test_set <- read.table('test_set', sep='\t', header=TRUE)
#datasets <- list(evex_data, highly_variable, low_variable, developmental_modalities, conserved_set)
datasets <- list(conserved_set)
#datasets <- list(evex_data, highly_variable, low_variable, developmental_modalities)
#datasets <- list(test_set)
db <- EnsDb.Mmusculus.v79
filtered_dbs <- hash()

cores = detectCores()

process_exon <- function(exon)
{
		#print(exon)
		coords <- strsplit(exon, '_')
		#print(coords)
		chr <- str_replace(coords[[1]][1], "chr", "")
		#chr <- unlist(strsplit(coords[[1]][1], split="chr", fixed=TRUE))[2]
		#chr <- coords[[1]][1]
		start_coord <- strtoi(coords[[1]][2])
		end_coord <- strtoi(coords[[1]][3])
		strand <- coords[[1]][4]

		range <- IRanges(start=start_coord, end=end_coord)
		genome <- GRanges(seqnames=chr, 
							ranges=range,
							strand=strand)
		#if (!(chr %in% names(filtered_dbs)))
		#{
		#	#aaa <- filter(EnsDb.Mmusculus.v79, filter = ~ seq_name == chr)
		#	filtered_dbs[[chr]] <- filter(EnsDb.Mmusculus.v79, filter = ~ seq_name == chr)
		#}

		#print(chr)
		#print(start_coord)
		#print(end_coord)
		#filtered_db <- filter(EnsDb.Mmusculus.v79, filter = ~ seq_name == chr)
		#genome_prot <- genomeToProtein(genome, db)
		#filtered_db <- filter(EnsDb.Mmusculus.v79, filter = ~ seq_name == chr)
		#genome_prot <- genomeToProtein(genome, filtered_db)
		genome_prot <- genomeToProtein(genome, EnsDb.Mmusculus.v79)
		#genome_prot <- genomeToProtein(genome, filtered_dbs[[chr]])
		#print('Proteins!!!')
		#print(genome_prot)
		#print(genome_prot[3])
		#print(names(genome_prot[[1]]))
		print('catching...')
		return(tryCatch(
			{
				#print(names(genome_prot[[1]]))
				prots <- proteins(EnsDb.Mmusculus.v79, filter = ProteinIdFilter(names(genome_prot[[1]])))
				#proteins(EnsDb.Mmusculus.v79, filter = ProteinIdFilter(names(genome_prot[[1]])))
				#print(proteins)
				#return(proteins)
				prots$Exon <- rep(exon, length(prots$tx_id))
				print('Returning')
				return(prots)
			},
			error = function(err)
			{
				#print(proteins)
				#print(names(proteins))
				if(!is.null(names(proteins)))
				{
					return(proteins)
				}
				tx_id <- c("NaN")
				protein_id <- c("NaN")
				protein_sequence <- c("NaN")
				df_nan <- data.frame(tx_id, protein_id, protein_sequence)
				df_nan$Exon <- rep(exon, length(df_nan$tx_id))
				print('Failing')
				return(df_nan)
			}))
		#print(prt)
		#print(length(prt$tx_id))
		#print('Before')
		#print(prt)
		#prt$Exon <- rep(exon, length(prt$tx_id))
		#print('After')
		#print(prt)
		#return(prt)
}


exon_df_list = list()
for (df in datasets)
{
	#cluster <- makeCluster(cores[1]-10)
	cluster <- makeCluster(cores[1]-20)
	registerDoParallel(cluster)

	df_sequences <- foreach(i=1:length(df$Exon), .combine=rbind, .errorhandling = 'remove', 
							.packages=c("ensembldb","EnsDb.Mmusculus.v79", "stringr")) %dopar%
	{
		#print(i)
		prt <- process_exon(df$Exon[i])
		#print('Proteins')
		#print(prt)
		#print('Merged proteins')
		#print(df_sequences)
		#print(warnings())
		prt
		#exon_df_list <- append(exon_df_list, list(prt))
		#print(length(exon_df_list))
	}
	df_sequences
	stopCluster(cluster)
	print('DF sequences:')
	print(df_sequences)
	exon_df_list <- append(exon_df_list, list(df_sequences))
	#if (length(exon_df_list) > 1)
	#{
		combined_df <- do.call("rbind", exon_df_list)
	#}
	write.table(combined_df, 'protein_seq_for_exons_null.tsv', row.names = FALSE, quote=FALSE, sep='\t')
}
