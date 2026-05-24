#!/R/

library(Rsubread)

setwd(getwd())

#======= ARGS INPUT======

#1 REFERENCE GENOME
#2 RAW DATA LOCATION

args <- commandArgs(trailingOnly = TRUE)

#======= INDECES ======

index="genome_index"
#buildindex(basename = index,
 #          reference = args[[1]])

#======= MAPPING  ======


fq1 <- list.files(args[[2]],
	  	  pattern = "1P.fq$", full.names = TRUE)
fq2 <- sub("1P.fq$", "2P.fq", fq1)
bam <- sub("1P.fq$", ".bam", basename(fq1))
threads <- 24


align(index = index,
	readfile1 = fq1, readfile2 = fq2, type = 0,
	output_file = bam, nthreads = threads)
