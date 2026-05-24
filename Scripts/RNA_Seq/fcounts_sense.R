#!/R/

library(Rsubread)

setwd(getwd())

#======= ARGS INPUT======

#1 ANNOTATION GENOME
#2 BAM LOCATION


args <- commandArgs(trailingOnly = TRUE)

#======= Calling bam files  ======

bam.files <- list.files(path = args[[2]], pattern = ".bam$", full.names = TRUE)

#======= Counting  ======

threads <- 24

fc <- featureCounts(annot.ext = args[[1]],
			files = bam.files,
                       isPairedEnd = TRUE,
		       strandSpecific = 2,
                       allowMultiOverlap = FALSE, # TRUE, # 
                       isGTFAnnotationFile = TRUE,
		       #GTF.featureType = "CDS,rRNA,tRNA",
		       GTF.featureType = "CDS,tRNA",
		       fracOverlap = 0.5 ,
                       #GTF.attrType = "transcript_id",
		       nthreads = threads)

save(fc, file = "feature_counts_meta_ss_2_sense.RData")
