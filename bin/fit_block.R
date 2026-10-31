#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ActiveDriverWGS);library(BSgenome.Hsapiens.UCSC.hg38);library(GenomicRanges);library(data.table);library(jsonlite)})
setDTthreads(1)
stopifnot(as.character(packageVersion('ActiveDriverWGS'))=='1.2.1',as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38'))=='1.4.4')
args <- commandArgs(TRUE);mut <- readRDS(args[1]);el <- fread(args[2]);ids <- readLines(args[3]);label <- args[4]
stopifnot(length(ids)>0L,!anyDuplicated(ids),all(ids %in% el$id),length(unique(mut$chr))==1L,
          setequal(unique(mut$chr),unique(el$chr)),max(table(mut$patient))<=90000L)
g <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
mut <- format_muts(mut,this_genome=g,filter_hyper_MB=30)
mg <- GRanges(mut$chr,IRanges(mut$pos1,mut$pos2),mcols=mut[,c('patient','tag')])
eg <- GRanges(el$chr,IRanges(el$start+1L,el$end),mcols=el$id)
parts <- split(eg,eg$mcols)
started <- proc.time()[3]
results <- lapply(ids,function(id) {
  # Official ADWGS_test with one complete element and all chromosome mutations.
  element <- parts[[id]]
  r <- ADWGS_test(id,element,GRanges(),mg,50000,g,FALSE)
  independent_donors <- length(unique(mg$mcols.patient[queryHits(findOverlaps(mg,element,ignore.strand=TRUE))]))
  stopifnot(independent_donors>0L,!is.na(r$element_muts_obs),r$element_muts_obs==independent_donors)
  r
})
result <- do.call(rbind,results)
stopifnot(nrow(result)==length(ids),setequal(result$id,ids),!anyDuplicated(result$id))
fwrite(result,paste0(label,'.tsv'),sep='\t')
write_json(list(block=label,rows=nrow(result),elapsed_seconds=proc.time()[3]-started,donors=length(unique(mut$patient)),
                genome='hg38',engine_version=as.character(packageVersion('ActiveDriverWGS')),
                hg38_package=as.character(packageVersion('BSgenome.Hsapiens.UCSC.hg38')),session=capture.output(sessionInfo())),
           paste0(label,'.qc.json'),auto_unbox=TRUE,pretty=TRUE)
