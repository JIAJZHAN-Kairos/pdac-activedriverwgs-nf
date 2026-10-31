#!/usr/bin/env Rscript
# Reconstruct the pinned engine's indel counting at element/background boundaries.
# This validates native output; it does not change the fitted model or its inputs.
native_indel_observed_count <- function(element,mg,genome,window=50000L) {
  whole <- ActiveDriverWGS:::.create_background(element,window,genome)
  local <- mg[queryHits(findOverlaps(mg,whole,ignore.strand=TRUE))]
  background <- setdiff(whole,element,ignore.strand=TRUE)
  snv <- local[local$mcols.tag!='indel>X']
  snv_donors <- unique(snv$mcols.patient[queryHits(findOverlaps(snv,element,ignore.strand=TRUE))])
  indel <- local[local$mcols.tag=='indel>X']
  fg <- indel[queryHits(findOverlaps(indel,element,ignore.strand=TRUE))]
  fg <- fg[!duplicated(fg$mcols.patient)]
  bg <- indel[queryHits(findOverlaps(indel,background,ignore.strand=TRUE))]
  combined <- c(fg,bg)
  key <- paste(seqnames(combined),start(combined),end(combined),combined$mcols.patient,sep='::')
  combined <- combined[!duplicated(key)]
  indices <- unique(queryHits(findOverlaps(combined,element,ignore.strand=TRUE)))
  indices <- indices[!combined$mcols.patient[indices] %in% snv_donors]
  length(snv_donors)+length(indices)
}
