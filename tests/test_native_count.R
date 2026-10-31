suppressPackageStartupMessages({library(ActiveDriverWGS);library(BSgenome.Hsapiens.UCSC.hg38);library(GenomicRanges)})
source('bin/validate_observed_counts.R')
g <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
element <- GRanges('chr1',IRanges(1000101L,1000200L),mcols='synthetic')
for(crosses_boundary in c(FALSE,TRUE)) {
  indels <- GRanges('chr1',IRanges(c(1000180L,1000190L),c(1000181L,if(crosses_boundary) 1000210L else 1000191L)),
                    mcols=data.frame(patient='synthetic_donor',tag='indel>X'))
  native <- ADWGS_test('synthetic',element,GRanges(),indels,50000L,g,FALSE)
  expected <- if(crosses_boundary) 2L else 1L
  stopifnot(native$element_muts_obs==expected,
            native_indel_observed_count(element,indels,g)==expected,
            length(unique(indels$mcols.patient))==1L)
}
cat('PASS: unchanged engine reproduces boundary indel excess while independent donor count stays one.\n')
