#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ActiveDriverWGS);library(data.table)})
setDTthreads(1)
a <- commandArgs(TRUE)
files <- list.files(a[1],pattern='\\.tsv$',full.names=TRUE)
tests <- rbindlist(c(lapply(files,fread),list(fread(a[4]))),use.names=TRUE,fill=TRUE)
universe <- fread(a[2]);membership <- fread(a[3]);expected <- as.integer(a[5])
expected_ids <- readLines(a[6])
stopifnot(nrow(tests)==expected,!anyDuplicated(tests$id),all(tests$id %in% universe$id),
          length(expected_ids)==expected,!anyDuplicated(expected_ids),setequal(tests$id,expected_ids))
# Use the pinned engine's formatting rules, including its P=0 floor.
tests <- as.data.table(ActiveDriverWGS:::.fix_all_results(as.data.frame(tests)))
tests[is.na(pp_element),pp_element:=1]
all <- merge(data.table(id=unique(universe$id)),tests,by='id',all.x=TRUE)
all[is.na(pp_element),pp_element:=1]
stopifnot(!anyDuplicated(all$id),all(is.finite(all$pp_element)),all(all$pp_element>0 & all$pp_element<=1))
bounds <- universe[,.(chr=chr[1],start=min(start),end=max(end),effective_bp=sum(end-start)),by=id]
all <- merge(all,bounds,by='id')
counts <- fread(a[7])
stopifnot(!anyDuplicated(counts$id),setequal(counts$id,expected_ids))
all <- merge(all,counts,by='id',all.x=TRUE)
for(col in c('mutated_donors','SNV_donors','indel_donors','HMF_donors','APGI_donors','internal_donors'))
  set(all,which(is.na(all[[col]])),col,0L)
stopifnot(all(all$element_muts_obs[!is.na(all$element_muts_obs)]==all$mutated_donors[!is.na(all$element_muts_obs)]),
          all(!is.na(all$element_muts_obs[all$mutated_donors>0L])))
dir.create('results',showWarnings=FALSE);fwrite(all,'results/union_raw_tests.tsv',sep='\t')
summary <- list()
for(family_name in unique(membership$family)) {
  ids <- unique(membership[family==family_name]$id);r <- all[id %in% ids]
  stopifnot(nrow(r)==length(ids))
  r[,fdr_element:=p.adjust(pp_element,'BH',n=length(ids))]
  r[,significant:=!is.na(element_enriched) & element_enriched & fdr_element<.05]
  annotation <- membership[family==family_name,.(classes=paste(sort(unique(class)),collapse=';'),
                  directly_annotated_genes=paste(sort(unique(gene[nzchar(gene)])),collapse=';'),
                  sources=paste(sort(unique(source)),collapse=';')),by=id]
  r <- merge(r,annotation,by='id');setorder(r,fdr_element,pp_element)
  r[,`:=`(analysis_donors=373L,observed_expected_ratio=element_muts_obs/element_muts_exp)]
  folder <- file.path('results',family_name);dir.create(folder)
  fwrite(r,file.path(folder,'all_element_tests.tsv'),sep='\t')
  fwrite(r[significant==TRUE],file.path(folder,'significant_elements.tsv'),sep='\t')
  pdf(file.path(folder,'diagnostics.pdf'),width=10,height=5)
  par(mfrow=c(1,2));p <- sort(r$pp_element)
  plot(-log10(ppoints(length(p))),-log10(p),pch=16,cex=.3,xlab='Expected -log10(P)',ylab='Observed -log10(P)',main=family_name)
  abline(0,1,col='red');hist(r$mutated_donors,breaks=40,xlab='Independent mutated donors',main='Region recurrence');dev.off()
  summary[[family_name]] <- data.table(family=family_name,elements=nrow(r),significant=sum(r$significant))
}
fwrite(rbindlist(summary),'results/scenario_summary.tsv',sep='\t')
writeLines('Model collection passed; downstream scientific and VCF evidence review remains required.','results/MODELS_COLLECTED')
