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
  folder <- file.path('results',family_name);dir.create(folder)
  fwrite(r,file.path(folder,'all_element_tests.tsv'),sep='\t')
  fwrite(r[significant==TRUE],file.path(folder,'significant_elements.tsv'),sep='\t')
  summary[[family_name]] <- data.table(family=family_name,elements=nrow(r),significant=sum(r$significant))
}
fwrite(rbindlist(summary),'results/scenario_summary.tsv',sep='\t')
writeLines('Model collection passed; downstream scientific and VCF evidence review remains required.','results/MODELS_COLLECTED')
