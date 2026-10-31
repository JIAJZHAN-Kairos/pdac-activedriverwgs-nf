suppressPackageStartupMessages(library(data.table))
script <- normalizePath('bin/collect_tests.R')
root <- tempfile('adwgs-collector-');dir.create(root);dir.create(file.path(root,'chunks'))
row <- function(id,p,obs,exp,enriched) data.table(id=id,pp_element=p,element_muts_obs=obs,
  element_muts_exp=exp,element_enriched=enriched,pp_site=NA_real_,site_muts_obs=NA_real_,
  site_muts_exp=NA_real_,site_enriched=NA)
fwrite(row('synthetic1',.001,3,1,TRUE),file.path(root,'chunks/block.tsv'),sep='\t')
fwrite(row('synthetic2',.9,1,2,FALSE),file.path(root,'seed.tsv'),sep='\t')
fwrite(data.table(chr='chr1',start=c(100,300,500),end=c(200,400,600),id=paste0('synthetic',1:3)),
       file.path(root,'universe.tsv'),sep='\t')
fwrite(data.table(id=c('synthetic1','synthetic2','synthetic3','synthetic1','synthetic3'),
                 family=c('A','A','A','B','B'),class='synthetic',gene='',source='test'),
       file.path(root,'membership.tsv'),sep='\t')
writeLines(c('synthetic1','synthetic2'),file.path(root,'expected.txt'))
counts <- data.table(id=c('synthetic1','synthetic2'),mutated_donors=c(3,1),SNV_donors=c(2,1),
                     indel_donors=c(2,0),HMF_donors=c(2,0),APGI_donors=c(1,0),internal_donors=c(0,1))
fwrite(counts,file.path(root,'counts.tsv'),sep='\t')
setwd(root)
run <- function() suppressWarnings(system2(file.path(R.home('bin'),'Rscript'),
   vapply(c(script,'chunks','universe.tsv','membership.tsv','seed.tsv','2','expected.txt','counts.tsv'),shQuote,''),
   stdout='test.log',stderr='test.log'))
stopifnot(run()==0L)
a <- fread('results/A/all_element_tests.tsv');b <- fread('results/B/all_element_tests.tsv')
stopifnot(nrow(a)==3L,a[id=='synthetic3']$pp_element==1,
          abs(a[id=='synthetic1']$fdr_element-.003)<1e-12,
          abs(b[id=='synthetic1']$fdr_element-.002)<1e-12,
          a[id=='synthetic1']$mutated_donors==3)
file.copy('chunks/block.tsv','chunks/duplicate.tsv');stopifnot(run()!=0L);unlink('chunks/duplicate.tsv')
file.rename('chunks/block.tsv','missing.tsv');stopifnot(run()!=0L);file.rename('missing.tsv','chunks/block.tsv')
counts[id=='synthetic1',mutated_donors:=4L];fwrite(counts,'counts.tsv',sep='\t');stopifnot(run()!=0L)
cat('PASS: silent elements and family-wide BH; duplicate/missing tests and donor mismatches rejected.\n')
