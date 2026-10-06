#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
setDTthreads(2L)
source("code/enhancer/co_occ_functions/enhancer_site_pair_functions.R")
source("code/enhancer/co_occ_functions/enhancer_site_pair_plots.R")
stopifnot(abs(rarefied_pair_probability(2,5,2)-0.7)<1e-12,
  rarefied_pair_probability(0,20,10)==0,
  rarefied_pair_probability(20,20,10)==1)
# Enumerate actual subsets to verify the rarefaction formula.
subsets <- combn(1:5,2)
stopifnot(abs(mean(colSums(subsets<=2)>0)-rarefied_pair_probability(2,5,2))<1e-12)
contingencies <- list(matrix(c(40,0,0,40),2),matrix(c(0,40,40,0),2),
                     matrix(c(20,20,20,20),2),matrix(c(0,2,3,8),2))
for (m in contingencies) stopifnot(identical(cobinding_fisher_p(m),fisher.test(m)$p.value),
                                  identical(cobinding_fisher_p(m),fisher.test(m)$p.value))
pairs <- data.table(enhancer_id=c("a","b","c","d"),enhancer_class=c("active","inactive","active","inactive"),
  site_1_center=100,site_2_center=c(200,300,400,500),site_distance_bp=c(100,200,300,400),
  sample="LPS_0",time_min=0L,n_tf_tf=c(40L,0L,20L,80L),n_tf_naked=c(0L,40L,20L,0L),
  n_naked_tf=c(0L,40L,20L,0L),n_naked_naked=c(40L,0L,20L,0L),total=80L,
  obs_exp=c(2,0,1,1),p_value=c(fisher.test(contingencies[[1]])$p.value,
    fisher.test(contingencies[[2]])$p.value,1,NA_real_),q_value=NA_real_,
  status=c("tested","tested","tested","not_testable"))
pairs <- annotate_cobinding_pairs(pairs)
stopifnot(identical(pairs$cobinding_class,c("enriched","depleted","no significant association","untestable")),
  identical(pairs$q_value[1:3],p.adjust(pairs$p_value[1:3],"BH")),
  all(abs(pairs[,occupancy_both+occupancy_site1_only+occupancy_site2_only+occupancy_neither]-1)<1e-12))
enh <- data.table(enhancer_id=c("a","b","c","d","zero"),chr="chr1",start0=0L,end0=500L,
                  enhancer_class=c("active","inactive","active","inactive","inactive"))
inventory <- data.table(enhancer_id=enh$enhancer_id,n_sites=c(2L,2L,2L,2L,1L),
                        n_possible_pairs=c(1,1,1,1,0))
counts <- data.table(enhancer_id=enh$enhancer_id,sample="LPS_0",time_min=0L,n_reads=80L)
summary <- summarize_enhancer_pairs(pairs,inventory,enh,counts)
stopifnot(nrow(summary)==5L,summary[enhancer_id=="zero",standardized_pairs]==0,
          summary[enhancer_id=="zero",n_cobound_pairs]==0,
          is.na(summary[enhancer_id=="zero",fraction_pairs_cobound]))
# One-sided tests must detect active > inactive, not the reversed comparison.
s <- data.table(enhancer_class=rep(c("inactive","active"),each=30),time_min=0L,
  n_reads=rep(10:39,2),length_bp=rep(300+(1:30)^2,2),n_possible_pairs=100)
s[, standardized_pairs := rep(seq(0.1,1,length.out=30),2)+as.integer(enhancer_class=="active")*3]
s[,`:=`(n_cobound_pairs=standardized_pairs*5,mean_pairs_per_fiber=standardized_pairs,
         fraction_pairs_cobound=standardized_pairs/20)]
comparison <- compare_enhancer_pair_counts(s)
stopifnot(nrow(comparison)==8L,all(comparison$p_value<0.01),all(comparison$estimate>0))
s[,enhancer_class:=ifelse(enhancer_class=="active","inactive","active")]
reverse <- compare_enhancer_pair_counts(s)
stopifnot(all(reverse$p_value>0.99),all(reverse$estimate<0))
# Empty input and constant outcomes must return explicit results, not fake p-values.
stopifnot(nrow(compare_enhancer_pair_counts(s[0]))==0L)
s[,`:=`(standardized_pairs=0,n_cobound_pairs=0,mean_pairs_per_fiber=0,fraction_pairs_cobound=0)]
stopifnot(all(compare_enhancer_pair_counts(s)$p_value==1))
# Render all requested panel types, including zero ratios and missing classes/times.
result <- list(pairs=pairs,enhancers=summary,parameters=list(min_pair_reads=10L))
tmp <- tempfile("cobinding_panels_");dir.create(tmp)
save_enhancer_cobinding_plots(result,compare_enhancer_pair_counts(summary),file.path(tmp,"plots"))
empty <- result; empty$pairs <- empty$pairs[0];empty$enhancers <- empty$enhancers[0]
save_enhancer_cobinding_plots(empty,compare_enhancer_pair_counts(empty$enhancers),file.path(tmp,"empty","plots"))
stopifnot(length(list.files(file.path(tmp,"plots"),pattern="pdf$"))==6L)
cat("PASS: exact Fisher/cache, BH classes, four states, rarefaction, zero-pair enhancers, one-sided comparisons and plot rendering\n")
