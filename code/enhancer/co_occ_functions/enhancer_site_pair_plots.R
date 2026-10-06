# Descriptive pair-level panels; class inference uses enhancer-level summaries.
save_enhancer_cobinding_plots <- function(result, comparisons, plot_dir) {
  stopifnot(requireNamespace("hexbin",quietly=TRUE),requireNamespace("ggplot2",quietly=TRUE))
  dir.create(plot_dir,recursive=TRUE,showWarnings=FALSE)
  class_colors <- c(active="#D55E00",inactive="#0072B2")
  binding_levels <- c("depleted","no significant association","enriched")
  binding_colors <- c(depleted="#0072B2",`no significant association`="grey65",enriched="#D55E00")
  all_pairs <- data.table::copy(result$pairs)
  enh <- data.table::copy(result$enhancers)
  times <- sort(unique(enh$time_min))
  if (!length(times)) times <- c(0,5,10,15)
  for (dt in list(all_pairs,enh)) {
    dt[,enhancer_class:=factor(enhancer_class,levels=c("active","inactive"))]
    dt[,timepoint:=factor(time_min,levels=times,labels=paste0(times," min"))]
  }
  pair_data <- all_pairs[status=="tested" & is.finite(obs_exp) & is.finite(p_value)]
  pair_data[,`:=`(neg_log10_p=-log10(pmax(p_value,.Machine$double.xmin)),
                  cobinding_class=factor(cobinding_class,levels=binding_levels))]
  scaffold <- data.table::CJ(enhancer_class=factor(c("active","inactive"),levels=c("active","inactive")),
    timepoint=factor(paste0(times," min"),levels=paste0(times," min")))
  theme_common <- ggplot2::theme_classic(base_size=10) +
    ggplot2::theme(legend.position="bottom",strip.background=ggplot2::element_rect(fill="grey95"))
  save <- function(p,name,width=13,height=7) {
    ggplot2::ggsave(file.path(plot_dir,name),p,width=width,height=height,limitsize=FALSE)
    invisible(p)
  }
  a <- ggplot2::ggplot(pair_data,ggplot2::aes(obs_exp,neg_log10_p)) +
    ggplot2::geom_blank(data=scaffold,ggplot2::aes(x=0,y=0),inherit.aes=FALSE) +
    ggplot2::geom_hex(bins=45,na.rm=TRUE) +
    ggplot2::geom_hline(yintercept=2,colour="red",linetype="dotted") +
    ggplot2::geom_vline(xintercept=1,colour="grey50",linetype="dashed") +
    ggplot2::scale_x_continuous(trans="log1p",breaks=c(0,0.5,1,2,5,10,50,100,500)) +
    ggplot2::scale_fill_viridis_c(trans="log10",name="Site pairs") +
    ggplot2::facet_grid(enhancer_class~timepoint,drop=FALSE) + theme_common +
    ggplot2::labs(title="A  Normalized co-binding and pair-association significance",
      subtitle="All tested pairs; red dotted line: raw p = 0.01; x-axis uses log1p spacing",
      x="Normalized co-binding (observed / independent expected)",y="-log10(two-sided Fisher p)")
  save(a,"05_A_cobinding_hexbin.pdf")
  plot_counts <- function(metric,title,ylabel) {
    ggplot2::ggplot(enh,ggplot2::aes(enhancer_class,.data[[metric]],fill=enhancer_class)) +
      ggplot2::geom_blank(data=scaffold,ggplot2::aes(x=enhancer_class,y=0),inherit.aes=FALSE) +
      ggplot2::geom_boxplot(outlier.size=0.3,outlier.alpha=0.15,na.rm=TRUE) +
      ggplot2::scale_y_continuous(trans="log1p",breaks=c(0,1,2,5,10,20,50,100,500,1000,5000)) +
      ggplot2::scale_fill_manual(values=class_colors,limits=names(class_colors),drop=FALSE,guide="none") +
      ggplot2::facet_wrap(~timepoint,nrow=1,drop=FALSE) + theme_common +
      ggplot2::labs(title=title,subtitle="One observation per enhancer, including zero-pair enhancers; log1p y spacing",
        x=NULL,y=ylabel)
  }
  depth <- result$parameters$min_pair_reads
  b <- plot_counts("standardized_pairs","B  Co-occupied pairs at active and inactive enhancers",
                   sprintf("Expected distinct pairs in %d fibers",depth))
  labels <- merge(
    comparisons[metric=="standardized_pairs" & method=="wilcoxon_active_greater"],
    comparisons[metric=="standardized_pairs" & method=="adjusted_HC3_active_greater",
                .(time_min, adjusted_q=q_value)], by="time_min",all.x=TRUE)
  if (nrow(labels)) {
    labels <- data.table::copy(labels)
    labels[,timepoint:=factor(paste0(time_min," min"),levels=paste0(times," min"))]
    labels[,label:=sprintf("Wilcoxon q = %s\nAdjusted q = %s",
      format.pval(q_value,digits=2,eps=1e-300),format.pval(adjusted_q,digits=2,eps=1e-300))]
    b <- b + ggplot2::geom_text(data=labels,ggplot2::aes(x=1.5,y=Inf,label=label),
      inherit.aes=FALSE,vjust=1.2,size=2.8)
  }
  save(b,"05_B_cobound_pairs_per_enhancer.pdf",height=5)
  save(plot_counts("n_cobound_pairs","B  Raw distinct co-occupied pairs per enhancer",
                   "Distinct site pairs with >=1 both-bound fiber"),"05_B_raw_cobound_pairs_per_enhancer.pdf",height=5)
  class_counts <- pair_data[,.(n_pairs=.N),by=.(enhancer_class,timepoint,cobinding_class)]
  counts <- ggplot2::ggplot(class_counts,ggplot2::aes(enhancer_class,n_pairs,fill=cobinding_class)) +
    ggplot2::geom_blank(data=scaffold,ggplot2::aes(x=enhancer_class,y=0),inherit.aes=FALSE) +
    ggplot2::geom_col() + ggplot2::facet_wrap(~timepoint,nrow=1,drop=FALSE) +
    ggplot2::scale_fill_manual(values=binding_colors,limits=binding_levels,drop=FALSE) + theme_common +
    ggplot2::labs(title="B  Total tested pairs by association class",x=NULL,y="Site pairs",
      fill=NULL,subtitle="Descriptive totals depend on the number of enhancers and possible pairs")
  save(counts,"05_B_association_class_counts.pdf",height=5)
  binding_scaffold <- scaffold[rep(seq_len(nrow(scaffold)),each=length(binding_levels))]
  binding_scaffold[,cobinding_class:=factor(rep(binding_levels,times=nrow(scaffold)),levels=binding_levels)]
  c <- ggplot2::ggplot(pair_data,ggplot2::aes(cobinding_class,site_distance_bp,fill=cobinding_class)) +
    ggplot2::geom_blank(data=binding_scaffold,ggplot2::aes(x=cobinding_class,y=0),inherit.aes=FALSE) +
    ggplot2::geom_boxplot(outlier.size=0.25,outlier.alpha=0.15,na.rm=TRUE) +
    ggplot2::scale_fill_manual(values=binding_colors,limits=binding_levels,drop=FALSE,guide="none") +
    ggplot2::scale_x_discrete(drop=FALSE,labels=c("depleted","no significant\nassociation","enriched")) +
    ggplot2::scale_y_continuous(trans="log1p") +
    ggplot2::facet_grid(enhancer_class~timepoint,drop=FALSE) + theme_common +
    ggplot2::labs(title="C  Distance between site centers by co-binding class",
      x=NULL,y="Site-center distance (bp; log1p spacing)")
  save(c,"05_C_distance_by_cobinding_class.pdf")
  occupancies <- data.table::melt(pair_data,
    id.vars=c("enhancer_class","timepoint","cobinding_class"),
    measure.vars=c("occupancy_both","occupancy_site1_only","occupancy_site2_only","occupancy_neither"),
    variable.name="state",value.name="occupancy")
  occupancies[,state:=factor(state,levels=c("occupancy_both","occupancy_site1_only","occupancy_site2_only","occupancy_neither"),
    labels=c("Both","Site 1\nonly","Site 2\nonly","Neither"))]
  d <- ggplot2::ggplot(occupancies,ggplot2::aes(state,occupancy,fill=enhancer_class)) +
    ggplot2::geom_blank(data=binding_scaffold,ggplot2::aes(x="Both",y=0),inherit.aes=FALSE) +
    ggplot2::geom_boxplot(outlier.size=0.2,outlier.alpha=0.1,na.rm=TRUE) +
    ggplot2::scale_fill_manual(values=class_colors,limits=names(class_colors),drop=FALSE) +
    ggplot2::scale_y_continuous(limits=c(0,1)) +
    ggplot2::facet_grid(cobinding_class~timepoint,drop=FALSE) + theme_common +
    ggplot2::labs(title="D  Binding-state occupancies in the three co-binding classes",
      x="Pair binding state",y="Fraction of spanning fibers",fill="Enhancer class")
  save(d,"05_D_binding_state_occupancies.pdf",height=10)
  table_dir <- file.path(dirname(plot_dir),"tables")
  dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
  data.table::fwrite(class_counts,file.path(table_dir,"05_cobinding_class_counts.tsv"),sep="\t")
  invisible(list(A=a,B=b,C=c,D=d))
}
