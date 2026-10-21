# Plotting functions for the LCL haplotype / allele skew report (LCL_hsca.Rmd)
#
# Taken from the code of Vollger et al., Cell 2026 (Somatic epimutations cap genetic determinism
# in the human diploid chromatin epigenome):
#   fire-figures  https://github.com/mrvollger/fire-figures @ 01e9345dd15ef9e203b8f7e816b441e3627995f8
#     Rscripts/utils.R          FONT_SIZE, my_grid, my_hgrid, my_vgrid, reverselog_trans,
#                               scientific_10, my_read_bed (verbatim)
#     Rscripts/hap-diffs.Rmd    hap1-vs-hap2 hex plot and volcano (wrapped in functions below)
#     Rscripts/Fire-peak-df.R   imprinted = overlap with the Akbari B-lymphocyte DMRs
#   FIRE  workflow/scripts/hap-diffs.R  (same two plots; the test the per-peak table comes from)

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(scales)
  library(cowplot)
  library(ggpubr)
  library(ggrepel)
  library(glue)
  library(valr)
})

# Akbari et al. 143 DMRs in 12 B-lymphocyte cell lines, from fire-figures data/ (hg38)
AKBARI_DMR = "/project/spott/cshan/annotations/lcl_dmr_coordinates_Akbari.bed.gz"

#################################################
# fire-figures Rscripts/utils.R (verbatim)
#################################################

FONT_SIZE=6
my_grid = function(...){
    theme_minimal_grid(font_size=FONT_SIZE, ...)
}

my_hgrid = function(...){
    theme_minimal_hgrid(font_size=FONT_SIZE, ...)
}

my_vgrid = function(...){
    theme_minimal_vgrid(font_size=FONT_SIZE, ...)
}

reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv,
              log_breaks(base = base),
              domain = c(1e-100, Inf))
}

scientific_10 <- function(x) {
    is_one = as.numeric(x) == 1
    text = gsub("e", " %*% 10^", scales::scientific_format()(x))
    #print(text)
    text = str_remove(text, "^1 %\\*% ") # remove leading one
    #print(text)
    text[is_one] = "10^0"
    rtn = parse(text=text)
    rtn
}

my_read_bed = function(...){
    df=fread(...)
    names = colnames(df)
    names[names == "start"] = "start.other"
    names[names == "end"] = "end.other"
    names[1:3] = c("chrom", "start", "end")
    colnames(df) = names
    # sort by chrom then start then end
    df = df[order(df$chrom, df$start, df$end)]
}

#################################################
# imprinted flag (fire-figures Rscripts/Fire-peak-df.R)
#################################################

# df needs chrom/start/end; TRUE where the interval overlaps an Akbari DMR by >=1 bp
flag_imprinted = function(df, dmr_file = AKBARI_DMR){
    imprinted = my_read_bed(dmr_file)
    hits = bed_intersect(as_tibble(df) %>% mutate(.row = row_number()), as_tibble(imprinted))
    seq_len(nrow(df)) %in% hits$.row.x
}

#################################################
# fire-figures Rscripts/hap-diffs.Rmd plots, as functions
#################################################

# Edits from hap-diffs.Rmd, each also marked below:
#   facet by `facet` (sample_name) instead of autosome (all tables here are autosomal)
#   axis labels passed in (hap1/hap2 or REF/ALT); HiPhase haplotypes are not paternal/maternal
#   the BH-equivalent p-value line (cor_p_threshold) is computed per facet, because BH is per sample

# hap1-vs-hap2 accessibility, one hexbin panel per sample
plot_hap1_vs_hap2 = function(tdf, x, y, x_label, y_label, facet = "sample_name"){
    p = tdf %>%
        ggplot(aes(x=.data[[x]], y=.data[[y]])) +
        stat_cor(size=2) +
        geom_hex(bins=75) +
        geom_abline(aes(intercept=0, slope=1), linetype="dashed")+
        scale_fill_distiller("", palette = "Spectral", trans="log10") +
        scale_x_continuous(x_label, labels=percent) +  # edit: was "Paternal accessibility"
        scale_y_continuous(y_label, labels=percent) +  # edit: was "Maternal accessibility"
        #annotation_logticks(sides="lb") +
        my_grid()
    if (!is.null(facet)) p = p + facet_wrap(vars(.data[[facet]]))  # edit: was ~autosome
    p
}

# volcano of the accessibility difference against the Fisher p-value; dashed lines at p = 0.05
# (blue) and at the largest p-value with an adjusted p <= 0.05 (red); red crosses = imprinted
plot_hap_volcano = function(tdf, x_label, p_adjust = "p_adjust", facet = "sample_name", p_threshold = 0.05){
    tdf = as.data.table(tdf)[!is.na(p_value)]
    tdf$p_adjust = tdf[[p_adjust]]
    tdf$facet = if (is.null(facet)) "all" else tdf[[facet]]
    # edit: one cor_p_threshold per facet instead of one per table
    thresholds = tdf[, .(cor_p_threshold = if (any(p_adjust <= p_threshold)) max(p_value[p_adjust <= p_threshold]) else NA_real_), by = facet]
    y_lim = ceiling(max(max(-log10(tdf$p_value)), max(-log10(tdf$p_adjust))))
    y_by = 1
    if(y_lim > 10){
        y_by = 2
    }
    # counts: nominal-only at the blue line, adjusted-significant at the red line
    labels = tdf[p_value < 0.05, .(count = .N), by = .(facet, adjusted = p_adjust <= p_threshold)] %>%
        merge(thresholds, by = "facet") %>%
        mutate(y = ifelse(adjusted, cor_p_threshold, p_threshold)) %>%
        filter(!is.na(y))
    n=comma(nrow(tdf))
    p = tdf %>%
        ggplot(aes(x=diff, y=p_value)) +
        geom_hex(bins=100) + scale_fill_distiller("", palette = "Spectral", trans="log10") +
        geom_hline(aes(yintercept=(p_threshold)), linetype="dashed", color="darkblue")+
        geom_hline(data = thresholds %>% filter(!is.na(cor_p_threshold)),
                   aes(yintercept=(cor_p_threshold)), linetype="dashed", color="darkred")+
        scale_x_continuous(x_label, labels=percent) +  # edit: was "Difference between paternal and maternal accessibility"
        scale_y_continuous(
            glue("p-value   (n = {n})"),
            #labels=comma,
            breaks=10**(-seq(0, y_lim, y_by)),
            minor_breaks=10**(-seq(0, y_lim, 0.1)),
            trans=reverselog_trans(10),
            labels=scientific_10,
        ) +
        geom_text_repel(
            data = labels,
            aes(
                label = paste0("n = ", comma(count)),
                x=0.750,
                y=y,
            ),
            direction="y",
            force=2,
            nudge_y=0.01,
            size=2.5,
        )+
        my_grid()
    if ("imprinted" %in% names(tdf)) {
        p = p + geom_point(
            data = tdf %>% filter(imprinted),
            color="darkred",
            shape=3,
            size=.5,
            alpha=0.65,
        )
    }
    if (!is.null(facet)) p = p + facet_wrap(~facet)  # edit: was ~autosome
    p
}

# PDF only (fire-figures my_ggsave also writes tmp.<ext> into the working directory)
save_pdf = function(p, file, ...){
    ggsave(file, p, bg='transparent', ...)
    message(file)
}
