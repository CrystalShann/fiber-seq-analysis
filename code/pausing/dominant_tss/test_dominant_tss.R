suppressPackageStartupMessages({library(data.table); library(GenomicRanges)})
source("code/pausing/dominant_tss/select_dominant_tss.R")
tx <- data.table(gene_id = c("A", "A", "B", "C", "D", "E"),
                 tx_id = paste0("tx", 1:6), gene_name = "test", gene_type = "test",
                 chrom = "chr1", strand = c("+", "+", "-", "+", "+", "+"),
                 tss = c(100L, 200L, 500L, 700L, 720L, 900L),
                 tes = c(400L, 450L, 300L, 800L, 820L, 1000L),
                 canonical = c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE))
cage <- data.table(chrom = "chr1", start0 = c(99L, 201L, 497L, 709L, 899L),
                   end = c(100L, 202L, 498L, 710L, 900L), ctss_id = letters[1:5],
                   score = c(5, 20, 8, 50, 100), strand = c("+", "+", "-", "+", "-"))
out <- select_dominant_tss(tx, cage, 20)
stopifnot(identical(out$gene_id, c("A", "B")),
          identical(out$tss, c(202L, 498L)), identical(out$tes, c(450L, 300L)),
          out[gene_id == "A", tx_id] == "tx2")
# A CTSS tied between two genes is omitted; opposite-strand signal cannot rescue E.
stopifnot(!any(c("C", "D", "E") %in% out$gene_id))
stopifnot(nrow(select_dominant_tss(tx, cage[0], 20)) == 0L)
signal <- GRanges("chr1", IRanges(c(120L, 125L, 470L, 480L), c(121L, 125L, 470L, 480L)),
                  score = c(3, 9, 5, 5))
stopifnot(identical(find_pause_peak(c(100L, 1000L), "+", signal), c(25L, NA_integer_)),
          identical(find_pause_peak(500L, "-", signal), 20L),
          identical(find_pause_peak(100L, "+", GRanges()), NA_integer_))
cat("Dominant TSS selection and strand-aware pause-peak tests passed.\n")
