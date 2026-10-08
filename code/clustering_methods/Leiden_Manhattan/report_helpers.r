# PDF previews shared by the Leiden and autocorrelation reports. Keep the LCL
# font mappings so embedding an existing figure renders the same PNG.

# ---------------------------------------------------------------------------
# Write a temporary xpdf config that maps the standard PDF fonts (Helvetica,
# Times, Courier, Symbol, ZapfDingbats) to the Liberation Sans / DejaVu font
# files on the cluster, so pdftopng renders previews with consistent fonts.
#
# Inputs:
#   none
# Output:
#   path of the temporary xpdfrc file; the caller deletes it
# ---------------------------------------------------------------------------
report_xpdf_config <- function() {
  xpdfrc <- tempfile("xpdfrc-")
  sans <- "/usr/share/fonts/liberation-sans/LiberationSans-"
  mono <- "/usr/share/fonts/dejavu/DejaVuSansMono"
  writeLines(c(
    paste0("fontFile Helvetica ", sans, "Regular.ttf"),
    paste0("fontFile Helvetica-Bold ", sans, "Bold.ttf"),
    paste0("fontFile Helvetica-Oblique ", sans, "Italic.ttf"),
    paste0("fontFile Helvetica-BoldOblique ", sans, "BoldItalic.ttf"),
    paste0("fontFile Times-Roman ", sans, "Regular.ttf"),
    paste0("fontFile Times-Bold ", sans, "Bold.ttf"),
    paste0("fontFile Times-Italic ", sans, "Italic.ttf"),
    paste0("fontFile Times-BoldItalic ", sans, "BoldItalic.ttf"),
    paste0("fontFile Courier ", mono, ".ttf"),
    paste0("fontFile Courier-Bold ", mono, "-Bold.ttf"),
    paste0("fontFile Courier-Oblique ", mono, "-Oblique.ttf"),
    paste0("fontFile Courier-BoldOblique ", mono, "-BoldOblique.ttf"),
    "fontFile Symbol /usr/share/fonts/dejavu/DejaVuSans.ttf",
    "fontFile ZapfDingbats /usr/share/fonts/dejavu/DejaVuSans.ttf"), xpdfrc)
  xpdfrc
}

# ---------------------------------------------------------------------------
# Embed the first page of a PDF figure in a knitted HTML report: render it to
# PNG with xpdf's pdftopng (fonts from report_xpdf_config()) and print an HTML
# <figure> with the PNG as a base64 data URI and `label` as its caption. Use
# in a chunk with results = "asis".
#
# Inputs:
#   path     - existing, non-empty PDF
#   label    - caption and alt text
#   dpi      - rendering resolution
#   pdftopng - path of the pdftopng executable
# Output:
#   invisible NULL; the HTML is printed with cat()
# ---------------------------------------------------------------------------
inline_png <- function(path, label, dpi = 120,
                       pdftopng = "/software/xpdf-4.05-el8-x86_64/bin/pdftopng") {
  stopifnot(file.exists(path), file.info(path)$size > 0)
  xpdfrc <- report_xpdf_config()
  prefix <- tempfile("fig-")
  png <- paste0(prefix, "-000001.png")
  on.exit(unlink(c(xpdfrc, png)), add = TRUE)
  status <- system2(pdftopng, c("-cfg", shQuote(xpdfrc), "-q", "-r", dpi, "-f", 1, "-l", 1,
    shQuote(path), shQuote(prefix)))
  stopifnot(status == 0, file.exists(png))
  cat(as.character(htmltools::tags$figure(
    htmltools::tags$img(
      src = base64enc::dataURI(file = png, mime = "image/png"),
      alt = label, loading = "lazy",
      style = "display:block;max-width:100%;height:auto;"),
    htmltools::tags$figcaption(label))), "\n\n")
  invisible(NULL)
}
