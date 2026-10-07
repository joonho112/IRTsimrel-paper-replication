# Shared paths, display conventions, and portable graphics devices.
# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(scales)
  library(patchwork)
})

source("code/00_setup.R")
MS_ROOT <- repo_root()
CB_ROOT <- repo_path("data", "precomputed")
FIG_DIR <- repo_path("output", "figures")
DATA_DIR <- repo_path("data", "precomputed", "display")
DERIVED_DIR <- repo_path("output", "derived")
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
cb <- function(...) file.path(CB_ROOT, ...)
ms_data <- function(name) file.path(DATA_DIR, name)
read_locked <- function(path) {
  stopifnot(file.exists(path))
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

PAL <- list(
  ink = "#1A1A1A", ink2 = "#4D4D4D", ink3 = "#767676",
  grid = "grey90", band = "#EFF3F8", band2 = "#F5EFE6",
  blue = "#0072B2", sky = "#56B4E9", green = "#009E73",
  orange = "#E69F00", vermillion = "#D55E00", purple = "#CC79A7",
  yellow = "#F0E442", grey = "#999999"
)

PAL_FUN <- c(rho_tilde = PAL$blue, rho_psd = PAL$green, w_bar = PAL$vermillion)
LAB_FUN <- c(rho_tilde = "widetilde(rho)", rho_psd = "rho[PW]", w_bar = "bar(w)")

PAL_SHAPE <- c(
  Normal = PAL$blue, Skewed = PAL$purple,
  Bimodal = PAL$vermillion, `Heavy-tailed` = PAL$green
)
SHAPE_LEVELS <- c("Normal", "Skewed", "Bimodal", "Heavy-tailed")
shape_lab <- function(x) {
  y <- c(
    normal = "Normal", skewed = "Skewed", skew_pos = "Skewed",
    skew = "Skewed", bimodal = "Bimodal",
    heavy = "Heavy-tailed", heavy_tail = "Heavy-tailed",
    heavy_tailed = "Heavy-tailed", `heavy-tailed` = "Heavy-tailed",
    t5 = "Heavy-tailed"
  )[tolower(x)]
  factor(unname(y), levels = SHAPE_LEVELS)
}

PAL_ALGO <- c(
  EQC = PAL$blue, `SAC` = PAL$vermillion,
  `SAC (superpopulation)` = PAL$purple
)

BASE_SIZE <- 8.3
BASE_FAMILY <- "sans"

theme_hs <- function(base_size = BASE_SIZE) {
  theme_minimal(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = PAL$grid, linewidth = 0.28),
      axis.title = element_text(size = base_size, color = PAL$ink),
      axis.text = element_text(size = base_size - 0.6, color = PAL$ink2),
      axis.ticks = element_blank(),
      strip.text = element_text(
        size = base_size + 0.2, face = "bold",
        color = PAL$ink, hjust = 0,
        margin = margin(b = 3)
      ),
      strip.placement = "outside",
      legend.title = element_text(size = base_size, color = PAL$ink),
      legend.text = element_text(size = base_size - 0.6, color = PAL$ink2),
      legend.key.height = unit(3.2, "mm"),
      legend.key.width = unit(3.2, "mm"),
      plot.margin = margin(2, 5, 2, 2),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(
        size = base_size + 1.2, face = "bold",
        family = BASE_FAMILY
      )
    )
}

theme_schematic <- function(base_size = BASE_SIZE) {
  theme_void(base_size = base_size, base_family = BASE_FAMILY) +
    theme(
      strip.text = element_text(
        size = base_size + 0.2, face = "bold",
        color = PAL$ink, hjust = 0,
        margin = margin(b = 3)
      ),
      plot.margin = margin(2, 2, 2, 2),
      plot.title = element_blank(),
      plot.tag = element_text(
        size = base_size + 1.2, face = "bold",
        family = BASE_FAMILY
      )
    )
}

annot <- function(x, y, label, size = 2.55, color = PAL$ink2, ...) {
  annotate("text",
    x = x, y = y, label = label, size = size,
    family = BASE_FAMILY, color = color, ...
  )
}

valid_pdf <- function(path) {
  file.exists(path) && file.info(path)$size > 100 &&
    identical(readBin(path, "raw", n = 5), charToRaw("%PDF-"))
}
.pdf_device <- NULL
choose_pdf_device <- function() {
  if (!is.null(.pdf_device)) {
    return(.pdf_device)
  }
  candidates <- list(
    Cairo = function(filename, width, height, ...) {
      grDevices::cairo_pdf(filename,
        width = width, height = height,
        family = BASE_FAMILY, ...
      )
    }
  )
  if (isTRUE(capabilities("aqua"))) {
    candidates$Quartz <- function(filename, width, height, ...) {
      grDevices::quartz(
        type = "pdf", file = filename, width = width,
        height = height, family = BASE_FAMILY, ...
      )
    }
  }
  candidates$Base <- function(filename, width, height, ...) {
    grDevices::pdf(
      file = filename, width = width, height = height,
      family = "Helvetica", useDingbats = FALSE, ...
    )
  }
  for (name in names(candidates)) {
    probe <- tempfile(fileext = ".pdf")
    before <- grDevices::dev.cur()
    ok <- suppressWarnings(tryCatch(
      {
        candidates[[name]](probe, 1, 1)
        if (grDevices::dev.cur() == before) stop("device did not open")
        graphics::par(mar = rep(0, 4))
        graphics::plot.new()
        graphics::text(.5, .5, "IRT alpha x -", family = BASE_FAMILY)
        grDevices::dev.off()
        valid_pdf(probe)
      },
      error = function(e) FALSE
    ))
    if (grDevices::dev.cur() != before) try(grDevices::dev.off(), silent = TRUE)
    unlink(probe)
    if (isTRUE(ok)) {
      .pdf_device <<- candidates[[name]]
      message("PDF device verified by temporary-file probe: ", name)
      return(.pdf_device)
    }
    message("PDF device unavailable after actual probe: ", name)
  }
  stop("No usable PDF device; retained figures have not been overwritten")
}

dev_portable_pdf <- function(filename, width, height, ...) {
  choose_pdf_device()(filename, width, height, ...)
}

save_fig <- function(plot, name, width_mm, height_mm) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  path <- file.path(FIG_DIR, paste0(name, ".pdf"))
  temporary_pdf <- tempfile(
    pattern = paste0(name, "-"), tmpdir = FIG_DIR,
    fileext = ".pdf"
  )
  on.exit(unlink(temporary_pdf), add = TRUE)
  ggsave(temporary_pdf, plot, width = w, height = h, device = dev_portable_pdf)
  if (!valid_pdf(temporary_pdf)) stop("PDF was not actually created: ", name)
  if (!file.rename(temporary_pdf, path)) stop("Cannot publish PDF: ", path)

  png_path <- file.path(FIG_DIR, paste0(name, ".png"))
  ggsave(png_path, plot, width = w, height = h, dpi = 300, bg = "white")
  message(sprintf("wrote %s (%g x %g mm)", name, width_mm, height_mm))
}

lab_nolead <- function(x) {
  x <- ifelse(abs(x) < 1e-9, 0, x)
  vapply(x, function(xi) {
    if (is.na(xi)) {
      return("")
    }
    s <- format(xi, drop0trailing = TRUE, scientific = FALSE, trim = TRUE)
    sub("^(-?)0\\.", "\\1.", s)
  }, character(1))
}

MM_FULL <- 165
MM_MID <- 130
MM_SINGLE <- 85

logistic <- function(x) 1 / (1 + exp(-x))
item_p <- function(theta, beta, lambda, g = 0) {
  g + (1 - g) * logistic(lambda * (theta - beta))
}
item_info <- function(theta, beta, lambda, g = 0) {
  p <- item_p(theta, beta, lambda, g)
  lambda^2 * (p - g)^2 * (1 - p) / ((1 - g)^2 * p)
}
test_info <- function(theta, beta, lambda, g = 0) {
  rowSums(mapply(
    function(b, l, gg) item_info(theta, b, l, gg),
    beta, lambda, if (length(g) == 1) rep(g, length(beta)) else g
  ))
}

fun_three <- function(J, w, v = 1) {
  c(
    rho_tilde = {
      A <- sum(w * J)
      v * A / (v * A + 1)
    },
    rho_psd = sum(w * (v * J) / (v * J + 1)),
    w_bar = {
      m <- sum(w / J)
      v / (v + m)
    }
  )
}
