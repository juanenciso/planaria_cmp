# ===========================
# plots_markers.R
# Marker-family plots (planaria)
# ===========================

suppressPackageStartupMessages({
  need <- c("readr","dplyr","tidyr","ggplot2","scales","forcats")
  miss <- need[!need %in% rownames(installed.packages())]
  if (length(miss) > 0) install.packages(miss, repos = "https://cloud.r-project.org")
  lapply(need, library, character.only = TRUE)
})

# Optional: composite layout
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork", repos = "https://cloud.r-project.org")
}
library(patchwork)

# ---------- I/O ----------
dir.create("figures", showWarnings = FALSE)
long_file <- "markers/family_TPM_long.tsv"   # columns: species, family, TPM_sum, (optional n)
comp_file <- "markers/compare_family_TPM.tsv" # columns: family, smed_TPM, smed_n, dug_TPM, dug_n, (optional log2fc)

stopifnot(file.exists(long_file), file.exists(comp_file))

# ---------- Helpers ----------
lab_fun <- if (utils::packageVersion("scales") >= "1.2.0") {
  scales::label_number(scale_cut = scales::cut_short_scale())
} else {
  scales::label_number_si()
}

fam_levels <- c("Vasa/DDX4","PL10/DDX3","Piwi","Argonaute/Ago","Nanos")

pal <- c(
  "Vasa/DDX4"     = "#4E79A7",
  "PL10/DDX3"     = "#F28E2B",
  "Piwi"          = "#59A14F",
  "Argonaute/Ago" = "#E15759",
  "Nanos"         = "#B07AA1"
)

# ---------- Read LONG table ----------
long <- readr::read_tsv(long_file, show_col_types = FALSE)

# normalize minimal names
nmL  <- tolower(names(long))
mapL <- c(species="species", family="family", tpm_sum="TPM_sum", sum_tpm="TPM_sum")
for (k in names(mapL)) {
  hit <- which(nmL == k)
  if (length(hit) == 1) names(long)[hit] <- mapL[[k]]
}
stopifnot(all(c("species","family","TPM_sum") %in% names(long)))

long <- long %>%
  dplyr::mutate(
    species = as.character(.data$species),
    family  = factor(.data$family, levels = fam_levels),
    TPM_sum = suppressWarnings(as.numeric(.data$TPM_sum))
  )

# ---------- Read COMP table ----------
comp_try <- readr::read_tsv(comp_file, show_col_types = FALSE)
bad_cols <- !("family" %in% names(comp_try)) ||
  any(!grepl("[A-Za-z]", names(comp_try)[-1]))

if (bad_cols) {
  comp <- readr::read_tsv(comp_file, col_names = FALSE, show_col_types = FALSE)
  names(comp) <- c("family","smed_TPM","smed_n","dug_TPM","dug_n","log2fc")
} else {
  comp <- comp_try
  nm  <- tolower(names(comp))
  map <- c(
    family="family",
    smed_tpm="smed_TPM", smed_sum_tpm="smed_TPM", smed_sum="smed_TPM",
    smed_n="smed_n",
    dug_tpm="dug_TPM",  dug_sum_tpm="dug_TPM",  dug_sum="dug_TPM",
    dug_n="dug_n",
    log2fc_smed_over_dug="log2fc", log2fc="log2fc"
  )
  for (k in names(map)) {
    hit <- which(nm == k)
    if (length(hit) == 1) names(comp)[hit] <- map[[k]]
  }
  if (!"log2fc" %in% names(comp)) {
    cand <- names(comp)[grepl("log2.*smed.*dug|log2fc", tolower(names(comp)))]
    if (length(cand) == 1) names(comp)[names(comp) == cand] <- "log2fc"
  }
}

stopifnot(all(c("family","smed_TPM","dug_TPM") %in% names(comp)))

comp <- comp %>%
  dplyr::mutate(
    family   = factor(.data$family, levels = fam_levels),
    smed_TPM = suppressWarnings(as.numeric(.data$smed_TPM)),
    dug_TPM  = suppressWarnings(as.numeric(.data$dug_TPM)),
    smed_n   = if ("smed_n" %in% names(.)) suppressWarnings(as.numeric(.data$smed_n)) else NA_real_,
    dug_n    = if ("dug_n"  %in% names(.)) suppressWarnings(as.numeric(.data$dug_n))  else NA_real_,
    log2fc   = suppressWarnings(as.numeric(.data$log2fc))
  ) %>%
  dplyr::mutate(
    log2fc = dplyr::coalesce(.data$log2fc, log2((.data$smed_TPM + 1e-6)/(.data$dug_TPM + 1e-6)))
  )

comp2 <- comp %>%
  dplyr::mutate(family = forcats::fct_rev(.data$family)) %>%
  dplyr::arrange(.data$family)

# ---------- Figures ----------
# 1) Total TPM by family
p_abs <- long %>%
  ggplot2::ggplot(ggplot2::aes(x = family, y = TPM_sum, fill = family)) +
  ggplot2::geom_col(width = 0.75, color = "grey20") +
  ggplot2::facet_wrap(~ species, nrow = 1) +
  ggplot2::scale_fill_manual(values = pal, guide = "none") +
  ggplot2::scale_y_continuous(labels = lab_fun) +
  ggplot2::labs(title = "Total TPM by marker family", x = NULL, y = "TPM (sum)") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.x = element_blank(),
                 axis.text.x = element_text(angle = 20, hjust = 1),
                 strip.text = element_text(face = "bold"))

ggplot2::ggsave("figures/family_TPM_sum.png", p_abs, width = 9, height = 4.2, dpi = 300)
ggplot2::ggsave("figures/family_TPM_sum.pdf", p_abs, width = 9, height = 4.2)

# 2) Percent composition within species
long_pct <- long %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(TPM_pct = 100 * TPM_sum / sum(TPM_sum, na.rm = TRUE)) %>%
  dplyr::ungroup()

p_pct <- long_pct %>%
  ggplot2::ggplot(ggplot2::aes(x = family, y = TPM_pct, fill = family)) +
  ggplot2::geom_col(width = 0.75, color = "grey20") +
  ggplot2::facet_wrap(~ species, nrow = 1) +
  ggplot2::scale_fill_manual(values = pal, guide = "none") +
  ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggplot2::labs(title = "Within-species composition by family",
                x = NULL, y = "% of total TPM") +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.x = element_blank(),
                 axis.text.x = element_text(angle = 20, hjust = 1),
                 strip.text = element_text(face = "bold"))

ggplot2::ggsave("figures/family_TPM_pct.png", p_pct, width = 9, height = 4.2, dpi = 300)
ggplot2::ggsave("figures/family_TPM_pct.pdf", p_pct, width = 9, height = 4.2)

# 3) Lollipop log2FC (Smed / Dug)
xr <- range(comp2$log2fc, na.rm = TRUE)
xmax <- max(abs(xr), na.rm = TRUE)
xlim_sym <- c(-xmax, xmax) + c(-0.15, 0.15)

p_fc <- ggplot2::ggplot(comp2, ggplot2::aes(y = family, x = log2fc, color = family)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = log2fc, yend = family), linewidth = 1) +
  ggplot2::geom_point(size = 3, stroke = 0) +
  ggplot2::scale_color_manual(values = pal, guide = "none") +
  ggplot2::scale_x_continuous(limits = xlim_sym, breaks = scales::pretty_breaks()) +
  ggplot2::labs(title = "log2FC (Schmidtea / Dugesia) by family",
                x = "log2 fold-change (TPMsum Smed / TPMsum Dug)",
                y = NULL) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(panel.grid.major.y = element_blank())

ggplot2::ggsave("figures/family_log2FC_lollipop.png", p_fc, width = 7.2, height = 4.2, dpi = 300)
ggplot2::ggsave("figures/family_log2FC_lollipop.pdf", p_fc, width = 7.2, height = 4.2)

# Composite page (A/B/C)
combo <- (p_abs / p_pct / p_fc) +
  patchwork::plot_annotation(
    title = "Germline/gene-silencing marker families across planarian species",
    tag_levels = "A",
    theme = ggplot2::theme(plot.title = element_text(face = "bold", size = 16))
  )

ggplot2::ggsave("figures/markers_family_composite_en.png", combo, width = 10, height = 12, dpi = 300)
ggplot2::ggsave("figures/markers_family_composite_en.pdf", combo, width = 10, height = 12)

message("✅ Figures saved in 'figures/'")
