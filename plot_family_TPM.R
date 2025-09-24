# plot_family_TPM.R
# Requiere: markers/family_TPM_long.tsv y markers/compare_family_TPM.tsv
# Salidas: figures/*.png y *.pdf

# ---------- paquetes ----------
need <- c("readr","dplyr","tidyr","ggplot2","scales","forcats")
inst <- need[!need %in% rownames(installed.packages())]
if (length(inst) > 0) install.packages(inst, repos = "https://cloud.r-project.org")

library(readr); library(dplyr); library(tidyr)
library(ggplot2); library(scales); library(forcats)

# ---------- I/O ----------
dir.create("figures", showWarnings = FALSE)
long_file  <- "markers/family_TPM_long.tsv"
comp_file  <- "markers/compare_family_TPM.tsv"

stopifnot(file.exists(long_file), file.exists(comp_file))

long <- read_tsv(long_file, show_col_types = FALSE)
comp <- read_tsv(comp_file, show_col_types = FALSE)

# Asegurar orden bonito de familias (ajusta si quieres otro orden)
fam_levels <- c("Vasa/DDX4","PL10/DDX3","Piwi","Argonaute/Ago","Nanos")
long$family <- factor(long$family, levels = fam_levels)
comp$family <- factor(comp$family, levels = fam_levels)

# Colores discretos consistentes (puedes cambiarlos)
pal <- c("Vasa/DDX4"="#4E79A7","PL10/DDX3"="#F28E2B","Piwi"="#59A14F",
         "Argonaute/Ago"="#E15759","Nanos"="#B07AA1")

# =============================
# 1) Barras comparativas (TPM suma)
# =============================
p_abs <- long %>%
  ggplot(aes(x = family, y = TPM_sum, fill = family)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.75, color = "grey20") +
  facet_wrap(~ species, nrow = 1) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = label_number_si()) +
  labs(title = "Suma de TPM por familia de marcadores",
       x = NULL, y = "TPM (suma)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave("figures/family_TPM_sum.png", p_abs, width = 9, height = 4.2, dpi = 300)
ggsave("figures/family_TPM_sum.pdf", p_abs, width = 9, height = 4.2)

# =============================
# 2) Porcentaje dentro de especie
# =============================
long_pct <- long %>%
  group_by(species) %>%
  mutate(TPM_pct = 100 * TPM_sum / sum(TPM_sum, na.rm = TRUE)) %>%
  ungroup()

p_pct <- long_pct %>%
  ggplot(aes(x = family, y = TPM_pct, fill = family)) +
  geom_col(width = 0.75, color = "grey20") +
  facet_wrap(~ species, nrow = 1) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "Composición porcentual por familia (dentro de cada especie)",
       x = NULL, y = "% del TPM total") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 20, hjust = 1),
        strip.text = element_text(face = "bold"))

ggsave("figures/family_TPM_pct.png", p_pct, width = 9, height = 4.2, dpi = 300)
ggsave("figures/family_TPM_pct.pdf", p_pct, width = 9, height = 4.2)

# =============================
# 3) Lollipop de log2FC (Smed/Dug)
# =============================
# comp: family, smed_TPM, smed_n, dug_TPM, dug_n, log2FC_smed_over_dug
comp2 <- comp %>%
  mutate(family = fct_rev(family))  # para que quede Piwi abajo si quieres

p_fc <- ggplot(comp2, aes(y = family, x = log2FC_smed_over_dug, color = family)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = 0, xend = log2FC_smed_over_dug, yend = family), linewidth = 1) +
  geom_point(size = 3, stroke = 0) +
  scale_color_manual(values = pal, guide = "none") +
  scale_x_continuous(breaks = pretty_breaks()) +
  labs(title = "log2FC(Schmidtea / Dugesia) por familia",
       x = "log2 fold-change (TPMsum Smed / TPMsum Dug)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.major.y = element_blank())

ggsave("figures/family_log2FC_lollipop.png", p_fc, width = 7.2, height = 4.2, dpi = 300)
ggsave("figures/family_log2FC_lollipop.pdf", p_fc, width = 7.2, height = 4.2)

cat("Listo.\nFiguras en ./figures :\n",
    " - family_TPM_sum.(png|pdf)\n",
    " - family_TPM_pct.(png|pdf)\n",
    " - family_log2FC_lollipop.(png|pdf)\n", sep = "")
