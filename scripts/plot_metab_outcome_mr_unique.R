#### DESCRIPTION #######################################
# Scatter plot: standard IVW vs unique-instrument IVW
# MR estimates for metabolites on outcomes.
# Faceted by outcome. Tests sensitivity of results to
# shared-instrument removal.
########################################################

#### SET INPUT #########################################
results_file  <- snakemake@input[["results_file"]]
unique_file   <- snakemake@input[["unique_file"]]
name_map_file <- snakemake@input[["name_map_file"]]
output_file   <- snakemake@output[["scatter"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
results_file  <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/metabolites_on_outcomes.tsv")
unique_file   <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/tables/mr_results/metabolites_on_outcomes_unique.tsv")
name_map_file <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "scripts/gwas_metab_name_map.xlsx")
output_file   <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output/figures/outcomes/metab_outcome_mr_unique_correlation.png")
# ─────────────────────────────────────────────────────────────────────────────
}

library(data.table)
library(ggplot2)
library(readxl)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
std  <- fread(results_file)
uniq <- fread(unique_file)

# ── 2. PREPARE STANDARD MR ───────────────────────────────────────────────────
# Wald ratio → IVW so single-SNP metabolites match
std[method == "Wald ratio", method := "Inverse variance weighted"]
# Non-Steiger IVW only
std <- std[
  steiger_filtering == FALSE & method == "Inverse variance weighted",
  .(exp_id, out_id, b_std = b, se_std = b_se, p_std = p, nsnp_std = nsnp)
]

# ── 3. PREPARE UNIQUE MR ─────────────────────────────────────────────────────
# Wald ratio → IVW
uniq[method == "Wald ratio", method := "Inverse variance weighted"]
uniq <- uniq[
  method == "Inverse variance weighted",
  .(exp_id, out_id, b_uniq = b, se_uniq = b_se, p_uniq = p,
    nsnp_uniq = nsnp, nsnp_removed = nsnp_removed_shared)
]

# ── 4. JOIN ───────────────────────────────────────────────────────────────────
dt <- merge(std, uniq, by = c("exp_id", "out_id"))

# Attach readable metabolite labels
map <- read_xlsx(name_map_file, sheet = 1) |> as.data.table()
dt[map, label := i.label, on = c("exp_id" = "gwas_id")]
dt[, label := fcoalesce(label, exp_id)]

# Significance tier based on FDR in the standard analysis (per outcome)
dt[, p_fdr_std := p.adjust(p_std, method = "fdr"), by = out_id]
dt[, sig := fcase(
  p_fdr_std < 0.05, "FDR p<0.05",
  p_std     < 0.05, "p<0.05",
  default           = "p≥0.05"
)]
dt[, sig := factor(sig, levels = c("FDR p<0.05", "p<0.05", "p≥0.05"))]

# Outcome facet labels
dt[, out_label := tools::toTitleCase(gsub("_", " ", out_id))]

# ── 5. PER-OUTCOME CORRELATION STATS ─────────────────────────────────────────
cor_dt <- dt[!is.na(b_std) & !is.na(b_uniq), {
  ct <- cor.test(b_std, b_uniq, method = "pearson")
  .(r = round(ct$estimate, 3), n = .N)
}, by = out_label]
cor_dt[, lab := sprintf("r = %.3f\nn = %d", r, n)]

# ── 6. PLOT ───────────────────────────────────────────────────────────────────
sig_cols <- c("FDR p<0.05" = "#B2182B", "p<0.05" = "#F4A582", "p≥0.05" = "grey60")

p <- ggplot(dt, aes(x = b_std, y = b_uniq)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = b_uniq - 1.96 * se_uniq, ymax = b_uniq + 1.96 * se_uniq, colour = sig),
    linewidth = 0.25, alpha = 0.4, width = 0
  ) +
  geom_errorbarh(
    aes(xmin = b_std - 1.96 * se_std, xmax = b_std + 1.96 * se_std, colour = sig),
    linewidth = 0.25, alpha = 0.4, height = 0
  ) +
  geom_point(aes(colour = sig), size = 1.5, alpha = 0.8) +
  geom_text(
    data = cor_dt, aes(label = lab),
    x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2,
    size = 2.8, colour = "grey20", inherit.aes = FALSE
  ) +
  scale_colour_manual(values = sig_cols, name = "Standard MR (FDR within outcome)") +
  facet_wrap(~out_label, scales = "free") +
  labs(
    x = "Standard IVW β (all instruments)",
    y = "Unique-instrument IVW β"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey20", colour = NA),
    strip.text        = element_text(colour = "white", face = "bold", size = 9),
    legend.position   = "bottom"
  )
p

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
n_out <- uniqueN(dt$out_label)
ggsave(output_file, p,
       width  = max(5, 4.5 * n_out),
       height = 5,
       dpi    = 300,
       bg     = "white",
       units  = "in")
cat("Saved:", output_file, "\n")
