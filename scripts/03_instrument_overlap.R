#### DESCRIPTION #######################################
# Assess SNP overlap between metabolite instruments and outcome
# instruments. Used to evaluate whether reverse MR results could
# be driven by shared genetic instruments (exclusion restriction
# violation).
########################################################

#### SET INPUT #########################################
metabolite_instrument_files <- snakemake@input[["metabolite_instrument_files"]]
outcome_instrument_files    <- snakemake@input[["outcome_instrument_files"]]
consistent_ids_file         <- snakemake@input[["consistent_ids_file"]]
cluster_file                <- snakemake@input[["cluster_file"]]
name_map_file               <- snakemake@input[["name_map_file"]]
overlap_table               <- snakemake@output[["overlap_table"]]
overlap_heatmap             <- snakemake@output[["overlap_heatmap"]]
outcome_snp_table           <- snakemake@output[["outcome_snp_table"]]
log_file                    <- snakemake@log[["log"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
repo_dir <- Sys.getenv("HF_METABOLITE_REPO2")
metabolite_instrument_files <- list.files(file.path(repo_dir, "output/tables/instruments/genome_wide"), pattern = "^GCST", full.names = TRUE)
outcome_instrument_files    <- file.path(repo_dir, "output/tables/instruments/genome_wide",
                                         c("heart_failure_instrument.tsv",
                                           "hip_knee_osteoarthritis_instrument.tsv",
                                           "kidney_cancer_instrument.tsv"))
consistent_ids_file <- file.path(repo_dir, "output/tables/replication/consistent_gwas_ids.txt")
cluster_file        <- file.path(repo_dir, "output/tables/instruments/cluster_membership.tsv")
name_map_file       <- file.path(repo_dir, "scripts/gwas_metab_name_map.xlsx")
overlap_table       <- file.path(repo_dir, "output/tables/instruments/metab_outcome_instrument_overlap.tsv")
overlap_heatmap     <- file.path(repo_dir, "output/figures/instruments/metab_outcome_instrument_overlap.png")
outcome_snp_table   <- file.path(repo_dir, "output/tables/instruments/outcome_snps_metab_overlap.tsv")
log_file            <- file.path(repo_dir, "output/logs/instrument_overlap.log")
# ─────────────────────────────────────────────────────────────────────────────
}

library(ggplot2)
library(readxl)
suppressPackageStartupMessages(library(data.table))

dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
on.exit(close(log_con))
tee <- function(...) { cat(...); cat(..., file = log_con, append = TRUE) }

thr <- 0.25  # overlap coefficient threshold — same as 02_instrument_clusters.R


# metadata ====
map            <- read_xlsx(name_map_file, sheet = 1) |> as.data.table()

# disambiguate labels that appear on more than one platform
platform_suffix <- c("nightingale" = "NMR", "metabolon" = "MS")
dup_labels      <- map[, .(n = uniqueN(gwas_id)), by = label][n > 1, label]
map[label %in% dup_labels,
    label := paste0(label, " (", platform_suffix[platform], ")")]

consistent_ids <- fread(consistent_ids_file, header = FALSE, col.names = "gwas_id")
consistent_ids[map, c("label", "pathway_group") := .(i.label, i.pathway_group), on = "gwas_id"]
clust <- fread(cluster_file, sep = "\t")


# load metabolite instruments (consistent metabolites only) ====
metab_instruments <- lapply(metabolite_instrument_files, function(fp) {
  id0 <- sub("(GCST[0-9]+)_instrument\\.tsv$", "\\1", basename(fp))
  if (!(id0 %in% consistent_ids$gwas_id)) return(NULL)
  d <- fread(fp)
  if (nrow(d) == 0 || !("rsid" %in% names(d))) return(NULL)
  d[!is.na(rsid), .(trait = id0, rsid)]
}) |> rbindlist(use.names = TRUE, fill = TRUE)

metab_instruments[map, label := i.label, on = c("trait" = "gwas_id")]
metab_instruments <- metab_instruments[!is.na(label)]

tee(sprintf(
  "Metabolite instruments loaded: %i metabolites, %i unique SNPs\n",
  metab_instruments[, uniqueN(trait)],
  metab_instruments[, uniqueN(rsid)]
))


# load outcome instruments (keep all columns for formatted table) ====
out_instruments <- lapply(outcome_instrument_files, function(fp) {
  out_id_val <- sub("_instrument\\.tsv$", "", basename(fp))
  d <- fread(fp)
  if (nrow(d) == 0 || !("rsid" %in% names(d))) return(NULL)
  d[!is.na(rsid), out_id := out_id_val]
  d[!is.na(rsid)]
}) |> rbindlist(use.names = TRUE, fill = TRUE)

tee(sprintf(
  "Outcome instruments loaded: %i outcomes, %i unique SNPs\n",
  out_instruments[, uniqueN(out_id)],
  out_instruments[, uniqueN(rsid)]
))


# build SNP sets per metabolite and per outcome ====
metab_snp_sets <- metab_instruments[, .(snps = list(unique(rsid))), by = .(trait, label)]
out_snp_sets   <- out_instruments[,   .(snps = list(unique(rsid))), by = out_id]


# pairwise overlap: all metabolites × all outcomes ====
pairs <- CJ(
  metab_idx = seq_len(nrow(metab_snp_sets)),
  out_idx   = seq_len(nrow(out_snp_sets))
)

overlap_dt <- pairs[, {
  m_snps   <- metab_snp_sets$snps[[metab_idx]]
  o_snps   <- out_snp_sets$snps[[out_idx]]
  shared   <- intersect(m_snps, o_snps)
  n_shared <- length(shared)
  n_metab  <- length(m_snps)
  n_out    <- length(o_snps)
  pct_metab <- if (n_metab > 0) n_shared / n_metab else 0
  pct_out   <- if (n_out   > 0) n_shared / n_out   else 0
  .(
    gwas_id          = metab_snp_sets$trait[metab_idx],
    label            = metab_snp_sets$label[metab_idx],
    out_id           = out_snp_sets$out_id[out_idx],
    n_metab_snps     = n_metab,
    n_out_snps       = n_out,
    n_shared         = n_shared,
    pct_metab_shared = round(pct_metab, 4),
    pct_out_shared   = round(pct_out,   4),
    overlap_coef     = round(max(pct_metab, pct_out), 4),
    shared_rsids     = if (n_shared > 0) paste(shared, collapse = ";") else NA_character_
  )
}, by = .(metab_idx, out_idx)][, c("metab_idx", "out_idx") := NULL]

overlap_dt[clust, cluster := i.cluster, on = "label"]
overlap_dt[is.na(cluster), cluster := 0L]


# write overlap table ====
dir.create(dirname(overlap_table), recursive = TRUE, showWarnings = FALSE)
fwrite(overlap_dt, overlap_table, sep = "\t")


# formatted outcome SNP table with metabolite co-usage annotation ====
out_snp_tab <- copy(out_instruments)

# (out_id, rsid, trait, label) triples for all SNPs shared with any metabolite
rsid_metab_map <- metab_instruments[, .(trait, label, rsid)]
shared_triples <- out_snp_tab[!is.na(rsid), unique(.SD), .SDcols = c("out_id","rsid")][
  rsid_metab_map, on = "rsid", nomatch = NULL, allow.cartesian = TRUE
]
shared_triples[
  overlap_dt,
  c("n_shared", "pct_out_shared", "pct_metab_shared") :=
    .(i.n_shared, i.pct_out_shared, i.pct_metab_shared),
  on = c("out_id", "trait" = "gwas_id")
]
shared_triples <- shared_triples[!is.na(n_shared)]

# per-(out_id, rsid): n_metab_overlap, metab_overlap string, max_metab_overlap
snp_summary <- shared_triples[, {
  ord  <- order(-n_shared, label)
  best <- which.max(n_shared)
  .(
    n_metab_overlap   = uniqueN(trait),
    metab_overlap     = paste(
      sprintf("%s|%s|%.1f%%|%.1f%%",
              label[ord], trait[ord],
              100 * pct_out_shared[ord],
              100 * pct_metab_shared[ord]),
      collapse = "; "),
    max_metab_overlap = sprintf("%i (%.1f%%|%.1f%%)",
                                n_shared[best],
                                100 * pct_out_shared[best],
                                100 * pct_metab_shared[best])
  )
}, by = .(out_id, rsid)]

out_snp_tab[snp_summary,
            c("n_metab_overlap","metab_overlap","max_metab_overlap") :=
              .(i.n_metab_overlap, i.metab_overlap, i.max_metab_overlap),
            on = c("out_id","rsid")]
out_snp_tab[is.na(n_metab_overlap),   n_metab_overlap   := 0L]
out_snp_tab[is.na(metab_overlap),     metab_overlap     := ""]
out_snp_tab[is.na(max_metab_overlap), max_metab_overlap := ""]

col_order <- c("out_id", "rsid",
               intersect(c("chr","bp","ea","oa","beta","se","p","eaf","n"), names(out_snp_tab)),
               "n_metab_overlap", "max_metab_overlap", "metab_overlap")
setcolorder(out_snp_tab, intersect(col_order, names(out_snp_tab)))
setorder(out_snp_tab, out_id, -n_metab_overlap, rsid)

fwrite(out_snp_tab, outcome_snp_table, sep = "\t")


# log summary ====
tee("\n--- Overlap summary per outcome ---\n")
for (oid in sort(out_snp_sets$out_id)) {
  sub     <- overlap_dt[out_id == oid]
  n_any   <- sub[n_shared > 0, .N]
  n_above <- sub[overlap_coef >= thr, .N]
  tee(sprintf(
    "\n  %s  (outcome instrument n=%i SNPs)\n    metabolites tested=%i  any overlap=%i  above threshold (%.2f): %i\n",
    oid,
    out_snp_sets[out_id == oid, lengths(snps)],
    nrow(sub), n_any, thr, n_above
  ))
  if (n_any > 0) {
    flagged <- sub[n_shared > 0][order(-overlap_coef)]
    for (i in seq_len(nrow(flagged))) {
      tee(sprintf(
        "    %-35s  n_shared=%i  pct_metab=%.1f%%  pct_out=%.1f%%  ov_coef=%.3f  SNPs: %s\n",
        flagged$label[i],
        flagged$n_shared[i],
        100 * flagged$pct_metab_shared[i],
        100 * flagged$pct_out_shared[i],
        flagged$overlap_coef[i],
        flagged$shared_rsids[i]
      ))
    }
  }
}


# bar chart ====
# x = metabolites ordered by cluster; bars = n_shared; reference line = n_out_snps;
# faceted vertically by outcome
label_order_dt <- overlap_dt[, .(mean_shared = mean(n_shared), cluster = cluster[1L]), by = label]
non_sing       <- sort(unique(label_order_dt[cluster > 0L, cluster]))
label_levels <- c(
  unlist(lapply(non_sing, function(cl)
    label_order_dt[cluster == cl][order(mean_shared), label]
  )),
  label_order_dt[cluster == 0L][order(mean_shared), label]
)

cluster_group_levels <- c(paste0("Cluster ", non_sing), "Singletons")
overlap_dt[, label         := factor(label,  levels = label_levels)]
overlap_dt[, out_id        := factor(out_id, levels = sort(unique(as.character(out_id))))]
overlap_dt[, cluster_group := factor(
  fifelse(cluster == 0L, "Singletons", paste0("Cluster ", cluster)),
  levels = cluster_group_levels
)]

# reference line data: total outcome instrument size per outcome
ref_dt <- out_snp_sets[, .(out_id, n_out_snps = lengths(snps))]
ref_dt[, out_id := factor(out_id, levels = levels(overlap_dt$out_id))]
ref_dt[, label  := sprintf("Total outcome\ninstrument\n(n=%i SNPs)", n_out_snps)]

cluster_palette <- c(
  "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#1B9E77", "#D95F02", "#7570B3"
)
cluster_cols <- setNames(
  c(rep_len(cluster_palette, length(non_sing)), "grey70"),
  cluster_group_levels
)

p <- ggplot(overlap_dt, aes(x = label, y = n_shared, fill = cluster_group)) +
  geom_col(width = 0.7) +
  geom_text(
    data = overlap_dt[n_shared > 0],
    aes(label = n_shared),
    vjust = -0.3, size = 2.2, colour = "grey20"
  ) +
  geom_hline(
    data     = ref_dt,
    aes(yintercept = n_out_snps),
    colour   = "grey30", linewidth = 0.4, linetype = "dashed"
  ) +
  geom_text(
    data     = ref_dt,
    aes(x = Inf, y = n_out_snps, label = label),
    hjust = 1, vjust = -0.2, size = 2.2, colour = "grey30",
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = cluster_cols, name = "Instrument cluster") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.25))) +
  facet_grid(out_id ~ ., scales = "free_y") +
  labs(x = "Metabolite (exposure instrument)", y = "Shared SNPs with outcome instrument") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 55, hjust = 1, size = 7),
    axis.text.y      = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.background = element_rect(fill = "grey20", colour = NA),
    strip.text       = element_text(colour = "white", face = "bold", size = 9),
    legend.position  = "right"
  )

n_metabs   <- uniqueN(overlap_dt$label)
n_outcomes <- uniqueN(overlap_dt$out_id)
fig_w <- max(6, 1.5 + n_metabs * 0.18)
fig_h <- max(4, n_outcomes * 2.5)

dir.create(dirname(overlap_heatmap), recursive = TRUE, showWarnings = FALSE)
ggsave(overlap_heatmap, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
tee(sprintf("Saved: %s\n", overlap_heatmap))
