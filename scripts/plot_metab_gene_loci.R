#### DESCRIPTION #######################################
# This script runs the locus plots for genes vs metabolite GWAS

#### SET INPUT #########################################
gwases       <- snakemake@params[["gwases"]]
genes        <- snakemake@params[["genes"]]
window_kb    <- as.integer(snakemake@params[["window_kb"]])
metab_dir    <- snakemake@input[["metab_dir"]]
cores        <- as.integer(snakemake@resources[["cpus_per_task"]])
locus_plot   <- snakemake@output[["locus_plot"]]
########################################################

# testing
if (FALSE) {
  cores <- 10
  window_kb <- 500
  locus_plot <- file.path(Sys.getenv("HF_METABOLITE_REPO2"), "output", "figures", "cis_mr", "test_gene_loci.png")
  metab_dir <- "/Users/xx20081/Documents/local_data/metabolite_gwas"
  gwases    <- list(
    valine =  list(
      path  = "metabolomic_karjalainen_2024/GCST90302122.fst",
      build = "Hg19"
    ),
    leucine = list(
      path  = "metabolomic_karjalainen_2024/GCST90301994.fst",
      build = "Hg19"
    ),
    kiv = list( # α-ketoisovalerate / 3-methyl-2-oxobutyrate
      path  = "metabolomic_chen_2023/GCST90199677.fst",
      build = "Hg38"
    )
  )
  genes  <- list(
    ppm1k = list(
      chr = 4,
      start = 89178772,
      end = 89205921,
      build = "Hg19"
    ),
    bckdk = list(
      chr = 16,
      start = 31117428,
      end = 31124110,
      build = "Hg19"
    )
  )

}

message("[DEBUG] gwases structure:"); print(str(gwases))
message("[DEBUG] genes structure:");  print(str(genes))
message("[DEBUG] metab_dir: ", metab_dir)
message("[DEBUG] locus_plot: ", locus_plot)

# requirements
message("[DEBUG] loading libraries...")
library(S7)
library(genepi.utils)
library(data.table)
library(fst)
library(ggplot2)
message("[DEBUG] libraries loaded")


# read metabolite GWAS data from fst files & lift if necessary
metab_gwases <- lapply(gwases, function(g) {
  fst_path <- file.path(metab_dir, g$path)
  message("[DEBUG] reading: ", fst_path)
  data <- read_fst(fst_path, as.data.table = TRUE)
  message("[DEBUG] read ", nrow(data), " rows; build=", g$build)
  data[, chromosome := as.character(chromosome)]
  if (g$build != "Hg19") {
    message("[DEBUG] lifting from ", g$build, " to Hg19")
    data <- genepi.utils::lift(data,
                              from    = g$build,
                              to      = "Hg19",
                              chr_col = "chromosome",
                              pos_col = "base_pair_location",
                              ea_col  = "effect_allele",
                              oa_col  = "other_allele",
                              remove_duplicates = FALSE)
    message("[DEBUG] lift done; ", nrow(data), " rows remain")
  }
  data
})
message("[DEBUG] all GWASes loaded: ", paste(names(metab_gwases), collapse=", "))


# extract the gene regions for each gwas, for each gene
regions <- list()

for (i in seq_along(genes)) {
  for (j in seq_along(metab_gwases)) {

    m_data <- metab_gwases[[j]]
    m_name <- names(metab_gwases)[j]

    g_name <- names(genes)[i]
    g_info <- genes[[i]]
    g_chr  <- as.character(g_info$chr)
    g_start<- g_info$start
    g_end  <- g_info$end
    win    <- window_kb * 1000

    message("[DEBUG] subsetting ", m_name, " @ ", g_name, " (chr", g_chr, ":", g_start, "-", g_end, ")")
    m_data <- m_data[chromosome == g_chr &
                       base_pair_location >= g_start - win &
                       base_pair_location <= g_start + win ] # see PMID: 35527238 for TSS-based windowing
    message("[DEBUG] ", nrow(m_data), " variants in window")

    m_data <-  m_data[, .(metabolite = factor(m_name, levels = names(metab_gwases)),
                          gene       = factor(g_name, levels = names(genes)),
                          in_gene    = base_pair_location >= g_start & base_pair_location <= g_end,
                          rsid       = variant_id,
                          chr        = chromosome,
                          bp         = as.integer(base_pair_location),
                          ea         = effect_allele,
                          oa         = other_allele,
                          eaf        = as.numeric(effect_allele_frequency),
                          beta       = as.numeric(beta),
                          se         = as.numeric(standard_error),
                          p          = as.numeric(p_value))]

    regions[[length(regions) + 1]] <- m_data
  }
}

regions_tbl <- rbindlist(regions)
message("[DEBUG] regions_tbl: ", nrow(regions_tbl), " rows, ", ncol(regions_tbl), " cols")

# build per-panel gene region rectangles for shading
gene_rects <- rbindlist(lapply(names(genes), function(g_name) {
  g_info <- genes[[g_name]]
  data.table(
    gene  = factor(g_name, levels = names(genes)),
    xmin  = g_info$start / 1e6,
    xmax  = g_info$end   / 1e6,
    ymin  = -Inf,
    ymax  = Inf
  )
}))

# significance thresholds (labelled once per facet column via geom_text)
sig_lines <- data.frame(
  yintercept = c(-log10(5e-8), -log10(5e-6)),
  label      = c("gws", "suggestive"),
  linetype   = c("dashed", "dotted"),
  color      = c("firebrick", "royalblue")
)
sig_labels <- c(
  gws        = expression("GWS " * (5 %*% 10^{-8})),
  suggestive = expression("Suggestive " * (5 %*% 10^{-6}))
)

# panel layout counts (used in facet_wrap and ggsave)
n_genes  <- length(genes)
n_gwases <- length(gwases)

# capitalise factor labels for display
regions_tbl[, gene_label       := factor(toupper(as.character(gene)),
                                          levels = toupper(names(genes)))]
fmt_metab <- function(x) ifelse(nchar(x) <= 3, toupper(x), tools::toTitleCase(x))
regions_tbl[, metabolite_label := factor(
  fmt_metab(as.character(metabolite)),
  levels = fmt_metab(names(metab_gwases))
)]
gene_rects[, gene_label := factor(toupper(as.character(gene)),
                                   levels = toupper(names(genes)))]

# create a combined panel factor: genes vary across cols, metabolites across rows
# order: for each metabolite (slow), cycle through genes (fast)
panel_levels <- unlist(lapply(levels(regions_tbl$metabolite_label), function(m) {
  paste0(levels(regions_tbl$gene_label), "\n", m)
}))
regions_tbl[, panel_label := factor(
  paste0(gene_label, "\n", metabolite_label),
  levels = panel_levels
)]

# expand gene_rects so each rect is matched to its specific metabolite panel
gene_rects <- merge(
  gene_rects[, .(gene_label, xmin, xmax, ymin, ymax)],
  unique(regions_tbl[, .(gene_label, panel_label)]),
  by = "gene_label"
)

p <- ggplot(regions_tbl, aes(x = bp / 1e6, y = -log10(p))) +
  # gene body shading
  geom_rect(
    data        = gene_rects,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill        = "steelblue", alpha = 0.10
  ) +
  # significance lines
  geom_hline(
    data        = sig_lines,
    aes(yintercept = yintercept, linetype = label, color = label),
    linewidth   = 0.6
  ) +
  scale_linetype_manual(
    name   = NULL,
    values = setNames(sig_lines$linetype, sig_lines$label),
    labels = sig_labels
  ) +
  scale_color_manual(
    name   = NULL,
    values = setNames(sig_lines$color, sig_lines$label),
    labels = sig_labels
  ) +
  # points: in-gene variants darker
  geom_point(
    aes(fill = in_gene),
    shape  = 21,
    size   = 1.4,
    stroke = 0.2,
    color  = "grey30",
    alpha  = 0.8
  ) +
  scale_fill_manual(
    values = c("TRUE" = "#2166ac", "FALSE" = "#d1d1d1"),
    labels = c("TRUE" = "In gene body", "FALSE" = "Flanking region"),
    name   = NULL
  ) +
  facet_wrap(
    ~ panel_label,
    ncol   = n_genes,
    nrow   = n_gwases,
    scales = "free"
  ) +
  labs(
    x = "Position (Mb)",
    y = expression(-log[10](italic(p)))
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background   = element_rect(fill = "grey92", colour = "grey70"),
    strip.text         = element_text(face = "bold", size = 10),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border       = element_rect(colour = "grey70"),
    axis.title         = element_text(size = 11),
    legend.position    = "bottom",
    legend.key.width   = unit(1.5, "cm"),
    legend.text        = element_text(size = 9),
    legend.box         = "horizontal"
  ) +
  guides(
    fill     = guide_legend(override.aes = list(size = 3), order = 1),
    color    = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )

# dynamic save dimensions: scale with number of genes (cols) and GWASes (rows)
panel_w  <- 3.5   # inches per gene column
panel_h  <- 2.5   # inches per GWAS row
margin_w <- 2.5   # left axis label + right margin
margin_h <- 2.0   # top + bottom margins + legend
fig_w    <- margin_w + n_genes  * panel_w
fig_h    <- margin_h + n_gwases * panel_h

message("[DEBUG] building plot...")
message("[DEBUG] saving to: ", locus_plot, " (", round(fig_w, 1), " x ", round(fig_h, 1), " in)")
dir.create(dirname(locus_plot), recursive = TRUE, showWarnings = FALSE)
ggsave(locus_plot, p, width = fig_w, height = fig_h, dpi = 300, bg = "white")
message("[DEBUG] done")


