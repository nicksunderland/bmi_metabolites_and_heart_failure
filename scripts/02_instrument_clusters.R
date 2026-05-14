#### DESCRIPTION #######################################
# Look at instrument redundancy
#
########################################################

#### SET INPUT #########################################
instrument_files    <- snakemake@input[["instrument_files"]]
consistent_ids_file <- snakemake@input[["consistent_ids_file"]]
name_map_file       <- snakemake@input[["name_map_file"]]
ld_block_file       <- snakemake@input[["ld_block_file"]]
overlap_fig         <- snakemake@output[["overlap_fig"]]
heatmap_fig         <- snakemake@output[["heatmap_fig"]]
ari_fig             <- snakemake@output[["ari_fig"]]
cluster_table       <- snakemake@output[["cluster_table"]]
cluster_bar_fig     <- snakemake@output[["cluster_bar_fig"]]
ld_blk_ari_fig      <- snakemake@output[["ld_blk_ari_fig"]]
ld_blk_overlap_fig  <- snakemake@output[["ld_blk_overlap_fig"]]
########################################################

if (FALSE) {
# ── TESTING BLOCK ────────────────────────────────────────────────────────────
repo_dir            <- Sys.getenv("HF_METABOLITE_REPO2")
instrument_files    <- list.files(file.path(repo_dir, "output/tables/instruments/genome_wide"), full.names = T)
consistent_ids_file <- file.path(repo_dir, "output/tables/replication/consistent_gwas_ids.txt")
name_map_file       <- file.path(repo_dir, "scripts/gwas_metab_name_map.xlsx")
  ld_block_file    <- file.path(repo_dir, "scripts/ld_blocks_hg19.tsv")
overlap_fig         <- file.path(repo_dir, "output/figures/instruments/overlap_graph.png")
heatmap_fig         <- file.path(repo_dir, "output/figures/instruments/overlap_heatmap.png")
ari_fig             <- file.path(repo_dir, "output/figures/instruments/cluster_ari_stability.png")
cluster_table       <- file.path(repo_dir, "output/tables/instruments/cluster_membership.tsv")
cluster_bar_fig     <- file.path(repo_dir, "output/figures/instruments/cluster_bar.png")
ld_blk_ari_fig      <- file.path(repo_dir, "output/figures/instruments/ld_block_ari_stability.png")
ld_blk_overlap_fig  <- file.path(repo_dir, "output/figures/instruments/ld_block_overlap_graph.png")
# ─────────────────────────────────────────────────────────────────────────────
}

library(ggplot2)
library(ggrepel)
library(patchwork)
library(readxl)
library(vegan)
library(igraph)
library(tidygraph)
library(ggraph)
library(mclust)
suppressPackageStartupMessages(library(data.table))


# set overlap threshold for clustering
thr <- 0.2  # overlap coefficient threshold: ≥33% of the smaller instrument shared




# meta data ====
map <- read_xlsx(name_map_file, sheet = 1) |> as.data.table()
consistent_ids <- fread(
  consistent_ids_file,
  header = FALSE,
  col.names = "gwas_id"
)
consistent_ids[
  map,
  c("label", "pathway_group") := .(i.label, i.pathway_group),
  on = "gwas_id"
]


# load LD data ====
if (!is.null(ld_block_file) && file.exists(ld_block_file)) {
  ld_blocks <- data.table::fread(ld_block_file)
} else {
  ld_blocks <- readRDS(system.file('extdata', 'LD.blocks.EUR.hg19.rds', package='mapgen'))
}

# load instruments ====
instruments <- lapply(instrument_files, function(fp) {
  id <- sub("(GCST[0-9]+)_instrument.tsv", "\\1", basename(fp))
  if (id %in% consistent_ids$gwas_id) {
    d <- fread(fp)
    if (nrow(d) == 0) {
      d <- data.table(trait = id)
    }
    d
  } else {
    NULL
  }
}) |> rbindlist(use.names = TRUE, fill = TRUE)
instruments[map, label := i.label, on = c("trait" = "gwas_id")]


# Build SNP sets per trait ====
snp_list <- instruments[
  !is.na(rsid),
  .(snps = list(unique(rsid))),
  by = trait
]


# Build binary SNP matrix ====
all_snps <- sort(unique(unlist(snp_list$snps)))
mat_snp <- matrix(
  0,
  nrow = nrow(snp_list),
  ncol = length(all_snps),
  dimnames = list(snp_list$trait, all_snps)
)
snp_index <- setNames(seq_along(all_snps), all_snps)
for (i in seq_len(nrow(snp_list))) {
  snps_i <- snp_list$snps[[i]]
  idx <- snp_index[snps_i]
  idx <- idx[!is.na(idx)]
  mat_snp[i, idx] <- 1
}


# SNP similarity (overlap coefficient = |A∩B| / min(|A|,|B|)) ====
int_mat <- mat_snp %*% t(mat_snp)          # intersection counts
size_vec <- rowSums(mat_snp)               # instrument sizes
# asymmetric overlap: intersection / size of each instrument
overlap_i <- sweep(int_mat, 1, size_vec, "/")   # divide by row size
overlap_j <- sweep(int_mat, 2, size_vec, "/")   # divide by col size
sim <- pmax(overlap_i, overlap_j)

diag(sim) <- 0


# Threshold graph ====
thr_seq <- quantile(
  sim[sim > 0],
  probs = seq(0, 0.99, by = 0.01),
  na.rm = TRUE
)

get_clusters <- function(thr) {
  g <- graph_from_adjacency_matrix(
    sim,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  g <- delete_edges(g, E(g)[weight < thr])
  if (ecount(g) == 0) {
    return(rep(1, vcount(g)))
  }
  cl <- cluster_louvain(g)$membership
  return(cl)
}
cluster_list <- lapply(thr_seq, get_clusters)
names(cluster_list) <- round(thr_seq, 3)

ari_vals <- numeric(length(cluster_list) - 1)

for (i in seq(2, length(cluster_list))) {
  ari_vals[i - 1] <- adjustedRandIndex(
    cluster_list[[i]],
    cluster_list[[i - 1]]
  )
}

p_ari <- ggplot(
  data.frame(threshold = thr_seq[-1], ARI = ari_vals),
  aes(x = threshold, y = ARI)
) +
  geom_line(colour = "grey40") +
  geom_point(size = 2, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Overlap coefficient threshold", y = "Cluster stability (ARI)") +
  theme_bw(base_size = 11)

dir.create(dirname(ari_fig), recursive = TRUE, showWarnings = FALSE)
ggsave(ari_fig, p_ari, width = 6, height = 4, dpi = 300, bg = "white")




g <- graph_from_adjacency_matrix(
  sim,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

g <- delete_edges(g, E(g)[weight < thr])
V(g)$label <- map$label[match(V(g)$name, map$gwas_id)]


# Singletons + clustering ====
deg <- degree(g)
is_singleton <- deg == 0
cl <- cluster_louvain(g)
clusters <- cl$membership
clusters[is_singleton] <- 0


# colours ====
cluster_ids <- sort(unique(clusters[clusters > 0]))

cluster_palette <- c(
  "#377EB8", "#E41A1C", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#66C2A5",
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854",
  "#FFD92F", "#1B9E77", "#D95F02", "#7570B3"
)
cols <- setNames(
  rep_len(cluster_palette, length(cluster_ids)),
  cluster_ids
)

V(g)$cluster <- clusters
V(g)$color <- ifelse(
  clusters == 0,
  "grey80",
  cols[as.character(clusters)]
)

# export cluster membership ====
cluster_export <- data.table(
  gwas_id = V(g)$name,
  label   = V(g)$label,
  cluster = as.integer(clusters),
  color   = V(g)$color
)
dir.create(dirname(cluster_table), recursive = TRUE, showWarnings = FALSE)
fwrite(cluster_export, cluster_table, sep = "\t")


# tidygraph ====
# label and color are already set as vertex attributes; as_tbl_graph carries them through
tg <- as_tbl_graph(g) %>%
  activate(nodes) %>%
  mutate(cluster = clusters)


# plot ====
p_graph <- ggraph(tg, layout = "fr", niter = 2000) +
  geom_edge_link(alpha = 0.2, colour = "grey60") +
  geom_node_point(aes(color = color), size = 3) +
  geom_node_text(aes(label = label), repel = TRUE, size = 2.8,
                 max.overlaps = 20) +
  scale_color_identity() +
  theme_graph(base_family = "sans")

dir.create(dirname(overlap_fig), recursive = TRUE, showWarnings = FALSE)
ggsave(overlap_fig, p_graph, width = 12, height = 12, dpi = 300, bg = "white")


# LD block clustering (parallel analysis) ====
ld_dt <- as.data.table(ld_blocks)
setnames(ld_dt, c("chr", "start", "end", "locus"))
ld_dt[, chr := as.character(chr)]

# locate chr/pos columns — prefer hg19-specific if present
chr_col <- "chr"
pos_col <- "bp"

inst_pos <- instruments[
  !is.na(rsid) & !is.na(get(chr_col)) & !is.na(get(pos_col)),
  .(trait, rsid,
    chr = sub("^chr", "", as.character(get(chr_col))),
    pos = as.integer(get(pos_col)))
][trait %in% snp_list$trait]

# assign each SNP to its LD block via non-equi join
inst_pos[ld_dt, locus := i.locus, on = .(chr, pos >= start, pos <= end)]

# build LD block sets per trait
blk_list <- inst_pos[
  !is.na(locus),
  .(blocks = list(unique(locus))),
  by = trait
]

# binary LD block matrix (traits × loci)
all_blocks <- sort(unique(unlist(blk_list$blocks)))
mat_blk <- matrix(
  0L,
  nrow = nrow(blk_list),
  ncol = length(all_blocks),
  dimnames = list(blk_list$trait, as.character(all_blocks))
)
blk_index <- setNames(seq_along(all_blocks), as.character(all_blocks))
for (i in seq_len(nrow(blk_list))) {
  idx <- blk_index[as.character(blk_list$blocks[[i]])]
  idx <- idx[!is.na(idx)]
  mat_blk[i, idx] <- 1L
}

# overlap coefficient (Szymkiewicz-Simpson) on LD blocks
int_blk  <- mat_blk %*% t(mat_blk)
size_blk <- rowSums(mat_blk)
ov_blk_i <- sweep(int_blk, 1, size_blk, "/")
ov_blk_j <- sweep(int_blk, 2, size_blk, "/")
sim_blk  <- pmax(ov_blk_i, ov_blk_j)
diag(sim_blk) <- 0

# ARI stability sweep over thresholds (LD blocks)
thr_seq_blk <- quantile(
  sim_blk[sim_blk > 0],
  probs = seq(0, 0.99, by = 0.01),
  na.rm = TRUE
)

get_clusters_blk <- function(thr) {
  g2 <- graph_from_adjacency_matrix(
    sim_blk, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  g2 <- delete_edges(g2, E(g2)[weight < thr])
  if (ecount(g2) == 0) return(rep(1L, vcount(g2)))
  cluster_louvain(g2)$membership
}

cluster_list_blk <- lapply(thr_seq_blk, get_clusters_blk)

ari_vals_blk <- numeric(length(cluster_list_blk) - 1)
for (i in seq(2, length(cluster_list_blk))) {
  ari_vals_blk[i - 1] <- adjustedRandIndex(
    cluster_list_blk[[i]],
    cluster_list_blk[[i - 1]]
  )
}

p_ari_blk <- ggplot(
  data.frame(threshold = thr_seq_blk[-1], ARI = ari_vals_blk),
  aes(x = threshold, y = ARI)
) +
  geom_line(colour = "grey40") +
  geom_point(size = 2, colour = "grey20") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Overlap coefficient threshold (LD blocks)",
       y = "Cluster stability (ARI)") +
  theme_bw(base_size = 11)

dir.create(dirname(ld_blk_ari_fig), recursive = TRUE, showWarnings = FALSE)
ggsave(ld_blk_ari_fig, p_ari_blk, width = 6, height = 4, dpi = 300, bg = "white")

# apply same threshold, derive LD block Louvain clusters
g_blk <- graph_from_adjacency_matrix(
  sim_blk, mode = "undirected", weighted = TRUE, diag = FALSE
)
g_blk       <- delete_edges(g_blk, E(g_blk)[weight < thr])
deg_blk     <- degree(g_blk)
cl_blk      <- cluster_louvain(g_blk)
clusters_blk <- cl_blk$membership
clusters_blk[deg_blk == 0] <- 0L

# ARI between SNP-based and LD block-based clusterings
snp_names <- V(g)$name
blk_names <- V(g_blk)$name
shared    <- intersect(snp_names, blk_names)
if (length(shared) > 0) {
  ari_cross <- adjustedRandIndex(
    clusters[match(shared, snp_names)],
    clusters_blk[match(shared, blk_names)]
  )
  message(sprintf("ARI (SNP-based vs LD block-based clusters): %.4f", ari_cross))
}

# LD block network graph plot
blk_cluster_ids <- sort(unique(clusters_blk[clusters_blk > 0]))
cols_blk <- setNames(
  rep_len(cluster_palette, length(blk_cluster_ids)),
  blk_cluster_ids
)
V(g_blk)$label   <- map$label[match(V(g_blk)$name, map$gwas_id)]
V(g_blk)$cluster <- clusters_blk
V(g_blk)$color   <- ifelse(
  clusters_blk == 0,
  "grey80",
  cols_blk[as.character(clusters_blk)]
)

tg_blk <- as_tbl_graph(g_blk) %>%
  activate(nodes) %>%
  mutate(cluster = clusters_blk)

p_graph_blk <- ggraph(tg_blk, layout = "fr", niter = 2000) +
  geom_edge_link(alpha = 0.2, colour = "grey60") +
  geom_node_point(aes(color = color), size = 3) +
  geom_node_text(aes(label = label), repel = TRUE, size = 2.8,
                 max.overlaps = 20) +
  scale_color_identity() +
  theme_graph(base_family = "sans")

ggsave(ld_blk_overlap_fig, p_graph_blk, width = 12, height = 12, dpi = 300, bg = "white")

# append cluster_blk and re-save cluster table
cluster_export[, cluster_blk := clusters_blk[match(gwas_id, blk_names)]]
fwrite(cluster_export, cluster_table, sep = "\t")


# heat map (ggplot2) ====
# Order: pathway_broad → non-singleton cluster → singleton → label
ht <- copy(cluster_export)
ht[consistent_ids, pathway_group := i.pathway_group, on = "gwas_id"]
ht[, pathway_broad := fcase(
  pathway_group == "Lipid metabolism",      "Lipid metabolism",
  pathway_group == "Lipoproteins",          "Lipoproteins",
  pathway_group == "Amino acid metabolism", "Amino acids",
  default = "Other"
)]
ht[, pathway_broad := factor(pathway_broad,
                             levels = c("Lipid metabolism", "Lipoproteins", "Amino acids", "Other"))]
ht[, is_singleton := (cluster == 0L)]
setorder(ht, pathway_broad, is_singleton, cluster, label)
ht[, row_idx := .I]
n_ht             <- nrow(ht)
heat_trait_order <- ht$gwas_id

sim_ht   <- sim[heat_trait_order, heat_trait_order]
sim_long <- melt(
  as.data.table(sim_ht, keep.rownames = "xid"),
  id.vars = "xid", variable.name = "yid", value.name = "sim_val"
)
sim_long[, xpos := match(xid, heat_trait_order)]
sim_long[, ypos := n_ht - match(yid, heat_trait_order) + 1L]

# group annotation positions (axis labels at group midpoints, separator lines between groups)
pw_pos <- ht[, .(pmin = min(row_idx) - 0.5, pmax = max(row_idx) + 0.5, mid = mean(row_idx)),
             by = pathway_broad]
setorder(pw_pos, pmin)
pw_pos[, ymid := n_ht + 1 - mid]
group_bounds <- pw_pos$pmax[-nrow(pw_pos)]

# cluster boxes: one box per (cluster × pathway_broad); label C1 if single group, C1a/C1b if multiple
clust_boxes <- ht[is_singleton == FALSE,
                  .(xmin = min(row_idx) - 0.5, xmax = max(row_idx) + 0.5,
                    ymin = n_ht - max(row_idx) + 0.5, ymax = n_ht - min(row_idx) + 1.5),
                  by = .(cluster, pathway_broad)]
setorder(clust_boxes, cluster, pathway_broad)
clust_boxes[, sub_label := {
  if (.N == 1L) paste0("C", cluster[1])
  else          paste0("C", cluster, letters[seq_len(.N)])
}, by = cluster]

p_hm <- ggplot(sim_long, aes(xpos, ypos, fill = sim_val)) +
  geom_raster() +
  geom_vline(xintercept = group_bounds, colour = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_hline(yintercept = n_ht + 1 - group_bounds, colour = "grey45", linewidth = 0.35, linetype = "dashed") +
  geom_rect(data = clust_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA, colour = "firebrick", linewidth = 0.6) +
  geom_text(data = clust_boxes,
            aes(x = xmin + 0.3, y = ymax - 0.3, label = sub_label),
            inherit.aes = FALSE,
            size = 1.8, colour = "firebrick", fontface = "bold", hjust = 0, vjust = 1) +
  scale_fill_gradientn(
    colours = c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"),
    values  = scales::rescale(c(0, 0.05, 0.2, 0.5, 1)),
    limits  = c(0, 1), name = "Overlap\ncoefficient"
  ) +
  scale_x_continuous(
    expand = c(0, 0), limits = c(0.5, n_ht + 0.5),
    breaks = pw_pos$mid, labels = pw_pos$pathway_broad
  ) +
  scale_y_continuous(
    expand = c(0, 0), limits = c(0.5, n_ht + 0.5),
    breaks = pw_pos$ymid, labels = pw_pos$pathway_broad
  ) +
  coord_fixed(clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8, colour = "grey20"),
    axis.text.y     = element_text(hjust  = 1, size = 8, colour = "grey20"),
    axis.ticks      = element_blank(),
    panel.grid      = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 8),
    legend.text     = element_text(size = 7),
    plot.margin     = margin(5, 5, 40, 80)
  )

dir.create(dirname(heatmap_fig), recursive = TRUE, showWarnings = FALSE)
ggsave(heatmap_fig, p_hm, width = 10, height = 10, dpi = 300, bg = "white")


# cluster bar chart ====
pg_levels  <- sort(unique(ht$pathway_group))
n_pg       <- length(pg_levels)
group_cols <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_pg),
  pg_levels
)
n_sing  <- sum(ht$is_singleton)
bar_dt  <- copy(ht)
bar_dt[, bar_label := ifelse(is_singleton,
                             paste0("Singletons\n(n=", n_sing, ")"),
                             paste0("C", cluster))]

clust_sizes <- bar_dt[is_singleton == FALSE, .N, by = cluster][order(-N)]
cluster_levels <- c(
  paste0("C", clust_sizes$cluster),
  paste0("Singletons\n(n=", n_sing, ")")
)
bar_dt[, bar_label := factor(bar_label, levels = cluster_levels)]

bar_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x        = element_text(size = 8, colour = "grey20"),
    axis.text.y        = element_text(size = 8, colour = "grey20"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank()
  )

n_clust <- nrow(clust_sizes)
max_bar <- max(clust_sizes$N)

# up to 5 metabolites per cluster (alphabetical selection)
annot_dt <- bar_dt[is_singleton == FALSE, {
  mets <- head(sort(label), 5)
  .(y_anch = .N,
    text   = ifelse(.N>5,
                    paste("e.g.\n", paste(mets, collapse = "\n")),
                    paste(mets, collapse = "\n")))
}, by = bar_label]

p_clusters <- ggplot(bar_dt[is_singleton == FALSE],
                     aes(x = bar_label, fill = pathway_group)) +
  geom_bar(width = 0.7, colour = "white", linewidth = 0.3) +
  ggrepel::geom_label_repel(
    data               = annot_dt,
    aes(x = bar_label, y = y_anch, label = text),
    inherit.aes        = FALSE,
    size               = 1.25,
    label.size         = 0.1,
    label.padding      = unit(0.10, "lines"),
    direction          = "y",
    nudge_y            = max_bar * 0.6,
    ylim               = c(max_bar * 1.01, NA),
    segment.colour     = "grey50",
    segment.size       = 0.2,
    min.segment.length = 0,
    fill               = "white",
    colour             = "grey20",
    lineheight         = 0.85
  ) +
  scale_fill_manual(values = group_cols, name = "Pathway group") +
  scale_x_discrete(limits = paste0("C", clust_sizes$cluster)) +
  scale_y_continuous(expand = expansion(mult = c(0, 1.2)), limits = c(0, NA)) +
  labs(x = "Cluster", y = "Number of metabolites") +
  bar_theme +
  theme(plot.margin = margin(5, 2, 5, 5), legend.position = "right")

p_singletons <- ggplot(bar_dt[is_singleton == TRUE],
                       aes(x = "Singletons", fill = pathway_group)) +
  geom_bar(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = group_cols, guide = "none") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  labs(x = NULL, y = paste0("n = ", n_sing)) +
  bar_theme +
  theme(
    axis.title.y = element_text(size = 8, colour = "grey40"),
    plot.margin  = margin(5, 5, 5, 2)
  )

p_bar <- p_clusters + p_singletons +
  plot_layout(widths = c(n_clust, 1))

ggsave(cluster_bar_fig, p_bar, width = 12, height = 4, dpi = 300, bg = "white")


