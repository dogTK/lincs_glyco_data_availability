#!/usr/bin/env Rscript
# Integrated ATC x WGCNA module figure (ComplexHeatmap reimplementation)
#
# R port of notebooks/analytics/15_integrated_atc_phh_module_figure.ipynb.
# Constants (MODULE_ORDER / MODULE_LABELS / MODULE_COLORS / ATC_ORDER /
# ATC_LABELS) and the helper logic (prepare_heatmap_matrix,
# prepare_delta_me_table, prepare_distribution_table, p_to_stars[_full]) are
# carried over verbatim so the underlying tables match the notebook output.
#
# Inputs (same as the notebook):
#   results/03_atc_pathway_gsea/atc_L1_module_gsea.csv
#   results/02_wgcna/module_eigengenes.csv
#   results/04_phh_similarity/phh_similarity_scores.csv
#
# Outputs:
#   results/figures/figure_integrated_atc_phh_module_complexheatmap/figure_integrated_atc_phh_module_complexheatmap.{png,pdf,svg}
#   results/figures/figure_integrated_atc_phh_module_complexheatmap_nes_matrix.tsv
#   results/figures/figure_integrated_atc_phh_module_complexheatmap/figure_integrated_atc_phh_module_complexheatmap_delta_me.tsv

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  if (requireNamespace("svglite", quietly = TRUE)) library(svglite)
})

# ============================================================
# Configuration (mirrors the Python notebook)
# ============================================================
resolve_base <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    script <- sub("--file=", "", file_arg)
    return(normalizePath(file.path(dirname(script), "..")))
  }
  getwd()
}
BASE_DIR <- resolve_base()
if (!dir.exists(file.path(BASE_DIR, "results"))) BASE_DIR <- getwd()

GSEA_FILE <- file.path(BASE_DIR, "results/03_atc_pathway_gsea/atc_L1_module_gsea.csv")
ME_FILE   <- file.path(BASE_DIR, "results/02_wgcna/module_eigengenes.csv")
SIM_FILE  <- file.path(BASE_DIR, "results/04_phh_similarity/phh_similarity_scores.csv")
STEM      <- "figure_integrated_atc_phh_module_complexheatmap"
OUTDIR    <- file.path(BASE_DIR, "results/figures", STEM)
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MODULE_ORDER <- c("ME2", "ME1", "ME18", "ME5", "ME3", "ME9", "ME7", "ME10", "ME14", "ME6")

# Multi-line display labels.
MODULE_LABELS <- c(
  ME2  = "ME2\nXenobiotic\nglucuronidation",
  ME1  = "ME1\nGlycan biosynthesis\nregulation",
  ME18 = "ME18\nGPI anchor\nbiosynthesis",
  ME5  = "ME5\nDolichol precursor\nbiosynthesis",
  ME3  = "ME3\nN-linked glycan\nlysosomal degradation",
  ME9  = "ME9\nKeratan sulfate\ndegradation",
  ME7  = "ME7\nUDP-GlcNAc\nsynthesis",
  ME10 = "ME10\nUDP glucose\nsynthesis",
  ME14 = "ME14\nGanglioside core\nstructure biosynthesis",
  ME6  = "ME6\nHyaluronic acid\nextracellular degradation"
)

MODULE_COLORS <- c(
  ME2  = "#2C5F8A",
  ME1  = "#6BAED6",
  ME18 = "#E67E22",
  ME5  = "#D9D9D9",
  ME3  = "#C0392B",
  ME9  = "#17A589",
  ME7  = "#7B241C",
  ME10 = "#5D4037",
  ME14 = "#A0785A",
  ME6  = "#C8C400"
)

ATC_ORDER <- c("N", "C", "D", "J", "L", "B", "P", "G", "A", "S", "M", "R")
ATC_LABELS <- c(
  N = "N (Nervous system)",
  C = "C (Cardiovascular)",
  D = "D (Dermatologicals)",
  J = "J (Antiinfectives)",
  L = "L (Antineoplastic)",
  B = "B (Blood)",
  P = "P (Antiparasitic)",
  G = "G (Genitourinary)",
  A = "A (Alimentary)",
  S = "S (Sensory organs)",
  M = "M (Musculoskeletal)",
  R = "R (Respiratory)"
)

COLOR_SIMILAR    <- "#D6604D"
COLOR_DISSIMILAR <- "#4393C3"
PHH_TOP_PCT      <- 0.10
PHH_BOTTOM_PCT   <- 0.10
NES_VMAX         <- 2.0
DIST_XLIM        <- c(-5.5, 5.5)

# ============================================================
# Helpers (ports of the notebook's Python helpers)
# ============================================================
p_to_stars <- function(q) {
  ifelse(is.na(q), "",
  ifelse(q < 0.001, "***",
  ifelse(q < 0.01,  "**",
  ifelse(q < 0.05,  "*", ""))))
}

p_to_stars_full <- function(q) {
  ifelse(is.na(q), "ns",
  ifelse(q < 0.001, "***",
  ifelse(q < 0.01,  "**",
  ifelse(q < 0.05,  "*", "ns"))))
}

prepare_heatmap_matrix <- function(gsea_df, module_order, atc_order) {
  gsea_filt <- gsea_df %>% filter(.data$atc_class %in% atc_order)
  nes_wide <- gsea_filt %>%
    select(.data$pathway, .data$atc_class, .data$nes) %>%
    pivot_wider(names_from = "atc_class", values_from = "nes")
  q_wide <- gsea_filt %>%
    select(.data$pathway, .data$atc_class, .data$q_value) %>%
    pivot_wider(names_from = "atc_class", values_from = "q_value")

  ridx <- match(module_order, nes_wide$pathway)
  nes_mat <- as.matrix(nes_wide[ridx, atc_order, drop = FALSE])
  q_mat   <- as.matrix(q_wide  [ridx, atc_order, drop = FALSE])
  rownames(nes_mat) <- module_order
  rownames(q_mat)   <- module_order
  list(nes = nes_mat, q = q_mat)
}

prepare_delta_me_table <- function(me_df, sim_df, module_order,
                                   top_pct = 0.10, bottom_pct = 0.10) {
  sim_sorted <- sim_df %>% arrange(desc(.data$similarity_spearman))
  n_total    <- nrow(sim_sorted)
  n_top      <- floor(n_total * top_pct)
  n_bottom   <- floor(n_total * bottom_pct)
  top_ids    <- sim_sorted$compound_id[seq_len(n_top)]
  bottom_ids <- sim_sorted$compound_id[(n_total - n_bottom + 1):n_total]

  me <- me_df %>%
    mutate(group = case_when(
      .data$pertname %in% top_ids    ~ "PHH-similar",
      .data$pertname %in% bottom_ids ~ "PHH-dissimilar",
      TRUE                            ~ "middle"
    ))

  rows <- lapply(module_order, function(mod) {
    if (!mod %in% colnames(me)) return(NULL)
    s <- me[me$group == "PHH-similar",    mod, drop = TRUE]
    d <- me[me$group == "PHH-dissimilar", mod, drop = TRUE]
    wt <- suppressWarnings(
      wilcox.test(s, d, alternative = "two.sided", exact = FALSE)
    )
    # Rank-biserial correlation: r = 2U / (n1 * n2) - 1, bounded in [-1, 1].
    # In R, wilcox.test(x, y)$statistic counts pairs with x_i > y_j, so this
    # becomes positive when PHH-similar (x) tends to rank above PHH-dissimilar.
    rbc <- 2 * as.numeric(wt$statistic) / (length(s) * length(d)) - 1
    data.frame(
      module          = mod,
      n_similar       = length(s),
      n_dissimilar    = length(d),
      mean_similar    = mean(s),
      mean_dissimilar = mean(d),
      delta_me        = mean(s) - mean(d),
      rank_biserial   = rbc,
      p_value         = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  res <- bind_rows(rows)
  res$q_value       <- p.adjust(res$p_value, method = "BH")
  res$effect_tier   <- effect_size_tier(res$rank_biserial)
  res$stars         <- p_to_stars_full(res$q_value)  # kept for provenance/TSV
  res
}

# Rank-biserial |r| tiers (Cohen-style): negligible / small / medium / large.
effect_size_tier <- function(r) {
  a <- abs(r)
  ifelse(is.na(r),  "",
  ifelse(a < 0.10,  "ns",
  ifelse(a < 0.30,  "+",
  ifelse(a < 0.50,  "++", "+++"))))
}

prepare_distribution_table <- function(me_df, sim_df, module_order,
                                       top_pct = 0.10, bottom_pct = 0.10) {
  sim_sorted <- sim_df %>% arrange(desc(.data$similarity_spearman))
  n_total    <- nrow(sim_sorted)
  top_ids    <- sim_sorted$compound_id[seq_len(floor(n_total * top_pct))]
  bottom_ids <- sim_sorted$compound_id[(n_total - floor(n_total * bottom_pct) + 1):n_total]

  me <- me_df %>%
    mutate(group = case_when(
      .data$pertname %in% top_ids    ~ "PHH-similar",
      .data$pertname %in% bottom_ids ~ "PHH-dissimilar",
      TRUE                            ~ NA_character_
    )) %>%
    filter(!is.na(.data$group))

  module_cols <- intersect(module_order, colnames(me))
  me %>%
    select("group", all_of(module_cols)) %>%
    pivot_longer(cols = -"group", names_to = "module", values_to = "me_value")
}

# Build a (n_modules x max_samples) matrix of ME values for a single PHH group,
# padded with NA so that anno_boxplot can plot one boxplot per row.
build_per_row_value_matrix <- function(dist_table, group_name, module_order) {
  values_list <- lapply(module_order, function(m) {
    dist_table %>%
      filter(.data$module == m, .data$group == group_name) %>%
      pull(.data$me_value)
  })
  max_n <- max(vapply(values_list, length, integer(1)))
  mat <- matrix(NA_real_, nrow = length(module_order), ncol = max_n)
  for (i in seq_along(values_list)) {
    v <- values_list[[i]]
    if (length(v)) mat[i, seq_along(v)] <- v
  }
  rownames(mat) <- module_order
  mat
}

# ============================================================
# Heatmap builder
# ============================================================
build_integrated_heatmap <- function(nes_mat, q_mat, delta_table, dist_table,
                                     module_order, atc_order,
                                     module_labels, atc_labels, module_colors) {
  nes_col <- colorRamp2(c(-NES_VMAX, 0, NES_VMAX),
                        c("#2166AC", "white", "#B2182B"))

  cell_fun <- function(j, i, x, y, width, height, fill) {
    qv <- q_mat[i, j]
    nv <- nes_mat[i, j]
    star <- p_to_stars(qv)
    if (nzchar(star)) {
      tcol <- if (!is.na(nv) && abs(nv) > 0.8) "white" else "black"
      grid.text(star, x, y,
                gp = gpar(fontsize = 11, col = tcol, fontface = "bold"))
    }
  }

  # Left: WGCNA module color strip
  left_anno <- rowAnnotation(
    Module = module_order,
    col = list(Module = module_colors),
    show_annotation_name = FALSE,
    show_legend          = FALSE,
    simple_anno_size     = unit(0.45, "cm")
  )

  # Right side annotations
  rbc_vec <- delta_table$rank_biserial[match(module_order, delta_table$module)]

  sim_mat <- build_per_row_value_matrix(dist_table, "PHH-similar",    module_order)
  dis_mat <- build_per_row_value_matrix(dist_table, "PHH-dissimilar", module_order)

  # Combined density-ridge annotation (PHH-similar + PHH-dissimilar per row).
  # Per-row densities are normalized to the same peak height so that rows remain
  # visually comparable regardless of absolute variance.
  paired_density_fun <- local({
    sim   <- sim_mat
    dis   <- dis_mat
    xlim  <- DIST_XLIM
    fill_s <- COLOR_SIMILAR
    fill_d <- COLOR_DISSIMILAR
    ticks  <- c(-4, -2, 0, 2, 4)

    # Pre-compute per-row densities + a global peak so all rows share a scale.
    compute_density <- function(v) {
      v <- v[!is.na(v)]
      if (length(v) < 2) return(NULL)
      d <- stats::density(v, from = xlim[1], to = xlim[2], n = 256,
                          bw = "nrd0")
      list(x = d$x, y = d$y)
    }
    n_rows <- nrow(sim)
    dens_s <- lapply(seq_len(n_rows), function(i) compute_density(sim[i, ]))
    dens_d <- lapply(seq_len(n_rows), function(i) compute_density(dis[i, ]))
    peak <- max(
      vapply(c(dens_s, dens_d),
             function(d) if (is.null(d)) 0 else max(d$y), numeric(1))
    )

    function(index) {
      n <- length(index)
      pushViewport(viewport(
        xscale = xlim, yscale = c(0.5, n + 0.5)
      ))
      # Row separators + vertical reference line
      for (k in seq_len(n - 1)) {
        grid.lines(x = unit(c(xlim[1], xlim[2]), "native"),
                   y = unit(c(k + 0.5, k + 0.5), "native"),
                   gp = gpar(col = "#E5E5E5", lwd = 0.4))
      }
      grid.lines(x = unit(c(0, 0), "native"),
                 y = unit(c(0.5, n + 0.5), "native"),
                 gp = gpar(col = "#555555", lty = 2, lwd = 0.6))

      # Height scaling: each density peak maps to `ridge_h` row units.
      ridge_h <- 0.82
      baseline_off <- 0.0  # both ridges share the same row baseline

      draw_ridge <- function(d, y_base, fill) {
        if (is.null(d)) return(invisible())
        y_scaled <- d$y / peak * ridge_h
        # Encode alpha into the hex color (works reliably across PNG/PDF/SVG
        # devices; relying on gpar(alpha=) can produce solid black on some
        # SVG backends).
        fill_rgba    <- grDevices::adjustcolor(fill, alpha.f = 0.55)
        outline_rgba <- grDevices::adjustcolor(fill, alpha.f = 0.95)
        xs <- c(d$x, rev(d$x))
        ys <- c(y_base + rep(0, length(d$x)), rev(y_base + y_scaled))
        grid.polygon(x = unit(xs, "native"),
                     y = unit(ys, "native"),
                     gp = gpar(fill = fill_rgba, col = NA))
        grid.lines(x = unit(d$x, "native"),
                   y = unit(y_base + y_scaled, "native"),
                   gp = gpar(col = outline_rgba, lwd = 0.9))
      }

      for (k in seq_along(index)) {
        i <- index[k]
        y <- n - k + 1 - 0.35  # baseline slightly below row center
        # Draw dissimilar first so PHH-similar (the focal group) sits on top
        draw_ridge(dens_d[[i]], y + baseline_off, fill_d)
        draw_ridge(dens_s[[i]], y + baseline_off, fill_s)
      }

      grid.xaxis(at = ticks, gp = gpar(fontsize = 8), main = TRUE)
      popViewport()
    }
  })

  rbc_text <- sprintf("%.2f", rbc_vec)
  right_anno <- rowAnnotation(
    `Module Eigengene (ME)` = AnnotationFunction(
      fun           = paired_density_fun,
      fun_name      = "paired_density_ridge",
      which         = "row",
      width         = unit(5.5, "cm"),
      n             = length(module_order),
      subsettable   = TRUE,
      var_import    = list()
    ),
    `r` = anno_text(
      rbc_text,
      gp       = gpar(fontsize = 10, fontface = "bold"),
      just     = "center",
      location = unit(0.5, "npc"),
      width    = max_text_width(rbc_text, gp = gpar(fontsize = 10, fontface = "bold")) +
                 unit(4, "mm")
    ),
    annotation_name_gp      = gpar(fontsize = 10, fontface = "bold"),
    annotation_name_side    = "top",
    annotation_name_rot     = 0,
    annotation_name_offset  = unit(2, "mm"),
    gap                     = unit(3, "mm")
  )

  ht <- Heatmap(
    nes_mat,
    name = "GSEA NES",
    col  = nes_col,

    cluster_rows    = FALSE,
    cluster_columns = FALSE,

    row_labels       = module_labels[module_order],
    row_names_side   = "left",
    row_names_gp     = gpar(fontsize = 9, lineheight = 1.05),
    row_names_max_width = max_text_width(module_labels[module_order],
                                         gp = gpar(fontsize = 9, lineheight = 1.05)) +
                         unit(2, "mm"),
    row_title        = "WGCNA Module (ME)",
    row_title_gp     = gpar(fontsize = 12, fontface = "bold"),
    row_title_side   = "left",

    column_labels           = atc_order,
    column_names_side       = "top",
    column_names_rot        = 0,
    column_names_centered   = TRUE,
    column_names_gp         = gpar(fontsize = 10),
    column_names_max_height = unit(5.0, "cm"),
    column_title            = "ATC Class (L1)",
    column_title_gp         = gpar(fontsize = 12, fontface = "bold"),
    column_title_side       = "bottom",

    rect_gp  = gpar(col = "white", lwd = 0.5),
    cell_fun = cell_fun,

    left_annotation  = left_anno,
    right_annotation = right_anno,

    heatmap_legend_param = list(
      title            = "GSEA NES",
      title_gp         = gpar(fontsize = 11, fontface = "bold"),
      labels_gp        = gpar(fontsize = 10),
      at               = c(-2, -1, 0, 1, 2),
      color_bar        = "continuous",
      legend_direction = "horizontal",
      legend_width     = unit(4, "cm"),
      title_position   = "topcenter"
    ),

    width  = unit(11, "cm"),
    height = unit(length(module_order) * 1.15, "cm")
  )

  list(ht = ht)
}

# ============================================================
# Extra legends + draw / save
# ============================================================
make_extra_legends <- function() {
  group_lgd <- Legend(
    labels    = c("PHH-similar (top 10%)", "PHH-dissimilar (bottom 10%)"),
    legend_gp = gpar(fill = c(COLOR_SIMILAR, COLOR_DISSIMILAR), col = "black"),
    title     = "PHH similarity group",
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9),
    border    = "black"
  )

  # Heatmap cell stars (GSEA q_value from atc_L1_module_gsea.csv)
  gsea_sig_lgd <- Legend(
    labels    = c("*    q < 0.05", "**   q < 0.01", "***  q < 0.001"),
    type      = "grid",
    legend_gp = gpar(fill = "white", col = "black"),
    border    = "black",
    title     = "Heatmap cell\n(GSEA FDR)",
    title_gp  = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )

  packLegend(group_lgd, gsea_sig_lgd,
             direction = "vertical", gap = unit(4, "mm"))
}

draw_and_decorate <- function(ht, extra_lgd, title) {
  draw(
    ht,
    heatmap_legend_side    = "bottom",
    annotation_legend_side = "right",
    annotation_legend_list = list(extra_lgd),
    merge_legends          = FALSE,
    padding                = unit(c(8, 8, 12, 8), "mm"),
    column_title           = title,
    column_title_gp        = gpar(fontsize = 13, fontface = "bold")
  )
}

save_outputs <- function(ht, extra_lgd, outdir, stem,
                         width_in = 12, height_in = 7) {
  title <- ""

  png(file.path(outdir, paste0(stem, ".png")),
      width = width_in * 300, height = height_in * 300, res = 300)
  draw_and_decorate(ht, extra_lgd, title); dev.off()
  message("PNG saved.")

  pdf_path <- file.path(outdir, paste0(stem, ".pdf"))
  if (.Platform$OS.type == "unix" && Sys.info()[["sysname"]] == "Darwin") {
    quartz(type = "pdf", file = pdf_path, width = width_in, height = height_in)
  } else {
    pdf(pdf_path, width = width_in, height = height_in)
  }
  draw_and_decorate(ht, extra_lgd, title); dev.off()
  message("PDF saved.")

  svg_path <- file.path(outdir, paste0(stem, ".svg"))
  tryCatch({
    if (requireNamespace("svglite", quietly = TRUE)) {
      svglite::svglite(svg_path, width = width_in, height = height_in)
    } else {
      svg(svg_path, width = width_in, height = height_in)
    }
    draw_and_decorate(ht, extra_lgd, title); dev.off()
    # svglite emits properties via inline `style='...'` which some vector
    # editors (notably Affinity Designer) mishandle for `fill-opacity`,
    # rendering translucent polygons as solid black. Flatten each style
    # block into direct XML attributes for maximal compatibility.
    flatten_svglite_styles(svg_path)
    message("SVG saved.")
  }, error = function(e) message("SVG skipped (", conditionMessage(e), ")"))
}

# Affinity Designer ignores `fill-opacity` when it lives inside svglite's
# inline `style='...'` attribute (rendering translucent fills as solid
# black). Every style block we touch is left intact for compliant renderers
# (browsers, Inkscape, rsvg-convert), but we additionally duplicate
# `fill-opacity` as a direct XML attribute on the same element so Affinity
# picks it up from the attribute parser.
flatten_svglite_styles <- function(svg_path) {
  txt <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  # Match any element style block that contains fill-opacity
  pat <- "style='([^']*fill-opacity:\\s*([0-9.]+)[^']*)'"
  txt2 <- gsub(pat,
               "fill-opacity='\\2' style='\\1'",
               txt, perl = TRUE)
  writeLines(txt2, svg_path)
}

# ============================================================
# Main
# ============================================================
if (sys.nframe() == 0) {
  message("Loading data...")
  gsea_df <- read_csv(GSEA_FILE, show_col_types = FALSE)
  me_df   <- read_csv(ME_FILE,   show_col_types = FALSE)
  sim_df  <- read_csv(SIM_FILE,  show_col_types = FALSE)
  message(sprintf("GSEA records: %d", nrow(gsea_df)))
  message(sprintf("Module eigengenes: %d x %d", nrow(me_df), ncol(me_df)))
  message(sprintf("PHH similarity: %d x %d", nrow(sim_df), ncol(sim_df)))

  message("Preparing intermediate tables...")
  hm <- prepare_heatmap_matrix(gsea_df, MODULE_ORDER, ATC_ORDER)
  delta_table <- prepare_delta_me_table(me_df, sim_df, MODULE_ORDER,
                                        PHH_TOP_PCT, PHH_BOTTOM_PCT)
  dist_table  <- prepare_distribution_table(me_df, sim_df, MODULE_ORDER,
                                            PHH_TOP_PCT, PHH_BOTTOM_PCT)
  message(sprintf("NES matrix: %d x %d", nrow(hm$nes), ncol(hm$nes)))
  message("Delta ME table:")
  print(delta_table[, c("module", "delta_me", "q_value", "stars")], row.names = FALSE)

  nes_out <- as.data.frame(hm$nes) %>% rownames_to_column("module")
  write_tsv(nes_out, file.path(OUTDIR, paste0(STEM, "_nes_matrix.tsv")))
  write_tsv(delta_table, file.path(OUTDIR, paste0(STEM, "_delta_me.tsv")))
  message("Source data saved.")

  message("Building heatmap...")
  built <- build_integrated_heatmap(hm$nes, hm$q, delta_table, dist_table,
                                    MODULE_ORDER, ATC_ORDER,
                                    MODULE_LABELS, ATC_LABELS, MODULE_COLORS)

  extra_lgd <- make_extra_legends()

  message("Saving outputs...")
  save_outputs(built$ht, extra_lgd, OUTDIR, STEM)
  message("Done.")
}
