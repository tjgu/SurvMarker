# ======================================================================
# PC scatter + Top-2 survival PCs
# Fully customizable, user-friendly ggplot functions
# ======================================================================

#' Internal: apply ggplot layers passed via ...
#' @noRd
.apply_gg_dots <- function(p, ...) {
  dots <- list(...)
  if (length(dots) == 0) return(p)
  for (layer in dots) p <- p + layer
  p
}

#' Internal: package-style theme builder (users can override via args or ...)
#' @noRd
.sg_theme <- function(
    theme_fn = ggplot2::theme_minimal,
    base_size = 14,
    title_size = 18,
    axis_title_size = 14,
    axis_text_size = 12,
    legend_title_size = 12,
    legend_text_size = 11,
    legend_position = "right",
    panel_grid = TRUE,
    panel_border = FALSE
) {
  th <- theme_fn(base_size = base_size) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(hjust = 0.5, size = title_size, face = "bold"),
      axis.title   = ggplot2::element_text(size = axis_title_size, face = "bold"),
      axis.text    = ggplot2::element_text(size = axis_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.text  = ggplot2::element_text(size = legend_text_size),
      legend.position = legend_position
    )

  if (!panel_grid) {
    th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  }
  if (panel_border) {
    th <- th + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA))
  }
  th
}

#' Internal: prepare metadata merge safely
#' @noRd
.merge_meta <- function(scores, meta = NULL, id_col = "PATIENT_ID") {
  scores <- as.data.frame(scores)
  scores$SampleID <- rownames(scores)

  if (is.null(meta)) return(scores)

  meta <- as.data.frame(meta)

  # If rownames are missing but id_col exists, use it
  if ((is.null(rownames(meta)) || all(is.na(rownames(meta))) || all(rownames(meta) == "")) &&
      id_col %in% colnames(meta)) {
    rownames(meta) <- as.character(meta[[id_col]])
  }

  meta$SampleID <- rownames(meta)
  merge(scores, meta, by = "SampleID", all.x = TRUE)
}

#' PC scatter plot (PCa vs PCb), fully customizable
#'
#' Creates a scatter plot of two selected PCs. Users can customize labels, fonts,
#' legends, colors, shapes, and add arbitrary ggplot layers via `...`.
#'
#' @param res Result object from `pca_based_weighted_score()`. Must contain `pc_scores`.
#' @param meta Optional metadata data.frame. If rownames are NULL/empty and `id_col` exists, it will be used as rownames.
#' @param id_col Metadata column to use as rownames when rownames are missing. Default "PATIENT_ID".
#' @param pc_x,pc_y PC indices for x- and y-axes.
#' @param color_by Optional column name in merged data to map to color.
#' @param shape_by Optional column name in merged data to map to shape.
#'
#' @param point_size Point size.
#' @param alpha Point transparency.
#' @param point_stroke Stroke for shapes 21–25 (ignored otherwise).
#'
#' @param title Plot title (character or expression()).
#' @param xlab,ylab Axis labels. If NULL, defaults to "PCk".
#' @param color_title,shape_title Legend titles. If NULL, uses `color_by` / `shape_by`.
#'
#' @param color_values Named vector for manual discrete colors (e.g., c("Low"="#1b9e77","High"="#d95f02")).
#' @param color_palette Vector of colors for continuous gradients (if `color_by` is numeric).
#' @param na_color Color for NA values in legends. Default "grey80".
#'
#' @param theme_fn ggplot theme function (e.g., `ggplot2::theme_minimal`, `ggplot2::theme_classic`).
#' @param base_size,title_size,axis_title_size,axis_text_size,legend_title_size,legend_text_size Font sizes.
#' @param legend_position Legend position ("right","left","top","bottom","none").
#' @param panel_grid Logical; show grid lines.
#' @param panel_border Logical; draw a border.
#'
#' @param ... Additional ggplot layers to add (e.g., `theme()`, `scale_color_*()`, `guides()`, `coord_*()`, etc.).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' # Basic:
#' # plot_pc12(res, pc_x=1, pc_y=2)
#' #
#' # With meta coloring and user theme:
#' # plot_pc12(res, meta, color_by="Risk",
#' #   theme_fn = ggplot2::theme_classic,
#' #   base_size = 16,
#' #   ggplot2::labs(subtitle = "TCGA-LAML"),
#' #   ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5)))
#' # )
plot_pc12 <- function(
    res,
    meta = NULL,
    id_col = "PATIENT_ID",
    pc_x = 1,
    pc_y = 2,
    color_by = NULL,
    shape_by = NULL,
    point_size = 3,
    alpha = 0.9,
    point_stroke = 0.6,
    title = NULL,
    xlab = NULL,
    ylab = NULL,
    color_title = NULL,
    shape_title = NULL,
    color_values = NULL,
    color_palette = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    na_color = "grey80",
    theme_fn = ggplot2::theme_minimal,
    base_size = 14,
    title_size = 18,
    axis_title_size = 14,
    axis_text_size = 12,
    legend_title_size = 12,
    legend_text_size = 11,
    legend_position = "right",
    panel_grid = TRUE,
    panel_border = FALSE,
    ...
) {

  if (is.null(res$pc_scores)) stop("`res$pc_scores` is missing.")

  scores <- .merge_meta(res$pc_scores, meta = meta, id_col = id_col)

  xname <- paste0("PC", pc_x)
  yname <- paste0("PC", pc_y)

  if (!xname %in% colnames(scores)) stop("PC not found in pc_scores: ", xname)
  if (!yname %in% colnames(scores)) stop("PC not found in pc_scores: ", yname)

  if (is.null(title)) title <- paste0(xname, " vs ", yname)
  if (is.null(xlab))  xlab  <- xname
  if (is.null(ylab))  ylab  <- yname
  if (is.null(color_title)) color_title <- color_by
  if (is.null(shape_title)) shape_title <- shape_by

  # Build aes using tidy evaluation (modern + robust)
  mapping <- ggplot2::aes(x = .data[[xname]], y = .data[[yname]])
  if (!is.null(color_by)) {
    if (!color_by %in% colnames(scores)) stop("color_by not found in merged data: ", color_by)
    mapping$colour <- rlang::sym(color_by)
  }
  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(scores)) stop("shape_by not found in merged data: ", shape_by)
    mapping$shape <- rlang::sym(shape_by)
  }

  p <- ggplot2::ggplot(scores, mapping) +
    ggplot2::geom_point(size = point_size, alpha = alpha, stroke = point_stroke) +
    ggplot2::labs(
      title = title,
      x = xlab,
      y = ylab,
      color = color_title,
      shape = shape_title
    ) +
    .sg_theme(
      theme_fn = theme_fn,
      base_size = base_size,
      title_size = title_size,
      axis_title_size = axis_title_size,
      axis_text_size = axis_text_size,
      legend_title_size = legend_title_size,
      legend_text_size = legend_text_size,
      legend_position = legend_position,
      panel_grid = panel_grid,
      panel_border = panel_border
    )

  # Auto-scale color: discrete vs continuous
  if (!is.null(color_by)) {
    v <- scores[[color_by]]

    if (is.numeric(v)) {
      p <- p + ggplot2::scale_color_gradientn(colors = color_palette, na.value = na_color)
    } else {
      if (!is.null(color_values)) {
        # If user provides some but not all levels, ggplot will recycle NA -> na.value
        p <- p + ggplot2::scale_color_manual(values = color_values, na.value = na_color)
      } else {
        # default discrete palette but still allow NA color
        p <- p + ggplot2::scale_color_discrete(na.translate = TRUE)
      }
    }
  }

  # Let user override ANYTHING at the end
  p <- .apply_gg_dots(p, ...)

  p
}


#' Plot the top-2 most survival-significant PCs (smallest padj), fully customizable
#'
#' Finds the two PCs with smallest adjusted p-values in `res$pc_table` and calls `plot_pc12()`.
#' All plot customization options are available and forwarded.
#'
#' @param res Result object from `pca_based_weighted_score()`. Must contain `pc_table` and `pc_scores`.
#' @inheritParams plot_pc12
#' @param title Plot title.
#'
#' @return A ggplot object.
#' @export
plot_top2_survival_pcs <- function(
    res,
    meta = NULL,
    id_col = "PATIENT_ID",
    color_by = NULL,
    shape_by = NULL,
    point_size = 3,
    alpha = 0.9,
    point_stroke = 0.6,
    title = "Top 2 Significant PCs",
    xlab = NULL,
    ylab = NULL,
    color_title = NULL,
    shape_title = NULL,
    color_values = NULL,
    color_palette = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    na_color = "grey80",
    theme_fn = ggplot2::theme_minimal,
    base_size = 14,
    title_size = 18,
    axis_title_size = 14,
    axis_text_size = 12,
    legend_title_size = 12,
    legend_text_size = 11,
    legend_position = "right",
    panel_grid = TRUE,
    panel_border = FALSE,
    ...
) {

  if (is.null(res$pc_table)) stop("`res$pc_table` is missing.")
  pc_tab <- res$pc_table

  if (!all(c("PC", "padj") %in% colnames(pc_tab))) {
    stop("`res$pc_table` must have columns: PC, padj")
  }

  pc_tab <- pc_tab[!is.na(pc_tab$padj), , drop = FALSE]
  pc_tab <- pc_tab[order(pc_tab$padj), , drop = FALSE]
  if (nrow(pc_tab) < 2) stop("Not enough PCs with non-NA padj to plot.")

  pc1 <- as.integer(sub("^PC", "", pc_tab$PC[1]))
  pc2 <- as.integer(sub("^PC", "", pc_tab$PC[2]))

  plot_pc12(
    res = res,
    meta = meta,
    id_col = id_col,
    pc_x = pc1,
    pc_y = pc2,
    color_by = color_by,
    shape_by = shape_by,
    point_size = point_size,
    alpha = alpha,
    point_stroke = point_stroke,
    title = title,
    xlab = xlab,
    ylab = ylab,
    color_title = color_title,
    shape_title = shape_title,
    color_values = color_values,
    color_palette = color_palette,
    na_color = na_color,
    theme_fn = theme_fn,
    base_size = base_size,
    title_size = title_size,
    axis_title_size = axis_title_size,
    axis_text_size = axis_text_size,
    legend_title_size = legend_title_size,
    legend_text_size = legend_text_size,
    legend_position = legend_position,
    panel_grid = panel_grid,
    panel_border = panel_border,
    ...
  )
}


