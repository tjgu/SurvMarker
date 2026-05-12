# R/plot_pca.R
# ============================================================
# PCA visualization helpers for the survival-guided package
# - Fully customizable via arguments + `...` (ggplot layers)
# - Default plots are clean (NO threshold lines)
# - Users can optionally add threshold annotations
# ============================================================

#' @noRd
.apply_gg_dots <- function(p, ...) {
  dots <- list(...)
  if (length(dots) == 0) return(p)
  for (layer in dots) p <- p + layer
  p
}

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

  if (!panel_grid) th <- th + ggplot2::theme(panel.grid = ggplot2::element_blank())
  if (panel_border) th <- th + ggplot2::theme(panel.border = ggplot2::element_rect(fill = NA))
  th
}

#' Scree plot of PCA (variance explained) with optional threshold annotation
#'
#' By default, this draws a clean scree plot (no vertical line).
#' Users can optionally add a vertical line at a chosen PC and display the
#' cumulative variance at that PC.
#'
#' @param res Result object from `pca_based_weighted_score()`. Must contain `pc_table`.
#' @param n_pcs Number of PCs to display (default 50).
#' @param title Plot title.
#' @param xlab,ylab Axis labels.
#'
#' @param show_threshold Logical; if TRUE, adds a vertical line at `threshold_pc` and a label.
#' @param threshold_pc Integer PC index for the vertical line (required if show_threshold=TRUE).
#' @param threshold_color Color for threshold line/label.
#' @param threshold_linetype Linetype for threshold line (e.g., "dashed").
#' @param threshold_linewidth Line width for threshold line.
#' @param label_y_frac Placement of the label as a fraction of y-range (0-1).
#' @param label_digits Digits for rounding the cumulative variance in the label.
#'
#' @param theme_fn ggplot theme function (e.g., `ggplot2::theme_minimal`, `ggplot2::theme_classic`).
#' @param base_size,title_size,axis_title_size,axis_text_size,legend_title_size,legend_text_size Font sizes.
#' @param legend_position Legend position.
#' @param panel_grid Logical; show grid lines.
#' @param panel_border Logical; draw a border.
#'
#' @param line_color,line_width Point/line color and width for the scree curve.
#' @param point_color,point_size Point color and size.
#'
#' @param ... Additional ggplot layers to add (e.g., `theme()`, `scale_x_continuous()`, etc.).
#'
#' @return A ggplot object.
#' @export
plot_scree <- function(
    res,
    n_pcs = 50,
    title = "Scree Plot of Principal Components",
    xlab = "Principal Component",
    ylab = "Proportion of Variance Explained",
    show_threshold = FALSE,
    threshold_pc = NULL,
    threshold_color = "red",
    threshold_linetype = "dashed",
    threshold_linewidth = 1,
    label_y_frac = 0.90,
    label_digits = 2,
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
    line_color = "black",
    line_width = 1,
    point_color = "black",
    point_size = 3,
    ...
) {

  if (is.null(res$pc_table)) stop("`res$pc_table` is missing.")
  df0 <- res$pc_table

  req_cols <- c("var_explained", "cum_var")
  if (!all(req_cols %in% colnames(df0))) {
    stop("`res$pc_table` must contain columns: ", paste(req_cols, collapse = ", "))
  }

  df0$PC_Index <- seq_len(nrow(df0))
  n_show <- min(n_pcs, nrow(df0))
  df <- df0[seq_len(n_show), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC_Index, y = var_explained)) +
    ggplot2::geom_line(color = line_color, linewidth = line_width) +
    ggplot2::geom_point(color = point_color, size = point_size) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
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

  if (show_threshold) {
    if (is.null(threshold_pc)) stop("threshold_pc must be provided when show_threshold = TRUE")
    if (!is.numeric(threshold_pc) || length(threshold_pc) != 1) stop("threshold_pc must be a single integer.")
    threshold_pc <- as.integer(threshold_pc)

    if (threshold_pc < 1 || threshold_pc > nrow(df0)) {
      stop("threshold_pc must be between 1 and ", nrow(df0))
    }

    # We allow threshold_pc beyond n_pcs shown, but the vline must be within visible x-range
    threshold_pc_plot <- min(threshold_pc, n_show)
    cum_var <- df0$cum_var[threshold_pc]

    y_max <- max(df$var_explained, na.rm = TRUE)
    y_lab <- y_max * label_y_frac

    p <- p +
      ggplot2::geom_vline(
        xintercept = threshold_pc_plot,
        linetype = threshold_linetype,
        color = threshold_color,
        linewidth = threshold_linewidth
      ) +
      ggplot2::annotate(
        "label",
        x = threshold_pc_plot,
        y = y_lab,
        label = paste0("PC", threshold_pc, " (Cumulative Variance = ", round(100 *cum_var, label_digits),"%)"),
        fill = "white",
        color = threshold_color,
        size = 4
      )
  }

  p <- .apply_gg_dots(p, ...)
  p
}

#' Cumulative variance explained plot with optional threshold annotation
#'
#' By default, this draws cumulative variance explained across PCs (no threshold lines).
#' Users can optionally add a horizontal line at `cum_threshold` and a vertical line at the
#' first PC reaching that threshold, with a label.
#'
#' @param res Result object from `pca_based_weighted_score()`. Must contain `pc_table`.
#' @param n_pcs Number of PCs to display (default 30).
#' @param title Plot title.
#' @param xlab,ylab Axis labels.
#'
#' @param show_threshold Logical; if TRUE, adds threshold lines + label.
#' @param cum_threshold Numeric in (0,1]; cumulative variance threshold (used only if show_threshold=TRUE).
#' @param threshold_color Color for threshold lines/label.
#' @param threshold_linetype Linetype for threshold lines.
#' @param threshold_linewidth Line width for threshold lines.
#' @param label_digits Digits for rounding the cumulative variance in the label.
#'
#' @param theme_fn ggplot theme function.
#' @param base_size,title_size,axis_title_size,axis_text_size,legend_title_size,legend_text_size Font sizes.
#' @param legend_position Legend position.
#' @param panel_grid Logical; show grid lines.
#' @param panel_border Logical; draw a border.
#'
#' @param line_color,line_width Point/line color and width for the curve.
#' @param point_color,point_size Point color and size.
#'
#' @param ... Additional ggplot layers to add.
#'
#' @return A ggplot object.
#' @export
plot_cumvar <- function(
    res,
    n_pcs = 30,
    title = "Cumulative Variance Explained",
    xlab = "Principal Component",
    ylab = "Cumulative Variance Explained",
    show_threshold = FALSE,
    cum_threshold = 0.50,
    threshold_color = "red",
    threshold_linetype = "dashed",
    threshold_linewidth = 1,
    label_digits = 2,
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
    line_color = "black",
    line_width = 1,
    point_color = "black",
    point_size = 3,
    ...
) {

  if (is.null(res$pc_table)) stop("`res$pc_table` is missing.")
  df0 <- res$pc_table

  if (!("cum_var" %in% colnames(df0))) stop("`res$pc_table` must contain column: cum_var")

  df0$PC_Index <- seq_len(nrow(df0))
  n_show <- min(n_pcs, nrow(df0))
  df <- df0[seq_len(n_show), , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = PC_Index, y = cum_var)) +
    ggplot2::geom_line(color = line_color, linewidth = line_width) +
    ggplot2::geom_point(color = point_color, size = point_size) +
    ggplot2::labs(title = title, x = xlab, y = ylab) +
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

  if (show_threshold) {
    if (!is.numeric(cum_threshold) || length(cum_threshold) != 1 || cum_threshold <= 0 || cum_threshold > 1) {
      stop("cum_threshold must be a single numeric value in (0,1].")
    }

    idx <- which(df0$cum_var >= cum_threshold)[1]
    if (is.na(idx)) idx <- nrow(df0)

    idx_plot <- min(idx, n_show)
    cum_at_idx <- df0$cum_var[idx]
    y_lab <- min(0.98, cum_at_idx + 0.05)

    p <- p +
      ggplot2::geom_hline(
        yintercept = cum_threshold,
        linetype = threshold_linetype,
        color = threshold_color,
        linewidth = threshold_linewidth
      ) +
      ggplot2::geom_vline(
        xintercept = idx_plot,
        linetype = threshold_linetype,
        color = threshold_color,
        linewidth = threshold_linewidth
      ) +
      ggplot2::annotate(
        "label",
        x = idx_plot,
        y = y_lab,
        label = paste0("PC", idx, " (Cumulative variance = ", round(100 *cum_at_idx, label_digits), "%)"),
        fill = "white",
        color = threshold_color,
        size = 4
      )
  }

  p <- .apply_gg_dots(p, ...)
  p
}
