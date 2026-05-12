# ============================================================
# Empirical null vs observed score visualization
# - Fully customizable via arguments + `...` (ggplot layers)
# - Sensible defaults (publication-friendly)
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

#' Plot empirical null distribution vs observed score for a single feature
#'
#' Draws a histogram of the empirical null scores for one feature and overlays
#' the observed score as a vertical line. Users can fully customize appearance
#' via arguments and by passing additional ggplot layers through `...`.
#'
#' @param res Result from `pca_based_weighted_score(store_null=TRUE)`.
#'            Must contain `feature_table` and `null_scores`.
#' @param feature feature name (must be in `res$feature_table$feature` and `rownames(res$null_scores)`).
#'
#' @param bins Integer number of bins for histogram (default 30).
#' @param binwidth Optional numeric bin width. If provided, `binwidth` overrides `bins`.
#'
#' @param title Optional title; if NULL, a default title including p_emp is used.
#' @param subtitle Optional subtitle.
#' @param caption Optional caption.
#'
#' @param xlab X-axis label (default expression).
#' @param ylab Y-axis label.
#'
#' @param hist_fill Histogram fill color.
#' @param hist_color Histogram border color.
#' @param hist_alpha Histogram transparency.
#'
#' @param obs_line_color Color of observed vertical line.
#' @param obs_linetype Linetype for observed vertical line.
#' @param obs_linewidth Line width for observed vertical line.
#' @param show_obs_label Logical; if TRUE, label the observed score line.
#' @param obs_label_fmt Format string for observed label (uses `sprintf` with one numeric argument).
#' @param obs_label_y_frac Placement of observed label as a fraction of y-range.
#'
#' @param show_density Logical; if TRUE, overlay a density curve for the null distribution.
#' @param density_color Density line color.
#' @param density_linewidth Density line width.
#'
#' @param theme_fn ggplot theme function (e.g., `ggplot2::theme_minimal`, `ggplot2::theme_classic`).
#' @param base_size,title_size,axis_title_size,axis_text_size,legend_title_size,legend_text_size Font sizes.
#' @param legend_position Legend position.
#' @param panel_grid Logical; show grid lines.
#' @param panel_border Logical; draw a border.
#'
#' @param ... Additional ggplot layers to add (e.g., `theme()`, `scale_x_continuous()`, etc.).
#'
#' @return A ggplot object.
#' @export
#'
#' @examples
#' # Basic:
#' # plot_null_vs_observed_feature(res, feature="TAL1")
#' #
#' # Custom colors + bins:
#' # plot_null_vs_observed_feature(res, "TAL1", bins=40, hist_fill="grey85", obs_line_color="red")
#' #
#' # Add any ggplot styling via ...
#' # plot_null_vs_observed_feature(res, "TAL1",
#' #   ggplot2::theme_classic(base_size=16),
#' #   ggplot2::coord_cartesian(xlim=c(-2,2))
#' # )
plot_null_vs_observed <- function(
    res,
    feature,
    bins = 30,
    binwidth = NULL,
    title = NULL,
    subtitle = NULL,
    caption = NULL,
    xlab = expression(s[j]^{(null)}),
    ylab = "Count",
    hist_fill = "lightblue",
    hist_color = "black",
    hist_alpha = 0.9,
    obs_line_color = "red",
    obs_linetype = "dashed",
    obs_linewidth = 1,
    show_obs_label = TRUE,
    obs_label_fmt = "Observed = %.3f",
    obs_label_y_frac = 0.92,
    show_density = FALSE,
    density_color = "black",
    density_linewidth = 1,
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
  
  if (is.null(res$null_scores)) {
    stop("res$null_scores is NULL. Run pca_based_weighted_score(..., store_null=TRUE).")
  }
  if (is.null(res$feature_table)) stop("res$feature_table is missing.")
  tab <- res$feature_table
  
  if (!("feature" %in% colnames(tab))) stop("res$feature_table must contain a `feature` column.")
  if (!feature %in% tab$feature) stop("feature not found in res$feature_table: ", feature)
  if (!feature %in% rownames(res$null_scores)) stop("feature not found in res$null_scores: ", feature)
  
  if (!("sj" %in% colnames(tab))) stop("res$feature_table must contain `sj` (observed score).")
  if (!("p_emp" %in% colnames(tab))) stop("res$feature_table must contain `p_emp` (empirical p-value).")
  
  sj_obs  <- tab$sj[match(feature, tab$feature)]
  sj_null <- as.numeric(res$null_scores[feature, ])
  
  if (is.null(title)) {
    p_emp <- tab$p_emp[match(feature, tab$feature)]
    title <- paste0(feature, " - Observed vs Empirical Null \n(p=", signif(p_emp, 3), ")")
  }
  
  df <- data.frame(sj_null = sj_null)
  
  # Histogram base
  if (!is.null(binwidth)) {
    hist_layer <- ggplot2::geom_histogram(
      binwidth = binwidth, color = hist_color, fill = hist_fill, alpha = hist_alpha
    )
  } else {
    hist_layer <- ggplot2::geom_histogram(
      bins = bins, color = hist_color, fill = hist_fill, alpha = hist_alpha
    )
  }
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = sj_null)) +
    hist_layer +
    ggplot2::geom_vline(
      xintercept = sj_obs,
      linetype = obs_linetype,
      linewidth = obs_linewidth,
      color = obs_line_color
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      caption = caption,
      x = xlab,
      y = ylab
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
  
  # Optional: density overlay (scaled to counts)
  if (show_density) {
    # scale density to counts: density * N * binwidth
    bw <- binwidth
    if (is.null(bw)) {
      rng <- range(sj_null, finite = TRUE)
      bw <- (rng[2] - rng[1]) / max(1, bins)
    }
    n <- sum(is.finite(sj_null))
    p <- p + ggplot2::geom_density(
      ggplot2::aes(y = ggplot2::after_stat(density) * n * bw),
      color = density_color,
      linewidth = density_linewidth,
      fill = NA
    )
  }
  
  # Optional: label the observed score
  if (isTRUE(show_obs_label)) {
    b <- ggplot2::ggplot_build(p)
    y_max <- 0
    if (length(b$data) > 0 && "count" %in% names(b$data[[1]])) {
      y_max <- max(b$data[[1]]$count, na.rm = TRUE)
    }
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
    
    p <- p + ggplot2::annotate(
      "label",
      x = sj_obs,
      y = y_max * obs_label_y_frac,
      label = sprintf(obs_label_fmt, sj_obs),
      fill = "white",
      color = obs_line_color,
      size = 3.8
    )
  }
  
  # Let user override anything at the end
  p <- .apply_gg_dots(p, ...)
  
  p
}
