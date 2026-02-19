# R/plot_overview.R
# ============================================================
# Overview panel: scree + cumvar + PC scatter + top-2 survival PCs
# - Fully customizable: theme, fonts, legend positions, palettes
# - Users can also pass extra ggplot layers to individual panels
# ============================================================

#' Overview panel of PCA and survival-guided results
#'
#' Generates a multi-panel overview of PCA and survival-guided results:
#' 1) Scree plot (variance explained)
#' 2) Cumulative variance explained
#' 3) PC1 vs PC2 scatter (optional coloring/shaping)
#' 4) Top-2 survival-significant PC scatter (optional coloring/shaping)
#'
#' This function is designed to be user-friendly and highly customizable.
#' You can control global styling (theme, fonts, legend position), palettes,
#' and also pass additional ggplot layers to specific panels via `*_add`.
#'
#' @param res Result object from `survival_guided_loading_score()`.
#' @param meta Optional metadata data.frame (rownames or `id_col` must match samples).
#' @param id_col Metadata column to use as rownames when meta rownames are missing.
#'
#' @param color_by Optional metadata column name used for coloring points in scatter plots.
#' @param shape_by Optional metadata column name used for shaping points in scatter plots.
#' @param color_values Optional named vector for manual discrete colors (scatter plots).
#' @param color_palette Optional vector of colors for continuous gradients (scatter plots when color_by is numeric).
#' @param na_color Color for NA values in scatter plots.
#'
#' @param n_pcs Number of PCs to display in scree/cumulative plots.
#' @param show_threshold Logical; if TRUE, show variance threshold annotations.
#' @param threshold_pc PC index for scree vertical line (used if show_threshold=TRUE).
#' @param cum_threshold Cumulative variance threshold (used if show_threshold=TRUE).
#'
#' @param main_title Optional overall title for the full panel (drawn above the grid).
#' @param ncol Number of columns in the arranged grid.
#'
#' @param theme_fn ggplot theme function applied to all panels (default `ggplot2::theme_minimal`).
#' @param base_size,title_size,axis_title_size,axis_text_size,legend_title_size,legend_text_size Global font sizes.
#' @param legend_position Global legend position for all panels (scatter plots respect this; scree/cumvar typically have no legend).
#' @param panel_grid Logical; show grid lines.
#' @param panel_border Logical; draw borders.
#'
#' @param scree_args Optional named list of extra arguments passed to `plot_scree()`.
#' @param cumvar_args Optional named list of extra arguments passed to `plot_cumvar()`.
#' @param scatter_args Optional named list of extra arguments passed to `plot_pc_scatter()` for PC1 vs PC2.
#' @param top2_args Optional named list of extra arguments passed to `plot_top2_survival_pcs()`.
#'
#' @param scree_add,cumvar_add,scatter_add,top2_add Optional lists of ggplot layers
#'        added to each panel (e.g., `list(ggplot2::coord_cartesian(...), ggplot2::theme(...))`).
#'
#' @return Invisibly returns the grob produced by `gridExtra::grid.arrange()`.
#' @export
#'
#' @examples
#' # Overview with ELN coloring + custom palette
#' # eln_cols <- c(Favorable="#1b9e77", Intermediate="#bdbdbd", Adverse="#d95f02")
#' # plot_overview(res, meta, color_by="ELN_2022_risk", color_values=eln_cols,
#' #               n_pcs=50, show_threshold=TRUE, threshold_pc=20, cum_threshold=0.8,
#' #               theme_fn=ggplot2::theme_classic, base_size=16, legend_position="bottom")
plot_overview <- function(
    res,
    meta = NULL,
    id_col = "PATIENT_ID",
    color_by = NULL,
    shape_by = NULL,
    color_values = NULL,
    color_palette = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"),
    na_color = "grey80",
    n_pcs = 30,
    show_threshold = FALSE,
    threshold_pc = NULL,
    cum_threshold = 0.50,
    main_title = NULL,
    ncol = 2,
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
    scree_args = list(),
    cumvar_args = list(),
    scatter_args = list(),
    top2_args = list(),
    scree_add = NULL,
    cumvar_add = NULL,
    scatter_add = NULL,
    top2_add = NULL
) {
  
  # ---------- Panel 1: Scree ----------
  base_scree <- list(
    res = res,
    n_pcs = n_pcs,
    title = "Scree Plot",
    show_threshold = show_threshold,
    threshold_pc = threshold_pc,
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
  p1 <- do.call(plot_scree, utils::modifyList(base_scree, scree_args))
  if (!is.null(scree_add)) for (layer in scree_add) p1 <- p1 + layer
  
  # ---------- Panel 2: CumVar ----------
  base_cum <- list(
    res = res,
    n_pcs = n_pcs,
    title = "Cumulative Variance",
    show_threshold = show_threshold,
    cum_threshold = cum_threshold,
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
  p2 <- do.call(plot_cumvar, utils::modifyList(base_cum, cumvar_args))
  if (!is.null(cumvar_add)) for (layer in cumvar_add) p2 <- p2 + layer
  
  # ---------- Panel 3: PC1 vs PC2 ----------
  base_scatter <- list(
    res = res,
    meta = meta,
    id_col = id_col,
    pc_x = 1,
    pc_y = 2,
    color_by = color_by,
    shape_by = shape_by,
    title = expression("PC1 vs PC2"),
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
    panel_border = panel_border
  )
  p3 <- do.call(plot_pc12, utils::modifyList(base_scatter, scatter_args))
  if (!is.null(scatter_add)) for (layer in scatter_add) p3 <- p3 + layer
  
  # ---------- Panel 4: Top-2 survival PCs ----------
  base_top2 <- list(
    res = res,
    meta = meta,
    id_col = id_col,
    color_by = color_by,
    shape_by = shape_by,
    title = "Top 2 Survival-Significant PCs",
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
    panel_border = panel_border
  )
  p4 <- do.call(plot_top2_survival_pcs, utils::modifyList(base_top2, top2_args))
  if (!is.null(top2_add)) for (layer in top2_add) p4 <- p4 + layer
  
  # ---------- Arrange ----------
  if (is.null(main_title)) {
    g <- gridExtra::grid.arrange(p1, p2, p3, p4, ncol = ncol)
    return(invisible(g))
  }
  
  # Add a global title above the grid (nice for reports)
  title_grob <- grid::textGrob(
    main_title,
    gp = grid::gpar(fontsize = title_size, fontface = "bold")
  )
  
  g <- gridExtra::grid.arrange(
    title_grob,
    gridExtra::arrangeGrob(p1, p2, p3, p4, ncol = ncol),
    ncol = 1,
    heights = grid::unit.c(grid::unit(1.2, "lines"), grid::unit(1, "null"))
  )
  
  invisible(g)
}


