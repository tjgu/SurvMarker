#' Run survival-guided selection across multiple PC cutoffs
#'
#' @param X Numeric matrix (n_samples x p_features).
#' @param time Survival time vector (length n_samples).
#' @param status Event indicator 0/1 (length n_samples).
#' @param covar Optional data.frame of covariates (n_samples rows).
#' @param pcs_to_run Integer vector of requested PC counts.
#' @param max_pcs Hard cap on PCs (passed to pca_based_weighted_score).
#' @param pc_fdr_cutoff PC-level FDR cutoff.
#' @param feature_fdr_cutoff feature-level FDR cutoff.
#' @param null_B Number of null resamples.
#' @param seed Random seed.
#' @param scale_pca Logical.
#' @param store_null Logical.
#' @param verbose Logical.
#'
#' @return A list with feature_sets, cumvar_map, k_used_map, settings.
#' @export
run_survival_pca_multi_pc <- function(
    X, time, status, covar = NULL,
    pcs_to_run = c(10, 20, 30, 50),
    max_pcs = 2000,
    pc_fdr_cutoff = 0.05,
    feature_fdr_cutoff = 0.05,
    null_B = 500,
    seed = 1,
    scale_pca = TRUE,
    store_null = FALSE,
    verbose = TRUE
) {
  pcs_to_run <- as.integer(unique(pcs_to_run))
  if (length(pcs_to_run) < 2) stop("pcs_to_run must contain at least 2 values.")

  feature_lists <- list()
  cumvar_map <- setNames(rep(NA_real_, length(pcs_to_run)), paste0("PC", pcs_to_run))
  k_used_map <- setNames(rep(NA_integer_, length(pcs_to_run)), paste0("PC", pcs_to_run))

  for (k in pcs_to_run) {
    if (verbose) message("Running n_pcs = ", k)

    res_k <- pca_based_weighted_score(
      X = X, time = time, status = status, covar = covar,
      n_pcs = k, max_pcs = max_pcs,
      pc_fdr_cutoff = pc_fdr_cutoff,
      feature_fdr_cutoff = feature_fdr_cutoff,
      null_B = null_B, seed = seed,
      scale_pca = scale_pca,
      store_null = store_null,
      verbose = verbose
    )

    key <- paste0("PC", k)
    feature_lists[[key]] <- res_k$selected_features

    if (!is.null(res_k$pc_table) && "cum_var" %in% colnames(res_k$pc_table)) {
      k_avail <- nrow(res_k$pc_table)
      k_used_map[[key]] <- k_avail
      cumvar_map[[key]] <- res_k$pc_table$cum_var[k_avail]
    }
  }

  list(
    feature_sets = feature_lists,
    cumvar_map = cumvar_map,
    k_used_map = k_used_map,
    settings = list(
      pcs_to_run = pcs_to_run, max_pcs = max_pcs,
      pc_fdr_cutoff = pc_fdr_cutoff, feature_fdr_cutoff = feature_fdr_cutoff,
      null_B = null_B, seed = seed, scale_pca = scale_pca
    )
  )
}

#' Plot Venn diagram of selected feature sets across PC cutoffs
#'
#' @param pcset_obj Output from run_survival_pca_multi_pc()
#' @param pick Integer vector of requested PC counts to plot (e.g., c(100,200,500)).
#'        Must contain between 2 and 5 values.
#'
#' @param label_mode What to show in each set label. One of:
#'        "both" (default), "n", "cum_var", "none".
#' @param show_used Logical; if TRUE, add "used=K" (actual PCs available) to labels.
#'
#' @param main Main title.
#' @param fill Fill colors. If NULL, uses defaults. If provided, must have length >= number of sets.
#' @param alpha Transparency.
#' @param cex Number size.
#' @param cat_cex Label size.
#' @param margin Venn margin.
#'
#' @return The grob returned by VennDiagram::venn.diagram()
#' @export
plot_venn <- function(
    pcset_obj,
    pick,
    label_mode = c("both", "n", "cum_var", "none"),
    show_used = FALSE,
    main = "Increasing PC Count Refines the Selected Feature Set",
    fill = NULL,
    alpha = 0.5,
    cex = 1.2,
    cat_cex = 1.2,
    margin = 0.08
) {
  if (is.null(pcset_obj$feature_sets)) {
    stop("pcset_obj must contain $feature_sets (use run_survival_pca_multi_pc()).")
  }

  label_mode <- match.arg(label_mode)
  pick <- as.integer(unique(pick))
  if (length(pick) < 2 || length(pick) > 5) {
    stop("pick must contain between 2 and 5 PC values.")
  }

  keys <- paste0("PC", pick)

  missing_keys <- setdiff(keys, names(pcset_obj$feature_sets))
  if (length(missing_keys) > 0) {
    stop(
      "These PC sets were not found in pcset_obj$feature_sets: ",
      paste(missing_keys, collapse = ", "),
      ". Run them first in pcs_to_run."
    )
  }

  feature_sets <- pcset_obj$feature_sets[keys]

  # ---- label helpers ----
  fmt_cum <- function(x) sprintf("%.1f%%", 100 * pmin(x, 1))

  n_features <- sapply(feature_sets, length)
  cumv <- if (!is.null(pcset_obj$cumvar_map)) pcset_obj$cumvar_map[keys] else rep(NA_real_, length(keys))
  used <- if (!is.null(pcset_obj$k_used_map)) pcset_obj$k_used_map[keys] else rep(NA_integer_, length(keys))

  # base title line per set
  lbl <- paste0("Top", pick, " PCs")

  # add requested label content
  add_line <- function(base, line) {
    if (is.na(line) || is.null(line) || length(line) == 0) return(base)
    paste0(base, "\n", line)
  }

  for (i in seq_along(lbl)) {
    if (label_mode == "both") {
      lbl[i] <- add_line(lbl[i], paste0("n=", n_features[i]))
      lbl[i] <- add_line(lbl[i], paste0("cum_var=", fmt_cum(cumv[i])))
    } else if (label_mode == "n") {
      lbl[i] <- add_line(lbl[i], paste0("n=", n_features[i]))
    } else if (label_mode == "cum_var") {
      lbl[i] <- add_line(lbl[i], paste0("cum_var=", fmt_cum(cumv[i])))
    } else if (label_mode == "none") {
      # keep only "TopXXX PCs"
    }

    if (isTRUE(show_used)) {
      lbl[i] <- add_line(lbl[i], paste0("used=", used[i]))
    }
  }

  names(feature_sets) <- lbl

  # ---- colors ----
  default_fill <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628")
  if (is.null(fill)) fill <- default_fill

  if (length(fill) < length(feature_sets)) {
    stop("Provide at least ", length(feature_sets), " colors in `fill` (or leave fill=NULL for defaults).")
  }

  venn.plot <- VennDiagram::venn.diagram(
    x = feature_sets,
    filename = NULL,
    col = "black",
    fill = fill[seq_along(feature_sets)],
    alpha = alpha,
    cex = cex,
    fontface = "bold",
    cat.cex = cat_cex,
    cat.fontface = "bold",
    margin = margin,
    main = main,
    main.cex = 1.8,
    main.fontface = "bold",
    fontfamily = "Arial",
    main.fontfamily = "Arial",
    cat.fontfamily = "Arial"
  )

  grid::grid.newpage()
  grid::grid.draw(venn.plot)

  invisible(venn.plot)
}

#' Plot trade-off between selected feature set size and variance explained across PC counts
#'
#' Draws a dual-axis curve showing how the number of selected features (left axis)
#' changes with the requested PC count, alongside the cumulative variance explained
#' (right axis).
#'
#' @param pcobj Output from \code{\link{run_survival_pca_multi_pc}}.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param ylab_left Left y-axis label (feature count).
#' @param ylab_right Right y-axis label (cumulative variance, %).
#' @param feature_color Line color for the selected feature count curve.
#' @param cumvar_color Line color for the cumulative variance curve.
#' @param feature_break_by Tick spacing for the feature-count axis (e.g., 100 or 200).
#' @param base_size Base font size for theme.
#'
#' @return A ggplot object.
#' @export
plot_feature_set_tradeoff <- function(
    pcobj,
    title = "Feature Set Size and Variance Explained Across \nPC Counts",
    xlab = "Number of PCs selected",
    ylab_left = "Number of selected features",
    ylab_right = "Cumulative variance explained (%)",
    feature_color = "#2C7BB6",
    cumvar_color = "#D7191C",
    feature_break_by = 200,
    base_size = 18
) {
  if (is.null(pcobj$feature_sets) || length(pcobj$feature_sets) == 0) {
    stop("pcobj$feature_sets is empty. Run run_survival_pca_multi_pc(..., pcs_to_run=...) first.")
  }
  if (is.null(pcobj$cumvar_map) || length(pcobj$cumvar_map) == 0) {
    stop("pcobj$cumvar_map is missing/empty. Ensure run_survival_pca_multi_pc() returns cumvar_map.")
  }
  
  keys <- intersect(names(pcobj$cumvar_map), names(pcobj$feature_sets))
  if (length(keys) < 2) stop("Need at least 2 PC sets with both feature sets and cumvar values.")
  
  pcs <- as.integer(sub("^PC", "", keys))
  
  df <- data.frame(
    pc_requested = pcs,
    n_features = as.integer(vapply(pcobj$feature_sets[keys], length, integer(1))),
    cumvar = as.numeric(pcobj$cumvar_map[keys]),
    stringsAsFactors = FALSE
  )
  
  df <- df[order(df$pc_requested), ]
  df$cumvar_pct <- 100 * pmin(df$cumvar, 1)
  
  # scale cumvar (%) onto feature axis (dual-axis display)
  scale_factor <- max(df$n_features, na.rm = TRUE) / max(df$cumvar_pct, na.rm = TRUE)
  df$cumvar_scaled <- df$cumvar_pct * scale_factor
  
  # feature axis breaks
  y_max <- max(df$n_features, na.rm = TRUE)
  y_lim <- max(feature_break_by, ceiling(y_max / feature_break_by) * feature_break_by)
  y_breaks <- seq(0, y_lim, by = feature_break_by)
  
  ggplot2::ggplot(df, ggplot2::aes(x = pc_requested)) +
    
    ggplot2::geom_line(
      ggplot2::aes(y = n_features, color = "Selected features"),
      linewidth = 1.5
    ) +
    
    ggplot2::geom_line(
      ggplot2::aes(y = cumvar_scaled, color = "Cumulative variance"),
      linewidth = 1.5
    ) +
    
    ggplot2::scale_color_manual(
      name = NULL,
      values = c(
        "Selected features"      = feature_color,
        "Cumulative variance" = cumvar_color
      )
    ) +
    
    ggplot2::scale_y_continuous(
      name = ylab_left,
      breaks = y_breaks,
      limits = c(0, y_lim),
      sec.axis = ggplot2::sec_axis(
        ~ . / scale_factor,
        name = ylab_right
      )
    ) +
    
    ggplot2::labs(title = title, x = xlab) +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(face = "bold")
    )
}

