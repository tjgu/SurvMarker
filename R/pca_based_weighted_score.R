#' Survival-based loading-based feature scoring with empirical null
#'
#' @param X Numeric matrix/data.frame (n_samples x p_features). Rows=samples, cols=features.
#' @param time Numeric survival times (same order as rows of X).
#' @param status Event indicator (0/1). 1=event, 0=censored.
#' @param covar Optional data.frame of covariates (n_samples rows).
#'
#' @param n_pcs Integer or NULL. If provided, use first n_pcs PCs (default 50).
#' @param cumvar_threshold Numeric in (0,1] or NULL. If provided, choose smallest K such that
#'        cumulative variance explained >= cumvar_threshold. If both n_pcs and cumvar_threshold
#'        are provided, cumvar_threshold takes precedence.
#' @param max_pcs Hard cap on number of PCs allowed (safety cap; default 50).
#'
#' @param pc_fdr_cutoff FDR cutoff for selecting survival-associated PCs (default 0.05).
#' @param feature_fdr_cutoff FDR cutoff for selecting features from empirical null (default 0.05).
#'
#' @param null_B Number of null resamples (default 500).
#' @param seed Random seed.
#' @param scale_pca Scale features in PCA.
#' @param use_abs_loadings Use |loading| in score definition.
#' @param store_null If TRUE, return full null score matrix (features x B).
#' @param verbose Logical.
#'
#' @return List with feature_table (feature, loading_* for each significant PC, sj, p_emp, fdr),
#'         pc_table, pca, pc_scores, null_scores (optional), selected_features, settings.
#' @export
pca_based_weighted_score <- function(
    X,
    time,
    status,
    covar = NULL,
    
    n_pcs = 50,
    cumvar_threshold = NULL,
    max_pcs = 50,
    
    pc_fdr_cutoff = 0.05,
    feature_fdr_cutoff = 0.05,
    
    null_B = 500,
    seed = 1,
    scale_pca = TRUE,
    use_abs_loadings = TRUE,
    store_null = TRUE,
    verbose = TRUE
) {
  
  # ---------------- checks ----------------
  if (is.data.frame(X)) X <- as.matrix(X)
  stopifnot(is.matrix(X))
  stopifnot(is.numeric(time), length(time) == nrow(X))
  stopifnot(length(status) == nrow(X))
  status <- as.numeric(status)
  
  if (anyNA(time) || anyNA(status)) stop("time/status contain NA; please remove/impute.")
  if (any(time <= 0)) stop("time must be > 0 for Cox; filter non-positive times.")
  if (!all(status %in% c(0, 1))) stop("status must be 0/1.")
  
  if (!is.null(covar)) {
    if (!is.data.frame(covar)) covar <- as.data.frame(covar)
    stopifnot(nrow(covar) == nrow(X))
  }
  
  if (!is.null(cumvar_threshold)) {
    if (!is.numeric(cumvar_threshold) || length(cumvar_threshold) != 1 ||
        cumvar_threshold <= 0 || cumvar_threshold > 1) {
      stop("cumvar_threshold must be a single number in (0,1].")
    }
  }
  
  if (!is.null(n_pcs)) {
    if (!is.numeric(n_pcs) || length(n_pcs) != 1 || n_pcs < 2) {
      stop("n_pcs must be NULL or a single integer >= 2.")
    }
  }
  
  if (!is.numeric(max_pcs) || length(max_pcs) != 1 || max_pcs < 2) {
    stop("max_pcs must be a single integer >= 2.")
  }
  
  if (!is.numeric(pc_fdr_cutoff) || pc_fdr_cutoff <= 0 || pc_fdr_cutoff > 1)
    stop("pc_fdr_cutoff must be in (0,1].")
  if (!is.numeric(feature_fdr_cutoff) || feature_fdr_cutoff <= 0 || feature_fdr_cutoff > 1)
    stop("feature_fdr_cutoff must be in (0,1].")
  
  if (null_B < 50) warning("null_B is small; consider >=500.")
  if (is.null(colnames(X))) colnames(X) <- paste0("feature_", seq_len(ncol(X)))
  
  # ---------------- PCA ----------------
  pca <- stats::prcomp(X, center = TRUE, scale. = scale_pca)
  eig <- (pca$sdev)^2
  p_total <- length(eig)
  cum_var_all <- cumsum(eig / sum(eig))
  
  # ---- choose K ----
  if (!is.null(cumvar_threshold)) {
    K_raw <- which(cum_var_all >= cumvar_threshold)[1]
  } else {
    K_raw <- if (is.null(n_pcs)) 50 else as.integer(n_pcs)
  }
  K <- min(K_raw, max_pcs, p_total)
  
  if (verbose) {
    msg_rule <- if (!is.null(cumvar_threshold)) {
      paste0("cumvar_threshold=", cumvar_threshold)
    } else {
      paste0("n_pcs=", K_raw)
    }
    message("PCA done. Using K=", K,
            " (cum var=", round(cum_var_all[K], 3),
            "; rule: ", msg_rule, "; max_pcs=", max_pcs, ").")
  }
  
  pc_scores <- pca$x[, seq_len(K), drop = FALSE]
  colnames(pc_scores) <- paste0("PC", seq_len(K))
  
  # ============================================================
  # Step 1: JOINT Cox model on ALL PCs (and covariates, if given)
  # ============================================================
  pc_data <- data.frame(time = time, status = status, pc_scores)
  
  if (!is.null(covar)) pc_data <- cbind(pc_data, covar)
  
  fit_all <- survival::coxph(
    survival::Surv(time, status) ~ .,
    data = pc_data,
    ties = "efron",
    control = survival::coxph.control(iter.max = 100)
  )
  
  sm_all <- summary(fit_all)$coefficients
  
  # Extract per-PC stats from the joint model
  pc_names <- colnames(pc_scores)
  pc_beta <- pc_hr <- pc_p <- rep(NA_real_, K)
  names(pc_beta) <- names(pc_hr) <- names(pc_p) <- pc_names
  
  for (k in seq_len(K)) {
    nm <- pc_names[k]
    if (nm %in% rownames(sm_all)) {
      pc_beta[k] <- sm_all[nm, "coef"]
      pc_hr[k]   <- sm_all[nm, "exp(coef)"]
      pc_p[k]    <- sm_all[nm, "Pr(>|z|)"]
    }
  }
  
  pc_padj <- stats::p.adjust(pc_p, method = "fdr")
  
  pc_table <- data.frame(
    PC = pc_names,
    eigenvalue = eig[seq_len(K)],
    var_explained = eig[seq_len(K)] / sum(eig),
    cum_var = cum_var_all[seq_len(K)],
    beta = pc_beta,
    HR = pc_hr,
    pval = pc_p,
    padj = pc_padj,
    stringsAsFactors = FALSE
  )
  
  sig_idx <- which(pc_padj <= pc_fdr_cutoff)
  nonsig_idx <- setdiff(seq_len(K), sig_idx)
  
  if (verbose) {
    message("Survival PCs (PC-level FDR <= ", pc_fdr_cutoff, "): ", length(sig_idx))
    if (length(sig_idx) > 0) message("  ", paste0("PC", sig_idx, collapse = ", "))
  }
  
  # Loadings (features x K) for later use (both score and null)
  loadings_all <- pca$rotation[, seq_len(K), drop = FALSE]
  
  # ---------------- No significant PCs: return early ----------------
  if (length(sig_idx) == 0) {
    feature_table <- data.frame(
      feature = colnames(X),
      sj = 0,
      p_emp = NA_real_,
      fdr = NA_real_,
      stringsAsFactors = FALSE
    )
    
    return(list(
      feature_table = feature_table,
      pc_table = pc_table,
      pca = pca,
      pc_scores = pc_scores,
      significant_pcs = character(0),
      significant_pc_indices = integer(0),
      pc_weights = numeric(0),
      null_scores = NULL,
      selected_features = character(0),
      settings = list(
        n_pcs = n_pcs, cumvar_threshold = cumvar_threshold, max_pcs = max_pcs,
        pc_fdr_cutoff = pc_fdr_cutoff, feature_fdr_cutoff = feature_fdr_cutoff,
        null_B = null_B, scale_pca = scale_pca, use_abs_loadings = use_abs_loadings,
        store_null = store_null, seed = seed
      )
    ))
  }
  
  # ---------------- Step 2: sj + significant-PC loadings ----------------
  sig_pc_names <- paste0("PC", sig_idx)
  
  sig_loadings <- loadings_all[, sig_idx, drop = FALSE]
  colnames(sig_loadings) <- paste0("loading_", sig_pc_names)
  
  wk <- eig[sig_idx] / sum(eig[sig_idx])
  names(wk) <- sig_pc_names
  
  L_for_score <- sig_loadings
  if (use_abs_loadings) L_for_score <- abs(L_for_score)
  
  sj <- as.numeric(L_for_score %*% wk)
  names(sj) <- colnames(X)
  
  # ---------------- Step 3: empirical null ----------------
  m <- length(sig_idx)
  if (length(nonsig_idx) < m) {
    stop("Not enough non-significant PCs for null sampling. Increase K (n_pcs/max_pcs) or relax pc_fdr_cutoff.")
  }
  
  set.seed(seed)
  
  if (store_null) {
    sj_null <- matrix(
      NA_real_,
      nrow = ncol(X),
      ncol = null_B,
      dimnames = list(colnames(X), paste0("b", seq_len(null_B)))
    )
    
    for (b in seq_len(null_B)) {
      pick <- sample(nonsig_idx, size = m, replace = FALSE)
      w_null <- eig[pick] / sum(eig[pick])
      Lnull <- loadings_all[, pick, drop = FALSE]
      if (use_abs_loadings) Lnull <- abs(Lnull)
      sj_null[, b] <- as.numeric(Lnull %*% w_null)
    }
    
    p_emp <- rowMeans(sj_null >= sj)
  } else {
    exceed <- rep(0, ncol(X))
    
    for (b in seq_len(null_B)) {
      pick <- sample(nonsig_idx, size = m, replace = FALSE)
      w_null <- eig[pick] / sum(eig[pick])
      Lnull <- loadings_all[, pick, drop = FALSE]
      if (use_abs_loadings) Lnull <- abs(Lnull)
      sjb <- as.numeric(Lnull %*% w_null)
      exceed <- exceed + as.integer(sjb >= sj)
    }
    
    p_emp <- exceed / null_B
    sj_null <- NULL
  }
  
  fdr <- stats::p.adjust(p_emp, method = "BH")
  selected <- names(fdr)[fdr <= feature_fdr_cutoff]
  
  feature_table <- data.frame(feature = colnames(X), stringsAsFactors = FALSE)
  feature_table <- cbind(feature_table, as.data.frame(sig_loadings))
  feature_table$sj <- sj
  feature_table$p_emp <- p_emp
  feature_table$fdr <- fdr
  feature_table <- feature_table[order(feature_table$fdr, feature_table$p_emp), ]
  
  if (verbose) message("Selected features (feature-level FDR <= ", feature_fdr_cutoff, "): ", length(selected))
  
  list(
    feature_table = feature_table,
    pc_table = pc_table,
    pca = pca,
    pc_scores = pc_scores,
    significant_pcs = sig_pc_names,
    significant_pc_indices = sig_idx,
    pc_weights = wk,
    null_scores = sj_null,
    selected_features = selected,
    settings = list(
      n_pcs = n_pcs, cumvar_threshold = cumvar_threshold, max_pcs = max_pcs,
      pc_fdr_cutoff = pc_fdr_cutoff, feature_fdr_cutoff = feature_fdr_cutoff,
      null_B = null_B, scale_pca = scale_pca, use_abs_loadings = use_abs_loadings,
      store_null = store_null, seed = seed
    )
  )
}
