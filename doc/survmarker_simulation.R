## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 5,
  fig.retina = 2
)

## ----libs, warning=FALSE------------------------------------------------------
library(SurvMarker)
library(survival)
library(ggplot2)

## ----simulate-----------------------------------------------------------------
set.seed(111)

n <- 200
p <- 6000
n_signal <- 25   # features driven by survival factor

# Make batch factor dominate variance 
z <- rnorm(n)

X <- matrix(rnorm(n * p), n, p)
X[, 1:n_signal] <- X[, 1:n_signal] + 1.5 * z

colnames(X) <- paste0("G", seq_len(p))
rownames(X) <- paste0("S", seq_len(n))

# Survival depends on z, not u
lp <- 0.9 * z
base_rate <- 0.08
T_event <- rexp(n, rate = base_rate * exp(scale(lp)))
C <- rexp(n, rate = 0.02)
time <- pmin(T_event, C)
status <- as.integer(T_event <= C)

# Simple covariates (optional)
AGE_CATEGORY <- sample(c("Young", "Old"), n, replace = TRUE)
SEX <- sample(c("F", "M"), n, replace = TRUE)

meta <- data.frame(
  AGE_CATEGORY = AGE_CATEGORY,
  SEX = SEX,
  stringsAsFactors = FALSE
)
rownames(meta) <- rownames(X)

# After z, AGE_CATEGORY, SEX, meta are created:

z_std <- as.numeric(scale(z))

# Adverse increases with survival factor z and age
p_adverse <- plogis(1.2*z_std + 0.8*(meta$AGE_CATEGORY == "Old") - 0.3)

# Favorable decreases with z and slightly with age
p_fav <- plogis(-1.0*z_std - 0.4*(meta$AGE_CATEGORY == "Old") + 0.2)

# Intermediate is what's left (with a floor to avoid negatives)
p_int <- 1 - p_adverse - p_fav
p_int[p_int < 0.05] <- 0.05

# Renormalize
s <- p_adverse + p_fav + p_int
p_adverse <- p_adverse/s; p_fav <- p_fav/s; p_int <- p_int/s

set.seed(111)
meta$ELN_2022_risk <- vapply(seq_len(nrow(meta)), function(i) {
  sample(c("Favorable","Intermediate","Adverse"), 1,
         prob = c(p_fav[i], p_int[i], p_adverse[i]))
}, character(1))

meta$ELN_2022_risk <- factor(meta$ELN_2022_risk,
                             levels = c("Favorable","Intermediate","Adverse"))


## ----run----------------------------------------------------------------------
# Ensure ELN is an ordered factor (recommended)
meta$ELN_2022_risk <- factor(
  meta$ELN_2022_risk,
  levels = c("Favorable", "Intermediate", "Adverse")
)

res <- pca_based_weighted_score(
  X = X,
  time = time,
  status = status,
  covar = meta[, c("ELN_2022_risk", "AGE_CATEGORY", "SEX"), drop = FALSE],
  n_pcs = 50,
  max_pcs = 50,
  pc_fdr_cutoff = 0.2,
  feature_fdr_cutoff = 0.05,
  null_B = 500,
  seed = 111,
  scale_pca = TRUE,
  store_null = TRUE,
  verbose = FALSE
)

head(res$feature_table, 10)


## ----pc12---------------------------------------------------------------------
# Define ELN-2022 color palette
eln_cols <- c(
  Favorable    = "#5E9F47",  # green
  Intermediate = "#bdbdbd",  # grey
  Adverse      = "#E35408"   # orange/red
)

# PC1–PC2 scatter colored by ELN-2022 risk
p_all_eln <- plot_pc12(
  res,
  meta,
  color_by = "ELN_2022_risk",
  color_values = eln_cols,
  theme_fn = ggplot2::theme_classic,
  base_size = 14,
  legend_position = "bottom",
  legend_text_size = 12
)

p_all_eln


## ----top2pcs------------------------------------------------------------------
p_sig_eln <- plot_top2_survival_pcs(
  res,
  meta = meta,
  color_by = "ELN_2022_risk",
  color_values = eln_cols,
  title = "Top 2 Survival-Significant PCs (ELN-2022 Risk)",
  theme_fn = ggplot2::theme_classic,
  base_size = 14,
  legend_position = "bottom",
  legend_text_size = 12
)

p_sig_eln

## ----scree--------------------------------------------------------------------
plot_scree(res, n_pcs = 100)

## ----cumvar-------------------------------------------------------------------
plot_cumvar(res, n_pcs = 100)

## ----null-vs-obs--------------------------------------------------------------
true_feature <- res$feature_table$feature[1]
null_feature <- res$feature_table$feature[res$feature_table$p_emp >0.1][5]

plot_null_vs_observed(
  res, feature = true_feature, bins = 40,
  theme_fn = ggplot2::theme_classic,
  base_size = 14
)

plot_null_vs_observed(
  res, feature = null_feature, bins = 40,
  theme_fn = ggplot2::theme_classic,
  base_size = 14
)

## ----multi-pc-----------------------------------------------------------------
pcobj <- run_survival_pca_multi_pc(
  X = X,
  time = time,
  status = status,
  covar = meta[, c("AGE_CATEGORY", "SEX"), drop = FALSE],
  pcs_to_run = c(10, 50, 20, 30),
  null_B = 300,
  store_null = TRUE,
  verbose = FALSE
)

venn_plot <- plot_venn(pcobj, pick = c(10, 50, 20, 30), cex = 1.0, cat_cex = 1.0,
          main = "Overlap of Features Selected Across Multiple Cutoffs")
venn_plot
tradeoff_plot <- plot_feature_set_tradeoff(pcobj)
tradeoff_plot

