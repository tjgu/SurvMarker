# SurvMarker

SurvMarker implements a PCA-based, weighted feature-scoring framework for survival analysis that identifies prognostically relevant molecular features by aggregating information across survival-associated principal components and assessing feature significance using an empirical null distribution.

---

## Installation

SurvMarker can be installed in **two ways**, depending on user preference.

---

### Option 1 — Install directly from GitHub (recommended)

This installs the latest development version directly from this repository.

```r
# Install remotes if needed
install.packages("remotes")

# Install SurvMarker
remotes::install_github("tjgu/SurvMarker")

# Load package
library(SurvMarker)
```
If you wish to build vignettes:

```r
remotes::install_github("tjgu/SurvMarker", build_vignettes = TRUE)
```

### Option 2 — Install from a source tarball (.tar.gz)

This option is useful for offline use or reproducible deployments.

1. Download the latest SurvMarker_*.tar.gz file from the Releases page.

2. Install locally in R:

```r
install.packages("SurvMarker_0.1.0.tar.gz", repos = NULL, type = "source")
library(SurvMarker)
```
## Quick Start

```r
library(SurvMarker)

res <- pca_based_weighted_score(
  X      = expr_matrix,
  time   = surv_time,
  status = surv_status
)

head(res$feature_table)
```

## Dependencies

SurvMarker depends on the following R packages:

- ggplot2
- survival
- VennDiagram

These will be installed automatically when installing SurvMarker from GitHub.

## Contacts
tgu at versiti.org

## Citations
SurvMarker: An R Package for Identifying Survival-Associated Molecular Features Using PCA-Based Weighted Scores. Dona Hasini Gammune, Tongjun Gu.
