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

### Option 2 — Build and install from source (offline / reproducible)

This option is useful for offline use or reproducible deployments.

#### Option 2A: Build using RStudio (recommended)

1. Download the latest Source code **(zip)** from the GitHub **Releases** page and unzip it.

2. Open the project:
   - Open `SurvMarker.Rproj` in **RStudio**.

3. Build and install:
   - In RStudio, click **Build → Install and Restart**

4. Load the package:
```r
library(SurvMarker)
```

#### Option 2B: Build manually from source

1. Download the latest Source code (zip) from the GitHub Releases page and unzip it.

2. Set your working directory to the parent directory containing the SurvMarker/ folder.

3. Build and install:

```r
install.packages("SurvMarker", repos = NULL, type = "source")
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
