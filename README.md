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

#### Option A: Build using RStudio (recommended)

1. Download the latest Source code **(SurvMarker-main.zip)** from the GitHub page and unzip it.

2. Open the project:
   - Open `SurvMarker.Rproj` in **RStudio**.

3. Build and install:
   - In RStudio, click **Build → Install and Restart**

4. Load the package:
```r
library(SurvMarker)
```

#### Option B: Build manually from source

1. Download the latest Source code (SurvMarker-main.zip) from the GitHub page and unzip it.

2. Go inside the SurvMarker-main folder

3. Build and install:

```r
list.files()#optional for checking files
install.packages(".", repos = NULL, type = "source")
library(SurvMarker)
```
or
```r
install.packages("devtools")  # if needed
devtools::install(".")
```

4. To build the vignettes, please first check whether Pandoc or Quarto is installed on your system. If neither is installed, please install one of them first. Once Pandoc or Quarto is available, you can build the vignettes using the following commands.
```r
install.packages("devtools")  # if needed
devtools::install(".", build_vignettes = TRUE)
browseVignettes("SurvMarker")
vignette(package="SurvMarker")
```
or
```r
install.packages(".", repos=NULL, type="source", vignette=T)
browseVignettes("SurvMarker")
vignette(package="SurvMarker")
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

## Contact
tgu at versiti.org

## Citation

If you use this code, please cite our manuscript:

> *SurvMarker: An R Package for Identifying Survival-Associated Molecular Features Using PCA-Based Weighted Scores*  
> Authors: [Dona Hasini Gammune, Tongjun Gu]  
> [BioRxiv], 2025.  

([Preprint/DOI link will be added here upon publication.](https://www.biorxiv.org/content/10.64898/2025.12.31.697184v1))
