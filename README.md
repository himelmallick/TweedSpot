# TweedSpot

**Spatially variable feature detection with Tweedie generalized additive models**

TweedSpot finds features (genes, proteins, peaks, metabolites) whose abundance
varies smoothly across a tissue. Each feature is fitted with a spatial
generalized additive model (GAM) in [`mgcv`](https://cran.r-project.org/package=mgcv):

```
E[y_i] = mu_i,   log(mu_i) = s(x_i, y_i) + covariates_i + log(size factor_i)
y_i ~ Tweedie(mu_i, phi, p)
```

The spatial smooth `s(x, y)` is tested against a flat surface, so a small
p-value means the feature is spatially structured after accounting for
sequencing depth and any covariates you supply. The Tweedie family covers
counts with many zeros (1 < p < 2) and continuous positive measurements, so
the same model works on raw counts, normalized intensities, and other
non-negative assays without transformation.

TweedSpot works on `SpatialExperiment` objects and writes its results back
into them, so it fits into standard Bioconductor workflows.

## Highlights

- **Models data on its own scale.** Raw counts are modeled directly with a
  size-factor offset; no log-transform or variance-stabilization step.
- **Adapts to each feature.** The Tweedie power `p` is estimated per feature,
  or chosen among Gaussian, Tweedie, Gamma, and inverse Gaussian by BIC/AIC
  with `tweedie_p = "auto"`.
- **Cell-type-specific testing.** With `group = "cell_type"`, each feature is
  tested for spatial variation *within* each cell type, so a marker that only
  reflects where that cell type sits is not called.
- **Modality-agnostic filtering.** `filter_features()` removes features the
  spatial model cannot support and low-quality locations, with optional
  mitochondrial and ribosomal gene removal for transcriptomics.
- **Publication-ready figures.** `plot_feature()` and `plot_top_features()`
  draw observed signal, fitted surface, and spatial effect, and save PDF, PNG,
  SVG, or TIFF files.
- **Scales to large data.** Sparse and `DelayedArray` assays, `mgcv::bam()`
  fitting, and parallel execution through `BiocParallel`.
- **Reproducible.** One `seed` argument controls all random steps and is
  reused by the plotting functions.

## Installation

From Bioconductor (once released):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TweedSpot")
```

Development version from GitHub:

```r
BiocManager::install("himelmallick/TweedSpot", dependencies = TRUE)
# or: remotes::install_github("himelmallick/TweedSpot", dependencies = TRUE)
```

TweedSpot requires R (>= 4.6.0). Optional packages, needed only for the
matching feature:

| Package | Needed for |
| --- | --- |
| `SpaNorm` | `size_factor_method = "spanorm"` |
| `ggplot2`, `patchwork` | `plot_feature()`, `plot_top_features()` |
| `svglite` | saving figures as SVG |
| `STexampleData` | example datasets used below and in the vignette |

## Quick start

```r
library(TweedSpot)
library(SpatialExperiment)
library(STexampleData)

# 10x Visium, human dorsolateral prefrontal cortex
spe <- Visium_humanDLPFC()
spe <- spe[, colData(spe)$in_tissue == 1]

# 1. Keep features the spatial model can support
spe <- filter_features(spe, smooth_k = 20, exclude_mito = TRUE, report = TRUE)

# 2. Fit one spatial GAM per feature
spe <- tweedspot(spe, smooth_k = 20, seed = 1)

# 3. Look at the strongest spatial features
top_spatial_genes(spe, n = 10)

# 4. Plot one of them
plot_feature(spe, "MBP")
```

`top_spatial_genes()` and `get_tweedspot_results()` return a plain data frame,
one row per feature, sorted by adjusted p-value:

```
  feature         feature_name   stat    pval        padj        edf   dev_expl
1 ENSG00000197971 MBP            2841.6  0           0           17.9  0.612
2 ENSG00000123560 PLP1           2310.4  0           0           17.6  0.548
3 ENSG00000131095 GFAP           1102.7  1.2e-220    3.1e-217    16.8  0.341
...
```

(Values are illustrative.)

## Workflow in detail

### 1. Input

TweedSpot expects a `SpatialExperiment` with

- an assay of non-negative values (`assay_name`, default `"counts"`), which
  may be a dense matrix, a sparse `dgCMatrix`, or a `DelayedArray`;
- two spatial coordinates in `spatialCoords()`. Coordinates are scaled by one
  common factor, so the shape of the tissue is kept;
- optionally, covariates in `colData()` and a precomputed
  `colData()$sizeFactor`.

### 2. Filtering with `filter_features()`

`filter_features()` removes features and locations that would only add noise
or fail to fit. Its thresholds are written in modality-neutral terms.

| Argument | Default | What it does |
| --- | --- | --- |
| `min_total` | `5` | Minimum total signal across all locations |
| `min_pct_locations` | `5` | Minimum percent of locations where the feature is detected |
| `min_locations` | `0` | Minimum number of locations where the feature is detected |
| `min_mean`, `min_var` | `0` | Minimum mean and variance |
| `model_aware` + `smooth_k` | `TRUE`, `NULL` | Require at least as many detected locations as the spatial basis has coefficients, so every kept feature can be fitted |
| `exclude`, `exclude_pattern` | `NULL` | Drop named features, for example negative-control probes (`"^NegPrb"`) |
| `exclude_mito`, `exclude_ribo` | `FALSE` | Transcriptomics only: drop mitochondrial / ribosomal genes |
| `min_total_per_location`, `min_features_per_location` | `0` | Location QC thresholds |
| `outlier_nmads` | `NULL` | Drop locations that are outliers by this many MADs |
| `max_excluded_fraction` | `NULL` | Drop locations dominated by excluded features (for example >20% mitochondrial) |
| `group`, `min_group_size` | `NULL`, `50` | Check model support within each cell type, for group-mode runs |
| `action` | `"drop"` | `"flag"` keeps everything and adds `filter_keep` / `filter_reason` columns |
| `report` | `FALSE` | Print a summary of what was removed and why |

```r
# Inspect before removing anything
flagged <- filter_features(spe, smooth_k = 20, exclude_mito = TRUE,
                           outlier_nmads = 3, action = "flag")
table(rowData(flagged)$filter_reason)
```

Use the same `smooth_k` here and in `tweedspot()`, so the support check matches
the model you will fit.

### 3. Normalization

Sequencing depth enters the model as an offset, `log(size factor)`, so counts
stay on their original scale.

| `size_factor_method` | Description |
| --- | --- |
| `"scran"` (default) | Pooling-based deconvolution size factors from `scran::computeSumFactors()` |
| `"library_size"` | Total signal per location, scaled to a median of one |
| `"spanorm"` | Spatially aware normalization with [SpaNorm](https://bioconductor.org/packages/SpaNorm), for spatial transcriptomics |

If `colData(spe)$sizeFactor` already exists, `"scran"` and `"library_size"` use
it as it is, so you can supply size factors from any method.

SpaNorm can enter in two ways: as gene- and location-specific size factors on
raw counts (`spanorm_output = "size_factor"`, default), or by modeling its
percentile-adjusted counts (`spanorm_output = "counts"`). The SpaNorm fit is
cached in `metadata(spe)$SpaNorm` and reused when the data and settings match.

> **Caution.** TweedSpot treats the SpaNorm estimate as known. When library
> size itself varies across the tissue, this can make the spatial test
> anti-conservative. Treat SpaNorm results as exploratory and check
> calibration on your data.

### 4. Model

```r
spe <- tweedspot(
  spe,
  assay_name = "counts",
  covariates = ~ batch,        # optional colData() variables
  tweedie_p = NULL,            # estimate p per feature (default)
  smooth_k = 20,               # spatial basis size
  use_bam = TRUE,              # faster fitting for large data
  size_factor_method = "scran",
  padj_method = "BY",
  seed = 1,
  BPPARAM = BiocParallel::MulticoreParam(4)
)
```

**Variance power.** `tweedie_p` sets the mean-variance relationship,
`Var(y) = phi * mu^p`:

| `tweedie_p` | Family | Typical use |
| --- | --- | --- |
| `NULL` (default) | `mgcv::tw()`, `p` estimated in (1, 2) | Counts with excess zeros |
| `"auto"` | Lowest BIC (or AIC, `p_criterion = "AIC"`) among Gaussian, Tweedie, Gamma, inverse Gaussian | Let each feature choose |
| `0` | Gaussian, log link | Continuous, roughly normal |
| `1` | Poisson | Counts without overdispersion |
| `1 < p < 2` | Fixed-power Tweedie | Counts with a chosen power |
| `2` | Gamma | Strictly positive continuous |
| `3` | Inverse Gaussian | Strictly positive, heavily skewed |

With `"auto"`, candidates are fitted only when they are valid for the feature
(Gamma and inverse Gaussian need all values positive), and the selected power
and family are reported per feature. Poisson is not an `"auto"` candidate,
because its fixed dispersion makes the test anti-conservative on overdispersed
data; use `tweedie_p = 1` to fit it on purpose. For integer data, the
estimated Tweedie power is bounded below at 1.1 to avoid the density spikes
near p = 1. All families use a log link.

**Two-part model.** `two_part = TRUE` fits a binomial model for detection and
a positive-part model for abundance where detected, then combines the two
spatial p-values with the Cauchy combination test. Use it when spatial
patterns may appear in *whether* a feature is detected as well as in *how
much*.

**Multiple testing.** p-values are adjusted across features with
`padj_method` (default Benjamini–Yekutieli, `"BY"`); any method accepted by
`p.adjust()` works.

### 5. Results

```r
res <- get_tweedspot_results(spe)                   # all features, by padj
res <- get_tweedspot_results(spe, sort_by = "dev_expl")
top <- top_spatial_genes(spe, n = 20, rank_by = "stat")
```

| Column | Meaning |
| --- | --- |
| `feature` | Row name of the feature |
| `feature_name` | Readable name from `rowData()` (for example a gene symbol), when available |
| `stat` | Test statistic for the spatial smooth (larger = stronger evidence) |
| `pval` | p-value for the spatial smooth |
| `padj` | Adjusted p-value (`padj_method`) |
| `edf` | Effective degrees of freedom of the smooth: about 1–2 for a simple gradient, larger for more complex patterns |
| `dev_expl` | Fraction of deviance explained by the full model |
| `p`, `family` | Variance power and model family used for the feature |

The same values are stored in `rowData(spe)` with a `tweedspot_` prefix
(`tweedspot_padj`, `tweedspot_edf`, and so on), so they travel with the
object. Features that could not be fitted have `NA` results.

### 6. Cell-type-specific testing

In single-cell-resolution data (Xenium, CosMx, MERFISH, Slide-seq, segmented
imaging), a cell-type marker looks spatially variable simply because that cell
type is clustered. Group mode separates this composition effect from spatial
variation within a cell type:

```r
spe <- filter_features(spe, group = "cell_type", smooth_k = 20)
spe <- tweedspot(spe, group = "cell_type", smooth_k = 20, seed = 1)

get_tweedspot_results(spe, group = "Astrocyte")     # one cell type
get_tweedspot_results(spe, group = TRUE)            # all feature x type pairs
top_spatial_genes(spe, n = 10, group = "Astrocyte")
plot_feature(spe, "Gfap", group = "Astrocyte")
```

How it works:

- Each cell type is fitted on its own locations with its own spatial smooth;
  size factors and coordinate scaling are shared across types.
- Cell types with fewer than `min_group_size` (default 50) locations are
  skipped, as are features detected in too few locations of a type.
- Covariates that are constant within a type are dropped for that type.
- Adjustment across feature × type pairs, `group_adjust`:
  - `"stagewise"` (default): first screen features (any type spatial?), then
    confirm which types, controlling the false discovery rate at `alpha`;
  - `"global"`: one adjustment over all pairs;
  - `"per_group"`: a separate adjustment within each type.

Group results include `group`, `n_locations`, `n_detected`, `tested`, and a
`reason` for untested pairs (`group_too_small`, `low_model_support`,
`fit_failed`). Per-feature summaries are added to `rowData(spe)` as
`tweedspot_group_padj` and `tweedspot_group_sig` (the types where the feature
is significant).

### 7. Figures

```r
# Observed signal | TweedSpot fit | spatial effect, for one feature
plot_feature(spe, "MBP", output_dir = "figures", formats = c("pdf", "png"))

# Spatial effects of the top features on a shared colour scale
plot_top_features(spe, n_top = 9, output_dir = "figures")
```

Both functions refit with the exact settings of the `tweedspot()` run, return
`ggplot`/`patchwork` objects you can modify further, and use colour-blind-safe
Mako and Rocket palettes. Set `output_dir` to save files in PDF, PNG, SVG,
or TIFF at `dpi` resolution; `width` and `height` are in inches.

### 8. Reproducibility

`seed` makes a run reproducible (scran pooling and SpaNorm subsampling use
random numbers). It is stored with the results and reused by the plotting
functions, and your own random-number state is restored after the call.
Results from an earlier run are cleared before each new `tweedspot()` call,
so results from different settings are never mixed.

### 9. Speed and memory

- `smooth_k` controls the size of the spatial basis. Smaller values are faster
  but can miss fine patterns; 20–50 works well for most tissues.
- `use_bam = TRUE` uses `mgcv::bam()` with discretization
  (`bam_discrete = TRUE`), which is much faster on tens of thousands of
  locations.
- `BPPARAM` fits features in parallel, for example
  `BiocParallel::MulticoreParam(8)`.
- Sparse and `DelayedArray` assays are read one feature at a time, without a
  dense copy of the whole matrix.

## Function reference

| Function | Purpose |
| --- | --- |
| `filter_features()` | Remove or flag low-information features and low-quality locations |
| `tweedspot()` | Fit spatial GAMs and test each feature for spatial variation |
| `get_tweedspot_results()` | Results table for all features, or per cell type |
| `top_spatial_genes()` | The top `n` features by `padj`, `pval`, `stat`, or `dev_expl` |
| `plot_feature()` | Observed, fitted, and spatial-effect panels for one feature |
| `plot_top_features()` | Spatial effects of the top `n_top` features |

## Where results are stored

| Location | Contents |
| --- | --- |
| `rowData(spe)$tweedspot_*` | Per-feature results (`stat`, `pval`, `padj`, `edf`, `dev_expl`, `p`, `family`) |
| `rowData(spe)$tweedspot_group_padj`, `tweedspot_group_sig` | Per-feature summaries from group mode |
| `metadata(spe)$tweedspot` | Settings of the last run, used by the plotting functions |
| `metadata(spe)$tweedspot_group` | Feature × group results table |
| `metadata(spe)$filter_features` | Settings and summary of filtering |
| `metadata(spe)$SpaNorm` | Cached SpaNorm fit |

## Learn more

The vignette walks through a full analysis, including a real cell-type
example on Slide-seqV2 mouse hippocampus data:

```r
browseVignettes("TweedSpot")
```

Bug reports and feature requests:
<https://github.com/himelmallick/TweedSpot/issues>

## License

MIT
