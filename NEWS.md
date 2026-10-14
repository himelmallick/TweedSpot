# TweedSpot NEWS

## TweedSpot 0.99.2

### New features

- Cell-type-specific testing with `tweedspot(group = "cell_type")`: each
  feature is tested for spatial variation within each group, so pure
  composition effects are not called. Group p-values are adjusted stage-wise
  by default (`group_adjust`); results are available through
  `get_tweedspot_results(group = )` and `top_spatial_genes(group = )`.
- `tweedie_p = "auto"` chooses the variance power per feature by BIC or AIC
  (`p_criterion`) among Gaussian, Tweedie, Gamma, and inverse Gaussian models.
- `size_factor_method = "spanorm"` adds optional SpaNorm normalization, as size
  factors or as normalized counts (`spanorm_output`).
- New `seed` argument to `tweedspot()` for reproducible results.
- New plotting functions `plot_feature()` (observed, fitted, and spatial
  effect) and `plot_top_features()` (top `n_top` features), which can save
  publication-ready files to `output_dir`.

### Changes

- `filter_genes()` is renamed `filter_features()` with modality-agnostic
  arguments (`min_total`, `min_pct_locations`, `min_locations`, `min_mean`,
  `min_var`, `drop_empty_locations`). It adds named exclusions, a model-aware
  minimum, location QC, and `action = "flag"`; `exclude_mito` and
  `exclude_ribo` remain for transcriptomics.
- `plot_spatial_gene()` and `plot_spatial_fit()` are replaced by
  `plot_feature()` and `plot_top_features()`.
- `get_tweedspot_results()` and `top_spatial_genes()` return a plain table
  with columns `feature`, `feature_name`, `stat`, `pval`, `padj`, `edf`,
  `dev_expl`, `p`, and `family`; `rowData()` keeps the `tweedspot_` prefix.
- The two-part model always combines its p-values with the Cauchy combination
  test; the `combine` argument is removed.
- All families use a log link, and spatial coordinates are scaled by one
  common factor so the tissue's shape is kept.
- Sparse and `DelayedArray` assays are supported without dense copies.

### Bug fixes

- `tweedspot()` reused one Tweedie family object across features, so each fit
  started from the previous feature's estimate and could get stuck at the
  lower bound, inflating false positives. Results from earlier versions with
  `family = "tw"` should be recomputed.
- Stage-wise group adjustment now controls the false discovery rate at its
  target level.
- Results from an earlier `tweedspot()` run are cleared before each new run,
  so ordinary and group-specific results cannot be mixed.
- Mitochondrial and ribosomal gene matching no longer removes genes such as
  `MTOR` or the `RPS6K` kinases.
- Clearer handling of unusual inputs (non-syntactic covariate names, missing
  values, constant or extra coordinates, duplicated row names).

### Dependencies

- `MatrixGenerics` and `withr` added to Imports; `SpaNorm`, `patchwork`, `svglite`,
  `Matrix`, and `DelayedArray` added to Suggests.

## TweedSpot 0.99.1

- Added SCRAN-based size factor estimation through `size_factor_method = "scran"`,
  while preserving total-count library-size scaling as
  `size_factor_method = "library_size"`.
- Added `tweedie_p` support for `mgcv`-compatible variance power settings:
  estimated Tweedie, Gaussian, Poisson, fixed Tweedie powers in `(1, 2]`, and
  inverse Gaussian for strictly positive responses.
- Removed Stouffer p-value combination from the two-part model interface.

## TweedSpot 0.99.0

- Initial Bioconductor-oriented development release.
- Added agnostic spatially variable gene detection with Tweedie GAMs.
- Added gene filtering, result extraction, ranking helpers, and spatial plotters.
