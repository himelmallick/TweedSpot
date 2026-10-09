# TweedSpot NEWS

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
