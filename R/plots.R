#' Plot one feature: observed signal, TweedSpot fit, and spatial effect
#'
#' `plot_feature()` refits the TweedSpot model for one feature and draws three
#' aligned panels: the observed normalized signal, the fitted mean from the
#' TweedSpot model on the same colour scale, and the fitted spatial effect as a
#' log2 fold change relative to the feature's average. The model is refitted
#' with the settings stored by [tweedspot()] in
#' `S4Vectors::metadata(input)$tweedspot`, so the figure shows exactly the model
#' that was tested; any of those settings can be overridden through `...`.
#'
#' @param input A [SpatialExperiment::SpatialExperiment], typically the output
#'   of [tweedspot()].
#' @param feature A single row index, row name, or value of
#'   `rowData(input)$gene_name` / `rowData(input)$feature_name`.
#' @param group For results of `tweedspot(group = ...)`, the group (for example
#'   the cell type) to plot. The model is then refitted on that group's
#'   locations only, as it was tested, and only those locations are drawn.
#'   Required for group results and not allowed otherwise. Default `NULL`.
#' @param transform Display transform for the observed and fitted panels:
#'   `"log1p"` (default) or `"none"`. The spatial-effect panel is always on the
#'   log2 scale.
#' @param flip_y Logical; reverse the y axis so that images with the origin at
#'   the top-left are drawn upright. Default `TRUE`.
#' @param point_size Size of the location markers. Default `NULL` chooses a
#'   size from the spacing between locations and the figure size so that spots
#'   tile without overlapping.
#' @param base_size Base font size in points. Default `8`, suited to a
#'   full-width journal figure.
#' @param base_family Font family. Default `""` (the device's sans-serif).
#' @param output_dir Optional directory in which to save the figure. It is
#'   created if it does not exist. Default `NULL` (do not save).
#' @param formats File formats to write when `output_dir` is given; any of
#'   `"pdf"`, `"png"`, `"svg"`, `"tiff"`. Default `"pdf"`.
#' @param width,height Figure size in inches, used for saving and for choosing
#'   `point_size`. Defaults `7.2` by `2.9`, a two-column journal width.
#' @param dpi Resolution for raster formats. Default `300`.
#' @param ... Model settings that override those stored by [tweedspot()]:
#'   `assay_name`, `covariates`, `family`, `tweedie_p`, `p_criterion`,
#'   `fit_method`, `use_bam`, `bam_discrete`, `smooth_k`,
#'   `size_factor_method`, `spanorm_output`, `spanorm_args`, `seed`. Unless the
#'   normalization is overridden, the refit uses the size factors the run used.
#'
#' @details For a two-part run (`two_part = TRUE`) the panels show the count
#'   (abundance) part of the model.
#'
#' @return A `patchwork` object. When `output_dir` is supplied, the paths of the
#'   written files are attached as `attr(, "files")`.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("patchwork", quietly = TRUE)) {
#'   set.seed(1)
#'   grid <- expand.grid(x = 1:12, y = 1:12)
#'   counts <- rbind(
#'     band = stats::rpois(144, ifelse(grid$y > 6, 12, 2)),
#'     flat = stats::rpois(144, 5)
#'   )
#'   colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#'   spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts),
#'     spatialCoords = as.matrix(grid)
#'   )
#'   spe <- tweedspot(spe, smooth_k = 10, size_factor_method = "library_size",
#'                    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#'   plot_feature(spe, "band")
#' }
#' @export
plot_feature <- function(input,
                         feature,
                         group = NULL,
                         transform = c("log1p", "none"),
                         flip_y = TRUE,
                         point_size = NULL,
                         base_size = 8,
                         base_family = "",
                         output_dir = NULL,
                         formats = "pdf",
                         width = 7.2,
                         height = 2.9,
                         dpi = 300,
                         ...) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  tweedspot_plot_require(c("ggplot2", "patchwork"), "plot_feature")
  transform <- match.arg(transform)
  tweedspot_plot_check_args(flip_y, point_size, base_size, width, height, dpi)
  formats <- tweedspot_plot_formats(formats)

  info <- tweedspot_resolve_gene(input, feature)
  settings <- tweedspot_plot_settings(input, list(...))
  locations <- tweedspot_plot_locations(input, group)
  fitted <- tweedspot_plot_fits(input, info$index, settings, locations)[[1L]]

  coords <- tweedspot_plot_xy(input)
  if (!is.null(locations)) {
    coords <- coords[locations, , drop = FALSE]
  }
  if (is.null(point_size)) {
    point_size <- tweedspot_auto_point_size(coords, width / 3 - 0.15, height - 1.0)
  }
  disp <- if (identical(transform, "log1p")) log1p else identity
  observed <- disp(fitted$observed)
  fit <- disp(fitted$fitted)
  upper <- stats::quantile(c(observed, fit), 0.995, names = FALSE)
  lim_expr <- c(min(0, observed, fit), if (upper > 0) upper else max(c(observed, fit), 1))
  lim_eff <- max(abs(fitted$effect), 1e-8) * c(-1, 1)
  expr_title <- if (identical(transform, "log1p")) {
    "Normalized signal, log(1 + x)"
  } else {
    "Normalized signal"
  }

  panel <- function(value, title, scale) {
    df <- data.frame(x = coords[, 1], y = coords[, 2], value = value)
    tweedspot_spatial_panel(df, scale, point_size, flip_y, base_size, base_family) +
      ggplot2::labs(title = title)
  }
  p_obs <- panel(observed, "Observed",
                 tweedspot_scale_sequential(expr_title, lim_expr))
  p_fit <- panel(fit, "TweedSpot fit",
                 tweedspot_scale_sequential(expr_title, lim_expr))
  p_eff <- panel(fitted$effect, "Spatial effect",
                 tweedspot_scale_diverging("log2 fold change vs. feature mean", lim_eff))

  # Collecting the guides merges the two identical expression scales into one
  # legend and puts both legends in one row, so all three panels stay aligned.
  out <- patchwork::wrap_plots(p_obs, p_fit, p_eff, nrow = 1L) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(
      title = if (is.null(group)) info$label else paste0(info$label, " in ", group),
      theme = tweedspot_pub_theme(base_size, base_family) +
        ggplot2::theme(plot.title = ggplot2::element_text(
          size = base_size * 1.35, face = "bold", hjust = 0,
          margin = ggplot2::margin(b = base_size * 0.6)
        ))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.spacing.x = ggplot2::unit(base_size * 4, "pt")
    )

  if (!is.null(output_dir)) {
    raw_name <- if (is.null(rownames(input))) paste0("feature", info$index) else rownames(input)[info$index]
    stem <- tweedspot_file_stem(raw_name)
    if (!identical(stem, raw_name)) {
      # The name had to be cleaned for a file name; add the row number so two
      # features (for example "A/B" and "A_B") can never overwrite each other.
      stem <- paste0(stem, "_row", info$index)
    }
    if (!is.null(group)) {
      stem <- paste0(stem, "_", tweedspot_file_stem(group))
    }
    files <- tweedspot_save_plot(
      out, output_dir, paste0(stem, "_tweedspot"),
      formats, width, height, dpi
    )
    attr(out, "files") <- files
  }
  out
}

#' Plot the spatial effects of the top spatially variable features
#'
#' `plot_top_features()` draws the fitted spatial effect of the most
#' significant features from a [tweedspot()] run as small multiples on one
#' shared diverging scale, so that the strength of different patterns can be
#' compared directly. Each panel is labelled with the feature name and, by
#' default, its adjusted p-value and deviance explained.
#'
#' @inheritParams plot_feature
#' @param n_top Number of top features to show, ranked by adjusted p-value and
#'   then by test statistic. Default `6`.
#' @param features Optional vector of features (row indices, row names, or
#'   feature names) to show instead of the top `n_top`, in the given order.
#' @param group For results of `tweedspot(group = ...)`, the group to plot:
#'   features are ranked by that group's adjusted p-values and refitted on that
#'   group's locations. Required for group results and not allowed otherwise.
#'   When given and `title` is `NULL`, the group is used as the title.
#' @param ncol Number of panel columns. Default `NULL` uses 3 for up to nine
#'   panels and 4 otherwise.
#' @param label_stats Logical; add the adjusted p-value and deviance explained
#'   under each feature name. Default `TRUE`.
#' @param title Optional figure title. Default `NULL` (none).
#' @param width Figure width in inches, used for saving and for choosing
#'   `point_size`. Default `7.2`, a two-column journal width.
#' @param height Figure height in inches. Default `NULL` derives it from the
#'   number of rows and the tissue's aspect ratio.
#'
#' @return A `ggplot` object. When `output_dir` is supplied, the paths of the
#'   written files are attached as `attr(, "files")`.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   grid <- expand.grid(x = 1:12, y = 1:12)
#'   counts <- rbind(
#'     band = stats::rpois(144, ifelse(grid$y > 6, 12, 2)),
#'     spot = stats::rpois(144, 2 + 10 * exp(-((grid$x - 4)^2 + (grid$y - 4)^2) / 6)),
#'     flat = stats::rpois(144, 5)
#'   )
#'   colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
#'   spe <- SpatialExperiment::SpatialExperiment(
#'     assays = list(counts = counts),
#'     spatialCoords = as.matrix(grid)
#'   )
#'   spe <- tweedspot(spe, smooth_k = 10, size_factor_method = "library_size",
#'                    BPPARAM = BiocParallel::SerialParam(), verbose = FALSE)
#'   plot_top_features(spe, n_top = 2)
#' }
#' @export
plot_top_features <- function(input,
                              n_top = 6,
                              features = NULL,
                              group = NULL,
                              ncol = NULL,
                              label_stats = TRUE,
                              title = NULL,
                              flip_y = TRUE,
                              point_size = NULL,
                              base_size = 8,
                              base_family = "",
                              output_dir = NULL,
                              formats = "pdf",
                              width = 7.2,
                              height = NULL,
                              dpi = 300,
                              ...) {
  stopifnot(methods::is(input, "SpatialExperiment"))
  tweedspot_plot_require("ggplot2", "plot_top_features")
  tweedspot_check_flag(label_stats, "label_stats")
  tweedspot_plot_check_args(flip_y, point_size, base_size, width, height, dpi)
  formats <- tweedspot_plot_formats(formats)
  if (!is.null(title) && (!is.character(title) || length(title) != 1L)) {
    stop("`title` must be NULL or a single character string.")
  }

  locations <- tweedspot_plot_locations(input, group)
  top <- if (is.null(group)) {
    tweedspot_top_features_ordinary(input, n_top, features)
  } else {
    tweedspot_top_features_group(input, group, n_top, features)
  }
  idx <- top$index
  n <- length(idx)
  if (is.null(ncol)) {
    ncol <- if (n <= 9L) min(n, 3L) else 4L
  }
  tweedspot_check_scalar(ncol, "ncol", min_value = 1)
  ncol <- as.integer(min(ncol, n))
  nrow <- as.integer(ceiling(n / ncol))

  coords <- tweedspot_plot_xy(input)
  if (!is.null(locations)) {
    coords <- coords[locations, , drop = FALSE]
  }
  aspect <- diff(range(coords[, 2])) / max(diff(range(coords[, 1])), 1e-8)
  panel_w <- width / ncol - 0.12
  strip_h <- if (label_stats) base_size * 2.6 / 72 else base_size * 1.5 / 72
  if (is.null(height)) {
    height <- nrow * (panel_w * aspect + strip_h + 0.1) + 0.75 +
      if (is.null(title)) 0 else 0.3
  }
  if (is.null(point_size)) {
    point_size <- tweedspot_auto_point_size(
      coords, panel_w, (height - 0.75) / nrow - strip_h - 0.1
    )
  }

  settings <- tweedspot_plot_settings(input, list(...))
  fits <- tweedspot_plot_fits(input, idx, settings, locations)
  labels <- vapply(seq_along(idx), function(j) {
    lab <- tweedspot_resolve_gene(input, idx[j])$label
    if (label_stats) {
      lab <- paste0(
        lab, "\n",
        "adj. p ", tweedspot_format_p(top$padj[j]),
        ", ", sprintf("%.0f%%", 100 * top$dev_expl[j]), " deviance"
      )
    }
    lab
  }, character(1))
  df <- do.call(rbind, lapply(seq_along(idx), function(j) {
    data.frame(x = coords[, 1], y = coords[, 2], value = fits[[j]]$effect,
               panel = labels[j])
  }))
  df$panel <- factor(df$panel, levels = labels)
  cap <- stats::quantile(abs(df$value), 0.99, names = FALSE)
  cap <- max(cap, 1e-8)

  out <- tweedspot_spatial_panel(
    df,
    tweedspot_scale_diverging("log2 fold change vs. feature mean", c(-cap, cap)),
    point_size, flip_y, base_size, base_family
  ) +
    ggplot2::facet_wrap(~panel, ncol = ncol) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(
        size = base_size, hjust = 0, lineheight = 1.1,
        margin = ggplot2::margin(b = base_size * 0.3)
      ),
      panel.spacing.x = ggplot2::unit(base_size * 0.9, "pt"),
      panel.spacing.y = ggplot2::unit(base_size * 0.9, "pt")
    )
  if (is.null(title) && !is.null(group)) {
    title <- as.character(group)
  }
  if (!is.null(title)) {
    out <- out + ggplot2::labs(title = title) +
      ggplot2::theme(plot.title = ggplot2::element_text(
        size = base_size * 1.35, face = "bold", hjust = 0,
        margin = ggplot2::margin(b = base_size * 0.6)
      ))
  }

  if (!is.null(output_dir)) {
    stem <- "tweedspot_top_features"
    if (!is.null(group)) {
      stem <- paste0(stem, "_", tweedspot_file_stem(group))
    }
    files <- tweedspot_save_plot(out, output_dir, stem, formats, width, height, dpi)
    attr(out, "files") <- files
  }
  out
}

# ---- internals ---------------------------------------------------------------

# Colours. Expression magnitude uses the perceptually uniform "Mako" map
# (light teal to near-black). The spatial effect uses a diverging map built from
# Mako (below the feature's average) and "Rocket" (above it) around a light
# neutral midpoint. Values come from hcl.colors(9, "Mako", rev = TRUE) and
# hcl.colors(9, "Rocket", rev = TRUE), written out so figures look identical
# across R versions; the near-black ends are dropped from the diverging map so
# both extremes keep a readable hue.
tweedspot_seq_colours <- function() {
  c("#E0F7E1", "#A0DDC2", "#48C2B4", "#00A5B3", "#007FA8", "#375788",
    "#423460", "#311A35", "#070707")
}

tweedspot_div_colours <- function() {
  c(
    "#375788", "#007FA8", "#00A5B3", "#48C2B4", "#A0DDC2",   # Mako: below average
    "#F6F6F6",                                              # neutral midpoint
    "#F1B889", "#ED8353", "#E6464B", "#C30E62", "#970064"    # Rocket: above average
  )
}

tweedspot_plot_require <- function(pkgs, fun) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(sprintf("`%s()` requires the %s package%s: install.packages(c(%s)).",
                 fun, paste0("`", missing, "`", collapse = " and "),
                 if (length(missing) > 1L) "s" else "",
                 paste0("\"", missing, "\"", collapse = ", ")))
  }
}

tweedspot_plot_check_args <- function(flip_y, point_size, base_size, width, height, dpi) {
  tweedspot_check_flag(flip_y, "flip_y")
  if (!is.null(point_size)) {
    tweedspot_check_scalar(point_size, "point_size", min_value = 0)
  }
  tweedspot_check_scalar(base_size, "base_size", min_value = 1)
  tweedspot_check_scalar(width, "width", min_value = 0.5)
  if (!is.null(height)) {
    tweedspot_check_scalar(height, "height", min_value = 0.5)
  }
  tweedspot_check_scalar(dpi, "dpi", min_value = 1)
}

tweedspot_plot_formats <- function(formats) {
  allowed <- c("pdf", "png", "svg", "tiff")
  if (!is.character(formats) || !length(formats) ||
      !all(tolower(formats) %in% allowed)) {
    stop("`formats` must be one or more of: ", paste(allowed, collapse = ", "), ".")
  }
  unique(tolower(formats))
}

tweedspot_plot_xy <- function(input) {
  coords <- SpatialExperiment::spatialCoords(input)
  if (ncol(coords) < 2L) {
    stop("`input` must contain at least two spatial coordinates per location.")
  }
  coords[, seq_len(2), drop = FALSE]
}

# Model settings for refitting: package defaults, overridden by the settings a
# tweedspot() run stored in the metadata, overridden by explicit arguments.
tweedspot_plot_settings <- function(input, overrides = list()) {
  settings <- list(
    assay_name = "counts",
    covariates = NULL,
    family = "tw",
    tweedie_p = NULL,
    p_criterion = "BIC",
    fit_method = "REML",
    use_bam = FALSE,
    bam_discrete = TRUE,
    smooth_k = NULL,
    size_factor_method = "scran",
    spanorm_output = "size_factor",
    spanorm_args = list(),
    seed = NULL
  )
  stored <- S4Vectors::metadata(input)$tweedspot
  if (is.list(stored)) {
    for (nm in intersect(names(stored), names(settings))) {
      settings[nm] <- list(stored[[nm]])
    }
  }
  if (length(overrides)) {
    bad <- setdiff(names(overrides), names(settings))
    if (is.null(names(overrides)) || any(!nzchar(names(overrides))) || length(bad)) {
      stop("Unknown model setting(s) in `...`: ",
           paste(if (length(bad)) bad else "<unnamed>", collapse = ", "),
           ". Allowed: ", paste(names(settings), collapse = ", "), ".")
    }
    for (nm in names(overrides)) {
      settings[nm] <- list(overrides[[nm]])
    }
  }
  # Reuse the size factors the run actually used, unless the normalization is
  # overridden here or the locations changed since the run.
  if (is.list(stored) && !is.null(stored$size_factors) &&
      !any(c("size_factor_method", "assay_name", "spanorm_output") %in% names(overrides)) &&
      identical(colnames(input), stored$locations)) {
    settings$size_factors <- stored$size_factors
  }
  settings$size_factor_method <- match.arg(settings$size_factor_method,
                                           tweedspot_size_factor_methods())
  settings$spanorm_output <- match.arg(settings$spanorm_output,
                                       tweedspot_spanorm_outputs())
  settings
}

# Refit features and return, per feature, the observed and fitted signal on the
# normalized scale (divided by the size factor) and the spatial smooth as a log2
# fold change relative to the feature's average.
tweedspot_plot_fits <- function(input, idx, settings, locations = NULL) {
  tweedspot_local_seed(settings$seed)
  fit_input <- tweedspot_prepare_input(
    input = input,
    assay_name = settings$assay_name,
    covariates = settings$covariates,
    size_factor_method = settings$size_factor_method,
    spanorm_args = settings$spanorm_args,
    spanorm_output = settings$spanorm_output,
    seed = settings$seed
  )
  if (!is.null(settings$size_factors) && !is.matrix(fit_input$libsz) &&
      length(settings$size_factors) == length(fit_input$libsz)) {
    fit_input$libsz <- settings$size_factors
  }
  # Resolve the family from the whole assay, exactly as tweedspot() did, so the
  # refit uses the same Tweedie bounds as the test.
  full_Y <- fit_input$Y
  if (!is.null(locations)) {
    # Group results: size factors and coordinates come from all locations, as in
    # tweedspot(group = ...); the model is fitted on the group's locations only.
    fit_input$Y <- fit_input$Y[, locations, drop = FALSE]
    fit_input$coords <- fit_input$coords[locations, , drop = FALSE]
    fit_input$libsz <- if (is.matrix(fit_input$libsz)) {
      fit_input$libsz[, locations, drop = FALSE]
    } else {
      fit_input$libsz[locations]
    }
    fit_input$covariates <- tweedspot_subset_covariates(fit_input$covariates, locations, "plot")
  }
  lapply(idx, function(i) {
    expr <- as.numeric(fit_input$Y[i, ])
    libsz <- fit_input$libsz
    if (is.matrix(libsz)) {
      libsz <- as.numeric(libsz[i, ])
    }
    family <- tweedspot_resolve_family(settings$family, settings$tweedie_p, Y = full_Y)
    fit <- tweedspot_fit_single_gene(
      expr = expr,
      coords = fit_input$coords,
      libsz = libsz,
      covariates = fit_input$covariates,
      family = family,
      fit_method = settings$fit_method,
      use_bam = settings$use_bam,
      bam_discrete = settings$bam_discrete,
      bam_nthreads = 1L,
      smooth_k = settings$smooth_k,
      p_criterion = settings$p_criterion
    )
    if (is.null(fit)) {
      stop(sprintf("The TweedSpot model failed to fit for feature `%s`.",
                   rownames(input)[i]))
    }
    terms <- stats::predict(fit, type = "terms")
    smooth_col <- grep("^s\\(x,\\s*y", colnames(terms))[1L]
    list(
      observed = expr / libsz,
      fitted = as.numeric(stats::predict(fit, type = "response")) / libsz,
      effect = as.numeric(terms[, smooth_col]) / log(2)
    )
  })
}

# Which locations a plot uses: all of them for an ordinary run, or one group's
# locations for a tweedspot(group = ...) run.
tweedspot_plot_locations <- function(input, group) {
  run_group <- S4Vectors::metadata(input)$tweedspot$group
  if (is.null(run_group)) {
    if (!is.null(group)) {
      stop("`group` applies only to results from `tweedspot(group = ...)`.")
    }
    return(NULL)
  }
  labels <- tweedspot_group_labels(input, run_group)
  if (is.null(group)) {
    stop("This object holds group-specific results from `tweedspot(group = \"", run_group,
         "\")`. Choose the group to plot with `group =`, one of: ",
         paste(levels(labels), collapse = ", "), ".")
  }
  if (length(group) != 1L || !as.character(group) %in% levels(labels)) {
    stop("`group` must be one of: ", paste(levels(labels), collapse = ", "), ".")
  }
  which(labels == as.character(group))
}

tweedspot_check_n_top <- function(n_top) {
  tweedspot_check_scalar(n_top, "n_top", min_value = 1)
  if (n_top != round(n_top)) {
    stop("`n_top` must be a whole number.")
  }
}

tweedspot_take_top <- function(ranked, n_top) {
  if (n_top > length(ranked)) {
    message("Only ", length(ranked), " feature(s) have results; showing all of them.")
  }
  ranked[seq_len(min(as.integer(n_top), length(ranked)))]
}

tweedspot_top_features_ordinary <- function(input, n_top, features) {
  rd <- SummarizedExperiment::rowData(input)
  if (!"tweedspot_padj" %in% colnames(rd)) {
    stop("No TweedSpot results in `rowData(input)`; run `tweedspot()` first.")
  }
  tweedspot_check_locations(input, S4Vectors::metadata(input)$tweedspot$locations)
  if (!is.null(features)) {
    if (!length(features)) {
      stop("`features` must contain at least one feature.")
    }
    idx <- vapply(features, function(f) tweedspot_resolve_gene(input, f)$index, integer(1))
  } else {
    tweedspot_check_n_top(n_top)
    ok <- which(!is.na(rd$tweedspot_padj))
    if (!length(ok)) {
      stop("No feature has a TweedSpot result to plot.")
    }
    stat <- rd$tweedspot_stat[ok]
    stat[is.na(stat)] <- -Inf
    idx <- tweedspot_take_top(ok[order(rd$tweedspot_padj[ok], rd$tweedspot_pval[ok], -stat)], n_top)
  }
  list(index = unname(idx), padj = rd$tweedspot_padj[idx], dev_expl = rd$tweedspot_dev_expl[idx])
}

tweedspot_top_features_group <- function(input, group, n_top, features) {
  tab <- get_tweedspot_results(input, group = as.character(group), sort_by = "padj")
  if (!is.null(features)) {
    if (!length(features)) {
      stop("`features` must contain at least one feature.")
    }
    idx <- vapply(features, function(f) tweedspot_resolve_gene(input, f)$index, integer(1))
    row <- match(idx, tab$feature_index)
    return(list(index = unname(idx), padj = tab$padj[row], dev_expl = tab$dev_expl[row]))
  }
  tweedspot_check_n_top(n_top)
  tab <- tab[tab$tested & !is.na(tab$padj), , drop = FALSE]
  if (!nrow(tab)) {
    stop("No feature was tested in group '", group, "'.")
  }
  rows <- tweedspot_take_top(seq_len(nrow(tab)), n_top)
  list(index = tab$feature_index[rows], padj = tab$padj[rows], dev_expl = tab$dev_expl[rows])
}

tweedspot_format_p <- function(p) {
  if (is.na(p)) {
    return("NA")
  }
  if (p < 1e-300) {
    return("< 1e-300")
  }
  if (p >= 0.001) {
    return(paste("=", formatC(p, format = "f", digits = 3)))
  }
  paste("=", formatC(p, format = "e", digits = 1))
}

# Marker size (ggplot2 size units) at which neighbouring locations just touch,
# from the median nearest-neighbour distance and the space available per panel.
tweedspot_auto_point_size <- function(coords, panel_w_in, panel_h_in) {
  n <- nrow(coords)
  if (n < 2L) {
    return(2)
  }
  # Deterministic subsample so the user's random-number stream is untouched.
  idx <- if (n > 1500L) unique(round(seq(1, n, length.out = 1500L))) else seq_len(n)
  sub <- coords[idx, , drop = FALSE]
  nn <- vapply(seq_len(nrow(sub)), function(j) {
    d <- sqrt((coords[, 1] - sub[j, 1])^2 + (coords[, 2] - sub[j, 2])^2)
    min(d[d > 0], Inf)
  }, numeric(1))
  spacing <- stats::median(nn[is.finite(nn)])
  xr <- max(diff(range(coords[, 1])), 1e-8)
  yr <- max(diff(range(coords[, 2])), 1e-8)
  mm_per_unit <- min(max(panel_w_in, 0.3) / xr, max(panel_h_in, 0.3) / yr) * 25.4
  if (!is.finite(spacing) || spacing <= 0) {
    return(1)
  }
  # A filled circle (shape 16) of size s has a diameter of about 0.75 * s mm.
  size <- spacing * mm_per_unit / 0.75 * 0.98
  min(max(size, 0.05), 6)
}

tweedspot_pub_theme <- function(base_size = 8, base_family = "") {
  ink <- "#1a1a19"
  ink2 <- "#52514e"
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = ink),
      plot.title = ggplot2::element_text(
        size = base_size * 1.1, face = "bold", hjust = 0,
        margin = ggplot2::margin(b = base_size * 0.4)
      ),
      plot.title.position = "panel",
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = base_size * 0.95, colour = ink2),
      legend.text = ggplot2::element_text(size = base_size * 0.9, colour = ink2),
      legend.key.height = ggplot2::unit(base_size * 0.55, "pt"),
      # ggplot2 draws a colour bar five keys long: about 2 inches at base_size 8.
      legend.key.width = ggplot2::unit(base_size * 3.6, "pt"),
      legend.margin = ggplot2::margin(t = base_size * 0.3),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(base_size * 0.5, base_size * 0.5,
                                    base_size * 0.5, base_size * 0.5)
    )
}

# Legend title above a horizontal colour bar, in a way that works across
# ggplot2 versions (the theme element was added in ggplot2 3.5.0).
tweedspot_colourbar <- function() {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    ggplot2::guide_colourbar(theme = ggplot2::theme(legend.title.position = "top"))
  } else {
    ggplot2::guide_colourbar(title.position = "top", title.hjust = 0)
  }
}

tweedspot_scale_sequential <- function(name, limits) {
  ggplot2::scale_colour_gradientn(
    colours = tweedspot_seq_colours(), limits = limits, name = name,
    oob = scales_squish, guide = tweedspot_colourbar()
  )
}

tweedspot_scale_diverging <- function(name, limits) {
  ggplot2::scale_colour_gradientn(
    colours = tweedspot_div_colours(), limits = limits, name = name,
    oob = scales_squish, guide = tweedspot_colourbar()
  )
}

# Equivalent of scales::squish() without adding a dependency.
scales_squish <- function(x, range = c(0, 1), only.finite = TRUE) {
  finite <- if (only.finite) is.finite(x) else TRUE
  x[finite & x < range[1]] <- range[1]
  x[finite & x > range[2]] <- range[2]
  x
}

tweedspot_spatial_panel <- function(df, scale, point_size, flip_y, base_size, base_family) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = value)) +
    ggplot2::geom_point(size = point_size, shape = 16, stroke = 0) +
    ggplot2::coord_fixed() +
    scale +
    tweedspot_pub_theme(base_size, base_family)
  if (flip_y) {
    p <- p + ggplot2::scale_y_reverse()
  }
  p
}

tweedspot_file_stem <- function(x) {
  if (!length(x) || is.na(x[1L])) {
    return("feature")
  }
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  if (!nzchar(x)) "feature" else x
}

tweedspot_save_plot <- function(plot, output_dir, stem, formats, width, height, dpi) {
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) ||
      !nzchar(output_dir)) {
    stop("`output_dir` must be a single directory path.")
  }
  if (!dir.exists(output_dir)) {
    ok <- dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!ok) {
      stop("Could not create `output_dir`: ", output_dir)
    }
  }
  if ("svg" %in% formats && !requireNamespace("svglite", quietly = TRUE)) {
    stop("Saving SVG files requires the `svglite` package: install.packages(\"svglite\").",
         call. = FALSE)
  }
  files <- file.path(output_dir, paste0(stem, ".", formats))
  existing <- file.exists(files)
  if (any(existing)) {
    message("Overwriting existing file(s): ", paste(basename(files[existing]), collapse = ", "))
  }
  for (j in seq_along(formats)) {
    device <- switch(
      formats[j],
      pdf = if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf,
      formats[j]
    )
    ggplot2::ggsave(files[j], plot, width = width, height = height, units = "in",
                    dpi = dpi, device = device, bg = "white")
  }
  normalizePath(files)
}
