test_that("filter_genes filters genes and zero-library spots", {
  spe <- mock_spe()

  filtered <- filter_genes(
    input = spe,
    assay_name = "counts",
    filter_genes_ncounts = 8,
    filter_genes_pcspots = 40,
    filter_genes_nspots = 3,
    filter_genes_mean = 1,
    drop_zero_spots = TRUE
  )

  expect_s4_class(filtered, "SpatialExperiment")
  expect_true(nrow(filtered) < nrow(spe))
  expect_true(ncol(filtered) < ncol(spe))
  expect_true(all(colSums(SummarizedExperiment::assay(filtered, "counts")) > 0))
})

test_that("filter_genes reports retained and removed features", {
  spe <- mock_spe()

  expect_message(
    filter_genes(
      input = spe,
      assay_name = "counts",
      filter_genes_ncounts = 3,
      report = TRUE
    ),
    "filter_genes retained"
  )
})

test_that("filter_genes can exclude mitochondrial and ribosomal genes", {
  spe <- mock_spe()

  filtered <- filter_genes(
    input = spe,
    assay_name = "counts",
    exclude_mito = TRUE,
    exclude_ribo = TRUE,
    gene_name_col = "gene_name"
  )

  expect_false(any(rownames(filtered) %in% c("gene1", "gene2")))
})
