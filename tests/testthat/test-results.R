test_that("get_tweedspot_results sorts by adjusted p-value and statistic", {
  spe <- mock_spe()
  rd <- SummarizedExperiment::rowData(spe)
  rd$tweedspot_stat <- c(2, 5, 4, 1)
  rd$tweedspot_pval <- c(0.02, 0.03, 0.03, 0.2)
  rd$tweedspot_padj <- c(0.05, 0.05, 0.05, 0.3)
  rd$tweedspot_edf <- c(2, 2, 2, 2)
  rd$tweedspot_dev_expl <- c(0.1, 0.4, 0.3, 0.05)
  SummarizedExperiment::rowData(spe) <- rd

  res <- get_tweedspot_results(spe, sort_by = "padj")

  expect_identical(res$feature[1:3], c("gene2", "gene3", "gene1"))
  expect_identical(colnames(res),
                   c("feature", "feature_name", "stat", "pval", "padj", "edf", "dev_expl"))
  expect_identical(res$feature_name[1], "RPLP0")
  expect_identical(get_tweedspot_results(spe, sort_by = "dev_expl")$feature[1], "gene2")
  expect_identical(get_tweedspot_results(spe, sort_by = "none")$feature,
                   paste0("gene", 1:4))
})

test_that("top_spatial_genes uses ranked TweedSpot results", {
  spe <- mock_spe()
  rd <- SummarizedExperiment::rowData(spe)
  rd$tweedspot_stat <- c(2, 5, 4, 1)
  rd$tweedspot_pval <- c(0.02, 0.03, 0.03, 0.2)
  rd$tweedspot_padj <- c(0.05, 0.05, 0.05, 0.3)
  rd$tweedspot_edf <- c(2, 2, 2, 2)
  rd$tweedspot_dev_expl <- c(0.1, 0.4, 0.3, 0.05)
  SummarizedExperiment::rowData(spe) <- rd

  top_res <- top_spatial_genes(spe, n = 2, rank_by = "padj")

  expect_equal(nrow(top_res), 2)
  expect_identical(top_res$feature, c("gene2", "gene3"))
})
