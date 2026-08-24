mock_spe <- function(counts = NULL) {
  if (is.null(counts)) {
    counts <- matrix(
      c(
        0, 0, 2, 3, 4, 0,
        0, 1, 0, 2, 1, 0,
        5, 4, 3, 2, 1, 0,
        1, 1, 1, 1, 1, 0
      ),
      nrow = 4,
      byrow = TRUE
    )
  }

  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("spot", seq_len(ncol(counts)))

  row_data <- S4Vectors::DataFrame(
    gene_id = rownames(counts),
    gene_name = c("MT-CO1", "RPLP0", "GENE3", "GENE4")[seq_len(nrow(counts))],
    feature_type = rep("Gene Expression", nrow(counts))
  )

  col_data <- S4Vectors::DataFrame(
    batch = rep(c("A", "B"), length.out = ncol(counts)),
    sample_id = paste0("sample", seq_len(ncol(counts)))
  )

  spatial_coords <- cbind(
    pxl_col_in_fullres = seq_len(ncol(counts)),
    pxl_row_in_fullres = seq_len(ncol(counts))
  )

  SpatialExperiment::SpatialExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = col_data,
    spatialCoords = spatial_coords
  )
}
