# =============================================================================
# gcf-package.R -- package-level documentation.
# =============================================================================

#' gcf: Generalized Covariate Field
#'
#' The generalized covariate field (GCF) method expands spatial covariates
#' into spatial-pattern and neighbourhood-distribution features and selects a
#' stable subset of them for geospatial prediction. The package implements the
#' two stages of the method:
#'
#' 1. **GCF variable generation** ([gcf_field()], with step functions
#'    [gcf_psi()], [gcf_zx()], and [gcf_reduce()]): one or multiple spatial
#'    variables with projected coordinates in, their GCF variables out.
#' 2. **Variable selection** ([gcf_select()], with the spatial-block helper
#'    [gcf_blocks()]): random forest importance combined with spatial-block
#'    stability resampling and group voting.
#'
#' The GCF method is prediction-oriented feature construction: the selected
#' variables feed any downstream regression learner (for example a random
#' forest), which stays outside this package.
#'
#' Two example datasets are included: [sim_grid] (simulation, 900 grid cells)
#' and [bio_grid] (south-western Australia plant-richness case study, 6229
#' grid cells).
#'
#' @references Song, Y. (2026). Generalized covariate field (GCF):
#'   spatial-pattern and neighbourhood-distribution feature expansion improves
#'   geospatial prediction. *International Journal of Geographical Information
#'   Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @keywords internal
"_PACKAGE"
