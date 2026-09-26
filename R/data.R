# =============================================================================
# data.R -- documentation of the two bundled datasets.
# =============================================================================

#' Simulated spatial grid (GCF paper simulation)
#'
#' The 30 x 30 simulation grid of the GCF paper: a response surface and three
#' spatially structured covariates on a 900-cell regular grid with unit
#' spacing. The paper builds the GCF variables with buffers 2, 4, 6 (with
#' `d_norm = 4`), quantile levels `seq(0, 1, 0.1)`, fine band \{2\} and broad
#' band \{6\}, and evaluates prediction of `y1` from the covariates by random
#' and spatial-block five-fold cross-validation (spatial blocks of side 6).
#'
#' @format A data frame with 900 rows and 6 columns:
#' \describe{
#'   \item{y1}{simulated response.}
#'   \item{x1, x2, x3}{simulated spatial covariates.}
#'   \item{x, y}{projected grid coordinates (column and row, 1--30, unit
#'     spacing).}
#' }
#'
#' @source Simulation data of Song, Y. (2026). Generalized covariate field
#'   (GCF): spatial-pattern and neighbourhood-distribution feature expansion
#'   improves geospatial prediction. *International Journal of Geographical
#'   Information Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @examples
#' data(sim_grid)
#' head(sim_grid)
#'
#' @seealso [gcf_field()], [gcf_select()], [bio_grid]
"sim_grid"

#' South-western Australia plant species richness grid (GCF paper case study)
#'
#' The case-study dataset of the GCF paper: vascular plant species richness
#' and twelve environmental covariates over the Southwest Australian
#' Floristic Region (SWAFR), aggregated to a 10-km grid of 6229 cells.
#' Species richness observations (4989 plot records) are aggregated into 958
#' grid cells; the remaining cells carry the covariates only. The rows with
#' `observed = TRUE` (equivalently `richness > 0`) are the observed sample
#' locations used to fit and validate prediction models; `richness` is
#' recorded as 0 where no sample exists (this is a missing observation, not a
#' true zero richness).
#'
#' The paper builds the GCF variables on the projected kilometre coordinates
#' `xkm`/`ykm` with buffers 20, 30, ..., 100 km (with `d_norm = 100`),
#' quantile levels `seq(0, 1, 0.05)`, fine band \{20, 30\} km and broad band
#' \{90, 100\} km, and selects variables with spatial blocks of side 132 km.
#'
#' @format A data frame with 6229 rows and 19 columns:
#' \describe{
#'   \item{GridID}{grid cell identifier.}
#'   \item{lon, lat}{cell centre longitude/latitude (degrees, GDA94).}
#'   \item{xkm, ykm}{projected cell centre coordinates in kilometres
#'     (EPSG:3577, GDA94 Australian Albers, equal-area); use these for GCF
#'     buffers and blocks.}
#'   \item{richness}{mean vascular plant species richness of the samples in
#'     the cell; 0 where the cell holds no sample.}
#'   \item{observed}{logical; `TRUE` for the 958 cells with observed
#'     richness.}
#'   \item{Elevation}{elevation (m).}
#'   \item{Slope}{terrain slope (degrees).}
#'   \item{Precipitation}{annual average precipitation (mm/year).}
#'   \item{Radiation}{annual average shortwave radiation (W/m^2).}
#'   \item{DistWater}{distance to the nearest water body (km).}
#'   \item{DistBuilt}{distance to artificial land or urban area (km).}
#'   \item{SoilN}{soil total nitrogen content (percent).}
#'   \item{SoilC}{soil organic carbon content (percent).}
#'   \item{SoilClay}{soil clay fraction (percent).}
#'   \item{SoilDepth}{depth to impermeable layer from surface (m).}
#'   \item{SoilpH}{average soil pH value.}
#'   \item{SoilBD}{soil bulk density (g/cm^3).}
#' }
#'
#' @source Case-study data of Song, Y. (2026). Generalized covariate field
#'   (GCF): spatial-pattern and neighbourhood-distribution feature expansion
#'   improves geospatial prediction. *International Journal of Geographical
#'   Information Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}.
#'   The richness grid and the covariates are the paper's case-study data,
#'   archived with the paper's reproducibility materials at
#'   \doi{10.6084/m9.figshare.29425739}; Table 3 of the paper lists the
#'   original source of each variable. The richness observations derive from
#'   Mokany, K. et al. (2022). Patterns and drivers of plant diversity across
#'   Australia. *Ecography*. \doi{10.1111/ecog.06426}. Covariates were
#'   derived from SRTM elevation, Digital Earth Australia, and CSIRO TERN
#'   Landscape soil products.
#'
#' @examples
#' data(bio_grid)
#' table(bio_grid$observed)
#' summary(bio_grid$richness[bio_grid$observed])
#'
#' @seealso [gcf_field()], [gcf_select()], [sim_grid]
"bio_grid"
