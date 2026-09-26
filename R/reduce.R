# =============================================================================
# reduce.R -- Step 3a of the GCF method: functional reduction. Collapses the
# collinear quantile/buffer sweeps of the psi and Zx layers into a compact
# set of interpretable functionals per (variable, scale band):
#   D layer: median, IQR, low tail (q0.10), high tail (q0.90), skew
#            (q0.90 + q0.10 - 2 q0.50) of the band-averaged quantile curve;
#   P layer: band-averaged buffered operators (single-scale gc and scalevar
#            kept as-is);
#   X layer: the raw covariates passed through unchanged.
# Arithmetic ported verbatim from dataR-v8
# Code/gcf/select_gcf.R::gcf_build_reduced_candidates.
# =============================================================================

# Verbatim port of gcf_build_reduced_candidates (dataR-v8 select_gcf.R).
# data: data frame holding raw covariate columns plus all psi and Zx feature
# columns; manifest: block ("X"/"P"/"D"), feature_name, base_variable,
# feature_category, buffer, prob.
gcf_build_reduced_candidates <- function(data, manifest,
                                         fine_band  = c(20, 30),
                                         broad_band = c(90, 100),
                                         vars = NULL,
                                         d_mode = c("functional",
                                                    "functional_allscale",
                                                    "qgrid", "full"),
                                         d_qgrid = c(0.05, 0.10, 0.25, 0.50,
                                                     0.75, 0.90, 0.95)) {
  d_mode <- match.arg(d_mode)
  man <- manifest
  if (is.null(vars)) vars <- unique(man$base_variable[man$block == "X"])

  bands <- list(fine = fine_band, broad = broad_band)
  cols  <- list()   # named numeric vectors
  meta  <- list()   # rows: feature, group, category, operator, scale

  add <- function(name, vec, group, category,
                  operator = NA_character_, scale = NA_character_) {
    cols[[name]]  <<- vec
    meta[[length(meta) + 1L]] <<- data.frame(
      feature = name, group = group, category = category,
      operator = operator, scale = scale, stringsAsFactors = FALSE)
  }

  # helper: band-averaged column for a manifest row subset
  band_avg <- function(colnames_in_band) {
    if (length(colnames_in_band) == 1L) return(data[[colnames_in_band]])
    rowMeans(as.matrix(data[, colnames_in_band, drop = FALSE]), na.rm = TRUE)
  }

  for (v in vars) {

    # ---- X: raw covariate (forced) ----
    add(v, data[[v]], group = v, category = "X")

    # ---- D: context-quantile features (representation set by d_mode) ----
    dman <- man[man$block == "D" & man$base_variable == v, ]
    if (d_mode %in% c("functional", "functional_allscale")) {
      # functionals of the quantile curve, per scale band
      d_bands <- if (d_mode == "functional") bands else
        stats::setNames(lapply(sort(unique(dman$buffer)), identity),
                        paste0("b", sort(unique(dman$buffer))))
      for (bn in names(d_bands)) {
        bset <- d_bands[[bn]]
        qcol <- function(p) {
          sub <- dman[dman$buffer %in% bset, ]
          avail <- unique(sub$prob)
          pn <- avail[which.min(abs(avail - p))]   # nearest available quantile
          band_avg(sub$feature_name[abs(sub$prob - pn) < 1e-9])
        }
        q10 <- qcol(0.10); q25 <- qcol(0.25); q50 <- qcol(0.50)
        q75 <- qcol(0.75); q90 <- qcol(0.90)
        pre <- paste0(v, "_D_", bn, "_")
        add(paste0(pre, "med"),    q50,                 v, "D", scale = bn)
        add(paste0(pre, "iqr"),    q75 - q25,           v, "D", scale = bn)
        add(paste0(pre, "lotail"), q10,                 v, "D", scale = bn)
        add(paste0(pre, "hitail"), q90,                 v, "D", scale = bn)
        add(paste0(pre, "skew"),   q90 + q10 - 2 * q50, v, "D", scale = bn)
      }
    } else if (d_mode == "qgrid") {
      # raw quantile values, all buffers x a coarse quantile grid (no functional)
      keep <- dman[vapply(dman$prob, function(p)
        any(abs(p - d_qgrid) < 1e-9), logical(1)), ]
      for (j in seq_len(nrow(keep))) {
        add(keep$feature_name[j], data[[keep$feature_name[j]]], v, "D",
            scale = paste0("b", keep$buffer[j]))
      }
    } else {  # "full": all D columns as-is, no reduction
      for (j in seq_len(nrow(dman))) {
        add(dman$feature_name[j], data[[dman$feature_name[j]]], v, "D",
            scale = paste0("b", dman$buffer[j]))
      }
    }

    # ---- P: pattern operators, band-averaged ----
    pman <- man[man$block == "P" & man$base_variable == v, ]
    ops  <- unique(pman$feature_category)
    for (op in ops) {
      opman <- pman[pman$feature_category == op, ]
      has_buffers <- any(!is.na(opman$buffer))
      if (has_buffers) {
        for (bn in names(bands)) {
          bset <- bands[[bn]]
          nm <- opman$feature_name[opman$buffer %in% bset]
          if (length(nm) == 0L) next
          add(paste0(v, "_P_", op, "_", bn), band_avg(nm), v, "P",
              operator = op, scale = bn)
        }
      } else {
        # single-scale operator (gc, scalevar)
        add(paste0(v, "_P_", op), data[[opman$feature_name[1]]], v, "P",
            operator = op, scale = "single")
      }
    }
  }

  meta_df <- do.call(rbind, meta)
  X <- do.call(cbind, cols)
  colnames(X) <- names(cols)
  rownames(meta_df) <- NULL
  stopifnot(identical(colnames(X), meta_df$feature))
  list(X = X, meta = meta_df)
}

# Assemble the modeling frame + manifest expected by
# gcf_build_reduced_candidates from a raw data frame and fitted psi/Zx layers.
gcf_assemble_manifest <- function(psi_map, zx_map, vars) {
  keep <- c("block", "feature_name", "base_variable", "feature_category",
            "buffer", "prob")
  rbind(
    data.frame(block = "X", feature_name = vars, base_variable = vars,
               feature_category = "raw", buffer = NA_real_, prob = NA_real_,
               stringsAsFactors = FALSE),
    transform(psi_map, block = feature_type)[, keep],
    transform(zx_map, block = feature_type)[, keep]
  )
}

#' Reduce psi and Zx layers to the GCF candidate variables
#'
#' Step 3a of the generalized covariate field (GCF) method. The buffer and
#' quantile sweeps produced by [gcf_psi()] and [gcf_zx()] are collinear;
#' `gcf_reduce()` collapses them into a compact set of interpretable
#' functionals over two scale bands (a leakage-free per-row transform):
#' \itemize{
#'   \item **D (context)**: per variable and band, from the band-averaged
#'     quantile curve, five functionals: `med` (q0.50), `iqr` (q0.75 -
#'     q0.25), `lotail` (q0.10), `hitail` (q0.90), and `skew` (q0.90 + q0.10
#'     - 2 q0.50).
#'   \item **P (pattern)**: per variable and operator, the band-averaged
#'     buffered operators; the single-scale operators (geocomplexity, scale
#'     variance) are kept as-is.
#'   \item **X**: the raw covariates passed through unchanged.
#' }
#'
#' @param data The data frame that was passed to [gcf_psi()] and [gcf_zx()]
#'   (supplies the raw covariate columns).
#' @param psi A `"gcf_psi"` object from [gcf_psi()].
#' @param zx A `"gcf_zx"` object from [gcf_zx()].
#' @param fine_band Numeric vector of buffer radii forming the fine scale
#'   band. Defaults to the smallest buffer used; the paper's case study uses
#'   `c(20, 30)` km.
#' @param broad_band Numeric vector of buffer radii forming the broad scale
#'   band. Defaults to the largest buffer used; the paper's case study uses
#'   `c(90, 100)` km.
#' @param d_mode Representation of the D layer: `"functional"` (default, the
#'   five functionals per band as in the paper), `"functional_allscale"`
#'   (functionals at every buffer), `"qgrid"` (raw quantile columns on a
#'   coarse grid), or `"full"` (all quantile columns unchanged).
#'
#' @return An object of class `"gcf_field"`; see [gcf_field()] for its
#'   structure.
#'
#' @seealso [gcf_field()], which runs [gcf_psi()], [gcf_zx()], and
#'   `gcf_reduce()` in one call.
#'
#' @examples
#' data(sim_grid)
#' sub <- sim_grid[sim_grid$x <= 10 & sim_grid$y <= 10, ]
#' psi <- gcf_psi(sub, coords = c("x", "y"), vars = c("x1", "x2"),
#'                buffers = c(2, 4), d_norm = 4)
#' zx <- gcf_zx(sub, coords = c("x", "y"), vars = c("x1", "x2"),
#'              buffers = c(2, 4), probs = seq(0, 1, 0.1))
#' field <- gcf_reduce(sub, psi, zx, fine_band = 2, broad_band = 4)
#' field
#'
#' @export
gcf_reduce <- function(data, psi, zx, fine_band = NULL, broad_band = NULL,
                       d_mode = "functional") {
  d_mode <- match.arg(d_mode, c("functional", "functional_allscale", "qgrid",
                                "full"))
  gcf_assert(inherits(psi, "gcf_psi"), "psi must be a gcf_psi object.")
  gcf_assert(inherits(zx, "gcf_zx"), "zx must be a gcf_zx object.")
  gcf_assert(identical(psi$vars, zx$vars),
             "psi and zx must be built from the same variables.")
  gcf_assert(nrow(psi$features) == nrow(zx$features),
             "psi and zx must be built from the same locations.")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  gcf_assert(nrow(data) == nrow(psi$features),
             "data must have the same rows as the psi/zx features.")
  missing <- setdiff(psi$vars, names(data))
  gcf_assert(!length(missing), "data is missing variable column(s): ",
             paste(missing, collapse = ", "))
  bufs <- sort(unique(c(psi$params$buffers, zx$params$buffers)))
  if (is.null(fine_band)) fine_band <- min(bufs)
  if (is.null(broad_band)) broad_band <- max(bufs)
  if (d_mode == "functional") {
    gcf_assert(any(zx$params$buffers %in% fine_band),
               "fine_band must contain at least one Zx buffer radius.")
    gcf_assert(any(zx$params$buffers %in% broad_band),
               "broad_band must contain at least one Zx buffer radius.")
  }
  modeling <- cbind(data[, psi$vars, drop = FALSE], psi$features, zx$features)
  man <- gcf_assemble_manifest(psi$map, zx$map, psi$vars)
  bad <- setdiff(man$feature_name, names(modeling))
  gcf_assert(!length(bad),
             "Feature name/manifest mismatch (buffer radii of mixed decimal ",
             "widths, e.g. c(1.5, 3), are not supported; use a consistent ",
             "buffer series): ", paste(utils::head(bad, 3), collapse = ", "))
  red <- gcf_build_reduced_candidates(modeling, man,
                                      fine_band = fine_band,
                                      broad_band = broad_band,
                                      d_mode = d_mode)
  out <- list(
    candidates = as.data.frame(red$X, stringsAsFactors = FALSE),
    meta = red$meta,
    vars = psi$vars,
    psi = psi,
    zx = zx,
    params = list(fine_band = fine_band, broad_band = broad_band,
                  d_mode = d_mode)
  )
  class(out) <- "gcf_field"
  out
}
