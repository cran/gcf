# =============================================================================
# psi.R -- Step 1 of the GCF method: spatial-pattern features (psi, the P
# layer). Per covariate the 11 spatial operators catalogued in the paper
# (LISA, local Geary's c, log local variance, rank-quantile entropy,
# geocomplexity, log scale-variance, local variogram exponent, and signed z-
# and MAD-outlier strengths) over the configured buffer radii. Never uses a
# response. Arithmetic ported verbatim from dataR-v8 Code/gcf/features_psi.R.
# =============================================================================

psi_fmt_scale <- function(x) {
  format(x, trim = TRUE, scientific = FALSE)
}

psi_multiscale_group_map <- function(vars, buffers, d_norm, gc_k = 23L,
                                     include_vario_exp = FALSE,
                                     vario_buffers = numeric(0)) {
  if (!isTRUE(include_vario_exp)) vario_buffers <- numeric(0)
  rows <- list()
  add_row <- function(v, tag, suffix, buffer = NA_real_) {
    rows[[length(rows) + 1L]] <<- data.frame(
      feature_name = paste0(v, "_", suffix),
      base_variable = v,
      group_id = paste0("P_", v),
      feature_category = tag,
      feature_type = "P",
      buffer = buffer,
      prob = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  for (v in vars) {
    for (b in buffers) {
      btxt <- psi_fmt_scale(b)
      add_row(v, "lisa", paste0("lisa_b", btxt, "_n", psi_fmt_scale(d_norm)), b)
    }
    for (b in buffers) {
      add_row(v, "geary", paste0("geary_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "lvar", paste0("lvar_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "qentropy", paste0("qentropy_b", psi_fmt_scale(b)), b)
    }
    add_row(v, "gc", paste0("gc_k", gc_k))
    add_row(v, "scalevar",
            paste0("scalevar_buffers", psi_fmt_scale(min(buffers)),
                   "_", psi_fmt_scale(max(buffers))))
    for (b in vario_buffers) {
      add_row(v, "vario_exp", paste0("vario_exp_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "poutlier_z", paste0("poutlier_z_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "noutlier_z", paste0("noutlier_z_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "poutlier_mad", paste0("poutlier_mad_b", psi_fmt_scale(b)), b)
    }
    for (b in buffers) {
      add_row(v, "noutlier_mad", paste0("noutlier_mad_b", psi_fmt_scale(b)), b)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Build neighbour indices from a distance matrix (verbatim).
gcf_neighbors_from_dist <- function(Dmat, radius, include_self = FALSE) {
  lapply(seq_len(nrow(Dmat)), function(i) {
    idx <- which(Dmat[i, ] <= radius & (include_self | Dmat[i, ] > 0))
    as.integer(idx)
  })
}

# Local normalization used internally for LISA (verbatim); not itself a
# pattern feature.
psi_local_norm <- function(x, nb_norm) {
  z <- numeric(length(x))
  for (i in seq_along(x)) {
    idx <- nb_norm[[i]]
    if (length(idx) >= 2L) {
      m <- mean(x[idx], na.rm = TRUE)
      s <- stats::sd(x[idx], na.rm = TRUE)
      z[i] <- if (is.finite(s) && s > 0) (x[i] - m) / s else 0
    }
  }
  z
}

psi_lisa <- function(z, nb_local) {
  vapply(seq_along(z), function(i) {
    idx <- nb_local[[i]]
    if (!length(idx)) return(0)
    zi <- mean(z[idx], na.rm = TRUE)
    out <- z[i] * zi
    if (is.finite(out)) out else 0
  }, numeric(1))
}

psi_geary <- function(x, nb_local, normalize = TRUE) {
  out <- vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (!length(idx)) return(0)
    mean((x[i] - x[idx])^2, na.rm = TRUE)
  }, numeric(1))
  if (isTRUE(normalize)) {
    v <- stats::var(x, na.rm = TRUE)
    if (is.finite(v) && v > 0) out <- out / v
  }
  out
}

psi_lvar <- function(x, nb_local) {
  log1p(vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (length(idx) < 2L) return(0)
    v <- stats::var(x[idx], na.rm = TRUE)
    if (is.finite(v)) v else 0
  }, numeric(1)))
}

psi_qentropy <- function(x, nb_local, bins = 10) {
  pct <- rank(x, ties.method = "average", na.last = "keep") / sum(!is.na(x))
  xb <- as.integer(cut(pct, breaks = seq(0, 1, length.out = bins + 1L),
                       include.lowest = TRUE, labels = FALSE))
  out <- vapply(seq_along(x), function(i) {
    idx <- nb_local[[i]]
    if (length(idx) < 2L) return(0)
    counts <- tabulate(xb[idx], nbins = bins)
    freq <- counts[counts > 0L] / length(idx)
    val <- -sum(freq * log(freq)) / log(bins)
    if (is.finite(val)) val else 0
  }, numeric(1))
  out
}

psi_scalevar <- function(x, coords, scales) {
  coords <- gcf_as_coords(coords)
  sv_mat <- matrix(NA_real_, nrow = length(x), ncol = length(scales))
  for (j in seq_along(scales)) {
    s <- scales[[j]]
    key <- paste(floor(coords$x / s), floor(coords$y / s), sep = "_")
    gm <- tapply(x, key, mean, na.rm = TRUE)
    sv_mat[, j] <- gm[key]
  }
  vals <- apply(sv_mat, 1, function(row) {
    v <- stats::var(row, na.rm = TRUE)
    if (is.finite(v)) v else 0
  })
  log1p(vals)
}

psi_outlier_z <- function(x, nb_local, theta = 2) {
  pos <- neg <- numeric(length(x))
  for (i in seq_along(x)) {
    idx <- nb_local[[i]]
    if (length(idx) >= 2L) {
      xi <- x[idx]
      s <- stats::sd(xi, na.rm = TRUE)
      if (is.finite(s) && s > 0) {
        z <- (xi - mean(xi, na.rm = TRUE)) / s
        pos[i] <- sum(abs(z[z > theta]), na.rm = TRUE)
        neg[i] <- sum(abs(z[z < -theta]), na.rm = TRUE)
      }
    }
  }
  data.frame(poutlier_z = pos, noutlier_z = neg)
}

psi_outlier_mad <- function(x, nb_local, theta = 2) {
  pos <- neg <- numeric(length(x))
  mad_floor <- 1e-3 * stats::mad(x, na.rm = TRUE) + 1e-9
  for (i in seq_along(x)) {
    idx <- nb_local[[i]]
    if (length(idx) >= 2L) {
      xi <- x[idx]
      ctr <- stats::median(xi, na.rm = TRUE)
      scl <- max(stats::mad(xi, na.rm = TRUE), mad_floor)
      z <- (xi - ctr) / scl
      pos[i] <- sum(abs(z[z > theta]), na.rm = TRUE)
      neg[i] <- sum(abs(z[z < -theta]), na.rm = TRUE)
    }
  }
  data.frame(poutlier_mad = pos, noutlier_mad = neg)
}

# Coordinate-only geocomplexity plan (verbatim).
psi_geocomplexity_plan <- function(coords, k = 23) {
  coords <- gcf_as_coords(coords)
  k <- min(as.integer(k), nrow(coords) - 1L)
  if (k < 1L) return(list(k = 0L))
  sf_base <- sf::st_as_sf(coords, coords = c("x", "y"), crs = sf::NA_crs_)
  nb <- spdep::knn2nb(spdep::knearneigh(sf_base, k = k))
  weights <- spdep::nb2mat(nb, style = "W")
  list(geometry = sf::st_geometry(sf_base), crs = sf::st_crs(sf_base),
       weights = weights, k = k)
}

# Geocomplexity vector from a saved plan (verbatim).
psi_geocomplexity <- function(plan, x, variable) {
  if (isTRUE(plan$k < 1L)) return(numeric(length(x)))
  sf_data <- sf::st_sf(stats::setNames(data.frame(x), variable),
                       geometry = plan$geometry, crs = plan$crs)
  gc <- geocomplexity::geocd_vector(sf_data, wt = plan$weights,
                                    method = "moran", normalize = TRUE)
  val <- gc[[paste0("GC_", variable)]]
  val[!is.finite(val)] <- 0
  val
}

psi_dist_pairs <- function(idx) {
  m <- length(idx)
  np <- m * (m - 1L) / 2L
  a <- integer(np)
  b <- integer(np)
  k <- 1L
  for (col in seq_len(m - 1L)) {
    for (row in seq.int(col + 1L, m)) {
      a[k] <- idx[row]
      b[k] <- idx[col]
      k <- k + 1L
    }
  }
  list(a = a, b = b)
}

# Deterministic local-variogram plan (verbatim); dependency-light backend, no
# external variogram package involved.
psi_vario_plan <- function(coords, nb_local, radius, nlags = 4L) {
  coords <- gcf_as_coords(coords)
  edges <- seq(0, radius, length.out = nlags + 1L)
  lapply(seq_len(nrow(coords)), function(i) {
    idx <- c(i, nb_local[[i]])
    m <- length(idx)
    if (m < 6L) return(NULL)
    d <- as.numeric(stats::dist(as.matrix(coords[idx, , drop = FALSE])))
    ok <- d > 0
    d <- d[ok]
    if (length(d) < nlags) return(NULL)
    pairs <- psi_dist_pairs(idx)
    pairs$a <- pairs$a[ok]
    pairs$b <- pairs$b[ok]
    lab <- cut(d, breaks = edges, include.lowest = TRUE, labels = FALSE)
    keep <- !is.na(lab)
    if (sum(keep) < nlags) return(NULL)
    d <- d[keep]
    lab <- lab[keep]
    pairs$a <- pairs$a[keep]
    pairs$b <- pairs$b[keep]
    hh <- tapply(d, lab, mean)
    wt <- tapply(d, lab, length)
    list(a = pairs$a, b = pairs$b, lab = lab, hh = hh, wt = wt)
  })
}

# Local variogram exponent from a saved plan (verbatim).
psi_vario_exp <- function(x, plan) {
  ex <- numeric(length(x))
  for (i in seq_along(plan)) {
    p <- plan[[i]]
    if (is.null(p)) next
    sv <- (x[p$a] - x[p$b])^2 / 2
    gh <- tapply(sv, p$lab, mean)
    hh <- p$hh[names(gh)]
    wt <- p$wt[names(gh)]
    keep <- is.finite(gh) & gh > 0 & is.finite(hh) & hh > 0
    if (sum(keep) >= 2L) {
      lg <- as.numeric(log(gh[keep] + 1e-9))
      lh <- as.numeric(log(hh[keep]))
      w <- as.numeric(wt[keep])
      W <- sum(w)
      mx <- sum(w * lh) / W
      my <- sum(w * lg) / W
      sxx <- sum(w * (lh - mx)^2)
      sxy <- sum(w * (lh - mx) * (lg - my))
      if (is.finite(sxx) && sxx > 0) ex[i] <- sxy / sxx
    }
  }
  ex
}

# Fit the psi plan: neighbour lists, geocomplexity plan, variogram plans, and
# feature/group metadata (verbatim arithmetic from psi_fit_recipe).
psi_fit_plan <- function(X_support, coords_support, psi_buffers, d_norm,
                         theta, bins, include_vario_exp, vario_buffers) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  coords_support <- gcf_as_coords(coords_support)
  psi_buffers <- sort(unique(as.numeric(psi_buffers)))
  gcf_assert(length(psi_buffers) > 0L && all(is.finite(psi_buffers)) &&
               all(psi_buffers > 0),
             "buffers must be positive numeric values.")
  d_norm <- as.numeric(d_norm)
  gcf_assert(length(d_norm) == 1L && is.finite(d_norm) && d_norm > 0,
             "d_norm must be a positive number.")
  scales <- psi_buffers
  theta <- as.numeric(theta)
  bins <- as.integer(bins)
  include_vario_exp <- isTRUE(include_vario_exp)
  vario_buffers <- if (include_vario_exp) {
    if (is.null(vario_buffers)) psi_buffers else as.numeric(vario_buffers)
  } else {
    numeric(0)
  }
  vario_buffers <- sort(unique(vario_buffers))
  if (include_vario_exp) {
    gcf_assert(length(vario_buffers) > 0L && all(is.finite(vario_buffers)) &&
                 all(vario_buffers > 0),
               "vario_buffers must be positive numeric values when include_vario_exp is TRUE.")
  }
  Dmat <- gcf_pairwise_dist(coords_support)
  neighbor_buffers <- sort(unique(c(psi_buffers, vario_buffers)))
  buffer_keys <- paste0("b", psi_fmt_scale(neighbor_buffers))
  nb_by_buffer <- stats::setNames(lapply(neighbor_buffers, function(b) {
    gcf_neighbors_from_dist(Dmat, b, include_self = FALSE)
  }), buffer_keys)
  nb_norm <- gcf_neighbors_from_dist(Dmat, d_norm, include_self = FALSE)
  gc_plan <- psi_geocomplexity_plan(coords_support, k = 23L)
  vario_plans <- if (include_vario_exp) {
    vario_keys <- paste0("b", psi_fmt_scale(vario_buffers))
    stats::setNames(lapply(seq_along(vario_buffers), function(i) {
      psi_vario_plan(coords_support, nb_by_buffer[[vario_keys[[i]]]],
                     radius = vario_buffers[[i]], nlags = 4L)
    }), vario_keys)
  } else {
    list()
  }
  group_map <- psi_multiscale_group_map(names(X_support), psi_buffers, d_norm,
                                        gc_k = 23L,
                                        include_vario_exp = include_vario_exp,
                                        vario_buffers = vario_buffers)
  list(
    vars = names(X_support),
    support_coords = coords_support,
    d_norm = d_norm,
    psi_buffers = psi_buffers,
    include_vario_exp = include_vario_exp,
    vario_buffers = vario_buffers,
    scales = scales,
    theta = theta,
    bins = bins,
    nb_by_buffer = nb_by_buffer,
    nb_norm = nb_norm,
    geocomplexity_plan = gc_plan,
    vario_plans = vario_plans,
    feature_names = group_map$feature_name,
    group_map = group_map
  )
}

# Apply the psi plan to its support rows (verbatim arithmetic from
# psi_apply_recipe).
psi_apply_plan <- function(plan, X_support) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  X_support <- X_support[, plan$vars, drop = FALSE]
  out <- matrix(NA_real_, nrow = nrow(X_support), ncol = length(plan$feature_names))
  colnames(out) <- plan$feature_names
  col <- 0L
  for (v in plan$vars) {
    x <- X_support[[v]]
    z <- psi_local_norm(x, plan$nb_norm)
    block <- list()
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      nb <- plan$nb_by_buffer[[bkey]]
      block[[paste0(v, "_lisa_b", btxt, "_n", psi_fmt_scale(plan$d_norm))]] <- psi_lisa(z, nb)
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_geary_b", btxt)]] <- psi_geary(x, plan$nb_by_buffer[[bkey]])
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_lvar_b", btxt)]] <- psi_lvar(x, plan$nb_by_buffer[[bkey]])
    }
    for (b in plan$psi_buffers) {
      bkey <- paste0("b", psi_fmt_scale(b))
      btxt <- psi_fmt_scale(b)
      block[[paste0(v, "_qentropy_b", btxt)]] <- psi_qentropy(x, plan$nb_by_buffer[[bkey]],
                                                              plan$bins)
    }
    block[[paste0(v, "_gc_k23")]] <- psi_geocomplexity(plan$geocomplexity_plan, x, v)
    block[[paste0(v, "_scalevar_buffers", psi_fmt_scale(min(plan$psi_buffers)),
                  "_", psi_fmt_scale(max(plan$psi_buffers)))]] <-
      psi_scalevar(x, plan$support_coords, plan$scales)
    if (isTRUE(plan$include_vario_exp)) {
      for (b in plan$vario_buffers) {
        bkey <- paste0("b", psi_fmt_scale(b))
        btxt <- psi_fmt_scale(b)
        block[[paste0(v, "_vario_exp_b", btxt)]] <-
          psi_vario_exp(x, plan$vario_plans[[bkey]])
      }
    }
    outlier_z_by_buffer <- lapply(plan$psi_buffers, function(b) {
      psi_outlier_z(x, plan$nb_by_buffer[[paste0("b", psi_fmt_scale(b))]], plan$theta)
    })
    outlier_mad_by_buffer <- lapply(plan$psi_buffers, function(b) {
      psi_outlier_mad(x, plan$nb_by_buffer[[paste0("b", psi_fmt_scale(b))]], plan$theta)
    })
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_poutlier_z_b", btxt)]] <- outlier_z_by_buffer[[i]]$poutlier_z
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_noutlier_z_b", btxt)]] <- outlier_z_by_buffer[[i]]$noutlier_z
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_poutlier_mad_b", btxt)]] <- outlier_mad_by_buffer[[i]]$poutlier_mad
    }
    for (i in seq_along(plan$psi_buffers)) {
      btxt <- psi_fmt_scale(plan$psi_buffers[[i]])
      block[[paste0(v, "_noutlier_mad_b", btxt)]] <- outlier_mad_by_buffer[[i]]$noutlier_mad
    }
    block <- as.data.frame(block, stringsAsFactors = FALSE, check.names = FALSE)
    rng <- seq.int(col + 1L, col + ncol(block))
    gcf_assert(identical(colnames(out)[rng], names(block)),
               "psi block column order mismatch for ", v)
    out[, rng] <- as.matrix(block)
    col <- col + ncol(block)
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Spatial-pattern features (psi) of spatial variables
#'
#' Step 1 of the generalized covariate field (GCF) method. For each input
#' variable, `gcf_psi()` computes the 11 spatial-pattern operators of the GCF
#' method over a series of buffer radii: local indicator of spatial
#' association (LISA) on locally normalized values, local Geary's c, log local
#' variance, rank-binned quantile entropy, geocomplexity (k = 23 nearest
#' neighbours, single scale), log scale-variance (single scale), local
#' variogram exponent, and positive/negative z-score and MAD outlier
#' strengths. The computation never uses a response variable.
#'
#' Non-finite covariate values are imputed by the column median and
#' zero-variance columns are dropped before feature construction.
#'
#' @param data A data frame with one row per location, containing the spatial
#'   variables (and optionally the coordinate columns).
#' @param coords Either a length-2 character vector naming the projected
#'   coordinate columns of `data` (for example `c("x", "y")`), or a
#'   two-column matrix or data frame of projected coordinates. Coordinate
#'   units must match `buffers` (for example kilometres).
#' @param vars Character vector of variable names to expand. Defaults to all
#'   columns of `data` except the coordinate columns; make sure the response
#'   and any non-numeric columns (which raise an error) are excluded when
#'   relying on the default.
#' @param buffers Numeric vector of positive buffer radii, in coordinate
#'   units (the paper uses 20--100 km for the case study and 2, 4, 6 grid
#'   units for the simulation).
#' @param d_norm Radius of the local normalization neighbourhood used by the
#'   LISA operator. Defaults to `max(buffers)`.
#' @param theta Outlier threshold for the z-score and MAD outlier strengths
#'   (default 2, as in the paper).
#' @param bins Number of rank bins for the quantile entropy operator
#'   (default 10, as in the paper).
#' @param include_vario_exp Logical; include the local variogram exponent
#'   operator (default `TRUE`, the full 11-operator catalogue).
#' @param vario_buffers Buffer radii for the local variogram exponent.
#'   Defaults to `buffers`.
#'
#' @return An object of class `"gcf_psi"`: a list with elements
#'   \describe{
#'     \item{features}{data frame of pattern features (one row per location).}
#'     \item{map}{data frame describing each feature column (base variable,
#'       group id, operator category, buffer).}
#'     \item{vars}{the retained variable names.}
#'     \item{params}{the parameters used.}
#'   }
#'
#' @seealso [gcf_field()] for the full GCF variable generation pipeline,
#'   [gcf_zx()] for the neighbourhood-distribution features.
#'
#' @references Song, Y. (2026). Generalized covariate field (GCF):
#'   spatial-pattern and neighbourhood-distribution feature expansion improves
#'   geospatial prediction. *International Journal of Geographical Information
#'   Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @examples
#' data(sim_grid)
#' sub <- sim_grid[sim_grid$x <= 10 & sim_grid$y <= 10, ]
#' psi <- gcf_psi(sub, coords = c("x", "y"), vars = c("x1", "x2"),
#'                buffers = c(2, 4), d_norm = 4)
#' psi
#' head(psi$features[, 1:4])
#'
#' @export
gcf_psi <- function(data, coords, vars = NULL, buffers,
                    d_norm = max(buffers), theta = 2, bins = 10,
                    include_vario_exp = TRUE, vario_buffers = NULL) {
  inp <- gcf_resolve_input(data, coords, vars)
  clean <- gcf_raw_clean(inp$data[, inp$vars, drop = FALSE])
  plan <- psi_fit_plan(clean$X, inp$coords, psi_buffers = buffers,
                       d_norm = d_norm, theta = theta, bins = bins,
                       include_vario_exp = include_vario_exp,
                       vario_buffers = vario_buffers)
  features <- psi_apply_plan(plan, clean$X)
  out <- list(
    features = features,
    map = plan$group_map,
    vars = plan$vars,
    dropped_vars = clean$dropped,
    params = list(buffers = plan$psi_buffers, d_norm = plan$d_norm,
                  theta = plan$theta, bins = plan$bins,
                  include_vario_exp = plan$include_vario_exp,
                  vario_buffers = plan$vario_buffers, gc_k = 23L)
  )
  class(out) <- "gcf_psi"
  out
}

#' @export
print.gcf_psi <- function(x, ...) {
  cat("GCF spatial-pattern features (psi)\n")
  cat("  locations: ", nrow(x$features), "\n", sep = "")
  cat("  variables: ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  operators: ", length(unique(x$map$feature_category)),
      " | buffers: ", paste(x$params$buffers, collapse = ", "), "\n", sep = "")
  cat("  features:  ", ncol(x$features), "\n", sep = "")
  invisible(x)
}
