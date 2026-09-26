# =============================================================================
# zx.R -- Step 2 of the GCF method: neighbourhood-distribution features (Zx,
# the D layer). For each location v, variable x and buffer b, summarise the
# distribution of x over the locations within b by quantile levels tau:
# Z_x(v; b, tau) = Q_tau({ x(u) : u in N(v, b) }). Never uses a response.
# Arithmetic ported verbatim from dataR-v8 Code/gcf/features_zx.R.
# =============================================================================

# Pairwise Euclidean distance matrix (verbatim).
gcf_pairwise_dist <- function(query_coords, support_coords = query_coords) {
  query_coords <- gcf_as_coords(query_coords)
  support_coords <- gcf_as_coords(support_coords)
  same <- nrow(query_coords) == nrow(support_coords) &&
    isTRUE(all.equal(query_coords, support_coords, check.attributes = FALSE))
  if (same) return(as.matrix(stats::dist(as.matrix(query_coords))))
  dx <- outer(query_coords$x, support_coords$x, "-")
  dy <- outer(query_coords$y, support_coords$y, "-")
  sqrt(dx * dx + dy * dy)
}

# Fit the Zx plan: buffers, quantile levels, feature names, and group map
# (verbatim arithmetic from zx_fit_recipe).
zx_fit_plan <- function(X_support, coords_support, buffers, probs) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  coords_support <- gcf_as_coords(coords_support)
  buffers <- sort(unique(as.numeric(buffers)))
  probs <- sort(unique(as.numeric(probs)))
  gcf_assert(length(buffers) > 0L && all(is.finite(buffers)) && all(buffers > 0),
             "buffers must be positive numeric values.")
  gcf_assert(length(probs) > 0L && all(is.finite(probs)) && all(probs >= 0) && all(probs <= 1),
             "probs must be in [0, 1].")
  feature_names <- unlist(lapply(names(X_support), function(v) {
    unlist(lapply(buffers, function(b) {
      paste0(v, "_b", format(b, trim = TRUE, scientific = FALSE),
             "_q", format(probs, trim = TRUE, scientific = FALSE))
    }), use.names = FALSE)
  }), use.names = FALSE)
  group_map <- do.call(rbind, lapply(names(X_support), function(v) {
    do.call(rbind, lapply(buffers, function(b) {
      data.frame(buffer = b, prob = probs, stringsAsFactors = FALSE)
    })) |>
      transform(
        feature_name = paste0(v, "_b", format(buffer, trim = TRUE, scientific = FALSE),
                              "_q", format(prob, trim = TRUE, scientific = FALSE)),
        base_variable = v,
        group_id = paste0("D_", v),
        feature_category = "context_quantile",
        feature_type = "D"
      ) |>
      subset(select = c("feature_name", "base_variable", "group_id",
                        "feature_category", "feature_type", "buffer", "prob"))
  }))
  rownames(group_map) <- NULL
  list(
    vars = names(X_support),
    support_coords = coords_support,
    buffers = buffers,
    probs = probs,
    feature_names = feature_names,
    group_map = group_map
  )
}

# Apply the Zx plan (verbatim arithmetic from zx_apply_recipe). Empty buffers
# deliberately produce NA.
zx_apply_plan <- function(plan, X_support, coords_query,
                          query_to_support_row = NULL) {
  X_support <- gcf_as_numeric_df(X_support, "X_support")
  X_support <- X_support[, plan$vars, drop = FALSE]
  coords_query <- gcf_as_coords(coords_query)
  if (!is.null(query_to_support_row)) {
    gcf_assert(length(query_to_support_row) == nrow(coords_query),
               "query_to_support_row length must match coords_query rows.")
    Dmat <- gcf_pairwise_dist(plan$support_coords[query_to_support_row, , drop = FALSE],
                              plan$support_coords)
  } else {
    Dmat <- gcf_pairwise_dist(coords_query, plan$support_coords)
  }
  out <- matrix(NA_real_, nrow = nrow(coords_query), ncol = length(plan$feature_names))
  colnames(out) <- plan$feature_names
  col_offset <- 0L
  for (v in plan$vars) {
    x <- X_support[[v]]
    for (b in plan$buffers) {
      idx_by_row <- lapply(seq_len(nrow(Dmat)), function(i) which(Dmat[i, ] <= b))
      for (q in plan$probs) {
        col_offset <- col_offset + 1L
        out[, col_offset] <- vapply(idx_by_row, function(idx) {
          if (!length(idx)) return(NA_real_)
          as.numeric(stats::quantile(x[idx], probs = q, na.rm = TRUE, names = FALSE, type = 7))
        }, numeric(1))
      }
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

#' Neighbourhood-distribution features (Zx) of spatial variables
#'
#' Step 2 of the generalized covariate field (GCF) method. For each location
#' \eqn{v}, variable \eqn{x} and buffer radius \eqn{b}, `gcf_zx()` summarises
#' the distribution of the variable over the locations within the buffer by
#' quantile levels \eqn{\tau}:
#' \deqn{Z_x(v; b, \tau) = Q_\tau(\{ x(u) : u \in N(v, b) \}).}
#' The buffer neighbourhood includes the location itself. The computation
#' never uses a response variable.
#'
#' Non-finite covariate values are imputed by the column median and
#' zero-variance columns are dropped before feature construction. A buffer
#' that contains no locations yields `NA`.
#'
#' @inheritParams gcf_psi
#' @param probs Numeric vector of quantile levels in \[0, 1\]. Defaults to
#'   `seq(0, 1, 0.05)` (21 levels, as in the paper's case study).
#'
#' @return An object of class `"gcf_zx"`: a list with elements
#'   \describe{
#'     \item{features}{data frame of context-quantile features (one row per
#'       location; columns named `<var>_b<buffer>_q<prob>`).}
#'     \item{map}{data frame describing each feature column (base variable,
#'       group id, buffer, quantile level).}
#'     \item{vars}{the retained variable names.}
#'     \item{params}{the parameters used.}
#'   }
#'
#' @seealso [gcf_field()] for the full GCF variable generation pipeline,
#'   [gcf_psi()] for the spatial-pattern features.
#'
#' @references Song, Y. (2026). Generalized covariate field (GCF):
#'   spatial-pattern and neighbourhood-distribution feature expansion improves
#'   geospatial prediction. *International Journal of Geographical Information
#'   Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @examples
#' data(sim_grid)
#' sub <- sim_grid[sim_grid$x <= 10 & sim_grid$y <= 10, ]
#' zx <- gcf_zx(sub, coords = c("x", "y"), vars = c("x1", "x2"),
#'              buffers = c(2, 4), probs = c(0.1, 0.5, 0.9))
#' zx
#' head(zx$features)
#'
#' @export
gcf_zx <- function(data, coords, vars = NULL, buffers,
                   probs = seq(0, 1, 0.05)) {
  inp <- gcf_resolve_input(data, coords, vars)
  clean <- gcf_raw_clean(inp$data[, inp$vars, drop = FALSE])
  plan <- zx_fit_plan(clean$X, inp$coords, buffers = buffers, probs = probs)
  features <- zx_apply_plan(plan, clean$X, inp$coords,
                            query_to_support_row = seq_len(nrow(clean$X)))
  out <- list(
    features = features,
    map = plan$group_map,
    vars = plan$vars,
    dropped_vars = clean$dropped,
    params = list(buffers = plan$buffers, probs = plan$probs)
  )
  class(out) <- "gcf_zx"
  out
}

#' @export
print.gcf_zx <- function(x, ...) {
  cat("GCF neighbourhood-distribution features (Zx)\n")
  cat("  locations: ", nrow(x$features), "\n", sep = "")
  cat("  variables: ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  buffers:   ", paste(x$params$buffers, collapse = ", "),
      " | quantile levels: ", length(x$params$probs), "\n", sep = "")
  cat("  features:  ", ncol(x$features), "\n", sep = "")
  invisible(x)
}
