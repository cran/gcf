# =============================================================================
# field.R -- the main GCF variable generator: raw spatial variables in, their
# generalized covariate field variables out (psi -> Zx -> functional
# reduction), plus print/summary methods for the result.
# =============================================================================

#' Generate generalized covariate field (GCF) variables
#'
#' The main generator of the GCF method: one or multiple spatial variables
#' with projected coordinates in, their GCF variables out. `gcf_field()` runs
#' the three feature-construction steps in one call:
#' \enumerate{
#'   \item [gcf_psi()] -- spatial-pattern features (11 operators over buffer
#'     radii);
#'   \item [gcf_zx()] -- neighbourhood-distribution features (buffer-wise
#'     quantiles);
#'   \item [gcf_reduce()] -- functional reduction of the buffer/quantile
#'     sweeps to a compact candidate set of X (raw), P (pattern), and D
#'     (context) variables.
#' }
#' The resulting candidate variables, together with the raw covariates, form
#' the GCF candidate field that [gcf_select()] screens for a stable subset.
#' No response variable is used at any point of the generation.
#'
#' @inheritParams gcf_psi
#' @inheritParams gcf_zx
#' @param fine_band,broad_band Buffer radii forming the fine and broad scale
#'   bands of the functional reduction. Default to `min(buffers)` and
#'   `max(buffers)`; the paper's case study uses `c(20, 30)` and
#'   `c(90, 100)` km with buffers 20--100 km.
#' @param d_mode Representation of the D layer in the reduction; see
#'   [gcf_reduce()]. Default `"functional"` (as in the paper).
#'
#' @return An object of class `"gcf_field"`: a list with elements
#'   \describe{
#'     \item{candidates}{data frame of candidate variables (one row per
#'       location): the raw covariates (category `"X"`) plus the GCF
#'       variables (categories `"P"` and `"D"`).}
#'     \item{meta}{data frame describing each candidate column: `feature`,
#'       `group` (base variable), `category` (`"X"`, `"P"`, `"D"`),
#'       `operator`, `scale`.}
#'     \item{vars}{the retained input variable names.}
#'     \item{psi}{the `"gcf_psi"` object (full pattern-feature sweep).}
#'     \item{zx}{the `"gcf_zx"` object (full context-quantile sweep).}
#'     \item{params}{the reduction parameters.}
#'   }
#'
#' @seealso [gcf_select()] for the follow-up variable selection;
#'   [gcf_psi()], [gcf_zx()], [gcf_reduce()] for the individual steps.
#'
#' @references Song, Y. (2026). Generalized covariate field (GCF):
#'   spatial-pattern and neighbourhood-distribution feature expansion improves
#'   geospatial prediction. *International Journal of Geographical Information
#'   Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @examples
#' # GCF variables of two simulated covariates on a 10 x 10 grid subset
#' data(sim_grid)
#' sub <- sim_grid[sim_grid$x <= 10 & sim_grid$y <= 10, ]
#' field <- gcf_field(sub, coords = c("x", "y"), vars = c("x1", "x2"),
#'                    buffers = c(2, 4), probs = seq(0, 1, 0.1), d_norm = 4)
#' field
#' head(field$candidates[, 1:6])
#'
#' \donttest{
#' # Full simulation grid, paper settings (buffers 2, 4, 6; 11 quantiles)
#' field <- gcf_field(sim_grid, coords = c("x", "y"),
#'                    vars = c("x1", "x2", "x3"),
#'                    buffers = c(2, 4, 6), probs = seq(0, 1, 0.1),
#'                    d_norm = 4, fine_band = 2, broad_band = 6)
#' summary(field)
#' }
#'
#' @export
gcf_field <- function(data, coords, vars = NULL, buffers,
                      probs = seq(0, 1, 0.05),
                      d_norm = max(buffers), theta = 2, bins = 10,
                      include_vario_exp = TRUE, vario_buffers = NULL,
                      fine_band = min(buffers), broad_band = max(buffers),
                      d_mode = "functional") {
  d_mode <- match.arg(d_mode, c("functional", "functional_allscale", "qgrid",
                                "full"))
  inp <- gcf_resolve_input(data, coords, vars)
  psi <- gcf_psi(inp$data, inp$coords, vars = inp$vars, buffers = buffers,
                 d_norm = d_norm, theta = theta, bins = bins,
                 include_vario_exp = include_vario_exp,
                 vario_buffers = vario_buffers)
  zx <- gcf_zx(inp$data, inp$coords, vars = inp$vars, buffers = buffers,
               probs = probs)
  gcf_reduce(inp$data, psi, zx, fine_band = fine_band,
             broad_band = broad_band, d_mode = d_mode)
}

#' @export
print.gcf_field <- function(x, ...) {
  tab <- table(factor(x$meta$category, levels = c("X", "P", "D")))
  cat("Generalized covariate field (GCF)\n")
  cat("  locations:  ", nrow(x$candidates), "\n", sep = "")
  cat("  variables:  ", length(x$vars), " (",
      paste(utils::head(x$vars, 5), collapse = ", "),
      if (length(x$vars) > 5) ", ..." else "", ")\n", sep = "")
  cat("  candidates: ", ncol(x$candidates),
      "  [X (raw) ", tab[["X"]], " | P (pattern) ", tab[["P"]],
      " | D (context) ", tab[["D"]], "]\n", sep = "")
  invisible(x)
}

#' Summarise a generalized covariate field
#'
#' @param object A `"gcf_field"` object from [gcf_field()] or [gcf_reduce()].
#' @param ... Unused.
#'
#' @return `object`, invisibly. Prints a per-variable and per-category
#'   breakdown of the candidate variables.
#'
#' @export
summary.gcf_field <- function(object, ...) {
  print(object)
  cat("\nCandidate variables per input variable and category:\n")
  tab <- table(object$meta$group, factor(object$meta$category,
                                         levels = c("X", "P", "D")))
  print(tab[object$vars, , drop = FALSE])
  cat("\nReduction: d_mode = ", object$params$d_mode,
      " | fine band {", paste(object$params$fine_band, collapse = ", "),
      "} | broad band {", paste(object$params$broad_band, collapse = ", "),
      "}\n", sep = "")
  cat("Full sweeps kept in $psi (", ncol(object$psi$features),
      " pattern features) and $zx (", ncol(object$zx$features),
      " context-quantile features).\n", sep = "")
  invisible(object)
}
