# =============================================================================
# select.R -- Step 3b of the GCF method: variable selection. Random forest
# impurity importance (top-K per subsample) combined with spatial-block
# stability resampling and medium (variable x category) group voting.
# Arithmetic ported verbatim from dataR-v8 Code/gcf/gcf_model.R::
# gcf_select_rfimp; gcf_blocks ported from Code/functions/validation.R.
# =============================================================================

# Verbatim core of gcf_select_rfimp (dataR-v8 gcf_model.R): X is the reduced
# candidate matrix, meta has feature/group/category, block_id covers all rows.
# The `sel` matrix is returned so the wrapper can report selection
# frequencies; everything the original computed is computed identically.
# The RNG is not seeded here: gcf_select() seeds (and afterwards restores)
# the R generator only when the user supplies `seed`, and the same seed is
# passed to ranger, so `seed = 1` reproduces the original engine exactly.
gcf_select_rfimp_core <- function(X, y, meta, train_rows, block_id,
                                  params = list()) {
  p <- utils::modifyList(list(B = 80L, subsample_frac = 0.7, pi = 0.6,
                              ktop = 20L, rf_trees = 200L), params)
  is_X <- meta$category == "X"; pen <- which(!is_X); x_cols <- meta$feature[is_X]
  Ztr <- X[train_rows, , drop = FALSE]; ytr <- y[train_rows]
  btr <- block_id[train_rows]; blocks <- unique(btr)
  kt <- floor(p$subsample_frac * length(blocks))
  sel <- matrix(FALSE, p$B, length(pen))
  for (b in seq_len(p$B)) {
    sub <- which(btr %in% sample(blocks, kt))
    rf <- ranger::ranger(x = Ztr[sub, , drop = FALSE], y = ytr[sub],
                         num.trees = p$rf_trees, importance = "impurity",
                         seed = p$seed, num.threads = 1)
    imp <- rf$variable.importance[meta$feature[pen]]; imp[is.na(imp)] <- 0
    sel[b, ] <- rank(-imp, ties.method = "first") <= p$ktop
  }
  key <- paste(meta$group[pen], meta$category[pen])
  reps <- unlist(lapply(unique(key), function(kk) {
    idx <- which(key == kk)
    if (mean(rowSums(sel[, idx, drop = FALSE]) > 0) < p$pi) return(NULL)
    meta$feature[pen][idx][which.max(colMeans(sel[, idx, drop = FALSE]))]
  }))
  list(selected = c(x_cols, reps), x_cols = x_cols, reps = reps,
       sel = sel, pen = pen, key = key, params = p)
}

#' Spatial blocks for stability selection and spatial cross-validation
#'
#' Assigns each location to a square spatial block of the given side length,
#' on projected coordinates. Block ids are used by [gcf_select()] for
#' spatial-block stability resampling, and can also define spatial
#' cross-validation folds.
#'
#' @param coords A two-column matrix or data frame of projected coordinates.
#' @param size Block side length, in coordinate units (the paper's case study
#'   uses 132 km, twice the residual variogram range; the simulation uses 6
#'   grid units).
#'
#' @return A character vector of block ids, one per location.
#'
#' @examples
#' data(sim_grid)
#' blocks <- gcf_blocks(sim_grid[, c("x", "y")], size = 6)
#' table(blocks)
#'
#' @export
gcf_blocks <- function(coords, size) {
  coords <- gcf_as_coords(coords)
  gcf_assert(length(size) == 1L && is.finite(size) && size > 0,
             "size must be a positive number.")
  cx <- coords[, 1]; cy <- coords[, 2]
  paste0("bx", floor((cx - min(cx)) / size), "_by", floor((cy - min(cy)) / size))
}

#' Select stable GCF variables
#'
#' Step 3b of the generalized covariate field (GCF) method: screens the
#' candidate variables of a [gcf_field()] for a stable subset. The raw
#' covariates (category `"X"`) are always kept; the derived P/D variables are
#' selected by the rf_imp kernel with spatial-block stability and medium
#' group voting:
#' \enumerate{
#'   \item **Per-subsample selector (rf_imp).** On each subsample, fit a
#'     random forest (\pkg{ranger}, `num_trees` trees) and keep the top
#'     `ktop` derived variables by impurity importance.
#'   \item **Spatial-block stability.** Repeat `B` times on subsamples drawn
#'     as `floor(subsample_frac * M)` of the `M` spatial blocks, recording
#'     which variables are kept each time.
#'   \item **Medium group voting.** Group the derived variables by (variable
#'     x category); a group qualifies if any member is kept in at least
#'     `pi_thr` of the subsamples, and contributes its most frequently kept
#'     member as representative.
#' }
#' Group voting de-dilutes collinear blocks: sibling variables split the
#' selection frequency so that none clears the threshold individually, while
#' the block as a whole fires reliably.
#'
#' The random forests always run single-threaded. When `seed` is supplied the
#' selection is fully reproducible: the seed governs both the block resampling
#' and the ranger forests, and the caller's random number generator state is
#' restored on exit. When `seed = NULL` (the default) no seed is set and the
#' procedure follows the current random number generator stream; set a seed
#' yourself before the call for reproducibility.
#'
#' @param x A `"gcf_field"` object from [gcf_field()] or [gcf_reduce()] (or a
#'   list with elements `candidates` and `meta` in the same format).
#' @param y Numeric response vector, one value per location of the field.
#'   The response is used for selection only; `gcf_field()` never sees it.
#' @param blocks Vector of spatial block ids, one per location; see
#'   [gcf_blocks()].
#' @param train Optional integer vector of training row indices to run the
#'   selection on (for example the training part of a cross-validation
#'   fold). Defaults to all rows.
#' @param B Number of stability resamples (default 80, as in the paper).
#' @param subsample_frac Fraction of spatial blocks drawn per resample
#'   (default 0.7).
#' @param pi_thr Group fire-frequency threshold (default 0.6).
#' @param ktop Number of top derived variables kept per subsample
#'   (default 20).
#' @param num_trees Number of trees of the random forest importance kernel
#'   (default 200).
#' @param seed Optional integer seed for the block resampling and the random
#'   forests. `NULL` (default) leaves the random number generator untouched;
#'   the paper uses `seed = 1`. When supplied, the generator state is
#'   restored on exit.
#'
#' @return An object of class `"gcf_selection"`: a list with elements
#'   \describe{
#'     \item{selected}{character vector of selected variable names (forced
#'       raw covariates first, then the group representatives).}
#'     \item{forced}{the raw covariate names (always kept).}
#'     \item{derived}{the selected derived (P/D) variable names.}
#'     \item{freq}{named numeric vector: per-variable selection frequency
#'       across the `B` resamples (derived variables, sorted decreasing).}
#'     \item{group_fire}{named numeric vector: per-group fire frequency.}
#'     \item{params}{the selection parameters.}
#'   }
#'
#' @seealso [gcf_field()], [gcf_blocks()].
#'
#' @references Song, Y. (2026). Generalized covariate field (GCF):
#'   spatial-pattern and neighbourhood-distribution feature expansion improves
#'   geospatial prediction. *International Journal of Geographical Information
#'   Science*, 40, 1--29. \doi{10.1080/13658816.2026.2729719}
#'
#' @examples
#' data(sim_grid)
#' sub <- sim_grid[sim_grid$x <= 12 & sim_grid$y <= 12, ]
#' field <- gcf_field(sub, coords = c("x", "y"), vars = c("x1", "x2", "x3"),
#'                    buffers = c(2, 4), probs = seq(0, 1, 0.1), d_norm = 4)
#' blocks <- gcf_blocks(sub[, c("x", "y")], size = 6)
#' # B reduced from the paper's 80 to keep the example fast
#' sel <- gcf_select(field, y = sub$y1, blocks = blocks, B = 5, seed = 1)
#' sel
#'
#' \donttest{
#' # Full simulation grid with the paper's field settings; B = 20 resamples
#' # here to keep the example short (the paper uses B = 80; see the vignette)
#' field <- gcf_field(sim_grid, coords = c("x", "y"),
#'                    vars = c("x1", "x2", "x3"),
#'                    buffers = c(2, 4, 6), probs = seq(0, 1, 0.1),
#'                    d_norm = 4, fine_band = 2, broad_band = 6)
#' blocks <- gcf_blocks(sim_grid[, c("x", "y")], size = 6)
#' sel <- gcf_select(field, y = sim_grid$y1, blocks = blocks, B = 20, seed = 1)
#' sel$selected
#' }
#'
#' @export
gcf_select <- function(x, y, blocks, train = NULL, B = 80,
                       subsample_frac = 0.7, pi_thr = 0.6, ktop = 20,
                       num_trees = 200, seed = NULL) {
  gcf_assert(is.list(x) && !is.null(x$candidates) && !is.null(x$meta),
             "x must be a gcf_field object (or a list with candidates and meta).")
  X <- as.matrix(x$candidates)
  gcf_assert(is.numeric(X), "candidate variables must be numeric.")
  meta <- x$meta
  y <- as.numeric(y)
  gcf_assert(length(y) == nrow(X), "y must have one value per location.")
  gcf_assert(all(is.finite(y)), "y must be finite.")
  gcf_assert(length(blocks) == nrow(X),
             "blocks must have one id per location; see gcf_blocks().")
  if (is.null(train)) train <- seq_len(nrow(X))
  if (!is.null(seed)) {
    gcf_assert(length(seed) == 1L && is.numeric(seed) && is.finite(seed),
               "seed must be a single finite number or NULL.")
    # Seed only on request, and hand the caller's RNG state back on exit.
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
              add = TRUE)
    } else {
      on.exit(rm(list = ".Random.seed", envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
  }
  core <- gcf_select_rfimp_core(
    X, y, meta, train_rows = train, block_id = blocks,
    params = list(B = as.integer(B), subsample_frac = subsample_frac,
                  pi = pi_thr, ktop = as.integer(ktop),
                  rf_trees = as.integer(num_trees), seed = seed))
  freq <- colMeans(core$sel)
  names(freq) <- meta$feature[core$pen]
  group_fire <- vapply(unique(core$key), function(kk) {
    idx <- which(core$key == kk)
    mean(rowSums(core$sel[, idx, drop = FALSE]) > 0)
  }, numeric(1))
  out <- list(
    selected = core$selected,
    forced = core$x_cols,
    derived = if (is.null(core$reps)) character(0) else core$reps,
    freq = sort(freq, decreasing = TRUE),
    group_fire = sort(group_fire, decreasing = TRUE),
    params = list(B = as.integer(B), subsample_frac = subsample_frac,
                  pi_thr = pi_thr, ktop = as.integer(ktop),
                  num_trees = as.integer(num_trees), seed = seed,
                  n_train = length(train)),
    n_candidates = ncol(X)
  )
  class(out) <- "gcf_selection"
  out
}

#' @export
print.gcf_selection <- function(x, ...) {
  cat("GCF variable selection (rf_imp + spatial-block stability + group voting)\n")
  cat("  candidates: ", x$n_candidates, " | resamples B = ", x$params$B,
      " | pi threshold = ", x$params$pi_thr, "\n", sep = "")
  cat("  selected:   ", length(x$selected), " (", length(x$forced),
      " forced raw + ", length(x$derived), " derived)\n", sep = "")
  cat("  forced:  ", paste(x$forced, collapse = ", "), "\n", sep = "")
  if (length(x$derived)) {
    cat("  derived: ", paste(x$derived, collapse = ", "), "\n", sep = "")
  } else {
    cat("  derived: (none passed the group vote)\n")
  }
  invisible(x)
}
