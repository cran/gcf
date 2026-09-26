# =============================================================================
# utils.R -- internal utilities: error helpers, input coercion (numeric data
# frames, coordinates), stage-1 raw covariate cleaning, and input resolution
# shared by all exported generators. Arithmetic ported verbatim from the
# original GCF engine (dataR-v8, Code/gcf/utils.R + preprocess.R).
# =============================================================================

gcf_stop <- function(...) {
  stop(paste0("[gcf] ", paste0(..., collapse = "")), call. = FALSE)
}

gcf_assert <- function(ok, ...) {
  if (!isTRUE(ok)) gcf_stop(...)
  invisible(TRUE)
}

# Validate user data as a numeric data frame with stable names. Column
# handling is verbatim from the original engine except that non-numeric
# columns (character, factor, logical) raise an informative error instead of
# being coerced silently; integer columns are accepted as numeric.
gcf_as_numeric_df <- function(X, arg_name = "X") {
  if (is.null(X)) return(NULL)
  X <- as.data.frame(X, stringsAsFactors = FALSE)
  if (!ncol(X)) gcf_stop(arg_name, " must have at least one column.")
  if (is.null(names(X)) || any(!nzchar(names(X)))) {
    gcf_stop(arg_name, " must have non-empty column names.")
  }
  if (anyDuplicated(names(X))) {
    gcf_stop(arg_name, " has duplicated column names: ",
             paste(unique(names(X)[duplicated(names(X))]), collapse = ", "))
  }
  bad <- names(X)[!vapply(X, is.numeric, logical(1))]
  if (length(bad)) {
    types <- vapply(X[bad], function(col) class(col)[1L], character(1))
    gcf_stop(arg_name, " must contain numeric columns only; non-numeric ",
             "column(s): ", paste0(bad, " (", types, ")", collapse = ", "),
             ". Convert them or exclude them via `vars`.")
  }
  X
}

# Convert coordinates to a two-column numeric data frame (verbatim).
gcf_as_coords <- function(coords) {
  coords <- as.data.frame(coords, stringsAsFactors = FALSE)
  gcf_assert(ncol(coords) >= 2, "coords must have at least two columns.")
  coords <- coords[, seq_len(2), drop = FALSE]
  names(coords) <- c("x", "y")
  coords$x <- suppressWarnings(as.numeric(coords$x))
  coords$y <- suppressWarnings(as.numeric(coords$y))
  gcf_assert(all(is.finite(coords$x)) && all(is.finite(coords$y)),
             "coords must be finite numeric values.")
  coords
}

# Stage-1 raw covariate cleaning (verbatim arithmetic from
# gcf_raw_x_clean_fit + gcf_raw_x_clean_apply): impute non-finite values by
# the column median (0 if the median is not finite) and drop zero-variance
# columns. Pattern and neighbourhood-distribution features are always built
# from this cleaned, raw-scale matrix.
gcf_raw_clean <- function(X) {
  X <- gcf_as_numeric_df(X, "data")
  impute_values <- vapply(X, function(col) {
    ok <- is.finite(col)
    med <- stats::median(col[ok], na.rm = TRUE)
    if (is.finite(med)) med else 0
  }, numeric(1))
  X_imp <- X
  for (nm in names(X_imp)) {
    bad <- !is.finite(X_imp[[nm]])
    if (any(bad)) X_imp[[nm]][bad] <- impute_values[[nm]]
  }
  zero_var <- vapply(X_imp, function(col) {
    vals <- unique(col[is.finite(col)])
    length(vals) <= 1L
  }, logical(1))
  retained <- names(X_imp)[!zero_var]
  dropped <- names(X_imp)[zero_var]
  gcf_assert(length(retained) > 0L,
             "Stage-1 cleaning dropped all covariate columns.")
  list(
    X = X_imp[, retained, drop = FALSE],
    retained = retained,
    dropped = dropped,
    impute_values = impute_values
  )
}

# Resolve the (data, coords, vars) user contract shared by gcf_psi, gcf_zx,
# and gcf_field: `coords` is either a length-2 character vector naming the
# coordinate columns of `data`, or a two-column matrix / data frame of
# projected coordinates; `vars` defaults to every other column of `data`.
gcf_resolve_input <- function(data, coords, vars = NULL) {
  gcf_assert(is.data.frame(data) || is.matrix(data),
             "data must be a data frame (or matrix).")
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  coord_cols <- character(0)
  if (is.character(coords)) {
    gcf_assert(length(coords) == 2L,
               "coords must name exactly two coordinate columns.")
    missing <- setdiff(coords, names(data))
    gcf_assert(!length(missing), "data is missing coordinate column(s): ",
               paste(missing, collapse = ", "))
    coord_cols <- coords
    coords_df <- gcf_as_coords(data[, coords, drop = FALSE])
  } else {
    coords_df <- gcf_as_coords(coords)
    gcf_assert(nrow(coords_df) == nrow(data),
               "coords must have one row per row of data.")
  }
  if (is.null(vars)) {
    vars <- setdiff(names(data), coord_cols)
  } else {
    vars <- as.character(vars)
    missing <- setdiff(vars, names(data))
    gcf_assert(!length(missing), "data is missing variable column(s): ",
               paste(missing, collapse = ", "))
  }
  gcf_assert(length(vars) > 0L, "At least one spatial variable is required.")
  list(data = data, coords = coords_df, vars = vars)
}
