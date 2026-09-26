# =============================================================================
# imports.R -- namespace imports.
# =============================================================================

#' @importFrom stats dist mad median quantile sd setNames var
#' @importFrom utils head modifyList
NULL

# Column names used non-standardly inside transform() calls (verbatim ports
# of the original engine's manifest construction).
utils::globalVariables(c("buffer", "prob", "feature_type"))
