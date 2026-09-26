## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4
)

## ----setup--------------------------------------------------------------------
library(gcf)
has_rf <- requireNamespace("randomForest", quietly = TRUE)

## ----data---------------------------------------------------------------------
data(sim_grid)
head(sim_grid)

## ----field--------------------------------------------------------------------
field <- gcf_field(sim_grid, coords = c("x", "y"),
                   vars = c("x1", "x2", "x3"),
                   buffers = c(2, 4, 6), probs = seq(0, 1, 0.1),
                   d_norm = 4, fine_band = 2, broad_band = 6)
summary(field)

## ----meta---------------------------------------------------------------------
head(field$meta, 10)

## ----select-------------------------------------------------------------------
blocks <- gcf_blocks(sim_grid[, c("x", "y")], size = 6)
sel <- gcf_select(field, y = sim_grid$y1, blocks = blocks, B = 10, seed = 1)
sel

## ----freq---------------------------------------------------------------------
round(head(sel$freq, 8), 2)

## ----folds--------------------------------------------------------------------
folds_random <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  sample(rep_len(seq_len(K), n))
}
folds_blockwise <- function(block_id, K = 5, seed = 36) {
  set.seed(seed)
  ub <- sample(unique(block_id))
  as.integer(stats::setNames(rep_len(seq_len(K), length(ub)), ub)[block_id])
}
fold_rd <- folds_random(nrow(sim_grid))
fold_sp <- folds_blockwise(blocks)

## ----cv, eval = has_rf--------------------------------------------------------
library(randomForest)

y <- sim_grid$y1
X <- field$candidates
x_cols <- field$meta$feature[field$meta$category == "X"]

cv_rf <- function(fold, cols_by_fold, ntree = 300) {
  ks <- sort(unique(fold))
  r2 <- rmse <- numeric(length(ks))
  for (i in seq_along(ks)) {
    te <- which(fold == ks[i]); tr <- which(fold != ks[i])
    cols <- cols_by_fold[[i]]
    set.seed(1)
    fit <- randomForest(X[tr, cols, drop = FALSE], y[tr], ntree = ntree)
    p <- as.numeric(predict(fit, X[te, cols, drop = FALSE]))
    r2[i] <- 1 - sum((y[te] - p)^2) / sum((y[te] - mean(y[te]))^2)
    rmse[i] <- sqrt(mean((y[te] - p)^2))
  }
  c(R2 = mean(r2), RMSE = mean(rmse))
}

# per-fold GCF selection on the training part of each fold (leakage-free)
select_by_fold <- function(fold) {
  lapply(sort(unique(fold)), function(k) {
    gcf_select(field, y = y, blocks = blocks,
               train = which(fold != k), B = 10, seed = 1)$selected
  })
}

results <- do.call(rbind, lapply(
  list(random = fold_rd, spatial = fold_sp), function(fold) {
    S <- select_by_fold(fold)
    base <- cv_rf(fold, rep(list(x_cols), 5))
    gcfv <- cv_rf(fold, S)
    data.frame(feature_set = c("base", "GCF"),
               R2 = c(base["R2"], gcfv["R2"]),
               RMSE = c(base["RMSE"], gcfv["RMSE"]))
  }))
results$partition <- rep(c("random", "spatial"), each = 2)
rownames(results) <- NULL
results[, c("partition", "feature_set", "R2", "RMSE")]

## ----improvement, eval = has_rf-----------------------------------------------
imp <- do.call(rbind, lapply(split(results, results$partition), function(d) {
  data.frame(partition = d$partition[1],
             R2_gain_pct = 100 * (d$R2[2] - d$R2[1]) / d$R2[1],
             RMSE_drop_pct = 100 * (d$RMSE[1] - d$RMSE[2]) / d$RMSE[1])
}))
round(imp[, -1], 1)

## ----case, eval = FALSE-------------------------------------------------------
# data(bio_grid)
# obs <- bio_grid[bio_grid$observed, ]
# covs <- c("Elevation", "Slope", "Precipitation", "Radiation", "DistWater",
#           "DistBuilt", "SoilN", "SoilC", "SoilClay", "SoilDepth", "SoilpH",
#           "SoilBD")
# field <- gcf_field(obs, coords = c("xkm", "ykm"), vars = covs,
#                    buffers = seq(20, 100, 10), probs = seq(0, 1, 0.05),
#                    d_norm = 100, fine_band = c(20, 30),
#                    broad_band = c(90, 100))
# blocks <- gcf_blocks(obs[, c("xkm", "ykm")], size = 132)
# sel <- gcf_select(field, y = obs$richness, blocks = blocks, B = 80)

