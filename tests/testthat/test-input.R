# Input handling: cleaning, defaults, error paths, and the bundled data.

test_that("non-finite values are median-imputed and zero-variance dropped", {
  fx <- gcf_test_fixture()
  fx$v3 <- 1                      # zero variance -> dropped
  fx$v1[3] <- NA                  # -> imputed by median before features
  zx <- gcf_zx(fx, coords = c("x", "y"), vars = c("v1", "v2", "v3"),
               buffers = 3, probs = 0.5)
  expect_identical(zx$vars, c("v1", "v2"))
  expect_identical(zx$dropped_vars, "v3")
  expect_false(anyNA(zx$features))
  # imputation happens on the variable before quantiles are taken
  fx2 <- gcf_test_fixture()
  fx2$v1[3] <- stats::median(fx$v1, na.rm = TRUE)
  zx2 <- gcf_zx(fx2, coords = c("x", "y"), vars = c("v1", "v2"),
                buffers = 3, probs = 0.5)
  expect_equal(zx$features[["v1_b3_q0.5"]], zx2$features[["v1_b3_q0.5"]],
               tolerance = 1e-15)
})

test_that("vars defaults to all non-coordinate columns", {
  fx <- gcf_test_fixture()[, c("x", "y", "v1", "v2")]
  zx <- gcf_zx(fx, coords = c("x", "y"), buffers = 3, probs = 0.5)
  expect_identical(zx$vars, c("v1", "v2"))
})

test_that("bad inputs raise clear errors", {
  fx <- gcf_test_fixture()
  expect_error(gcf_zx(fx, coords = c("x", "y", "v1"), buffers = 3),
               "exactly two")
  expect_error(gcf_zx(fx, coords = c("x", "nope"), buffers = 3),
               "missing coordinate")
  expect_error(gcf_zx(fx, coords = c("x", "y"), vars = "nope", buffers = 3),
               "missing variable")
  expect_error(gcf_zx(fx, coords = c("x", "y"), vars = "v1", buffers = -2),
               "positive")
  expect_error(gcf_zx(fx, coords = c("x", "y"), vars = "v1", buffers = 2,
                      probs = 2), "0, 1")
  expect_error(suppressWarnings(
    gcf_psi(fx, coords = c("x", "y"), vars = "v1", buffers = 2, d_norm = -1)),
    "d_norm")
  psi <- suppressWarnings(
    gcf_psi(fx, coords = c("x", "y"), vars = "v1", buffers = 2, d_norm = 3))
  zx <- gcf_zx(fx, coords = c("x", "y"), vars = c("v1", "v2"), buffers = 2)
  expect_error(gcf_reduce(fx, psi, zx), "same variables")
})

test_that("bundled datasets are intact", {
  data("sim_grid", package = "gcf", envir = environment())
  expect_identical(dim(sim_grid), c(900L, 6L))
  expect_identical(names(sim_grid), c("y1", "x1", "x2", "x3", "x", "y"))
  data("bio_grid", package = "gcf", envir = environment())
  expect_identical(nrow(bio_grid), 6229L)
  expect_identical(sum(bio_grid$observed), 958L)
  expect_true(all(bio_grid$richness[bio_grid$observed] > 0))
  expect_true(all(bio_grid$richness[!bio_grid$observed] == 0))
})

test_that("non-numeric variable columns raise an informative error", {
  fx <- gcf_test_fixture()
  fx$id <- paste0("s", seq_len(nrow(fx)))
  expect_error(gcf_zx(fx, coords = c("x", "y"), buffers = 3, probs = 0.5),
               "non-numeric column\\(s\\): id \\(character\\)")
  fx$id <- NULL
  fx$flag <- fx$v1 > 0
  expect_error(gcf_zx(fx, coords = c("x", "y"), vars = c("v1", "flag"),
                      buffers = 3, probs = 0.5), "flag \\(logical\\)")
  fx$flag <- NULL
  fx$grp <- factor(rep(c("a", "b"), length.out = nrow(fx)))
  expect_error(gcf_zx(fx, coords = c("x", "y"), vars = "grp", buffers = 3),
               "grp \\(factor\\)")
  # integer columns are numeric and pass
  fx$grp <- NULL
  fx$cnt <- seq_len(nrow(fx))
  expect_s3_class(gcf_zx(fx, coords = c("x", "y"), vars = "cnt", buffers = 3,
                         probs = 0.5), "gcf_zx")
})

test_that("d_mode is validated up front", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  expect_error(gcf_reduce(fx, field$psi, field$zx, d_mode = "bogus"),
               "should be one of")
  expect_error(gcf_field(fx, coords = c("x", "y"), vars = "v1", buffers = 2,
                         d_mode = "bogus"), "should be one of")
  full <- gcf_reduce(fx, field$psi, field$zx, d_mode = "full")
  expect_identical(full$params$d_mode, "full")
})
