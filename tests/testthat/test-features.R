# Numerical agreement with the original GCF engine (paper code, dataR-v8).
# The reference values in helper-reference.R were computed on the fixture of
# helper-fixture.R with the original engine; the package must reproduce them
# to machine precision.

test_that("psi features match the original engine", {
  fx <- gcf_test_fixture()
  psi <- suppressWarnings(
    gcf_psi(fx, coords = c("x", "y"), vars = c("v1", "v2"),
            buffers = c(2, 3), d_norm = 3)
  )
  ref <- gcf_test_reference()
  expect_identical(names(psi$features), ref$psi_names)
  expect_equal(unname(as.numeric(psi$features[1, ])), ref$psi_row1,
               tolerance = 1e-12)
  expect_equal(unname(colMeans(as.matrix(psi$features))), ref$psi_colmeans,
               tolerance = 1e-12)
  expect_identical(nrow(psi$features), nrow(fx))
  expect_identical(psi$map$feature_name, names(psi$features))
})

test_that("Zx features match the original engine", {
  fx <- gcf_test_fixture()
  zx <- gcf_zx(fx, coords = c("x", "y"), vars = c("v1", "v2"),
               buffers = c(2, 3), probs = c(0, 0.25, 0.5, 0.75, 1))
  ref <- gcf_test_reference()
  expect_identical(names(zx$features), ref$zx_names)
  expect_equal(unname(as.numeric(zx$features[1, ])), ref$zx_row1,
               tolerance = 1e-12)
  expect_equal(unname(colMeans(as.matrix(zx$features))), ref$zx_colmeans,
               tolerance = 1e-12)
  expect_identical(zx$map$feature_name, names(zx$features))
})

test_that("reduced candidate field matches the original engine", {
  field <- gcf_test_field()
  ref <- gcf_test_reference()
  expect_identical(names(field$candidates), ref$red_names)
  expect_equal(unname(as.numeric(field$candidates[1, ])), ref$red_row1,
               tolerance = 1e-12)
  expect_equal(unname(colMeans(as.matrix(field$candidates))), ref$red_colmeans,
               tolerance = 1e-12)
  expect_identical(field$meta$category, ref$red_categories)
})

test_that("gcf_field equals the composition of its step functions", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  psi <- suppressWarnings(
    gcf_psi(fx, coords = c("x", "y"), vars = c("v1", "v2"),
            buffers = c(2, 3), d_norm = 3)
  )
  zx <- gcf_zx(fx, coords = c("x", "y"), vars = c("v1", "v2"),
               buffers = c(2, 3), probs = c(0, 0.25, 0.5, 0.75, 1))
  field2 <- gcf_reduce(fx, psi, zx, fine_band = 2, broad_band = 3)
  expect_identical(field$candidates, field2$candidates)
  expect_identical(field$meta, field2$meta)
})

test_that("coords can be given as columns names or as a matrix", {
  fx <- gcf_test_fixture()
  zx1 <- gcf_zx(fx, coords = c("x", "y"), vars = "v1",
                buffers = 2, probs = 0.5)
  zx2 <- gcf_zx(fx["v1"], coords = as.matrix(fx[, c("x", "y")]),
                buffers = 2, probs = 0.5)
  expect_identical(zx1$features, zx2$features)
})
