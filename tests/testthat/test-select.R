# Selection: agreement with the original engine and structural behaviour.

test_that("gcf_select reproduces the original rf_imp selection", {
  # Bit-identical reproduction of the reference selection depends on ranger
  # internals (tie handling of impurity sums); run locally and in CI only.
  skip_on_cran()
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  sel <- gcf_select(field, y = fx$resp, blocks = blocks,
                    B = 5, ktop = 8, num_trees = 100, seed = 1)
  ref <- gcf_test_reference()
  # identical feature set, in identical order (forced covariates first)
  expect_identical(sel$selected, ref$selected)
})

test_that("gcf_select is deterministic for a fixed seed", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  s1 <- gcf_select(field, y = fx$resp, blocks = blocks, B = 3, num_trees = 50,
                   seed = 1)
  s2 <- gcf_select(field, y = fx$resp, blocks = blocks, B = 3, num_trees = 50,
                   seed = 1)
  expect_identical(s1$selected, s2$selected)
  expect_identical(s1$freq, s2$freq)
  # seed = NULL follows the caller's RNG stream
  set.seed(5)
  s3 <- gcf_select(field, y = fx$resp, blocks = blocks, B = 3, num_trees = 50)
  set.seed(5)
  s4 <- gcf_select(field, y = fx$resp, blocks = blocks, B = 3, num_trees = 50)
  expect_identical(s3$freq, s4$freq)
  expect_null(s3$params$seed)
})

test_that("gcf_select restores the caller's RNG state when seed is given", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  set.seed(123)
  before <- runif(3)
  set.seed(123)
  invisible(gcf_select(field, y = fx$resp, blocks = blocks, B = 2,
                       num_trees = 20, seed = 7))
  after <- runif(3)
  expect_identical(before, after)
})

test_that("selection frequencies and group fires are well formed", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  sel <- gcf_select(field, y = fx$resp, blocks = blocks,
                    B = 5, ktop = 8, num_trees = 100, seed = 1)
  expect_s3_class(sel, "gcf_selection")
  expect_identical(sel$forced, c("v1", "v2"))
  expect_identical(sel$selected, c(sel$forced, sel$derived))
  n_derived <- sum(field$meta$category != "X")
  expect_length(sel$freq, n_derived)
  expect_true(all(sel$freq >= 0 & sel$freq <= 1))
  expect_true(all(sel$group_fire >= 0 & sel$group_fire <= 1))
  # every selected derived variable belongs to a group that fired >= pi_thr
  expect_true(all(sel$group_fire[
    paste(field$meta$group, field$meta$category)[
      match(sel$derived, field$meta$feature)]] >= sel$params$pi_thr))
})

test_that("gcf_blocks tiles the coordinates", {
  fx <- gcf_test_fixture()
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  expect_length(blocks, nrow(fx))
  expect_identical(blocks[1], "bx0_by0")
  expect_true(all(grepl("^bx[0-9]+_by[0-9]+$", blocks)))
  expect_error(gcf_blocks(fx[, c("x", "y")], size = -1), "positive")
})

test_that("gcf_select validates its inputs", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  expect_error(gcf_select(field, y = fx$resp[-1], blocks = blocks),
               "one value per location")
  expect_error(gcf_select(field, y = fx$resp, blocks = blocks[-1]),
               "one id per location")
  expect_error(gcf_select(list(), y = fx$resp, blocks = blocks),
               "gcf_field")
  expect_error(gcf_select(field, y = fx$resp, blocks = blocks, seed = "a"),
               "seed")
})

test_that("print methods run", {
  fx <- gcf_test_fixture()
  field <- gcf_test_field(fx)
  blocks <- gcf_blocks(fx[, c("x", "y")], size = 3)
  sel <- gcf_select(field, y = fx$resp, blocks = blocks, B = 3, num_trees = 50,
                    seed = 1)
  expect_output(print(field), "Generalized covariate field")
  expect_output(summary(field), "Candidate variables per input variable")
  expect_output(print(field$psi), "psi")
  expect_output(print(field$zx), "Zx")
  expect_output(print(sel), "selected")
})
