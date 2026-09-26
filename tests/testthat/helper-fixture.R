# Deterministic 7x7 fixture shared by all tests. The reference values in
# helper-reference.R were computed on exactly this fixture with the original
# GCF engine (paper code, dataR-v8); the fixture code must not change.
gcf_test_fixture <- function() {
  set.seed(42)
  fx <- expand.grid(x = as.numeric(1:7), y = as.numeric(1:7))
  fx$v1 <- sin(fx$x / 2) + cos(fx$y / 3) + rnorm(nrow(fx), sd = 0.3)
  fx$v2 <- fx$x * 0.1 + rnorm(nrow(fx))
  fx$resp <- 2 * fx$v1 - fx$v2 + rnorm(nrow(fx), sd = 0.2)
  fx
}

# Field on the fixture with the settings the reference was computed under.
# (suppressWarnings: spdep warns because the geocomplexity k = 23 exceeds a
# third of the 49 fixture locations.)
gcf_test_field <- function(fx = gcf_test_fixture()) {
  suppressWarnings(
    gcf_field(fx, coords = c("x", "y"), vars = c("v1", "v2"),
              buffers = c(2, 3), probs = c(0, 0.25, 0.5, 0.75, 1),
              d_norm = 3, fine_band = 2, broad_band = 3)
  )
}
