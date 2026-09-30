# gcf 0.1.1

* Fixed platform-dependent geocomplexity features (`*_gc_k23` / `*_P_gc`).
  `spdep::knearneigh()` switches to a kd-tree search when the optional
  dbscan package is installed, and the two searches break ties at the k-th
  neighbour distance differently. On gridded data such ties are common, so
  the neighbour sets (and hence the geocomplexity values) depended on the
  machine; this made the reference tests fail on some CRAN check flavours.
  The geocomplexity neighbours now always use spdep's brute-force search,
  whose ties go to the lower row index, which is also what the paper's
  results were computed with. Data without distance ties (for example the
  `bio_grid` case study) are unaffected.
* New test pinning the tie-breaking rule of the geocomplexity neighbours.
* Now requires spdep >= 1.1-7 (the first version with the `use_kd_tree`
  argument of `knearneigh()`).

# gcf 0.1.0

* Initial CRAN release.
* GCF variable generation: `gcf_field()`, with step functions `gcf_psi()`
  (11 spatial-pattern operators over buffer radii), `gcf_zx()` (buffer-wise
  neighbourhood quantiles) and `gcf_reduce()` (functional reduction to the
  X/P/D candidate variables).
* Variable selection: `gcf_select()` (random forest impurity importance with
  spatial-block stability resampling and group voting) and `gcf_blocks()`.
  `gcf_select()` seeds the random number generator only when `seed` is
  supplied (`seed = 1` reproduces the paper) and restores the caller's
  generator state on exit.
* Datasets: `sim_grid` (paper simulation, 900 grid cells) and `bio_grid`
  (south-western Australia plant species richness, 6229 grid cells).
* Vignette: "The GCF workflow: from spatial variables to better predictions".
* Implements Song (2026) <doi:10.1080/13658816.2026.2729719>; feature
  generation and selection reproduce the paper's reference implementation
  to machine precision.
