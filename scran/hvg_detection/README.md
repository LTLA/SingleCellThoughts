# Description of the HVG machinery in _scran_

This directory contains some files describing the HVG detection machinery in _scran_.

- `description.tex` focuses on the theoretical basis behind `trendVar` and `decomposeVar`.
- `simulations/sim_*.R` contains simulation scripts to compare the performance of `decomposeVar` with `technicalCV2` and `improvedCV2` for detecting HVGs.
   Each simulation script describes a different type of variability, though this is not particularly important for gene-wise testing.
- `moresims/simple_tests.R` assesses the type I error rate for the tests in `decomposeVar`, `technicalCV2` and `improvedCV2`.

