# Development notes

- [Rename to vamp and archive compatibility validation](rename-to-vamp.md)
- [Paper comparison: additional epsilon measurements and figures](dysco-paper-comparison.md)
- [Matched calibrated I/V and rate comparison against Dysco](compare-dysco.md)
- [Calibrated archive profiling, optimization and standard I/V validation](calibrated-performance.md)

- [Fixed epsilon D: expanded bit ladder and calibrated image comparisons](epsd-dense.md)
- [Sixteen-percent payload: G/S/D/R combinations, block sizes and calibrated images](fixed16.md)
- [Dimensionless G/S/D/R error budgets and reproducible global sampling](normalized-budgets.md)
- [Ten-percent visibility payload: six-width selection and repeated calibrated images](fixed10.md)
- [Earlier 3/4/6/8/12/16-bit ladder and three-bit conformance](precision-ladder-3bit.md)
- [Selected archive: image difference versus reference residual-map RMS](fixed25-residual-rms.md)
- [Fixed 25% visibility-payload budget: parameters and repeated calibrated images](fixed25.md)
- [Absolute versus relative error limits: matched-rate calibrated imaging](absolute_relative_error.md)
- [Gamma selection from maximum/std: matched blocks and calibrated CLEAN](gamma_max_std.md)
- [Calibrated dataset: deeper CLEAN, elevation cut and compression tradeoffs](test_calibrated_dataset.md)
- [Calibrated-only archives, four-bit weights and the 15-slot imaging test](calibrated-imaging.md)
- [Median offset and adaptive power-law quantization](power-quantization.md)
- [Full-range unsigned quantization and zero-bias validation](uint-midrise.md)
- [Two-level MXFP microscaling experiment](microscaling.md)
- [Profiling and the under-15-second encoder](performance.md)
- [Methods paper and reproducible figures](../docs/paper/vamp.md)
- [Initial development manuscript, preserved 25 October 2026](../docs/paper/work_note_20261025_init_work.md)
- [Design and experiment specification](design.md)
- [Initial source inspection findings](initial-findings.md)
- [Initial measured archive and CASA results](experiment-results.md)
- [Geometry/antenna/frequency grouping and WSClean validation](grouping-wsclean.md)
- [Opt-in destructive flag policy (`burn_flag`)](burn-flag.md)
- [Earlier 4/8/12/16-bit ladder and grouping settings](precision-ladder.md)
- [Repository publishing and workspace arrangement](publishing.md)

Machine-readable measurements are in `../reports/`. Large generated products
are local under `../results/` and excluded from Git. The user-facing archive
format contract remains in `../docs/format.md`.
