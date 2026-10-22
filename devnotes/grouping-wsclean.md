# Grouping and WSClean experiment, 2026-10-22

This experiment tests the user's proposal to group numbers by baseline length,
orientation, antenna endpoint and frequency channel before quantization. It uses
the same copied 62,128-row, 192-channel, four-product integration as the initial
prototype. Results are pending while the runs complete.

## Matched compression configurations

`configs/benchmark-grouping.json` fixes adaptive 2/4/6/8-bit candidates, exact
fallback, stored float32 steps, quantile candidates and numerical ceilings.
The ceilings are engineering examples, not validated thermal-noise/science budgets.
Except for the single-baseline control, baseline groups have at most 32 rows.
The main ordering comparisons use 64-channel segments. Amplitude sorting also
tests 16, 32 and 192 channels, without changing archive channel resolution.

- `class-original`: auto/cross and no-valid classes are contiguous, original
  row order is retained within each class. This controls the earlier interleaved
  autocorrelation boundary confound.
- `length`: ascending norm of original metre UVW, preserving ordered pairs.
- `orientation`: atan2(V,U) modulo 2π, sorted within auto/cross classes and split
  at 24 fixed angular-bin edges (15°). This is projected UV bearing, not geographic
  antenna bearing. There is no π folding or antenna swapping.
- `length-orientation`: factor-of-two norm bins referenced to 1 metre, crossed
  with those angular bins; original row order inside each joint bin.
- `antenna`: original row order within a fixed ANTENNA1 or ANTENNA2 endpoint;
  ordered endpoint identity is never canonicalized.
- `amplitude` and `buckets`: the previously documented per-product statistics.

No-valid samples, zero projected UV and nonfinite geometry receive deterministic
separate groups as appropriate. Partial groups are recorded in size histograms;
all actual permutations and boundaries are stored and counted. Additional bin
boundaries are part of each experiment's cost and are not hidden.

Independent lanes run concurrently on a shared host, with one BLAS/OpenMP thread
per codec process. Timing is diagnostic only; this is not isolated throughput.

## WSClean validation

The user requested 2048×2048 images, 4-arcmin pixels and multiscale CLEAN.
The available local container is `astronrd/linc:5.2.1`, resolving to
`sha256:03ffbf70c829e220707b9843e96765f5886005b4b05a26344a0b5bb4d9ece8c6`.
It contains WSClean 3.7, commit `v3.7-79-g2d5c1ed8`. Docker networking is disabled
for imaging, and only a dedicated output directory is mounted writable.

The experiment uses DATA explicitly, Stokes I, Briggs 0, `-mgain 0.8`, automatic
mask threshold 5, stopping threshold 1, at most 10,000 CLEAN components, 8 threads
and an 8-GB WSClean memory limit. The input MS is copied before imaging and its
file hashes are checked afterward. Commands, image identity, source hashes and
output FITS inventory are saved for each run. No beam correction or calibration
is added, and no archive values or UVW are regenerated.

The initial reference reached the iteration cap. The fitted restoring beam varied
slightly between repeated lossless runs, motivating a common fixed restoring beam
for the controlled comparison. Independent automatic masks and nonlinear CLEAN
paths remain part of the end-to-end experiment.

The fixed-beam, eight-thread repeat gave a relative dirty-image RMS difference of
4.39668×10⁻⁷ (0.00004397%) but a restored-image RMS difference of 0.00681772
(0.68177%). That repeatability floor must accompany the quantized comparison;
the restored-image difference cannot automatically be attributed to compression.
Single-thread reference/repeat controls are being tested separately.

Comparisons require matching FITS axes, coordinates, units and restoring beam.
They measure dirty, PSF, model, residual and restored images; celestial-invalid
SIN corners are excluded. Reports include full valid-sky errors, bright-reference
pixels, a fixed aperture around the absolute reference peak, peak shifts and CLEAN
iteration counts. The peak is not assumed to be the Sun, aperture pixel sums are
not claimed as calibrated flux, and image RMS is not a thermal-noise estimate.

Reference documentation: [WSClean multiscale CLEAN](https://wsclean.readthedocs.io/en/latest/multiscale_cleaning.html),
[basic CLEAN controls](https://wsclean.readthedocs.io/en/v3.5/basic_cleaning.html).

## Reproduction

Use Python with `pip install -e '.[validation]'` and an available Docker daemon.
The WSClean runner resolves an already-installed image to its immutable ID; it
does not pull a container or enable container networking.

```bash
visarchive benchmark INPUT.ms --matrix configs/benchmark-grouping.json --output results/grouping
python scripts/run_wsclean.py LOSSLESS.ms --output results/reference-fit
python scripts/run_wsclean.py LOSSLESS.ms --output results/reference \
  --beam-from results/reference-fit/image-image.fits
python scripts/validate_grouping.py --matrix configs/benchmark-grouping.json \
  --benchmarks results/grouping --reference results/reference \
  --beam-from results/reference-fit/image-image.fits \
  --output results/images --reports reports
```

Create the imaging/report parent directories beforehand. Each output must be new.
`validate_grouping.py` may run while benchmark entries finish; it only images an
entry after source/export verification passes and the recorded configuration
matches the matrix. Use `--experiments` to select entries, or `--threads 1` with
`run_wsclean.py` for the separate single-thread controls.
