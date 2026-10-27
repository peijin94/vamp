# Three-epsilon Stokes V follow-up and paper revisions

Measured 27 October 2026. This follows the request to replace the old combined
I/V image figures with V counterparts of the three-epsilon I plots, use reference
residual RMS, shorten colourbar labels, hide repeated axes, and replace the speed
figure with a timing table.

Six new Stokes V CLEANs use the existing verified exports for epsilon D 0.01,
0.04 and 0.18 (first and repeat for each). They use the established pinned
WSClean container, 2048² pixels, 4 arcmin, CORRECTED_DATA, eight threads, multiscale,
scale bias .8, niter 150000, mgain .95, auto-threshold .5, auto-mask 3, local RMS,
horizon mask 10deg, Briggs 0, mem 40, reorder and no model update. Each runs on a
private copy with the V reference restoring beam. No archive was re-encoded and
no visibility threshold was changed.

All six reach the final CLEAN threshold. Their input MS hashes match those used
for the original I tests and remain unchanged;
reference reuse validates the original source hashes and imaging parameters.
Recorded reference/candidate FITS hashes and the fixed mask are rechecked by the
report builder. The comparison mask contains 2,245,726 pixels above 10° throughout
the observation. The V reference residual RMS is 0.3844394444764314 Jy/beam.

| epsilon D | Difference RMS mJy/beam, first / repeat | Difference/reference residual RMS %, first / repeat |
| --- | ---: | ---: |
| 0.01 | 4.1801 / 4.1802 | 1.0873249 / 1.0873441 |
| 0.04 | 8.0278 / 8.0278 | 2.0881795 / 2.0881794 |
| 0.18 | 38.6601 / 38.6601 | 10.0562314 / 10.0562317 |

These are new V realizations, distinct from the prior performance-study pair.
The original I results and the old I/V performance report remain unchanged.
The comparisons quantify agreement with the calibrated reference, not absolute
polarization accuracy or a universal acceptance threshold.
[Measured report](../reports/calibrated-epsd-v.json).

The paper retains the I full/zoom comparisons and replaces the old combined
panels with matching V full/zoom comparisons. Each is 2×3: recovered images above
signed differences, one epsilon per column. Difference pixels are
`100 × (S_quant − S_ref) / RMS(reference residual)`, for S=I or V, with a fixed
−30% to +30% display range and extended colourbars. The I reference denominator
is 1.031244811054661 Jy/beam. Denominators are uncentred RMS over the full fixed
mask, shared across epsilon values and zooms; the earlier std display convention
is superseded. V and differences use bwr; I intensity uses viridis. Colourbar
labels are compact mathematical expressions. Top-row x axes and y axes in
columns 2–3 are hidden. Each polarization uses its own reference-peak crop.

## Timing table replacing Figure 10

The table decomposes the existing **run 2**, which is the median of the three
unprofiled optimized full-MS timings: **72.68508086586371 s**. It does not use the
new, overlapping CLEAN runs as performance data.

| Non-overlapping phase | Seconds |
| --- | ---: |
| Global reference sampling | 6.140628 |
| Remaining bundle work: discovery, parallel inspection/encoding/source checks, finalization | 46.964300 |
| Validate input bundle for compaction | 6.116741 |
| Copy/share arrays into final archive | 4.848464 |
| Validate final archive | 5.308617 |
| Other compaction/finalization overhead | 3.051499 |
| Process startup/shutdown and remaining orchestration | 0.254831 |
| **Total** | **72.685081** |

The last two rows and the remaining-bundle row are differences of measured wall
intervals, not independent function profiles. Within the parallel stage, summed
work across the 15 partitions is 23.466356 s for inspection, 223.877001 s for
encoding/writing, and 73.777756 s for source/integrity checks. These overlap and
cannot be added to the wall-time breakdown. Finer additive encoding-stage times
were not measured. The [timing audit](../reports/calibrated-runtime-breakdown.json)
records input hashes, exact values and derivations.

## Reproduction and checks

```bash
export PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python scripts/run_error_limits.py --root results/calibrated-epsd-dense \
  --plan configs/calibrated-epsd-study.json --source data/working/lwa-15slot-flat.ms \
  --stokes V --image-only --jobs 2
python scripts/run_error_limits.py --root results/calibrated-epsd-dense \
  --plan configs/calibrated-epsd-study.json --source data/working/lwa-15slot-flat.ms \
  --stokes V --repeat-only --jobs 2
PYTHONPATH=src:scripts python scripts/report_paper_followup.py
PYTHONPATH=src:scripts python scripts/update_calibrated_paper.py
python scripts/compile_paper.py
```

The imaging commands require the previously verified exports and original config
files. Existing completed stages are validated/reused; incomplete outputs are
not accepted. Report collection checks status, settings, container, source and
FITS hashes, final CLEAN termination and residual normalization. The additive
phase table is checked against the original recorded median; parallel totals
are obtained from its verified partition records. Figures were regenerated from
the audited FITS and checked visually. PDF compilation, local links, page bounds
and Python syntax were checked; codec implementation and prior measurement
reports were not changed.
