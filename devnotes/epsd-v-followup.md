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
columns 2–3 are hidden. I uses its reference-peak crop; the V crop was subsequently
moved to high elevation as described below.

## Initial V display limits and high-elevation zoom (superseded below)

The later 27 October request fixes V intensity limits to **−30 to +30 Jy/beam**
in both the full-field and zoom figures. The symmetric-log normalization is
retained, with its linear threshold at 0.3 Jy/beam; both intensity and difference
panels use bwr. Differences remain linear at **−30% to +30%**, normalized by the
unchanged full-mask V reference residual RMS.

The V zoom is centred on the largest absolute reference pixel whose elevation
exceeds 60° at all 15 original integration times. The complete 128×128 crop is
also checked to exceed 45° at every epoch. Selection uses only the V reference
and the original mask's retained time, array location and IERS information.
It does not use reconstructed images or their errors.

- Centre: zero-based `(y, x) = (1204, 1231)`; reference V = 2.015388 Jy/beam.
- Crop: `y[1140:1268], x[1167:1295]`, with exclusive stop indices.
- Minimum centre elevation: **68.967343°** across the integration times.
- Minimum elevation over every crop pixel and time: **62.633077°**.

[Selection audit](../reports/calibrated-v-zoom.json) records the bounds,
reference/mask/IERS hashes and geometric inputs. Regenerate it with
`PYTHONPATH=src:scripts python scripts/select_v_zoom.py` using the retained
reference FITS, mask provenance and local IERS table. The plotting script checks
the selection's reference and mask hashes before using its fixed crop.

These are display changes. Scientific error measurements, the fixed above-10°
comparison domain, I figures, archives and source MS remain unchanged. Checks
cover every panel's colour limits, colormap and axis visibility, zoom dimensions
and bounds, recorded FITS errors, figure appearance and PDF compilation.

## Source-centred V zoom and ±120 Jy/beam limits

The next request replaces the diffuse high-elevation field with a compact V
feature at a radio source and widens V intensity limits to **−120 to +120
Jy/beam**. Both full-field and zoomed V images initially used those limits, bwr
and the existing symmetric-log linear threshold of 0.3 Jy/beam; the subsequent
linear-scale revision below supersedes that normalization. Differences remain
linear at ±30%, with the same full-mask reference residual RMS.

Searching reference V peaks above 45° at every integration selects `(y,x) =
(947,449)`: **+6.091340 Jy/beam**, coincident with the **488.629883 Jy/beam** I
peak. Its ICRS position is RA 83.692929°, Dec +22.017468°, consistent within one
restoring beam with [Taurus A / the Crab Nebula](https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Crab+Nebula).
The catalogue association is made after selecting the reference peak; candidate
quantized images do not influence selection. The two stronger absolute V peaks
found in the full mask are at only 19.95° and 13.71° elevation and are excluded.

The reference repeat gives +6.091304 Jy/beam at the same pixel; the V reference
residual there is +0.177968 Jy/beam. In an 8–16-pixel annulus, median reference V
is 0.003094 Jy/beam and scaled MAD is 0.628021 Jy/beam. This local fluctuation
statistic is not a thermal-noise estimate. The compact feature is measured V at
a known radio-source position; this test does not establish whether it is
intrinsic circular polarization or instrumental leakage.

The tighter **64×64** crop uses `y[915:979], x[417:481]` (zero-based, exclusive
stops). Its centre stays above **50.125659°**, and every crop pixel remains above
**46.710021°** at all 15 integration times. This replaces the previous 128×128
crop above 60°, which mostly showed diffuse fluctuations. The existing
[selection audit](../reports/calibrated-v-zoom.json), now schema version 2,
includes the catalogue position, angular separation, beam, supporting I/repeat
V/residual FITS hashes and values. The generator checks the catalogue match and
image grids. Scratch candidate maps remain ignored under `plt_tmp/`.

The full-mask science metrics and I figures are unchanged. V intensity and
difference bounds, bwr colormaps, hidden axes and new crop bounds are checked
while regenerating; dashed circles locate the reference peak in the recovered
V zoom panels. The paper and PDF captions describe the source association
and its polarization limitation.

## Linear V intensity scale (initial limits)

The subsequent request removes the symmetric-log stretch from V intensity.
Both full-field and Taurus A zoom panels now use a strictly linear mapping from
**−120 to +120 Jy/beam**, with colourbar ticks at −120, −60, 0, 60 and 120.
The bwr colormap, source crop, reference-peak circles and linear ±30% difference
panels are retained. I intensity keeps its existing display scale.

The figures and paper/PDF captions are regenerated from the existing FITS.
Checks verify the V normalization is linear, limits and ticks match the request,
and crop/axes remain correct. No imaging, compression or scientific metric is
recomputed; all measurement reports and I figures remain unchanged.

The subsequent display adjustment narrows the linear V intensity range to
**−30 to +30 Jy/beam** in both figures, with evenly spaced ticks at −30, −15,
0, 15 and 30. The zoom, bwr colormap, ±30% difference scale and measured data
remain unchanged.

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
