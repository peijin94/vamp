# Paper update: calibrated I/V, performance and storage

26 October 2026. This revision incorporates the completed measurements into
[the paper](../docs/paper/vamp.md) and its [PDF](../docs/paper/vamp.pdf).
No new quantization, CLEAN or timing experiment was run for this manuscript update.

The initial three-epsilon I results remain separate from the later I/V tests at
εD=0.18. Figures 7–8 now show I/V reference, reconstruction and signed difference,
with full-field and reference-peak views. I flux uses viridis; V and differences
use bwr. Each difference is divided by its matching scalar reference residual
RMS. Figure 9 retains its requested residual-normalized image-error axis.
Figure 10 reports the three repeated complete-workflow times, with median/range
and the unmet 15-second target. Implementation details remain in
[the performance note](calibrated-performance.md).

Size accounting now includes the original physical MS, logical CORRECTED_DATA,
the same-policy lossless shared-file control, all three quantized files and the
new optimized file. The lossless control's SHA-256 was reread and verified, and
its archived configuration confirms calibrated-only preservation, retained flags
and four-bit weights. The measured 94.393 MB structural compaction saving is
separate from quantization. The optimized file is 242,769,176 bytes versus
985,235,551 bytes for this retained-data lossless control, a 4.058× reduction.
Dropping original DATA/MODEL_DATA and flagged values is explicitly disclosed.

The physical categories sum to the complete file. The precision table counts all
178,928,640 complex samples and 29,224 blocks. Exact zero blocks contain 41,730,048
samples but occupy only 1,727,348 serialized bytes. Precompression component,
code, header and exception counts come from the original optimized run-1 bundle
manifest. Factor bytes are a subset of headers. Compressed codes and exceptions
share bodies, so final physical bytes are not assigned to those subcategories.
The audit is recorded in [calibrated-paper-storage.json](../reports/calibrated-paper-storage.json).

Validation of this revision:

- Regenerated figures from frozen reports and actual recorded FITS; I/V image
  difference RMS agrees with the saved measurements before plotting.
- Checked input report hashes and exact physical-category/sample/precision sums.
- Compiled with Pandoc 3.9 and Typst 0.15.0; reviewed the 16-page PDF's I/V panels,
  tables and timing figure and checked that text stays inside page boundaries.
- Retained the initial work note unchanged (SHA-256
  `22c079ca6a557d873c801a08c71a53a007606d5e315f41ff3c5843609fbf3b89`).
- Scientific measurement reports and codec implementation remain unchanged. The
  paper cites the performance study's 124 tests; this documentation update does
  not claim a new local codec test run.

`update_calibrated_paper.py` refreshes current tables and figures without running
experiments. `--skip-images` works from committed reports and figure files.
The historical `report_epsd_dense.py` now updates only its report and development
note, preventing it from replacing the current manuscript with the older I-only
text. SVG identifiers and vector-figure dates are deterministic.

## Grouping concept figure

The follow-up revision adds the existing `grouping` schematic to Section 2.1 as
Figure 1a and retains the candidate-selection flowchart as Figure 1b. Its caption
identifies the synthetic 2×4 example and the actual 128×48 calibrated-data blocks,
with stable whole-spectrum ordering and the inverse permutation. Later figure
numbers and measured results are unchanged. The rebuilt PDF has 17 pages; both
vector diagrams and captions were visually checked, local links resolve, and
text remains within page boundaries.

## Figures 5–6: recovered images and residual-std percentages

Both figures now use two rows and three columns, ordered by epsilon-D 0.01, 0.04
and 0.18. The top row shows each first-run recovered I image; the bottom row shows
`100 * (I_recovered - I_reference) / std(reference_residual)` in percent. The
population standard deviation (`ddof=0`) is 1.018715249293199 Jy/beam, evaluated
once over the full 2,245,726-pixel above-10° mask. The denominator and both display
scales are shared across cases and the full-field/128×128-pixel zoom figures.
I uses viridis; signed percentage differences use bwr.

The formula is now explicit in the paper. Existing tabulated error metrics and
Figure 9 retain their documented reference-residual RMS normalization. The
underlying imaging experiments and measurement reports are unchanged.

`image_comparisons()` in `build_epsd_figures.py` verifies the saved mask and image
hashes, matched grids/beams, reference residual standard deviation and original
image-difference RMS before drawing the panels. The generic historical plotting
routine can skip images so paper regeneration uses the new layout without
changing older development-note figures. Direct regeneration, Python syntax
checks, Markdown links and PDF compilation passed. Figures 5 and 6 were visually
checked on pages 11 and 12 of the rebuilt 18-page PDF; all PDF text remains within
page boundaries.

The subsequent display-range adjustment fixes the difference colour limits in
Figures 5–6 at −30% to +30%, with 10-percentage-point ticks and extended colourbar
ends. Values outside the range use endpoint colours; numerical differences,
normalization, measurements and top-row intensity scales are unchanged.
