# Paper update: calibrated I/V, performance and storage

26 October 2026. This revision incorporates the completed measurements into
[the paper](../docs/paper/varquantMS.md) and its [PDF](../docs/paper/varquantMS.pdf).
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
