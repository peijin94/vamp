# vamp

**Visibility-Adaptive Mixed Precision** is a research prototype for archiving
radio-interferometric Measurement Sets. It compresses complex visibilities in
locally scaled blocks and reconstructs a new MS for existing imaging tools.

## How it works

- Inspect the MS and group baselines and frequency channels within each
  configuration and correlation product. Ordering can use amplitude, baseline
  geometry or antenna identity; stored mappings restore the original order.
- Quantize each block with uniform or median-centered power-spaced levels.
  Adaptive mode tests 3/4/5/6/8/10/12/14/16 bits per real and imaginary
  component and selects the smallest encoded candidate that meets the configured
  error checks, with exact float32 fallback.
- By default, preserve original UVW in metres, metadata, flags, weights and
  flagged sample values. The default lossless mode preserves visibility values
  exactly; quantized modes are lossy. There is no time or channel averaging.

## Calibrated-data comparison

The following measurements use one calibrated OVRO-LWA observation with
15 integrations and 48 channels. They compare vamp with Dysco using the same
unflagged calibrated visibilities, retained flags and row weights, and WSClean
Stokes I and V settings. Both paths discard flagged visibility values; vamp
exports zeros and Dysco receives NaN markers. The vamp cases use amplitude-ordered
blocks of up to 128 baselines × 48 channels per product; blue point labels show
the dimensionless εD error setting relative to local median deviation.

![Measured vamp and Dysco rate–distortion comparison](docs/paper/figs/dysco-comparison-rate.png)

The horizontal axis is visibility storage: vamp block payloads, permutations
and indexes versus the Dysco DATA-manager file. The I/V axes show image-difference
RMS divided by the corresponding reference residual RMS; the third axis shows
relative error in the paired cross-baseline I visibilities. Error axes are
reversed, so lower error appears higher. Lines connect measured settings.

| Nearly matched visibility storage | Visibility MB | Complete output MB | I / residual RMS | V / residual RMS |
| --- | ---: | ---: | ---: | ---: |
| vamp, εD = 0.18 | 148.702 | 242.766 | 5.584% | 10.056% |
| Dysco RF, 3-bit, 1.5σ | 149.119 | 318.774 | 12.049% | 27.154% |

The I/V values are first CLEAN runs over the region above 10°; their denominator
is the reference residual RMS, not source flux or a verified thermal-noise level.
Complete sizes include different HDF5 archive and MS container costs. These
results describe this observation, not a general scientific-fidelity guarantee.
See the [paper](docs/paper/vamp.pdf) and [comparison notes](devnotes/compare-dysco.md)
for the full table, repeats, image differences and limitations.

## Quick start

Python 3.10+ is required.

```bash
python -m pip install -e .
vamp inspect INPUT.ms --report inspection.json
vamp encode INPUT.ms archive.h5 --config configs/lossless.json
vamp decode archive.h5 RESTORED.ms
vamp verify archive.h5 --source INPUT.ms --restored RESTORED.ms --report verification.json
```

Existing archives remain readable after the rename. The former `visarchive`
command and `varquanms` Python namespace remain compatibility aliases; new
integrations should use `vamp`.

Use `configs/engineering-adaptive.json` for an explicitly lossy experiment.
The [normalized calibrated preset](configs/calibrated-normalized.json) uses
dimensionless global-std, block-std, median-deviation and RMS error controls;
see [the policy definition](devnotes/normalized-budgets.md).
The [performance preset](configs/performance-adaptive.json) uses 24 CPU workers
and optional compiled kernels; see [measured timings and setup](devnotes/performance.md).
Optional `encode --burn-flag` discards MAIN flags and exports them as all-false;
no reflagging is performed. See [the destructive flag policy](devnotes/burn-flag.md).
Optional `encode --keep-calib-only` retains `CORRECTED_DATA`, removes raw/model
visibility columns, and zeros flagged values while preserving flags. Weights use
4-bit codebooks by default, with checked float32 fallback; see the
[calibrated profile and 15-slot test](devnotes/calibrated-imaging.md).
Export recreates `DATA` from decoded calibrated values for MS compatibility;
the discarded raw `DATA` is not recovered.
The current quantized profile supports one integration/configuration with
channels inside complex64 cells; the experimental partition-bundle script handles
contiguous integrations. Unsupported layouts fail explicitly.

```bash
python -m pip install -e '.[validation]'
python -m unittest discover -s tests -v
vamp benchmark INPUT.ms --matrix configs/benchmark-grouping.json --output results/my-run
```

## Documentation

- Methods paper, diagrams, imaging results and timing: [PDF](docs/paper/vamp.pdf) · [Markdown](docs/paper/vamp.md)
- [Archive and codec format](docs/format.md)
- [Development notes and experiment results](devnotes/README.md)
- [Initial MS inspection](reports/initial-inspection.txt)

Public repository: [peijin94/vamp](https://github.com/peijin94/vamp).
Observation payloads, reconstructed MSs and generated FITS files are excluded
from Git.
