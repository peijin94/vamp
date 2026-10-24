# varquanMS

Adaptive 3/4/6/8/12/16-bit visibility compression for CASA Measurement Sets, with exact
float32 fallback and reconstruction to a new MS. Experiments group samples by
baseline geometry, antenna, amplitude and frequency segment.

By default, the research prototype preserves original UVW, metadata, flags, weights and
protected sample bits exactly. Quantized visibilities are lossy. Default encoding
is lossless; science-ready tolerances have not been established.

## Quick start

Python 3.10+ is required.

```bash
python -m pip install -e .
visarchive inspect INPUT.ms --report inspection.json
visarchive encode INPUT.ms archive.h5 --config configs/lossless.json
visarchive decode archive.h5 RESTORED.ms
visarchive verify archive.h5 --source INPUT.ms --restored RESTORED.ms --report verification.json
```

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
visarchive benchmark INPUT.ms --matrix configs/benchmark-grouping.json --output results/my-run
```

## Documentation

- Methods paper, diagrams, imaging results and timing: [PDF](docs/paper/varquantMS.pdf) · [Markdown](docs/paper/varquantMS.md)
- [Archive and codec format](docs/format.md)
- [Development notes and experiment results](devnotes/README.md)
- [Initial MS inspection](reports/initial-inspection.txt)

Public repository: [peijin94/varquanMS](https://github.com/peijin94/varquanMS).
Observation payloads, reconstructed MSs and generated FITS files are excluded
from Git.
