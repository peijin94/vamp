# varquanMS

Adaptive 2/4/6/8-bit visibility compression for CASA Measurement Sets, with exact
float32 fallback and reconstruction to a new MS. Experiments group samples by
baseline geometry, antenna, amplitude and frequency segment.

The research prototype preserves original UVW, metadata, flags, weights and
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
The current quantized profile supports one integration/configuration with
channels inside complex64 cells; unsupported layouts fail explicitly.

```bash
python -m pip install -e '.[validation]'
python -m unittest discover -s tests -v
visarchive benchmark INPUT.ms --matrix configs/benchmark-grouping.json --output results/my-run
```

## Documentation

- [Archive and codec format](docs/format.md)
- [Development notes and experiment results](devnotes/README.md)
- [Initial MS inspection](reports/initial-inspection.txt)

Public repository: [peijin94/varquanMS](https://github.com/peijin94/varquanMS).
Observation payloads, reconstructed MSs and generated FITS files are excluded
from Git.
