# Dysco comparison in the paper: additional epsilon measurements

The existing matched comparison is documented in [compare-dysco.md](compare-dysco.md)
and [issue #3](https://github.com/peijin94/varquanMS/issues/3). Paper section 3.7
adds its four figures and the complete comparison table. The rate figure now
includes independently measured `epsilon_D = 0.09` and `0.3` points alongside
the existing 0.01, 0.04 and 0.18 cases. It connects measured settings within
each method and reverses the log error axes, so lower errors appear higher.
The other three figures are copies of the matched-rate issue figures.

## Experiment and results

Both new cases use the same 15-integration calibrated OVRO-LWA input,
128-baseline × 48-channel amplitude-sorted blocks, independent product ordering,
power-midrise codec, gamma policy, precision ladder and error guards as the
0.18 comparison. Only epsilon and execution parallelism changed. The source
contains 178,928,640 complex64 values (1,431,429,120 logical bytes). The
calibrated-only policy retains flags, UVW, metadata and weights; flagged
visibility values are deliberately purged. Row weights were bitwise exact for
this input. This is not a claim that all four-bit weight reconstructions are
exact for other inputs.

| epsilon D | Visibility storage, MB | Complete archive, MB | Paired-I visibility error | I / residual, first / repeat | V / residual, first / repeat |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.09 | 188.904 | 282.826 | 4.658% | 2.800% / 2.797% | 5.764% / 5.764% |
| 0.3 | 124.769 | 218.750 | 13.353% | 20.098% / 20.059% | 12.964% / 12.964% |

MB means 10^6 bytes. Visibility storage comprises the block payload plus
permutations and visibility indexes: 182,365,484 + 6,185,493 + 352,749 bytes
at 0.09, and 118,230,585 + 6,185,493 + 352,749 bytes at 0.3. These amounts
are 13.197% and 8.716% of the original complex64 logical bytes, respectively.
The complete archive includes additional metadata, coordinates, quality
statistics, flags, weights and container overhead. The Dysco DATA manager
and complete output MS use different representations, so the storage
boundaries are described separately in the paper and comparison report.

The paired-I figure uses RMS error relative to the unflagged cross-baseline
`(XX+YY)/2` visibilities. Each image number is instead
`100 × RMS(candidate restored − reference restored) / RMS(reference residual)`
on the same 2,245,726 pixels above 10° for all integrations. Reference residual
RMS is 1.031244811 Jy/beam for I and 0.384439444 Jy/beam for V. It is not a
flux-relative percentage or a demonstrated thermal-noise level. All eight new
CLEANs reached their final threshold. The corresponding I dirty-image
differences were 1.865% for 0.09 and 5.332% for 0.3, while the first
restored-I differences were 2.800% and 20.098%; nonlinear CLEAN can amplify or
rearrange differences, and these results do not establish causality.

The first image differences are slightly larger than their repeats. A separate
uncompressed-reference repeat differed by 0.421% for I and 0.0111% for V;
this is a repeat control, not a compression-error estimate. Two repeated
images on one observation do not establish general scientific fidelity.
Calibration solutions, spectral features, leakage separation and temporal
behavior remain untested by these runs.

## Integrity and reproducibility

The new [machine-readable report](../reports/dysco-paper-extra.json) records
the actual per-product and binned visibility diagnostics, sample precision
histograms, storage categories, image results and artifact SHA-256 values.
Its collection script checks the original input snapshot, mask checksum,
scientific configuration against the original 0.18 case, complete archive
inventories, expected sample count, restored-MS hashes and reference/imaging
identities. UVW, FLAG, FLAG_ROW, WEIGHT and SIGMA were checked bitwise after
each export, and the quantized column matched the archive decoder. The source
remained unchanged. All eight WSClean inputs were private copies and used the
same calibrated template and reference beams as the established comparison.

The configurations are
[`calibrated-epsd-009-parallel.json`](../configs/calibrated-epsd-009-parallel.json)
and [`calibrated-epsd-030-parallel.json`](../configs/calibrated-epsd-030-parallel.json).
With the recorded local input, images and archives under
`results/dysco-paper-extra/`, recreate the report at a *new* report path,
then render and update the paper:

```bash
PYTHONPATH=src:scripts python scripts/build_dysco_paper_comparison.py collect \
  --extra-report results/dysco-paper-extra/recollected.json
PYTHONPATH=src:scripts python scripts/build_dysco_paper_comparison.py render
PYTHONPATH=src:scripts python scripts/update_dysco_paper.py
python scripts/compile_paper.py
```

The normal renderer uses the committed reports and issue figures; it does not
need the local MS or FITS outputs. Collection performs the local data audit but
does not rerun encoding or imaging. The archived plan used four partition workers
and six workers per partition; runs shared the host, so their elapsed times
are not isolated performance benchmarks.
