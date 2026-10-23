# Opt-in destructive flag policy

`burn_flag` defaults to `false`. Enable it in an encode configuration with
`"burn_flag": true`, or use either CLI spelling:

```bash
visarchive encode INPUT.ms archive.h5 --burn-flag
visarchive encode INPUT.ms archive.h5 --config experiment.json --burn_flag
```

The switch overrides the configuration to `true`; omitting it leaves the
configuration unchanged. Non-Boolean configuration values fail validation.
It is supported with lossless, fixed-bit and adaptive visibility modes.

## Exact scope

The original **MAIN** `FLAG`, `FLAG_CATEGORY` and `FLAG_ROW` values are not stored.
Their descriptors, keywords, shapes and defined/undefined cell presence remain.
Uniform defined flag columns require only a shape declaration in the manifest,
with no HDF5 flag payload. Irregular/undefined cells use bounded shape-only JSON
chunks. Defined empty cells remain distinct from undefined cells.

MS export creates the original columns: every defined flag cell is all-false.
This keeps the Measurement Set structure valid; it does **not** regenerate flags.
The output carries the `VARQUANMS_FLAG_POLICY` table keyword, and encode/decode/
verify reports carry the destructive policy and warning. The archive records the
policy in configuration, column encodings and provenance. No source MS is edited.

Weights, sigmas, UVW, other metadata and subtable flags remain exactly preserved.
This switch does not discard any visibility values. Quantization still uses the
source `FLAG`/`FLAG_ROW` to exclude flagged samples from scale statistics and
preserve their values bitwise. Error metrics still use the source's unflagged
finite domain. Consequently, toggling only this switch does not change decoded
visibility values or encode-time error statistics. This is not a secure-erasure
feature: exception positions and aggregate statistics may reveal flag-related
information, although the original flag masks are not archived.

Even `mode: "lossless"` with `burn_flag: true` is **not an exact MS round trip**:
the mode describes visibility storage, not recovery of discarded flags.
Verification checks flag shapes/presence and all-false exported values, and
labels the original flag values intentionally discarded rather than bitwise
verified. All other declared checks still apply. Original flags cannot be
recovered from this profile; reflagging requires a separate workflow and may
produce different decisions, particularly for manual flags or external context.
Reflag before scientific use when the source contained flags.

## Size expectations and validation

The initial 59 MHz MS has zero true values in all three MAIN flag columns.
Its 95.491 MB of logical flags must **not** be described as 95.491 MB of archive
savings: the default archive already deduplicates/compresses those constants.
Compare measured whole-file sizes; manifest and HDF5 overhead count too.

Tests cover explicit CLI/config opt-in, strict Boolean validation, zero and
nonzero flags, irregular shapes, undefined/empty category cells, no flag payload
for regular arrays, exact visibility/protected-value preservation, unchanged
weights/subtable flags/source files, policy/encoding mismatch rejection and
rejection of nonzero flags introduced into a reconstructed MS. Existing
preserve-by-default round-trip tests remain mandatory.

## Measured full-MS check (2026-10-23)

All 51 unit/integration tests passed. The paired matrix
[`configs/benchmark-burn-flag.json`](../configs/benchmark-burn-flag.json) also
encoded, exported and verified the complete `20260904_190008_59MHz.ms` in
lossless-visibility mode under both policies:

| Policy | Whole archive bytes | MAIN flag dataset bytes |
| --- | ---: | ---: |
| Preserve (default) | 331,220,418 | 55 |
| `burn_flag: true` | 331,207,410 | 0 |

Measured whole-file saving: **13,008 bytes (0.00393%)**. Only 55 bytes were flag
dataset payloads; the rest is the net change in HDF5/manifest overhead. Weights,
sigmas and subtable flags still occupy their normal categories. This result
does not predict savings for nonconstant flags or other archive configurations.

Both exported MSs reopened with python-casacore and passed policy-aware
verification. All 65 source files retained the SHA-256 values recorded before
this feature was implemented. Full machine-readable results, software hashes,
byte categories and checks are in
[`reports/burn-flag-validation.json`](../reports/burn-flag-validation.json).
Raw generated archives/MSs/logs are under ignored `results/burn-flag-20261023/`.
No additional CASA imaging/calibration workflow was run for this change;
previous WSClean comparisons do not validate discarding nonzero flags.
