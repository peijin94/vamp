# Dimensionless distribution-referenced error budgets

The `distribution-v1` encoder policy replaces fixed absolute tolerances with
four configurable dimensionless coefficients: `epsilon_G`, `epsilon_S`,
`epsilon_D`, and `epsilon_R`. A null coefficient disables that constraint;
zero requests exact reconstruction. At least one must be enabled. The policy
rejects mixtures with legacy absolute, relative, per-product or noise budgets.
Historical configurations and archives remain readable; the lossless default
is unchanged. [Example configuration](../configs/calibrated-normalized.json).

For each block, using unflagged finite original values, define

$$S_b=\sqrt{\frac{\operatorname{Var}(\Re V)+\operatorname{Var}(\Im V)}2},$$

$$d_R=\operatorname{median}|\Re V-\operatorname{median}(\Re V)|,\quad
d_I=\operatorname{median}|\Im V-\operatorname{median}(\Im V)|,\quad
D_b=\sqrt{(d_R^2+d_I^2)/2},$$

$$R_b=\sqrt{\operatorname{mean}|V|^2/2},\qquad
T_b=\min(\epsilon_G G_{p,a},\epsilon_S S_b,\epsilon_D D_b,\epsilon_R R_b).$$

The minimum includes only enabled terms. Variance uses component means,
population normalization (`ddof=0`) and float64 arithmetic. MAD uses exact
float64 component medians, with no Gaussian conversion. RMS includes the
coherent centre. These references are fixed before candidate selection and
outlier exceptions; they are signal-distribution statistics, not thermal noise.
The codec still rounds its stored centre/range to float32 before evaluating
the actual decoded values.

Acceptance requires component MSE ≤T². Default dimensionless guards are maximum
complex error ≤10T, complex bias ≤T and maximum baseline component MSE ≤2T².
An optional p99 guard is `p99_factor*T`. These remain pooled real/imaginary
error limits; separate component MSEs are measured but not individually bounded.
Every bit-width, range and gamma candidate is measured after reconstruction.
Failed candidates use exact float32 fallback. There is no fixed absolute floor.

Zero MAD does not imply constant data: a mostly-zero block with nonzero outliers
can have zero MAD. If that reference is enabled, only exact reconstructions can
pass; exceptions or compressed float32 may be chosen. Empty/all-flagged blocks
have zero error budgets, with values handled according to the selected flag
preservation policy. Flags and nonfinite values are excluded from statistics;
genuine unflagged zeros remain included.

## Global sampling

`G` is the pooled within-component std on a reproducible sample of the complete
input, separately for each visibility column, physical configuration,
correlation product and auto/cross class. Four products in one configuration
therefore have eight global std references. The configuration identity contains
DATA_DESC_ID, FIELD_ID, ARRAY_ID, PROCESSOR_ID, STATE_ID and feed IDs, where
present. Time, scan and observation IDs do not split these references; original
MS identities and reconstruction partitions are still retained exactly.

The sampler reads bounded row batches. A SHA-256-derived stream key includes
the configured seed, configuration, column, correlation name and baseline class.
SplitMix64 maps original `row*nchannel+channel` indices into uint64 priorities;
priorities below `floor(fraction*2^64)` are selected. Fraction 1 selects all.
This is deterministic Bernoulli sampling, so the actual fraction can differ
from the requested one. It is independent of quantization blocks, permutations
and encoder worker count. It is not decoder randomness.

If fewer than `minimum_samples` were selected, the sampler uses the lowest
priorities, up to the number of eligible samples. The default minimum is 128.
No absolute error tolerance is introduced by this count safeguard. Sampling
fractions 2%, 3% and 5% with the same seed are nested except for minimum-count
fallback. Sample counts, means, separate component stds, pooled std, method and
sample SHA-256 are recorded. Means/M2 are merged in float64; batch-size changes
can change final rounding, while sample identities remain the same.

For partition bundles, references are measured once from the complete source
and reused in every time partition. The full-source export validator resamples
and compares the references. The standalone decoder validates reference
provenance and needs no source access or sampling to reconstruct visibilities.

The 15-slot calibrated source has separate OBSERVATION_IDs for its integrations
but otherwise one configuration. Pooling those slots gives the requested eight
references. A fixed 2% sample, nested 3%/5% samples and an independent 2% seed
were checked; the upcoming 16%-payload experiment reports their sensitivity.

## Configuration and scope

`normalized_bounds.reference` is `distribution-v1`. Alongside the four epsilon
values, it accepts `sampling_fraction`, `sampling_seed`, `minimum_samples`,
`max_error_factor`, `bias_factor`, `baseline_mse_factor` and optional
`p99_factor`. Unknown fields and nonfinite/negative coefficients fail. Active
epsilons may be zero; sampling fractions must be in (0,1].

One overall multiplier can scale a recorded coefficient recipe during a byte
search. The recipe ratios are part of the scientific policy and must also be
reported. Reference sampling percentage alone does not guarantee statistical
accuracy. Candidate comparisons use actual serialized bytes, and final files
must meet the physical payload ceiling before imaging.

The implementation adds encoder policy/provenance and keeps the packed codec
and decoder arithmetic unchanged. Validation covers gain scaling, zero MAD,
native/portable agreement, metadata and source-bit preservation, sampled class
separation, minimum-count fallback, sampling reproducibility, partition sharing
and malformed reference rejection. The local suite has 110 passing tests.

For the many-small-block benchmark, `timeslot_archive.py encode` also accepts
`--partition-workers 1..8`. Values above one require `workers=1` inside each
partition. Global references are still sampled once; independent contexts are
processed separately and assembled in original order only after every source
check passes. Runtime parallelism is recorded in the bundle's `execution`
field and does not change the scientific configuration. A fixture compares
all stored arrays and both exported MSs bitwise between serial and parallel
execution. This option changes throughput, not quantization criteria.
