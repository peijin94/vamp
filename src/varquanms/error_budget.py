"""Encoder-only error budgets; no decoder format change.

Distribution-v1 uses sampled global std G, block std S, raw component MAD D
and uncentered component RMS R. These signal references are not thermal noise.
It excludes legacy absolute/relative budgets, which remain supported separately.
Zero-reference blocks require exact values for each enabled constraint; there
is no implicit denominator floor or division by individual near-zero values.
"""
import math

import numpy as np


RELATIVE_FIELDS = {
    'max_component_rmse': ('max_component_mse', 2),
    'max_baseline_component_rmse': ('max_baseline_component_mse', 2),
    'max_complex_error': ('max_complex_error', 1),
    'max_complex_bias': ('max_complex_bias', 1),
    'max_p99_error': ('max_p99_error', 1),
}

EPSILONS = ('epsilon_G', 'epsilon_S', 'epsilon_D', 'epsilon_R')
NORMALIZED_DEFAULTS = dict(reference='distribution-v1', epsilon_G=None, epsilon_S=None,
    epsilon_D=None, epsilon_R=None, sampling_fraction=.02, sampling_seed=0,
    minimum_samples=128, max_error_factor=10., bias_factor=1., baseline_mse_factor=2.,
    p99_factor=None)


def normalized_policy(**overrides):
    result = dict(NORMALIZED_DEFAULTS, **overrides)
    validate_normalized_bounds(result)
    return result


def validate_normalized_bounds(spec):
    if not isinstance(spec, dict) or set(spec)-set(NORMALIZED_DEFAULTS) or spec.get('reference') != 'distribution-v1':
        raise ValueError('Unknown normalized error policy/options')
    values = dict(NORMALIZED_DEFAULTS, **spec)
    if all(values[k] is None for k in EPSILONS):
        raise ValueError('At least one epsilon_G/S/D/R must be enabled')
    for key in (*EPSILONS, 'max_error_factor', 'bias_factor', 'baseline_mse_factor', 'p99_factor'):
        value = values[key]
        if value is None and key in (*EPSILONS, 'p99_factor'):
            continue
        if type(value) not in (float, int) or not math.isfinite(value) or value < 0:
            raise ValueError('Normalized error coefficients must be finite and nonnegative: '+key)
    fraction = values['sampling_fraction']
    if type(fraction) not in (float, int) or not math.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError('sampling_fraction must be in (0,1]')
    for name, low, high in [('sampling_seed', 0, 2**64-1), ('minimum_samples', 2, 65536)]:
        if type(values[name]) is not int or not low <= values[name] <= high:
            raise ValueError('Invalid '+name)


def block_references(original, flags):
    """Unflagged finite source values; ddof=0, float64, raw component MAD.

    These references are fixed before candidate/outlier selection. The MAD
    centre is the exact float64 component median, not the rounded codec centre.
    """
    a = np.asarray(original)
    if a.shape != flags.shape:
        raise ValueError('Reference flags shape mismatch')
    valid = ~flags & np.isfinite(a)
    re, im = a.real[valid].astype(np.float64), a.imag[valid].astype(np.float64)
    if not len(re):
        return dict(valid_samples=0, S=0., D=0., R=0., median_real=0., median_imag=0.)
    mr, mi = float(np.median(re)), float(np.median(im))
    dr, di = np.median(np.abs(re-mr)), np.median(np.abs(im-mi))
    return dict(valid_samples=len(re), S=float(np.sqrt((np.var(re)+np.var(im))/2)),
                D=float(math.hypot(dr,di)/math.sqrt(2)),
                R=float(np.sqrt(np.mean(re*re+im*im)/2)), median_real=mr, median_imag=mi)


def normalized_limits(spec, references, global_std=None):
    """Resolve only data-referenced tolerances; zero spread means zero budget."""
    policy = dict(NORMALIZED_DEFAULTS, **spec)
    scales = dict(G=global_std, **{k:references[k] for k in ('S','D','R')})
    if references['valid_samples'] == 0:
        scales['G'] = 0.
    candidates = []
    for key in EPSILONS:
        epsilon = policy[key]
        if epsilon is None:
            continue
        scale = scales[key[-1]]
        if scale is None or not math.isfinite(scale) or scale < 0:
            raise ValueError('Missing or invalid distribution reference '+key[-1])
        candidates.append(epsilon*scale)
    if not candidates:
        raise ValueError('No enabled epsilon')
    threshold = min(candidates)
    try:
        mse = threshold**2
        result = dict(max_component_mse=mse, max_complex_error=policy['max_error_factor']*threshold,
            max_complex_bias=policy['bias_factor']*threshold,
            max_baseline_component_mse=policy['baseline_mse_factor']*mse,
            max_p99_error=None if policy['p99_factor'] is None else policy['p99_factor']*threshold)
    except OverflowError as exc:
        raise ValueError('Normalized error budget overflows') from exc
    if any(v is not None and not math.isfinite(v) for v in result.values()):
        raise ValueError('Normalized error budget overflows')
    return result


def configured_bounds(original, flags, config, absolute, global_std=None):
    if config.normalized_bounds is not None:
        return normalized_limits(config.normalized_bounds, block_references(original, flags), global_std)
    return resolve_bounds(original, flags, absolute, config.relative_bounds)


def validate_relative_bounds(spec):
    if not isinstance(spec, dict) or spec.get('reference') != 'block-rms-v1':
        raise ValueError('relative_bounds requires reference=block-rms-v1')
    if set(spec) - {'reference', *RELATIVE_FIELDS} or len(spec) == 1:
        raise ValueError('Unknown or empty relative error limits')
    for name, value in spec.items():
        if name != 'reference' and (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
            raise ValueError('Relative limits must be finite nonnegative fractions')


def component_signal_power(original, flags):
    a = np.asarray(original)
    if a.shape != flags.shape:
        raise ValueError('Relative error reference requires matching flags')
    valid = ~flags & np.isfinite(a)
    re = a.real[valid].astype(np.float64)
    im = a.imag[valid].astype(np.float64)
    return float(np.mean((re*re + im*im)/2)) if len(re) else 0.


def resolve_bounds(original, flags, absolute, relative):
    """Convert configured fractions into strict measured block limits."""
    result = dict(absolute)
    if relative is None:
        return result
    scale = math.sqrt(component_signal_power(original, flags))
    limits = dict(relative)
    if 'max_baseline_component_rmse' not in limits and 'max_component_rmse' in limits:
        limits['max_baseline_component_rmse'] = limits['max_component_rmse']
    for name, (target, exponent) in RELATIVE_FIELDS.items():
        if name not in limits:
            continue
        try:
            limit = (limits[name] * scale)**exponent
        except OverflowError as exc:
            raise ValueError('Relative error limit overflows') from exc
        if not math.isfinite(limit):
            raise ValueError('Relative error limit overflows')
        result[target] = limit if result[target] is None else min(result[target], limit)
    return result
