"""Encoder-only block-relative error limits; no decoder format change.

Reference S = sqrt(mean(|V|**2)/2) over unflagged finite complex samples.
S includes a coherent center: it is signal RMS, never a thermal-noise estimate.
Relative and absolute limits intersect. Zero-reference blocks require exact
values for every configured relative metric; there is no implicit denominator
floor or division by individual near-zero visibilities.
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
