"""Strict, versioned archive experiment configuration. Defaults are lossless."""

from dataclasses import asdict, dataclass, field
import json
import math
from pathlib import Path

DEFAULT_BITS = (3, 4, 6, 8, 12, 16)
# Two-bit archives and explicit historical configurations remain supported.
SUPPORTED_BITS = (2, *DEFAULT_BITS)
CODEBOOK_IDS = {"uniform-symmetric-v1": 0, "uniform-midrise-v1": 1, "power-midrise-v1": 2}


@dataclass
class Config:
    mode: str = "lossless"
    columns: list = field(default_factory=lambda: ["DATA"])
    deduplicate: bool = True
    burn_flag: bool = False
    keep_calib_only: bool = False
    weight_bits: int = 4
    weight_max_relative_error: float = .01
    compression: str = "gzip"
    compression_level: int = 1
    io_chunk_rows: int = 1024
    baseline_group: int = 1
    channel_segment: int = 192
    ordering: str = "original"
    shared_order: bool = False
    bucket_reference: float | None = None
    orientation_bins: int = 24
    length_reference_m: float = 1.0
    antenna_endpoint: int = 1
    bits: list = field(default_factory=lambda: list(DEFAULT_BITS))
    codebook: str = "uniform-symmetric-v1"
    rounding: str = "nearest-even"
    seed: int = 0
    scale_quantiles: list = field(default_factory=lambda: [1.0, 0.995])
    power_center: str = "median"
    power_gammas: list = field(default_factory=lambda: [1.0, 0.5])
    power_gamma_policy: str = "search"
    power_ratio_reference: float = 4.0
    max_component_mse: float | None = None
    max_complex_error: float | None = None
    max_complex_bias: float | None = None
    max_baseline_component_mse: float | None = None
    max_p99_error: float | None = None
    relative_bounds: dict | None = None
    product_bounds: dict = field(default_factory=dict)
    noise: dict | None = None
    max_visibility_bytes: int = 2 * 1024**3
    workers: int = 1

    @classmethod
    def load(cls, path=None, **overrides):
        data = {} if path is None else json.loads(Path(path).read_text())
        if not isinstance(data, dict):
            raise ValueError("Configuration must be a JSON object")
        unknown = (set(data) | set(overrides)) - set(cls.__dataclass_fields__)
        if unknown:
            raise ValueError(f"Unknown configuration options: {', '.join(sorted(unknown))}")
        data.update(overrides)
        if data.get("keep_calib_only") and "columns" not in data:
            data["columns"] = ["CORRECTED_DATA"]
        result = cls(**data)
        result.validate()
        return result

    def validate(self):
        if self.mode not in ("lossless", "fixed", "adaptive"):
            raise ValueError("mode must be lossless, fixed or adaptive")
        if not isinstance(self.columns, list) or not self.columns or any(c not in ("DATA", "CORRECTED_DATA", "MODEL_DATA") for c in self.columns) or len(set(self.columns)) != len(self.columns):
            raise ValueError("columns must be a unique nonempty list of complex visibility columns")
        for name in ("deduplicate", "shared_order", "burn_flag", "keep_calib_only"):
            if type(getattr(self, name)) is not bool:
                raise ValueError(f"{name} must be boolean")
        if self.keep_calib_only and (self.burn_flag or self.columns != ["CORRECTED_DATA"]):
            raise ValueError("keep_calib_only requires columns=[CORRECTED_DATA] and preserved flags (burn_flag=false)")
        if type(self.weight_bits) is not int or self.weight_bits not in (4, 32):
            raise ValueError("weight_bits must be 4 or 32 (exact)")
        if type(self.weight_max_relative_error) not in (float, int) or not math.isfinite(self.weight_max_relative_error) or not 0 <= self.weight_max_relative_error < 1:
            raise ValueError("weight_max_relative_error must be finite in [0,1)")
        for name, lower, upper in (("io_chunk_rows", 1, 65536), ("baseline_group", 1, 4096), ("channel_segment", 1, 65536), ("compression_level", 0, 9), ("seed", 0, 2**64 - 1), ("max_visibility_bytes", 1, 64 * 1024**3)):
            value = getattr(self, name)
            if type(value) is not int or not lower <= value <= upper:
                raise ValueError(f"{name} must be an integer in [{lower}, {upper}]")
        if self.baseline_group * self.channel_segment > 262144:
            raise ValueError("Quantization blocks are limited to 262144 complex samples")
        if type(self.workers) is not int or not 1 <= self.workers <= 64:
            raise ValueError("workers must be an integer in [1,64]")
        if self.compression not in ("gzip", "none"):
            raise ValueError("compression must be gzip or none")
        if self.ordering not in ("original", "class-original", "length", "orientation", "length-orientation", "antenna", "amplitude", "buckets"):
            raise ValueError("Unsupported baseline ordering")
        if type(self.orientation_bins) is not int or not 1 <= self.orientation_bins <= 3600:
            raise ValueError("orientation_bins must be an integer in [1,3600]")
        if type(self.length_reference_m) not in (int, float) or not math.isfinite(self.length_reference_m) or self.length_reference_m <= 0:
            raise ValueError("length_reference_m must be positive and finite")
        if type(self.antenna_endpoint) is not int or self.antenna_endpoint not in (1, 2):
            raise ValueError("antenna_endpoint must be 1 or 2 (ordered MS antenna columns)")
        if self.codebook not in CODEBOOK_IDS:
            raise ValueError("Unsupported codebook")
        if self.power_center not in ("median", "zero"):
            raise ValueError("power_center must be median or zero (control)")
        if not isinstance(self.power_gammas, list) or not self.power_gammas or any(type(g) not in (int, float) or g not in (1., .5, .25) for g in self.power_gammas) or len(set(self.power_gammas)) != len(self.power_gammas):
            raise ValueError("power_gammas must be a unique nonempty subset of [1, 0.5, 0.25]")
        if self.power_gamma_policy not in ("search", "max-std-v1", "max-std-inlier-v1"):
            raise ValueError("Unsupported power_gamma_policy")
        if type(self.power_ratio_reference) not in (int, float) or not math.isfinite(self.power_ratio_reference) or self.power_ratio_reference <= 0:
            raise ValueError("power_ratio_reference must be positive and finite")
        if self.power_gamma_policy != "search" and self.codebook != "power-midrise-v1":
            raise ValueError("Ratio-based gamma selection requires power-midrise-v1")
        if self.codebook == "power-midrise-v1" and self.rounding != "nearest-even":
            raise ValueError("power-midrise-v1 currently supports nearest-even in reconstructed-value space only")
        if self.rounding not in ("nearest-even", "stochastic-pcg64-v1"):
            raise ValueError("Unsupported rounding; subtractive dithering is not implemented")
        if not isinstance(self.bits, list) or not self.bits or any(type(b) is not int or b not in SUPPORTED_BITS for b in self.bits) or len(set(self.bits)) != len(self.bits):
            raise ValueError("bits must be a unique nonempty list selected from 3/4/6/8/12/16 (legacy 2 also supported)")
        if self.mode == "fixed" and len(self.bits) != 1:
            raise ValueError("fixed mode requires exactly one component bit width")
        if not isinstance(self.scale_quantiles, list) or not self.scale_quantiles or any(type(q) not in (int, float) or not math.isfinite(q) or not 0 < q <= 1 for q in self.scale_quantiles):
            raise ValueError("scale_quantiles must contain numbers in (0,1]")
        if len(self.scale_quantiles) != len(set(self.scale_quantiles)):
            raise ValueError("scale_quantiles must be unique")
        if self.bucket_reference is not None and (type(self.bucket_reference) not in (int, float) or not math.isfinite(self.bucket_reference) or self.bucket_reference <= 0):
            raise ValueError("bucket_reference must be positive and finite")
        bounds = {name: getattr(self, name) for name in ("max_component_mse", "max_complex_error", "max_complex_bias", "max_baseline_component_mse", "max_p99_error")}
        if not isinstance(self.product_bounds, dict):
            raise ValueError("product_bounds must be an object keyed by correlation name")
        for product, values in self.product_bounds.items():
            if not isinstance(product, str) or not isinstance(values, dict) or set(values) - set(bounds):
                raise ValueError("Unknown per-product error bound")
            self._validate_bounds({**bounds, **values})
        self._validate_bounds(bounds)
        if self.relative_bounds is not None:
            from .error_budget import validate_relative_bounds
            validate_relative_bounds(self.relative_bounds)
        if self.noise is not None:
            if not isinstance(self.noise, dict) or set(self.noise) != {"component_variance", "epsilon", "provenance", "validated"}:
                raise ValueError("noise requires component_variance, epsilon, provenance and validated")
            if self.noise["validated"] is not True or not isinstance(self.noise["provenance"], str) or not self.noise["provenance"].strip():
                raise ValueError("An explicitly validated noise model with provenance is required")
            for key in ("component_variance", "epsilon"):
                value = self.noise[key]
                if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
                    raise ValueError(f"noise.{key} must be positive and finite")
            if not math.isfinite(self.noise["component_variance"] * self.noise["epsilon"]):
                raise ValueError("Combined thermal-noise budget is not finite")
        if self.mode != "lossless":
            relative = self.relative_bounds or {}
            if self.max_component_mse is None and self.noise is None and 'max_component_rmse' not in relative:
                raise ValueError("Lossy encoding requires an absolute MSE bound, relative RMSE bound or a validated supplied noise model")
            if self.max_complex_error is None and 'max_complex_error' not in relative:
                raise ValueError("Lossy encoding also requires max_complex_error")
            for product in self.product_bounds:
                effective = self.bounds(product)
                if (effective["max_component_mse"] is None and 'max_component_rmse' not in relative or
                    effective["max_complex_error"] is None and 'max_complex_error' not in relative):
                    raise ValueError("Per-product overrides cannot remove required MSE/maximum-error bounds")

    @staticmethod
    def _validate_bounds(bounds):
        for key, value in bounds.items():
            if value is not None and (type(value) not in (int, float) or not math.isfinite(value) or value < 0):
                raise ValueError(f"{key} must be a finite nonnegative number or null")

    def bounds(self, product=None):
        names = ("max_component_mse", "max_complex_error", "max_complex_bias", "max_baseline_component_mse", "max_p99_error")
        result = {name: getattr(self, name) for name in names}
        if product is not None:
            result.update(self.product_bounds.get(product, {}))
        if self.noise is not None:
            thermal = self.noise["epsilon"] * self.noise["component_variance"]
            result["max_component_mse"] = min(thermal, result["max_component_mse"]) if result["max_component_mse"] is not None else thermal
        if result["max_baseline_component_mse"] is None:
            result["max_baseline_component_mse"] = result["max_component_mse"]
        return result

    def to_dict(self):
        return asdict(self)
