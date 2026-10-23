"""Explicit calibrated-only projection; the input table is never modified."""
from copy import deepcopy

import numpy as np

from .metadata import equal, unpack

OMITTED = ("DATA", "FLOAT_DATA", "MODEL_DATA")
KEYWORD = "VARQUANMS_CALIBRATED_ONLY"
WEIGHT_KEYWORD = "VARQUANMS_WEIGHT_POLICY"
DATA_KEYWORD = "VARQUANMS_DERIVED_DATA"


def data_marker():
    return {"source_column": "CORRECTED_DATA", "values": "bitwise copy of decoded calibrated visibilities",
            "original_raw_DATA_restored": False}


def materialize_ms_v2_data(path):
    """Add a declared calibrated DATA copy to a newly created export only."""
    from casacore.tables import table
    with table(str(path), readonly=False, ack=False) as t:
        if KEYWORD not in t.getkeywords() or 'DATA' in t.colnames():
            raise ValueError('DATA compatibility export requires a calibrated-only table without DATA')
        t.addcols({'DATA': t.getcoldesc('CORRECTED_DATA')})
        for start in range(0, t.nrows(), 1024):
            count = min(1024, t.nrows()-start)
            t.putcol('DATA', t.getcol('CORRECTED_DATA', start, count), startrow=start, nrow=count)
        t.putkeyword(DATA_KEYWORD, data_marker())


def marker():
    return {"profile": "calibrated-only-v1", "omitted_if_present": list(OMITTED),
            "visibility_column": "CORRECTED_DATA", "flagged_value_reconstruction": "complex positive zero",
            "flag_mask": "FLAG OR FLAG_ROW", "original_flag_values": "preserved",
            "weights": "explicit VARQUANMS_WEIGHT_POLICY; SIGMA/SIGMA_SPECTRUM remain lossless",
            "original_flagged_visibility_values": "discarded"}


def descriptor(original, configuration=None):
    configuration = configuration or {}
    result = deepcopy(original)
    for name in OMITTED:
        result.pop(name, None)
    private = result.get("_private_keywords_", {})
    for name in list(private):
        value = private[name]
        if name.startswith("Hypercolumn_") and any(c in OMITTED for k in ("data", "coord", "id") for c in value[k]):
            del private[name]
    # getdesc may expose both private records and constructor declarations.
    hypercolumns = result.get("_define_hypercolumn_", {})
    for name in list(hypercolumns):
        if any(c in OMITTED for k in ("HCdatanames", "HCcoordnames", "HCidnames") for c in hypercolumns[name].get(k, [])):
            del hypercolumns[name]
    keywords = result.setdefault("_keywords_", {})
    if KEYWORD in keywords and not equal(keywords[KEYWORD], marker()):
        raise ValueError("Conflicting calibrated-only provenance keyword")
    keywords[KEYWORD] = marker()
    policy = {"bits": configuration.get("weight_bits", 4),
              "max_relative_error": configuration.get("weight_max_relative_error", .01),
              "columns": ["WEIGHT", "WEIGHT_SPECTRUM"], "float32_fallback": True,
              "zero_and_special_weights": "exact", "profile_scope": "calibrated-only"}
    if WEIGHT_KEYWORD in keywords and not equal(keywords[WEIGHT_KEYWORD], policy):
        raise ValueError("Conflicting weight precision provenance keyword")
    keywords[WEIGHT_KEYWORD] = policy
    return result


class CalibratedView:
    """Read-through view exposing retained columns and zeroed flagged values."""
    def __init__(self, source, configuration=None):
        self.source = source
        self.configuration = configuration
        if "CORRECTED_DATA" not in source.colnames():
            raise ValueError("keep_calib_only requires CORRECTED_DATA")

    def __getattr__(self, name):
        return getattr(self.source, name)

    def colnames(self):
        return [c for c in self.source.colnames() if c not in OMITTED]

    def getdesc(self):
        return descriptor(self.source.getdesc(), self.configuration)

    def getcol(self, column, startrow=0, nrow=-1, rowincr=1):
        if column in OMITTED:
            raise ValueError("Column omitted by calibrated-only profile: " + column)
        values = self.source.getcol(column, startrow, nrow, rowincr)
        if column == "CORRECTED_DATA":
            flags = self.source.getcol("FLAG", startrow, nrow, rowincr)
            flags |= self.source.getcol("FLAG_ROW", startrow, nrow, rowincr)[:, None, None]
            if flags.shape != values.shape:
                raise ValueError("FLAG and CORRECTED_DATA shape mismatch")
            values[flags] = 0j
        return values

    def getcell(self, column, row):
        if column in OMITTED:
            raise ValueError("Column omitted by calibrated-only profile: " + column)
        value = self.source.getcell(column, row)
        if column == "CORRECTED_DATA":
            flags = self.source.getcell("FLAG", row) | self.source.getcell("FLAG_ROW", row)
            if flags.shape != value.shape:
                raise ValueError("FLAG and CORRECTED_DATA shape mismatch")
            value[flags] = 0j
        return value


def validate_manifest(manifest):
    enabled = manifest["config"].get("keep_calib_only", False)
    if type(enabled) is not bool:
        raise ValueError("keep_calib_only must be boolean")
    for spec in manifest['tables']:
        for name, column in spec['columns'].items():
            if column['encoding'] == 'weight-codebook4-v1' and (not enabled or spec['path'] != '.' or name not in ('WEIGHT', 'WEIGHT_SPECTRUM')):
                raise ValueError('Four-bit weight encoding is outside its declared calibrated-only policy')
    if not enabled:
        return
    if manifest["config"].get("burn_flag") or manifest["config"]["columns"] != ["CORRECTED_DATA"]:
        raise ValueError("Calibrated-only policy conflicts with flags or selected columns")
    main = next(t for t in manifest["tables"] if t["path"] == ".")
    names = main["column_names"]
    if "CORRECTED_DATA" not in names or any(c in names or c in main["columns"] for c in OMITTED):
        raise ValueError("Calibrated-only column inventory mismatch")
    if not equal(unpack(main["descriptor"]).get("_keywords_", {}).get(KEYWORD), marker()):
        raise ValueError("Calibrated-only archive lacks its explicit provenance marker")
    original = unpack(main["descriptor"])
    expected = descriptor(original, manifest["config"])
    if not equal(original.get("_keywords_", {}).get(WEIGHT_KEYWORD), expected["_keywords_"][WEIGHT_KEYWORD]):
        raise ValueError("Missing weight precision provenance")
    for name in ("WEIGHT", "WEIGHT_SPECTRUM"):
        if name in main["columns"] and ((main["columns"][name]["encoding"] == "weight-codebook4-v1") != (manifest["config"].get("weight_bits", 4) == 4)):
            raise ValueError("Weight encoding disagrees with calibrated-only precision policy")
