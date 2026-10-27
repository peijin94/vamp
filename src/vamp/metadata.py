"""Non-executable typed metadata and relocatable casacore table references."""

import base64
import json
from pathlib import Path

import numpy as np


def portable(array):
    a = np.asarray(array)
    if a.dtype.kind not in "biufc":
        raise ValueError(f"Unsupported numeric payload dtype: {a.dtype}")
    return np.ascontiguousarray(a.astype(a.dtype.newbyteorder("<"), copy=False))


def pack(value):
    if isinstance(value, dict):
        return ["record", [[str(k), pack(v)] for k, v in value.items()]]
    if isinstance(value, np.generic):
        return ["numpy-scalar", pack(np.asarray(value))]
    if isinstance(value, np.ndarray):
        a = np.asarray(value)
        if a.dtype.kind in "biufc":
            raw = portable(a).tobytes()
            return ["array", a.dtype.newbyteorder("<").str, list(a.shape), base64.b64encode(raw).decode("ascii")]
        if a.dtype.kind in "US":
            return ["string-array", a.dtype.str, list(a.shape), a.tolist()]
        raise ValueError(f"Unsupported metadata array dtype: {a.dtype}")
    if isinstance(value, (list, tuple)):
        return ["list", [pack(v) for v in value]]
    if isinstance(value, bool):
        return ["bool", value]
    if isinstance(value, int):
        return ["int", str(value)]
    if isinstance(value, float):
        return ["float", np.asarray(value, dtype="<f8").tobytes().hex()]
    if isinstance(value, complex):
        return ["complex", np.asarray(value, dtype="<c16").tobytes().hex()]
    if isinstance(value, str):
        return ["string", value]
    if value is None:
        return ["null"]
    raise ValueError(f"Unsupported metadata object: {type(value).__name__}")


def unpack(node):
    if not isinstance(node, list) or not node or not isinstance(node[0], str):
        raise ValueError("Invalid typed metadata")
    tag = node[0]
    if tag == "numpy-scalar" and len(node) == 2:
        value = unpack(node[1])
        if not isinstance(value, np.ndarray) or value.ndim != 0:
            raise ValueError("Invalid numpy scalar")
        return value[()]
    if tag in ("float", "complex") and len(node) == 2:
        dtype = "<f8" if tag == "float" else "<c16"
        value = np.frombuffer(bytes.fromhex(node[1]), dtype=dtype)
        if value.size != 1:
            raise ValueError("Invalid numeric scalar")
        return value[0].item()
    if tag == "record" and len(node) == 2:
        result = {}
        for key, value in node[1]:
            if not isinstance(key, str) or key in result:
                raise ValueError("Invalid or duplicate metadata key")
            result[key] = unpack(value)
        return result
    if tag == "array" and len(node) == 4:
        dtype = np.dtype(node[1])
        shape = tuple(node[2])
        if dtype.kind not in "biufc" or any(type(n) is not int or n < 0 for n in shape):
            raise ValueError("Invalid metadata dtype/shape")
        raw = base64.b64decode(node[3], validate=True)
        if len(raw) != int(np.prod(shape, dtype=object)) * dtype.itemsize:
            raise ValueError("Metadata array length mismatch")
        return np.frombuffer(raw, dtype=dtype).copy().reshape(shape)
    if tag == "string-array" and len(node) == 4:
        dtype = np.dtype(node[1])
        if dtype.kind not in "US":
            raise ValueError("Invalid string metadata dtype")
        return np.asarray(node[3], dtype=dtype).reshape(node[2])
    if tag == "list" and len(node) == 2:
        return [unpack(v) for v in node[1]]
    if tag == "int" and len(node) == 2:
        return int(node[1])
    if tag == "bool" and len(node) == 2 and type(node[1]) is bool:
        return node[1]
    if tag == "string" and len(node) == 2 and isinstance(node[1], str):
        return node[1]
    if tag == "null" and len(node) == 1:
        return None
    raise ValueError(f"Invalid metadata tag: {tag}")


def dumps(value):
    return json.dumps(value, separators=(",", ":"), allow_nan=False).encode("utf-8")


def references(value):
    if isinstance(value, str) and value.startswith("Table: "):
        yield value[7:]
    elif isinstance(value, dict):
        for v in value.values():
            yield from references(v)
    elif isinstance(value, (list, tuple)):
        for v in value:
            yield from references(v)


def relocate(value, root, restore=False):
    if isinstance(value, str) and value.startswith("Table: "):
        name = value[7:]
        path = Path(name)
        target = (root / path).resolve() if not path.is_absolute() else path.resolve()
        if not target.is_relative_to(root):
            raise ValueError(f"External table reference is unsupported: {value}")
        if restore:
            return "Table: " + str(target)
        return "Table: ./" + target.relative_to(root).as_posix()
    if isinstance(value, dict):
        return {k: relocate(v, root, restore) for k, v in value.items()}
    if isinstance(value, list):
        return [relocate(v, root, restore) for v in value]
    if isinstance(value, tuple):
        return tuple(relocate(v, root, restore) for v in value)
    return value


def logical_descriptor(description):
    """Storage-manager assignment is explicitly outside exact export equality."""
    result = {}
    for key, value in description.items():
        if isinstance(value, dict) and "valueType" in value:
            result[key] = {k: v for k, v in value.items()
                           if k not in ("dataManagerType", "dataManagerGroup")}
        else:
            result[key] = value
    return result


def casacore_keywords(value):
    """casacore returns string vectors as lists, including untyped empty lists."""
    if isinstance(value, list) and not value:
        return np.asarray([], dtype="U1")
    if isinstance(value, dict):
        return {k: casacore_keywords(v) for k, v in value.items()}
    return value


def equal(a, b):
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(equal(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(equal(x, y) for x, y in zip(a, b))
    if isinstance(a, (np.ndarray, np.generic)) or isinstance(b, (np.ndarray, np.generic)):
        x, y = np.asarray(a), np.asarray(b)
        if x.shape != y.shape or x.dtype.kind != y.dtype.kind or x.dtype.itemsize != y.dtype.itemsize:
            return False
        if x.dtype.kind in "biufc":
            return portable(x).tobytes() == portable(y).tobytes()
        return np.array_equal(x, y)
    if isinstance(a, (float, complex)) or isinstance(b, (float, complex)):
        return equal(np.asarray(a), np.asarray(b))
    return type(a) is type(b) and a == b
