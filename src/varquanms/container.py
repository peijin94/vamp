"""HDF5 payloads, explicit inventory and independently checksummed I/O chunks."""

from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np

from .metadata import dumps, portable
from .flags import burned_columns
from .calibrated import validate_manifest

FORMAT = "varquanms.hdf5"
SCHEMA_VERSION = 1
TARGET_CHUNK_BYTES = 4 * 1024**2
MAX_MANIFEST_BYTES = 32 * 1024**2


def new_output(output, source):
    output = Path(output).resolve()
    source = Path(source).resolve()
    if output.exists() or output == source or output.is_relative_to(source):
        raise ValueError(f"Output must be new and outside the source: {output}")
    if not output.parent.is_dir():
        raise ValueError(f"Output parent does not exist: {output.parent}")
    return output


def dataset_chunks(shape, itemsize):
    if not shape or any(n == 0 for n in shape):
        return None
    row_bytes = max(1, int(np.prod(shape[1:], dtype=object)) * itemsize)
    if row_bytes > 64 * 1024**2:
        raise ValueError("A single array cell exceeds the 64 MiB prototype limit")
    return (min(shape[0], max(1, TARGET_CHUNK_BYTES // row_bytes)), *shape[1:])


def slabs(dataset):
    step = dataset.chunks[0] if dataset.chunks else max(1, dataset.shape[0])
    for start in range(0, dataset.shape[0], step):
        yield slice(start, min(dataset.shape[0], start + step))


class Writer:
    def __init__(self, output, config):
        self.output = Path(output)
        fd, path = tempfile.mkstemp(prefix=self.output.name + ".", suffix=".partial", dir=self.output.parent)
        os.close(fd)
        self.partial = Path(path)
        # Keep at least two target chunks cached. With HDF5's smaller default
        # cache, row-slab writes repeatedly recompressed an entire 1–4 MiB
        # chunk whenever it was larger than the cache.
        self.h5 = h5py.File(self.partial, "w", rdcc_nbytes=2 * TARGET_CHUNK_BYTES)
        self.h5.attrs.update(format=FORMAT, schema_version=SCHEMA_VERSION, complete=False)
        self.config = config
        self.roles = {}

    def create(self, name, shape, dtype, role, compress=True):
        dtype = np.dtype(dtype).newbyteorder("<")
        chunks = dataset_chunks(shape, dtype.itemsize)
        options = {}
        if chunks is not None:
            options = {"chunks": chunks, "fletcher32": True}
            if compress and self.config.compression == "gzip":
                options.update(compression="gzip", compression_opts=self.config.compression_level, shuffle=True)
        self.h5.require_group(str(Path(name).parent))
        dataset = self.h5.create_dataset(name, shape=shape, dtype=dtype, **options)
        self.roles[name] = role
        return dataset

    def put(self, name, values, role, compress=True):
        values = portable(values)
        dataset = self.create(name, values.shape, values.dtype, role, compress)
        dataset[...] = values
        return name

    def put_json(self, name, value, role):
        return self.put(name, np.frombuffer(dumps(value), dtype=np.uint8), role)

    def create_stream(self, name, role):
        self.h5.require_group(str(Path(name).parent))
        dataset = self.h5.create_dataset(name, shape=(0,), maxshape=(None,), dtype="u1",
                                         chunks=(1024 * 1024,), fletcher32=True)
        self.roles[name] = role
        return dataset

    def finish(self, manifest):
        self.h5.flush()
        inventory = []
        categories = defaultdict(int)
        for name, role in self.roles.items():
            dataset = self.h5[name]
            digest = hashlib.sha256()
            for selection in slabs(dataset):
                digest.update(portable(dataset[selection]).tobytes())
            size = dataset.id.get_storage_size()
            categories[role] += size
            inventory.append({"path": name, "shape": list(dataset.shape), "dtype": dataset.dtype.str,
                              "sha256": digest.hexdigest(), "stored_bytes": size,
                              "logical_bytes": dataset.size * dataset.dtype.itemsize, "role": role})
        manifest.update(schema="varquanms.archive/1", inventory=inventory,
                        dataset_stored_bytes_by_category=dict(categories))
        encoded = dumps(manifest)
        if len(encoded) > MAX_MANIFEST_BYTES:
            raise ValueError("Archive manifest exceeds the 32 MiB prototype limit")
        self.h5.create_dataset("manifest", data=np.frombuffer(encoded, dtype=np.uint8))
        self.h5.attrs["manifest_sha256"] = hashlib.sha256(encoded).hexdigest()
        self.h5.attrs["complete"] = True
        self.h5.flush()
        self.h5.close()
        with self.partial.open("rb") as stream:
            os.fsync(stream.fileno())
        os.link(self.partial, self.output)
        self.partial.unlink()
        return manifest

    def close(self):
        if self.h5:
            self.h5.close()


class Reader:
    def __init__(self, path, verify=True):
        self.path = Path(path)
        self.h5 = h5py.File(self.path, "r")
        try:
            if self.h5.attrs.get("format") != FORMAT or self.h5.attrs.get("schema_version") != SCHEMA_VERSION or self.h5.attrs.get("complete") != True:
                raise ValueError("Unsupported or incomplete archive")
            if "manifest" not in self.h5 or self.h5["manifest"].ndim != 1 or self.h5["manifest"].size > MAX_MANIFEST_BYTES:
                raise ValueError("Missing or oversized manifest")
            encoded = self.h5["manifest"][...].tobytes()
            if hashlib.sha256(encoded).hexdigest() != self.h5.attrs.get("manifest_sha256"):
                raise ValueError("Manifest checksum mismatch")
            self.manifest = json.loads(encoded)
            if self.manifest.get("schema") != "varquanms.archive/1":
                raise ValueError("Unsupported scientific schema")
            burned_columns(self.manifest)
            validate_manifest(self.manifest)
            self.entries = {entry["path"]: entry for entry in self.manifest["inventory"]}
            if len(self.entries) != len(self.manifest["inventory"]):
                raise ValueError("Duplicate dataset inventory entries")
            paths = set()
            def visit(group):
                for name in group:
                    if not isinstance(group.get(name, getlink=True), h5py.HardLink):
                        raise ValueError("External and soft HDF5 links are unsupported")
                    obj = group[name]
                    if isinstance(obj, h5py.Dataset):
                        paths.add(obj.name.lstrip("/"))
                    else:
                        visit(obj)
            visit(self.h5)
            if paths != set(self.entries) | {"manifest"}:
                raise ValueError("Dataset inventory mismatch")
            if verify:
                self.verify()
        except BaseException:
            self.h5.close()
            raise

    def verify(self):
        for name, entry in self.entries.items():
            dataset = self.h5[name]
            if list(dataset.shape) != entry["shape"] or dataset.dtype.str != entry["dtype"]:
                raise ValueError(f"Dataset shape/dtype mismatch: {name}")
            digest = hashlib.sha256()
            for selection in slabs(dataset):
                digest.update(portable(dataset[selection]).tobytes())
            if digest.hexdigest() != entry["sha256"]:
                raise ValueError(f"Dataset checksum mismatch: {name}")

    def json(self, name):
        if name not in self.entries:
            raise ValueError(f"Dataset is absent from inventory: {name}")
        return json.loads(self.h5[name][...].tobytes())

    def close(self):
        self.h5.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
