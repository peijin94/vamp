"""Lossless HDF5 packing of ordered partitions with shared immutable arrays.

Each partition retains its original scientific manifest and relative datasets.
Only byte-identical arrays of the same role, dtype and shape share hard links.
No external links, source MS, requantization or downloaded decoder is needed.
"""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

import h5py
import numpy as np

from .container import MAX_MANIFEST_BYTES, Reader, new_output, slabs
from .imaging import file_sha256
from .metadata import dumps, portable

FORMAT = 'varquanms.shared-partitions'
SCHEMA = FORMAT + '/1'


def _json_bytes(h5, name, value):
    h5.create_dataset(name, data=np.frombuffer(value, dtype='u1'),
                      compression='gzip', compression_opts=1, fletcher32=True)


def _same_arrays(a, b):
    if a.shape != b.shape or a.dtype != b.dtype:
        return False
    return all(portable(a[s]).tobytes() == portable(b[s]).tobytes() for s in slabs(a))


def _paths(h5):
    """Reject escaping links and group cycles before following any payloads."""
    datasets, groups, seen = set(), set(), set()
    def visit(group, prefix=''):
        address = h5py.h5o.get_info(group.id).addr
        if address in seen:
            raise ValueError('Aliased or cyclic HDF5 groups are unsupported')
        seen.add(address)
        for name in group:
            if not isinstance(group.get(name, getlink=True), h5py.HardLink):
                raise ValueError('External and soft HDF5 links are unsupported')
            path = prefix + name
            obj = group[name]
            if isinstance(obj, h5py.Dataset):
                if obj.is_virtual or obj.external:
                    raise ValueError('Virtual or externally stored datasets are unsupported')
                datasets.add(path)
            else:
                groups.add(path)
                visit(obj, path + '/')
    visit(h5)
    return datasets, groups


def load_shared(path):
    """Validate every partition and all row mappings before permitting export."""
    from .partitioned import _signature
    with h5py.File(path, 'r') as h5:
        if (h5.attrs.get('format') != FORMAT or h5.attrs.get('schema_version') != 1
                or h5.attrs.get('complete') != True):
            raise ValueError('Unsupported or incomplete shared archive')
        paths, groups = _paths(h5)
        if ('manifest' not in paths or h5['manifest'].ndim != 1
                or h5['manifest'].dtype != np.dtype('u1') or h5['manifest'].size > MAX_MANIFEST_BYTES):
            raise ValueError('Missing or oversized shared manifest')
        encoded = h5['manifest'][...].tobytes()
        if hashlib.sha256(encoded).hexdigest() != h5.attrs.get('manifest_sha256'):
            raise ValueError('Shared manifest checksum mismatch')
        manifest = json.loads(encoded)
        if manifest.get('schema') != SCHEMA or manifest.get('status') != 'complete':
            raise ValueError('Unsupported or incomplete shared manifest')
        parts = manifest['partitions']
        if not 1 <= len(parts) <= 4096:
            raise ValueError('Invalid shared partition count')
        expected = {'manifest'}
        cursor, signature = 0, None
        for index, part in enumerate(parts):
            start, stop = part['start'], part['stop']
            if type(start) is not int or type(stop) is not int or start != cursor or stop <= start:
                raise ValueError('Shared partition mapping has gaps, overlaps or invalid bounds')
            cursor = stop
            group, verification = f'partitions/p{index:04d}', f'verification/v{index:04d}'
            if part.get('group') != group or part.get('verification_dataset') != verification:
                raise ValueError('Invalid shared member path')
            if group not in groups or verification not in paths:
                raise ValueError('Missing shared partition or verification')
            raw = h5[verification]
            if raw.ndim != 1 or raw.dtype != np.dtype('u1') or raw.size > MAX_MANIFEST_BYTES:
                raise ValueError('Invalid shared verification dataset')
            raw = raw[...].tobytes()
            if hashlib.sha256(raw).hexdigest() != part['verification_sha256']:
                raise ValueError('Shared verification checksum mismatch')
            if json.loads(raw).get('status') != 'passed':
                raise ValueError('Shared partition has no passed source verification')
            expected.add(verification)
            with Reader(path, group=group) as reader:
                if (reader.manifest['source']['main_row_count'] != stop-start
                        or reader.manifest['config'] != manifest['config']):
                    raise ValueError('Shared partition row count/config mismatch')
                current = _signature(reader.manifest)
                if signature is not None and current != signature:
                    raise ValueError('Shared partition schemas/subtables differ')
                signature = current
                expected.update(group + '/' + name for name in reader.entries)
                expected.add(group + '/manifest')
        if cursor != manifest['row_count']:
            raise ValueError('Shared archive row count mismatch')
        if paths != expected:
            raise ValueError('Shared dataset inventory mismatch')
        expected_groups = {str(parent) for name in expected for parent in Path(name).parents if str(parent) != '.'}
        if groups != expected_groups:
            raise ValueError('Shared group inventory mismatch')
    return manifest


def compact_bundle(source, output):
    """Share exact arrays in a new file; preserve original bundle and codec bits."""
    from .partitioned import _load, bundle_summary
    started = time.perf_counter()
    source = Path(source).resolve()
    if not source.is_dir():
        raise ValueError('Compaction currently requires an ordered directory bundle')
    output = new_output(output, source)
    _, original = _load(source)
    source_validation_seconds = time.perf_counter()-started
    manifest = deepcopy(original)
    manifest['schema'] = SCHEMA
    manifest['scope'] = 'Ordered partitions in one HDF5 file; exact shared arrays; unchanged scientific manifests and codec bits'
    manifest['storage_compaction'] = dict(
        source_manifest_sha256=file_sha256(source/'manifest.json'),
        source_bundle_bytes=bundle_summary(source, original)['archive_bytes'],
        scientific_manifests_unchanged=True, requantized=False,
        equality='Same role, shape and dtype; verified SHA256 plus bitwise comparison before hard linking',
        source_validation_seconds=source_validation_seconds)
    fd, partial = tempfile.mkstemp(prefix=output.name+'.', suffix='.partial', dir=output.parent)
    os.close(fd)
    partial = Path(partial)
    seen, counts, savings = {}, Counter(), Counter()
    copying = time.perf_counter()
    with h5py.File(partial, 'w') as h5:
        h5.attrs.update(format=FORMAT, schema_version=1, complete=False)
        h5.create_group('partitions'); h5.create_group('verification')
        for index, part in enumerate(manifest['partitions']):
            group = h5.create_group(f'partitions/p{index:04d}')
            part['group'] = group.name.lstrip('/')
            part['verification_dataset'] = f'verification/v{index:04d}'
            _json_bytes(h5, part['verification_dataset'], (source/part['verification']).read_bytes())
            with Reader(source/part['archive'], verify=False) as reader:
                for key, value in reader.h5.attrs.items():
                    group.attrs[key] = value
                _json_bytes(group, 'manifest', reader.h5['manifest'][...].tobytes())
                for name, entry in reader.entries.items():
                    key = (entry['role'], entry['dtype'], tuple(entry['shape']), entry['sha256'])
                    parent = group.require_group(str(Path(name).parent))
                    if key in seen:
                        previous = h5[seen[key]]
                        if not _same_arrays(reader.h5[name], previous):
                            raise ValueError('Hash-matched arrays differ; refusing to share')
                        parent[Path(name).name] = previous
                        savings[entry['role']] += entry['stored_bytes']
                        counts['shared_array_references'] += 1
                    else:
                        reader.h5.copy(reader.h5[name], parent, name=Path(name).name)
                        seen[key] = group.name + '/' + name
                        counts['unique_arrays'] += 1
                    counts['array_references'] += 1
        manifest['storage_compaction'].update(counts=dict(counts),
            duplicate_stored_bytes_by_category=dict(savings),
            copy_and_share_seconds=time.perf_counter()-copying)
        encoded = dumps(manifest)
        if len(encoded) > MAX_MANIFEST_BYTES:
            raise ValueError('Shared manifest exceeds size limit')
        _json_bytes(h5, 'manifest', encoded)
        h5.attrs.update(manifest_sha256=hashlib.sha256(encoded).hexdigest(), complete=True)
    checking = time.perf_counter()
    checked = load_shared(partial)
    output_validation_seconds = time.perf_counter()-checking
    with partial.open('rb') as stream:
        os.fsync(stream.fileno())
    os.link(partial, output)
    partial.unlink()
    result = shared_summary(output, checked)
    result.update(elapsed_seconds=time.perf_counter()-started,
                  output_validation_seconds=output_validation_seconds,
                  all_partition_inventories_verified=True)
    return result


def shared_summary(path, manifest=None):
    path = Path(path)
    manifest = manifest or load_shared(path)
    categories, seen = Counter(), set()
    with h5py.File(path, 'r') as h5:
        for part in manifest['partitions']:
            group = h5[part['group']]
            spec = json.loads(group['manifest'][...].tobytes())
            for entry in spec['inventory']:
                ds = group[entry['path']]
                address = h5py.h5o.get_info(ds.id).addr
                if address not in seen:
                    seen.add(address)
                    categories[entry['role']] += ds.id.get_storage_size()
    total = path.stat().st_size
    categories['container_manifests_verification_and_links'] = total-sum(categories.values())
    return dict(schema=SCHEMA, archive=str(path), archive_bytes=total,
                partitions=len(manifest['partitions']), row_count=manifest['row_count'],
                source_physical_bytes=manifest['source_physical_bytes'],
                stored_bytes_by_category=dict(categories), storage_compaction=manifest['storage_compaction'],
                whole_archive_ratio=manifest['source_physical_bytes']/total)
