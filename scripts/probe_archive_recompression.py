"""Measure lossless storage recompression; do not change quantized values."""
import argparse
from copy import deepcopy
import json
from pathlib import Path
import time

from varquanms.config import Config
from varquanms.container import Reader, Writer, slabs
from varquanms.imaging import file_sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    records = []
    for level in (1, 9):
        output = args.output/f'gzip{level}.h5'
        started = time.perf_counter()
        with Reader(args.archive) as reader:
            manifest = deepcopy(reader.manifest)
            original_hashes = {e['path']:e['sha256'] for e in manifest['inventory']}
            config = Config.load(None, **{**manifest['config'], 'compression_level':level, 'compression':'gzip'})
            writer = Writer(output, config)
            try:
                for name, entry in reader.entries.items():
                    source = reader.h5[name]
                    target = writer.create(name, source.shape, source.dtype, entry['role'])
                    for selection in slabs(source):
                        target[selection] = source[selection]
                manifest['storage_repack'] = dict(source_file_sha256=file_sha256(args.archive),
                    container_gzip_level=level, original_encoding_config_unchanged=True,
                    description='Lossless recompression of all logical arrays including existing compressed block streams; codec values unchanged')
                manifest = writer.finish(manifest)
            finally:
                writer.close()
        with Reader(output) as check:
            if any(e['sha256'] != original_hashes[e['path']] for e in check.manifest['inventory']):
                raise ValueError('Recompression changed logical values')
        record = dict(gzip_level=level, source_bytes=args.archive.stat().st_size, output_bytes=output.stat().st_size,
                      elapsed_seconds=time.perf_counter()-started, all_logical_array_checksums_unchanged=True,
                      categories=manifest['dataset_stored_bytes_by_category'])
        records.append(record)
        print(json.dumps(record), flush=True)
    (args.output/'report.json').write_text(json.dumps(records, indent=2)+'\n')


if __name__ == '__main__':
    main()
