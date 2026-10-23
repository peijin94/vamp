"""Packed 4-bit float32 weight codebooks with measured relative-error fallback."""
import json
import numpy as np

from .codec import pack_codes, unpack_codes

COLUMNS = ("WEIGHT", "WEIGHT_SPECTRUM")
ENCODING = "weight-codebook4-v1"


def encode_chunk(values, max_relative_error):
    a = np.ascontiguousarray(values, dtype='<f4').ravel()
    if not a.size:
        return None
    distinct, inverse = np.unique(a.view('<u4'), return_inverse=True)
    if len(distinct) <= 16:
        # Bit-pattern dictionary also preserves NaN payloads and signed zeros.
        return distinct.view('<f4'), pack_codes(inverse, 4), 0., 0., False
    if not np.all(np.isfinite(a)) or np.any(a < 0) or np.any((a == 0) & np.signbit(a)):
        return None
    book = np.linspace(float(a.min()), float(a.max()), 16, dtype=np.float64).astype('<f4')
    high = np.clip(np.searchsorted(book, a), 1, 15)
    low = high - 1
    dl, dh = abs(a.astype(np.float64)-book[low]), abs(a.astype(np.float64)-book[high])
    codes = np.where((dh < dl) | ((dh == dl) & (low % 2 == 1)), high, low)
    reconstructed = book[codes]
    error = abs(reconstructed.astype(np.float64)-a)
    positive = a > 0
    relative = float(np.max(error[positive]/a[positive], initial=0))
    if relative > max_relative_error or np.any(reconstructed[a == 0] != 0) or np.any(reconstructed[positive] <= 0):
        return None
    return book, pack_codes(codes, 4), float(error.max(initial=0)), relative, bool(np.any(error))


def write_weights(writer, t, column, prefix, config, info):
    shape = tuple(json.loads(next(iter(info['cell_shapes']))))
    spec = {'encoding': ENCODING, 'dtype': '<f4', 'row_count': t.nrows(), 'cell_shape': list(shape),
            'chunks': [], 'quality': {'four_bit_chunks': 0, 'float32_chunks': 0, 'lossy_chunks': 0,
                                     'max_absolute_error': 0., 'max_relative_error': 0., 'packed_code_bytes': 0}}
    for start in range(0, t.nrows(), config.io_chunk_rows):
        stop = min(t.nrows(), start+config.io_chunk_rows)
        a = np.asarray(t.getcol(column, start, stop-start), dtype='<f4')
        candidate = encode_chunk(a, config.weight_max_relative_error)
        record = {'start': start, 'stop': stop}
        if candidate is None:
            record.update(kind='float32', values=writer.put(prefix+f'/raw_{start}', a, 'flags_weights'))
            spec['quality']['float32_chunks'] += 1
        else:
            book, payload, absolute, relative, lossy = candidate
            record.update(kind='codebook4', count=a.size,
                          codebook=writer.put(prefix+f'/book_{start}', book, 'flags_weights'),
                          codes=writer.put(prefix+f'/codes_{start}', np.frombuffer(payload, dtype='u1'), 'flags_weights'))
            quality = spec['quality']
            quality['four_bit_chunks'] += 1; quality['lossy_chunks'] += int(lossy)
            quality['packed_code_bytes'] += len(payload)
            quality['max_absolute_error'] = max(quality['max_absolute_error'], absolute)
            quality['max_relative_error'] = max(quality['max_relative_error'], relative)
        spec['chunks'].append(record)
    spec['quality']['original_bits_exact'] = spec['quality']['lossy_chunks'] == 0
    return spec


def read_weights(reader, spec, start, stop):
    shape = spec['cell_shape']
    if spec['dtype'] != '<f4' or any(type(n) is not int or n < 0 for n in shape):
        raise ValueError('Invalid weight dtype/shape')
    output = np.empty((stop-start, *shape), dtype='<f4')
    cursor = start
    for chunk in spec['chunks']:
        lo, hi = max(start, chunk['start']), min(stop, chunk['stop'])
        if lo >= hi: continue
        if lo != cursor or chunk['stop'] <= chunk['start']:
            raise ValueError('Weight chunks have gaps/overlaps')
        count = (chunk['stop']-chunk['start'])*int(np.prod(shape, dtype=object))
        if count > 16*1024**2: raise ValueError('Weight chunk exceeds 64 MiB')
        if chunk['kind'] == 'float32':
            a = reader.h5[chunk['values']][...]
            if a.dtype != np.dtype('<f4') or a.shape != (chunk['stop']-chunk['start'], *shape):
                raise ValueError('Invalid raw weight chunk')
        elif chunk['kind'] == 'codebook4':
            book = reader.h5[chunk['codebook']][...]
            if chunk['count'] != count or book.dtype != np.dtype('<f4') or book.ndim != 1 or not 1 <= len(book) <= 16:
                raise ValueError('Invalid weight codebook/count')
            codes = unpack_codes(reader.h5[chunk['codes']][...].tobytes(), 4, count)
            if np.any(codes >= len(book)): raise ValueError('Weight code outside codebook')
            a = book[codes].reshape(chunk['stop']-chunk['start'], *shape)
        else: raise ValueError('Unknown weight chunk kind')
        output[lo-start:hi-start] = a[lo-chunk['start']:hi-chunk['start']]
        cursor = hi
    if cursor != stop: raise ValueError('Weight chunks have missing rows')
    return output


def verify_weights(reader, spec, t, column, configuration):
    maximum, absolute = 0., 0.
    exact = True
    for chunk in spec['chunks']:
        start, stop = chunk['start'], chunk['stop']
        original = np.asarray(t.getcol(column, start, stop-start))
        decoded = read_weights(reader, spec, start, stop)
        protected = ~np.isfinite(original) | (original <= 0)
        if original[protected].tobytes() != decoded[protected].tobytes():
            raise ValueError('Special/zero/negative weight bits changed')
        valid = ~protected
        error = abs(decoded[valid].astype(np.float64)-original[valid])
        relative = float(np.max(error/original[valid], initial=0))
        if not np.isfinite(decoded[valid]).all() or np.any(decoded[valid] <= 0) or relative > configuration['weight_max_relative_error']:
            raise ValueError('Weight reconstruction exceeds its relative-error bound')
        maximum = max(maximum, relative)
        absolute = max(absolute, float(error.max(initial=0)))
        exact &= original.tobytes() == decoded.tobytes()
    expected = spec['quality']
    if maximum != expected['max_relative_error'] or absolute != expected['max_absolute_error'] or exact != expected['original_bits_exact']:
        raise ValueError('Weight errors differ from encode-time accounting')
    return dict(expected)
