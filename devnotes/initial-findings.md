# Initial inspection findings

These are observations of the input, not codec benchmark results.

The copied input is `20260904_190008_59MHz.ms`, the lexically first MS tar in
the requested 59 MHz directory. The source tar is 409,661,440 bytes. Its SHA-256
is `3899b6b7111ec9a9ed8b69b9cee4bdf488b197c3f2379708b9d971ffec43409e`.
The local tar hash and all 65 extracted regular-file hashes match the source.
The observation payload remains outside the Git inventory.

MAIN contains 62,128 rows: the full upper triangle of 352 ordered antennas,
including 352 autocorrelation rows. All rows use one DATA_DESC_ID, one field
and one scan. No duplicate distinguishing row keys or reversed antenna pairs
were found. There is one 192-channel spectral window and one polarization
configuration. DATA is complex64, with Python cell shape `[192,4]` and product
order `XX,XY,YX,YY`. The logical DATA size is 381,714,432 bytes. FLOAT_DATA,
CORRECTED_DATA and MODEL_DATA are absent.

TIME is `5295265208.539846` and TIME_CENTROID is `5295265213.555357`, both float64
in seconds with UTC reference. The difference is approximately 5.0155 seconds;
these are separate original values and must remain separate. The stored interval
and exposure, along with their descriptors, are in the JSON inventory. No
coordinate regeneration or time conversion was performed.

Channel frequencies span 59,335,937.5 through 63,905,761.71875 Hz in original
order. Every width is 23,925.78125 Hz. MEAS_FREQ_REF is the original integer 5;
the full frequency measure keyword mapping is retained in the spectral-window
column keywords. No frequencies have been relabeled from the directory name.

UVW uses float64 metre values with ITRF reference. Of 62,128 vectors, 61,777 are
bitwise distinct. The 352 autocorrelation rows share the same geometry vector.
UVW occupies 1,491,072 logical bytes. A naive dictionary of 61,777 vectors plus
one uint32 per row would occupy 1,731,160 logical bytes, so that specific
dictionary representation would increase the uncompressed logical size. Other
exact representations and secondary compression remain unmeasured. There is no
192-fold UVW repetition across channel records in this input.

FLAG, FLAG_ROW, FLAG_CATEGORY, WEIGHT and SIGMA are present. All selected DATA
samples are unflagged and finite, with no negative-zero components detected.
This observation does not establish that DATA is RFI-free or that weights and
sigmas form a validated thermal-noise model. No such model has been accepted.

All MAIN columns share StandardStMan. The extracted MS files total 409,594,838
bytes, with 409,800,704 bytes allocated in file blocks on the current filesystem.
MAIN logical arrays total 485,406,064 bytes and subtables total 91,133 logical
value bytes. Logical bool arrays count one byte per NumPy value, whereas
casacore packs booleans. Shared storage files prevent exact per-column physical
attribution. These figures therefore do not establish achievable archive ratios.

The inspector's ten synthetic tests passed, covering both channel organizations,
mixed fields, duplicate rows, antenna ordering, heterogeneous shapes, all-zero
and all-flagged data, special-value bits, Stokes order, output protections and
unknown options. Inspecting a fixture changed none of its file hashes. Real-data
copy verification was repeated after inspection and passed for every file.

Next milestone: a conservative lossless archive and deterministic MS export,
then exact comparison of all declared supported values and metadata. No
low-bit quantization, Dysco comparison, scientific products, CASA reopening or
imaging/calibration workflow has been tested yet.
