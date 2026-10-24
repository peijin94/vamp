# Quantization difference relative to the reference residual RMS

For the [selected 25% payload profile](fixed25.md), the restored-image difference
has **0.482–0.508% of the reference residual-map RMS**, approximately 197–208
times smaller. This compares the existing two CLEAN runs; no new imaging was run.

All measurements use the same 2,245,726 pixels above 10° at all integrations,
with matching grids, restoring beams and Jy/beam units. RMS means
`sqrt(mean(x**2))`, without subtracting the mean, clipping or masking sources.

| Quantity | Run 1 | Run 2 |
| --- | ---: | ---: |
| Reference residual RMS (Jy/beam) | 1.031244811 | 1.031244811 |
| Quantized − reference restored-image RMS (mJy/beam) | 5.237975 | 4.966552 |
| Restored difference / reference residual RMS | 0.507927% | 0.481607% |
| Quantized residual-map RMS (Jy/beam) | 1.031052210 | 1.031277269 |
| Change in residual-map RMS | −0.018677% | +0.003147% |
| Quantized − reference residual-map RMS (mJy/beam) | 9.900162 | 10.308579 |

The primary ratio is

$$\frac{\sqrt{\langle(I_{\rm quant}-I_{\rm ref})^2\rangle}}
{\sqrt{\langle R_{\rm ref}^2\rangle}},$$

where I denotes the restored image and R the CLEAN residual map. It differs
from subtracting two residual maps, and from the earlier 0.262–0.276% metric
whose denominator was the **reference restored-image** RMS.

The reference residual mean is 0.160266 Jy/beam and its mean-subtracted standard
deviation is 1.018715 Jy/beam. Residual sky and deconvolution structure contribute
to this full-mask RMS; it is not a validated thermal-noise estimate. Restored
differences also include independent CLEAN paths, so no independent-error or
quadrature-addition assumption is made.

[Machine-readable measurements](../reports/calibrated-fixed25-residual-rms.json)
include input FITS hashes, the mask, definitions, RMS, mean and standard deviation.
The FITS hashes and primary RMS values match the existing experiment report.
