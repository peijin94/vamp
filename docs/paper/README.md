# Paper

- [PDF](varquantMS.pdf)
- [Markdown source](varquantMS.md)
- [Figures](figs/)
- [Preserved initial development manuscript](work_note_20261025_init_work.md)

Build from the repository root:

```bash
python -m pip install -r docs/paper/requirements.txt
python scripts/compile_paper.py
```

Pandoc converts the Markdown and equations; Typst produces the PDF. The PDF uses
vector diagrams/charts where available, retains the measured image panels, and
links repository references to their public GitHub pages. The build uses local
figures and fonts and does not download external templates or codebooks.

The current revision includes the three-epsilon Stokes I study, the additional
I/V validation at epsilon-D 0.18, complete size accounting and repeated full-MS
conversion timings. To refresh its measured tables and figures from saved reports:

```bash
PYTHONPATH=src:scripts python scripts/update_calibrated_paper.py
python scripts/compile_paper.py
```

Local FITS inputs are needed for image regeneration. With only the repository,
use `--skip-images` to retain the committed image panels. The underlying imaging
and timing experiments are not rerun by these commands. Source measurements are
in `reports/calibrated-epsd-dense.json`, `calibrated-performance.json`,
`calibrated-shared-storage.json` and `calibrated-paper-storage.json`.
