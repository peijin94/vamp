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

The current revision includes three-epsilon Stokes I and V comparisons, the
earlier I/V performance-validation runs, complete size accounting and a phase
breakdown of the median full-MS conversion, plus the matched Dysco comparison
and two additional epsilon measurements. To refresh its measured tables and
figures from saved reports:

```bash
PYTHONPATH=src:scripts python scripts/update_calibrated_paper.py
PYTHONPATH=src:scripts python scripts/build_dysco_paper_comparison.py render
PYTHONPATH=src:scripts python scripts/update_dysco_paper.py
python scripts/compile_paper.py
```

Local FITS inputs are needed for image regeneration. With only the repository,
use `--skip-images` to retain the committed image panels. The underlying imaging
and timing experiments are not rerun by these commands. Source measurements are
in `reports/calibrated-epsd-dense.json`, `calibrated-performance.json`,
`calibrated-shared-storage.json`, `calibrated-paper-storage.json`,
`calibrated-epsd-v.json` and `calibrated-runtime-breakdown.json`.
The comparison figures and table use `dysco-comparison.json` and
`dysco-paper-extra.json`.

With the recorded local run outputs, `PYTHONPATH=src:scripts python
scripts/report_paper_followup.py` audits the six new V runs and the existing
median timing run. It does not perform new imaging or benchmarks.
