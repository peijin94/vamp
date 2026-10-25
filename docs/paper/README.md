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
