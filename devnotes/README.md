# Development notes

- [Design and experiment specification](design.md)
- [Initial source inspection findings](initial-findings.md)
- [Initial measured archive and CASA results](experiment-results.md)
- [Geometry/antenna/frequency grouping and WSClean validation](grouping-wsclean.md)
- [Opt-in destructive flag policy (`burn_flag`)](burn-flag.md)
- [Default 4/8/12/16-bit precision ladder and grouping settings](precision-ladder.md)
- [Repository publishing and workspace arrangement](publishing.md)

Machine-readable measurements are in `../reports/`. Large generated products
are local under `../results/` and excluded from Git. The user-facing archive
format contract remains in `../docs/format.md`.
