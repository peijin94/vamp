# Development notes

- [Design and experiment specification](design.md)
- [Initial source inspection findings](initial-findings.md)
- [Initial measured archive and CASA results](experiment-results.md)
- [Geometry/antenna/frequency grouping and WSClean validation](grouping-wsclean.md)
- [Opt-in destructive flag policy (`burn_flag`)](burn-flag.md)
- [Repository publishing and workspace arrangement](publishing.md)

Machine-readable measurements are in `../reports/`. Large generated products
are local under `../results/` and excluded from Git. The user-facing archive
format contract remains in `../docs/format.md`.
