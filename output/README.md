<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Results from this checkout

Commands write new figures, tables, verification reports, and simulation results here. These files are excluded from version control. The retained paper inputs remain under `data/precomputed/`.

Each study has `runs/<study>/<mode>/results.csv`, `design.csv`, a session record, and resumable `checkpoints/`. Selected validation studies also write `comparison.csv` against the retained output. `summarize` adds `summary.csv` for a new run.

Move an existing study output directory before changing its design or numerical implementation. This preserves the old run and avoids mixing incompatible checkpoints.
