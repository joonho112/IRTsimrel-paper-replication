<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Data and package inputs

`precomputed/results/` contains the retained study outputs. `precomputed/tables/` contains derived summaries used by manuscript displays. `precomputed/display/` contains the final plotting inputs. These files are read-only inputs to the quick route. New computations belong under `output/`.

The retained files were assembled from the versioned study outputs and checked against the manuscript's inputs. Numerical values, identifiers, and seeds are preserved. Historical names such as `rho_psd`, `w_bar`, `s2`, and `s6` remain in column names and join keys; the guide translates them into the manuscript's terminology.

## IRW difficulty pool

`data/irw/irw_diff_pool.csv` is included so that all difficulty-pool conditions can run without an additional download. It contains 6,143 item estimates from 145 instruments, exported from `irw::diff_long` at commit `cc96e459448ee9d268eb233881e949c1151cbdb6`. It does not contain individual response records.

The IRW package authors compiled the underlying data. Their original MIT copyright and permission notice is preserved in `data/irw/LICENSE-MIT.txt`; source details are in the adjacent README. The replication package's authorship does not replace that attribution.

Study scripts verify the CSV against `config/external_inputs.csv`. Run `Rscript data/fetch_irw_pool.R` to check the bundled input. To independently recover it from the public source archive, use `--download`; the recovered copy goes to `data/external/`, which remains excluded from version control.

The IRW supplies empirical difficulties, not a joint empirical distribution of every item parameter. Discrimination and guessing follow the stated parametric distributions.

## Authorship and file integrity

All repository-authored files are attributed to JoonHo Lee (jlee296@ua.edu). Text source files carry an author header. `manifest/files.csv` records the author and SHA-256 for every distributed artifact, including CSV, compressed data, and package archives. For the bundled IRW input and its original license, the manifest additionally records the upstream source, copyright holder, and license; the `author` field identifies the replication-package author. This keeps machine-readable data unchanged while providing file-level attribution. Citations and third-party attribution remain intact.

See `docs/data-dictionary.md` for important fields and `docs/third-party-sources.md` for external sources.
