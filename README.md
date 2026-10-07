<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Reliability-targeted IRT simulation: replication package

**JoonHo Lee** · [jlee296@ua.edu](mailto:jlee296@ua.edu)

This repository accompanies *Reliability-Targeted Simulation of Item Response Data: Solving the Inverse Design Problem*. It reproduces the validation studies and displays used in manuscript v6.3. We ask a practical question: can we specify an information level before generating responses, instead of letting an item-generation recipe determine that level implicitly?

Start with the worked example, then reproduce the solver and fresh-draw checks. The [reader's guide](book/index.qmd) explains what each command computes and how to interpret its output. The [reproduction map](manifest/reproduction-map.csv) connects the paper's exhibits to their inputs and scripts.

## Start here

Open `IRTsimrel-paper-replication.Rproj`, or open a terminal in this directory. Use **R 4.3 or newer**; this release was checked with **R 4.6.0**. Python 3 is needed only to rebuild the LaTeX tables. Quarto is optional for reading the guide as a website.

```sh
Rscript code/install_dependencies.R
Rscript run_all.R verify
Rscript run_all.R example
```

The example targets .85 on a 20-item 2PL form. It should return a multiplier near **1.277857**, a solver-node index near **.850000**, and an independent evaluation near **.849322**. That small difference is part of the example: a numerical solution on one set of integration nodes need not equal an evaluation on fresh draws.

Next, rebuild the displays from the distributed results:

```sh
Rscript run_all.R quick
quarto render book
```

Figures appear in `output/figures/`, recalculated tables in `output/tables/`, and the guide in `book/_book/index.html`. The quick command rebuilds 18 numerical figures and 12 LaTeX tables. Figure 3 is a conceptual diagram; its TeX source is included in `verification/expected/`. Other approved table fragments are included as references, not described as recalculated results.

## Run a small validation study

```sh
Rscript run_all.R calibration --workers 2
Rscript run_all.R summarize --study calibration
Rscript run_all.R nodes --workers 2
Rscript run_all.R summarize --study nodes
```

The first command recalibrates 16 parametric conditions across all four item models and all four latent shapes. The node study compares 500 with 20,000 integration nodes and evaluates each returned form on 200,000 fresh draws. These runs retain the paper's numerical settings and seeds. They select fewer conditions; they do not relax the solver tolerance.

All commands default to `--mode smoke --pool parametric --workers 1`. Windows users should keep `--workers 1`; parallel execution uses forked R processes on macOS and Linux.

## Reproduce the full validation designs

The public IRW difficulty pool is included with its source and license. Verify it, then run:

```sh
Rscript data/fetch_irw_pool.R
Rscript run_all.R feasibility --mode full --pool both --workers 4
Rscript run_all.R calibration --mode full --pool both --workers 4
Rscript run_all.R algorithms --mode full --pool both --workers 4
Rscript run_all.R superpopulation --mode full --pool both --workers 4
Rscript run_all.R nodes --mode full --workers 4
Rscript run_all.R steps --mode full --workers 4
Rscript run_all.R sampling --mode full --pool both --workers 4
Rscript run_all.R recovery --mode full --workers 4
```

Full runs can take hours. Read [the validation guide](book/03-validation.qmd) and [the recorded compute budget](data/precomputed/tables/T_compute_budget.csv) before starting. Each study writes one checkpoint per condition and can resume after interruption. A changed design, implementation, package, or R version invalidates old checkpoints rather than silently reusing them.

`algorithms` runs EQC, warm-start and cold-start stochastic calibration, and the inverse-information target where defined. `superpopulation` uses a separate package version and targets the transform of expected information across forms. These are different estimands; they should not be pooled into one success rate.

Additional commands cover `recipes`, `treatment`, `treatment-b`, `dif`, and `scores`. Their scope and limitations are explained in [the supplementary guide](book/04-applications.qmd). In particular, `treatment-b` returns initial model fits; the paper's later boundary/retry analysis is distributed as retained evidence.

## What is included

| Directory | Contents |
|---|---|
| `code/01_simulation/` | Designs and executable study entry points |
| `code/lib/` | Information calculations, calibration adapters, model fits, and checkpoint runner |
| `code/02_analysis/` | Figure builders, table calculations, and summaries of new runs |
| `config/` | Published factors, numerical settings, recipes, and external-input checksum |
| `data/irw/` | Public item-difficulty pool, source, and original MIT notice |
| `data/software/` | Frozen IRTsimrel 0.3.0 and 0.3.1 source packages |
| `data/precomputed/` | Retained numerical evidence used by the paper |
| `manifest/` | File checksums, authorship, and paper-to-code mapping |
| `verification/` | Expected values and final display fragments |
| `book/` | Step-by-step Quarto guide |
| `output/` | New results; excluded from version control |

The full response-level simulations are **not** run by `quick` or `verify`. The release checks and their actual coverage are recorded in [verification/release-checks.md](verification/release-checks.md). See [reproduction scope](docs/reproduction-scope.md) for retained-only analyses, [data access](data/README.md), and [software environment](environment/README.md).

## Citation and license

Cite the accompanying paper when using this work; repository citation metadata is in `CITATION.cff`. No DOI is assigned in this package. Code and original derived artifacts are released under the MIT license in `LICENSE`. External sources retain their own attribution and terms; see [third-party sources](docs/third-party-sources.md).

Author and maintainer: **JoonHo Lee (jlee296@ua.edu)**.
