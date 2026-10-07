<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# External sources and attribution

Repository author and maintainer: JoonHo Lee (jlee296@ua.edu).

- **IRTsimrel:** source versions 0.3.0 and 0.3.1 are distributed in `data/software/` with their original package metadata and license. They are installed locally to preserve the study's numerical implementation.
- **Item Response Warehouse:** [public R package repository](https://github.com/itemresponsewarehouse/Rpkg), pinned at commit `cc96e459448ee9d268eb233881e949c1151cbdb6`. The package includes the `diff_long` difficulty estimates as a CSV in `data/irw/`. The original MIT copyright and permission notice is reproduced there unchanged. Individual response records and the external package code are not included.
- **Gilbert, Kim, and Miratrix:** the item-level treatment-effect design is discussed in the accompanying manuscript. The external replication deposit is [Dataverse DOI 10.7910/DVN/QARRYT](https://doi.org/10.7910/DVN/QARRYT), with CC BY-NC-SA 4.0 terms. This repository does not redistribute its code or raw deposit. The retained fidelity report records this study's numerical comparison with that source.
- **R packages:** their authorship and licenses are in their installed package metadata. The repository's MIT license does not replace dependency licenses.

The original public IRT simulation repository is [reliability-targeted-irt-simulation](https://github.com/joonho112/reliability-targeted-irt-simulation). The current paper uses later numerical implementations and a larger validation design; this package should not be run with the old repository's defaults.
