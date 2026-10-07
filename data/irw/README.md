<!-- Author: JoonHo Lee (jlee296@ua.edu) -->

# Public IRW difficulty input

This directory contains the item-difficulty estimates used by the replication studies: 6,143 items from 145 instruments. The CSV is an export of `diff_long` from `irw` 1.0.0. It contains item-level estimates and standard errors, not individual response records.

The data were compiled by the IRW package authors. Their copyright and MIT license are reproduced unchanged in `LICENSE-MIT.txt`. JoonHo Lee (jlee296@ua.edu) prepared this replication package and the CSV export; the underlying IRW data retain their original attribution.

Source: [IRW package at commit cc96e459448ee9d268eb233881e949c1151cbdb6](https://github.com/itemresponsewarehouse/Rpkg/tree/cc96e459448ee9d268eb233881e949c1151cbdb6).

The expected CSV SHA-256 is `07dce30787985ac19f3da6b576ebd52242688afd6f959108496466547797fc07`, also recorded in `config/external_inputs.csv`. Study scripts verify it before use. No account or download is required to use the bundled file. `Rscript data/fetch_irw_pool.R --download` can independently recover the input from the pinned public archive.
