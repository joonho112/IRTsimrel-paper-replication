# Author: JoonHo Lee (jlee296@ua.edu)
# License: MIT

# Recreate paper displays from the retained results; no simulations run here.
require_packages(c("dplyr", "tidyr", "readr", "ggplot2", "scales", "patchwork"))
paths <- list.files(repo_path("code", "02_analysis", "figures"),
  pattern = "\\.R$", full.names = TRUE
)
paths <- paths[basename(paths) != "fig_common.R"]
for (path in paths) {
  message("Figure script: ", basename(path))
  status <- system2(file.path(R.home("bin"), "Rscript"), shQuote(path))
  if (status != 0L) stop("Figure script failed: ", basename(path))
}
python <- Sys.which("python3")
if (!nzchar(python)) stop("Python 3 is required for the table scripts.")
for (name in c("tables.py", "tables_v5.py")) {
  path <- repo_path("code", "02_analysis", "tables", name)
  status <- system2(python, shQuote(path))
  if (status != 0L) stop("Table script failed: ", name)
}
record_session(repo_path("output", "display-session-info.txt"))
message("Displays written to output/figures and output/tables.")
