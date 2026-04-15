# =============================================================================
# _paths.R — portable path / helper stub for public release
# Replace values below via environment variables or local edit as needed.
# No individual-level data are distributed with this repository.
# =============================================================================
# DATA_DIR: directory containing the seven cohort-level harmonised files (not shipped)
# TAB_DIR : directory containing intermediate cohort-specific result CSVs (not shipped)
# FIG_DIR : output directory for figures
DATA_DIR <- Sys.getenv("STUDY24_DATA_DIR", "data")
TAB_DIR  <- Sys.getenv("STUDY24_TAB_DIR",  "aggregates")
FIG_DIR  <- Sys.getenv("STUDY24_FIG_DIR",  "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

log_step <- function(msg, tag = "study24") {
  cat(sprintf("[%s] %s  %s\n", format(Sys.time(), "%H:%M:%S"), tag, msg))
}
