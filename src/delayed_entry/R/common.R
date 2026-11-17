# Shared by the run's R scripts: the run's settings, folders and files.
#
# Every clinical and statistical value comes from run_spec.json, which the Python side writes from the locked
# settings (delayed_entry.run). Nothing here is a study's value. Ported from the CRC BPC project's R/paths.R and
# the constants of its R/passes.R, R/bias.R and R/stability.R, with each fixed value replaced by the spec's.
#
#   source("common.R"); load_spec(run_dir)
CLOUD_MARKERS <- c("onedrive", "dropbox", "google drive", "googledrive", "icloud", "box sync")
PASSES <- c("naive", "corrected")
BAR_TOLERANCE <- 1e-9   # floating point: 0.9 - 1 is -0.0999999999999999778, still "at the bar"

# Return the folder, or stop if it looks cloud-synced: row-level files stay on local disk.
check_not_cloud <- function(p) {
  low <- tolower(gsub("\\\\", "/", p))
  hit <- CLOUD_MARKERS[vapply(CLOUD_MARKERS, grepl, logical(1), x = low, fixed = TRUE)]
  if (length(hit)) stop(p, " looks cloud-synced ('", hit[1], "'). Row-level data must stay on local disk only.")
  p
}

# A grouping: the column that holds the groups, the groups with the reference first, and the reference's name.
grouping_from <- function(g) list(column = g$column, levels = unlist(g$levels), reference = g$reference)

# A JSON object of names to texts as a named character vector (empty if none).
named_text <- function(x) {
  v <- unlist(x)
  if (is.null(v)) character() else v
}

# Read the spec and set the values every script uses.
load_spec <- function(run_dir) {
  path <- file.path(run_dir, "run_spec.json")
  s <- jsonlite::read_json(path, simplifyVector = TRUE)
  s$secondary <- jsonlite::read_json(path, simplifyVector = FALSE)$secondary   # one list per analysis
  SPEC <<- s
  DAYS_PER_YEAR <<- s$days_per_year
  DAYS_PER_MONTH <<- DAYS_PER_YEAR / 12
  YEARS <<- s$survival_at_years
  RMST_DAYS <<- s$horizon_years * DAYS_PER_YEAR
  MEDIAN_CAP_DAYS <<- s$median_limit_years * DAYS_PER_YEAR
  ALL_PATIENTS <<- s$all_patients
  BOOT_RESAMPLES <<- as.integer(s$resamples)
  BOOT_SEED <<- as.integer(s$seed)
  # The interval's two percentiles, rounded so that a 95% level gives exactly 2.5% and 97.5%, as written.
  PROBS <<- round(c(1 - s$interval_level, 1 + s$interval_level) / 2, 10)
  Z <<- stats::qnorm(PROBS[2])
  BARS <<- unlist(s$bars)
  STABLE_RISK_SET <<- as.integer(s$stable_risk_set)
  GENOTYPE <<- grouping_from(s$genotype)
  LABELS <<- s$labels
  invisible(s)
}

# A string written as SQL text.
sql_text <- function(s) paste0("'", gsub("'", "''", s, fixed = TRUE), "'")

connect <- function() DBI::dbConnect(duckdb::duckdb(shared_home = FALSE))

read_table <- function(con, run_dir, file) {
  DBI::dbGetQuery(con, sprintf("SELECT * FROM read_parquet(%s)", sql_text(file.path(run_dir, file))))
}

# Save a data frame as a Parquet file: exact doubles, and the column types kept.
write_table <- function(con, x, run_dir, file) {
  DBI::dbWriteTable(con, "result_table", x, overwrite = TRUE)
  path <- normalizePath(file.path(run_dir, file), winslash = "/", mustWork = FALSE)
  invisible(DBI::dbExecute(con, sprintf("COPY result_table TO %s (FORMAT parquet)", sql_text(path))))
}

# Run expr, keeping its warnings (a model fitted to a small resample can warn) instead of printing each one.
with_warnings <- function(expr) {
  warned <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warned)
}

# Package versions, for the run report.
versions <- function(pkgs) {
  c(R = paste(R.version$major, R.version$minor, sep = "."),
    vapply(pkgs, function(p) as.character(utils::packageVersion(p)), character(1)))
}
