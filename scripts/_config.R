#!/usr/bin/env Rscript
# ============================================================================
# Single configuration loader for every R script in this repository.
# Reads config/config.ini so that no script contains a machine-specific path.
#
#   source("scripts/_config.R")
#   CFG$results_dir     # absolute path to results/
#   CFG$samples         # c("MOCK_1", ..., "PR8_3")
#   CFG$is_mock         # c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
# ============================================================================

.resolve_root <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) return(normalizePath(file.path(dirname(f[1]), ".."), mustWork = FALSE))
  # fall back: assume the working directory is the repository root or scripts/
  wd <- normalizePath(getwd(), mustWork = FALSE)
  if (basename(wd) == "scripts") return(dirname(wd))
  wd
}

.ROOT <- .resolve_root()
.CFG_PATH <- file.path(.ROOT, "config", "config.ini")
if (!file.exists(.CFG_PATH)) stop("STOP: configuration file not found: ", .CFG_PATH)

.read_ini <- function(path, section) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[!grepl("^\\s*#", ln) & !grepl("^\\s*;", ln)]
  sec <- NA_character_
  out <- list()
  for (l in ln) {
    if (grepl("^\\s*\\[", l)) { sec <- sub("^\\s*\\[([^]]+)\\].*$", "\\1", l); next }
    if (is.na(sec) || sec != section) next
    if (!grepl("=", l)) next
    k <- trimws(sub("=.*$", "", l)); v <- trimws(sub("^[^=]*=", "", l))
    v <- trimws(sub("\\s+#.*$", "", v))   # drop a trailing inline comment
    if (nzchar(k)) out[[k]] <- v
  }
  out
}

.sp <- function(x) trimws(strsplit(x, ",")[[1]])
.abs <- function(root, p) if (grepl("^(/|[A-Za-z]:)", p)) p else file.path(root, p)

p <- .read_ini(.CFG_PATH, "paths")
s <- .read_ini(.CFG_PATH, "samples")
a <- .read_ini(.CFG_PATH, "analysis")

data_root <- .abs(.ROOT, if (is.null(p$data_root)) "data" else p$data_root)

CFG <- list(
  root          = .ROOT,
  data_root     = data_root,
  raw_dir       = .abs(data_root, if (is.null(p$raw_dir)) "raw" else p$raw_dir),
  processed_dir = .abs(data_root, if (is.null(p$processed_dir)) "processed" else p$processed_dir),
  interim_dir   = .abs(data_root, if (is.null(p$interim_dir)) "interim" else p$interim_dir),
  results_dir   = .abs(.ROOT, if (is.null(p$results_dir)) "results" else p$results_dir),
  figdata_dir   = .abs(.ROOT, if (is.null(p$figdata_dir)) "figdata" else p$figdata_dir),
  figures_dir   = .abs(.ROOT, if (is.null(p$figures_dir)) "figures" else p$figures_dir),
  tables_dir    = .abs(.ROOT, if (is.null(p$tables_dir)) "tables" else p$tables_dir)
)

CFG$samples <- .sp(s$names)
CFG$groups  <- .sp(s$groups)
CFG$reference <- if (is.null(s$reference)) "MOCK" else s$reference
CFG$case      <- if (is.null(s$case)) "PR8" else s$case
CFG$is_mock   <- CFG$groups == CFG$reference

CFG$fdr_cut        <- as.numeric(if (is.null(a$fdr_cut)) 0.05 else a$fdr_cut)
CFG$log2fc_cut     <- as.numeric(if (is.null(a$log2fc_cut)) 1.0 else a$log2fc_cut)
CFG$pctmod_pp_cut  <- as.numeric(if (is.null(a$pctmod_pp_cut)) 5.0 else a$pctmod_pp_cut)
CFG$log2fc_flat    <- as.numeric(if (is.null(a$log2fc_flat)) 0.5 else a$log2fc_flat)
CFG$pctmod_pp_flat <- as.numeric(if (is.null(a$pctmod_pp_flat)) 10.0 else a$pctmod_pp_flat)

CFG$result   <- function(name, create = FALSE) {
  if (create) dir.create(CFG$results_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(CFG$results_dir, name)
}
CFG$interim  <- function(name, create = FALSE) {
  if (create) dir.create(CFG$interim_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(CFG$interim_dir, name)
}
CFG$figdata  <- function(name, create = FALSE) {
  if (create) dir.create(CFG$figdata_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(CFG$figdata_dir, name)
}
CFG$figure   <- function(name, create = FALSE) {
  if (create) dir.create(CFG$figures_dir, showWarnings = FALSE, recursive = TRUE)
  file.path(CFG$figures_dir, name)
}
CFG$report <- function() {
  cat(sprintf("repository root : %s\ndata root       : %s\nsamples         : %s\n",
              CFG$root, CFG$data_root,
              paste(sprintf("%s(%s)", CFG$samples, CFG$groups), collapse = ", ")))
  invisible(CFG)
}

# Alias: the workflow scripts refer to `cfg`, the list above is `CFG`.
# R is case-sensitive, so both names refer to the same configuration.
cfg <- CFG
