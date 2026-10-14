#!/usr/bin/env Rscript
# Regenerate the R3F golden outputs in data/golden/r3f/.
#
# Environment of the golden run (exact versions, as recorded in
# data/PROVENANCE.md; micromamba env "r-4.5.3"):
#   R             4.5.3
#   ape           5.8.1     read.tree / root / drop.tip / write.tree
#   phangorn      2.11.1
#   R.utils       2.13.0
#   FNN           1.1.4.1   KL.divergence(), used by ddbd(measure = "KL")
#   RColorBrewer  1.1-3     palette used inside the sourced R3F functions
# These versions are recorded, not enforced: nothing here aborts when a newer
# package is installed, so a regeneration on a different stack can differ in
# the last digits and would then need its own provenance note.  Nothing pins
# them either -- this directory carries no `renv.lock`, and neither does R3F at
# the commit named below (it ships renv/activate.R and renv/settings.dcf only),
# so the list above is the whole record of the R stack.
#
# Where R3F comes from, in this order:
#   1. the second positional argument, or the R3F_SOURCE_DIR environment
#      variable: the directory holding R3F's `R/` *sources*, which is what the
#      shipped goldens were produced from (files are sourced, not installed);
#   2. failing that, an installed `R3F` package on the library path;
#   3. otherwise the script stops and prints the clone of the exact commit the
#      article cites:
#        git clone https://github.com/cathyqqtao/R3F
#        git -C R3F checkout 6c30d09b650c83f9b34e10387ac57c4ea60a1624
#        Rscript reproduce/golden_r3f.R . "$PWD/R3F/R"
#   There is deliberately no default path into the author's private working
#   directory: a reader must be able to run this script without that directory
#   existing.
#
# Usage:
#   Rscript reproduce/golden_r3f.R [project_root] [R3F_source_dir]
#   Rscript reproduce/golden_r3f.R --help
#   Rscript reproduce/golden_r3f.R [project_root] [R3F_source_dir] --dry-run

args <- commandArgs(trailingOnly = TRUE)
flags <- args[grepl("^-", args)]
positional <- args[!grepl("^-", args)]

if (any(flags %in% c("--help", "-h"))) {
  cat(paste0(
    "Usage: Rscript reproduce/golden_r3f.R [project_root] [R3F_source_dir] ",
    "[--dry-run]\n\n",
    "  project_root     repository to write data/golden/r3f/ into; defaults to\n",
    "                   the parent of the directory holding this script\n",
    "  R3F_source_dir   directory of R3F's R/ sources (rrf_rates.R, corrtest.R,\n",
    "                   ddbd.R, ...).  Falls back to $R3F_SOURCE_DIR, then to an\n",
    "                   installed R3F package; with neither, the script stops and\n",
    "                   prints how to clone commit\n",
    "                   6c30d09b650c83f9b34e10387ac57c4ea60a1624 of\n",
    "                   https://github.com/cathyqqtao/R3F\n",
    "  --dry-run        resolve everything and report what would be used, without\n",
    "                   running R3F or touching data/golden/\n"
  ))
  quit(status = 0)
}
dry_run <- any(flags == "--dry-run")

root <- if (length(positional) >= 1) {
  positional[[1]]
} else {
  dirname(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = 1))))
}
root <- normalizePath(root, mustWork = TRUE)

r3f_src <- if (length(positional) >= 2) {
  positional[[2]]
} else {
  Sys.getenv("R3F_SOURCE_DIR", "")
}
r3f_installed <- FALSE
if (nzchar(r3f_src)) {
  r3f_src <- normalizePath(r3f_src, mustWork = TRUE)
} else if (requireNamespace("R3F", quietly = TRUE)) {
  r3f_installed <- TRUE
} else {
  stop(
    "no R3F found.  Pass the directory of R3F's R/ sources as the second ",
    "argument, set R3F_SOURCE_DIR, install R3F, or clone the commit the ",
    "article cites:\n",
    "  git clone https://github.com/cathyqqtao/R3F\n",
    "  git -C R3F checkout 6c30d09b650c83f9b34e10387ac57c4ea60a1624\n",
    "  Rscript reproduce/golden_r3f.R ", root, ' "$PWD/R3F/R"',
    call. = FALSE
  )
}

if (dry_run) {
  cat("project root:      ", root, "\n")
  cat("R3F source:        ",
      if (r3f_installed) "installed R3F package" else r3f_src, "\n")
  cat("golden output dir: ", file.path(root, "data", "golden", "r3f"), "\n")
  cat("input tree:        ",
      file.path(root, "data", "examples", "example.nwk"), "\n")
  cat("dry run: nothing was executed and no file was written\n")
  quit(status = 0)
}

suppressMessages({
  library(ape); library(phangorn); library(R.utils); library(FNN)
})
if (r3f_installed) {
  suppressMessages(library(R3F))
} else {
  invisible(lapply(file.path(r3f_src, c("rrf_rates.R", "rrf_times.R",
                                        "rrf_rates_times.R", "corrtest.R",
                                        "ddbd.R", "tree2table.R")), source))
}

golden <- file.path(root, "data", "golden", "r3f")
dir.create(golden, showWarnings = FALSE, recursive = TRUE)
setwd(golden)
tree.name <- file.path(root, "data", "examples", "example.nwk")
og <- c("Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus")

invisible(suppressWarnings(rrf_rates(tree.name = tree.name, type = "NEWICK", outgroup = og, filename = "example")))
invisible(suppressWarnings(rrf_times(tree.name = tree.name, type = "NEWICK", outgroup = og, filename = "example")))
invisible(suppressWarnings(rrf_rates_times(tree.name = tree.name, type = "NEWICK", outgroup = og, filename = "example", plot = FALSE)))
invisible(suppressWarnings(corrtest(tree.name = tree.name, type = "NEWICK", outgroup = og, sister.resample = 0, filename = "example_sr0")))
set.seed(42)
invisible(suppressWarnings(corrtest(tree.name = tree.name, type = "NEWICK", outgroup = og, sister.resample = 50, filename = "example_sr50")))
invisible(suppressWarnings(ddbd(tree.name = tree.name, type = "NEWICK", outgroup = og, sampling.frac = 0,
                                anchor.node = 272, anchor.time = 1.85, measure = "SSE",
                                filename = "example_anchor", plot = FALSE)))
invisible(suppressWarnings(ddbd(tree.name = tree.name, type = "NEWICK", outgroup = og, sampling.frac = 0,
                                anchor.node = 0, anchor.time = 1, measure = "SSE",
                                filename = "example_noanchor", plot = FALSE)))
invisible(tree2table(tree.name = tree.name, type = "NEWICK", filename = "example_tree2table"))

t <- read.tree(tree.name)
t3 <- drop.tip(root(t, og, resolve.root = TRUE), og)
write.tree(t3, file = "example_rooted_dropped.nwk")

cat("golden outputs written to", golden, "\n")
