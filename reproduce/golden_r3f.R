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
# the last digits and would then need its own provenance note.
#
# The R3F package sources are expected at
# ../OpenRelTime-参考软件/R3F/R; pass the directory as the second argument to
# use a different checkout.
#
# Usage:
#   Rscript reproduce/golden_r3f.R <project_root> <R3F_source_dir>

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else file.path(dirname(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = 1)))), "")
root <- normalizePath(root, mustWork = TRUE)
r3f <- if (length(args) >= 2) args[[2]] else file.path(dirname(root), "OpenRelTime-参考软件", "R3F", "R")
r3f <- normalizePath(r3f, mustWork = TRUE)

suppressMessages({
  library(ape); library(phangorn); library(R.utils); library(FNN)
})
invisible(lapply(file.path(r3f, c("rrf_rates.R", "rrf_times.R", "rrf_rates_times.R",
                                  "corrtest.R", "ddbd.R", "tree2table.R")), source))

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
