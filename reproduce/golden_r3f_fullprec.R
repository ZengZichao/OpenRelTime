#!/usr/bin/env Rscript
# Capture the R3F CorrScore at full precision.
#
# Why this exists: R3F's corrtest() writes the score with format(score,
# digits = 5) (corrtest.R:570-583), so data/golden/r3f/example_sr0_corrtest.txt
# records "score = 0.9996" and example_sr50_corrtest.txt records 0.99959. Those
# five significant digits are all the shipped golden file carries, which is too
# coarse to serve as an error budget for a faithful-implementation claim.
#
# This script re-sources the R3F functions with that one formatting constant
# widened to digits = 17, so the same code path prints the unrounded double.
# Nothing in the numerical path is touched, only the write() formatting.
#
# Usage:  Rscript reproduce/golden_r3f_fullprec.R <R3F_source_dir> <tree.nwk>
# Verified with the micromamba env r-4.5.3 (R 4.5.3, ape 5.8.1, phangorn 2.11.1,
# R.utils 2.13.0, FNN 1.1.4.1, RColorBrewer 1.1-3).

args <- commandArgs(trailingOnly = TRUE)
r3f.src <- args[[1]]
tree.name <- args[[2]]
og <- c("Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus")

work <- file.path(tempdir(), "r3f_fullprec")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
files <- list.files(r3f.src, pattern = "\\.R$", full.names = TRUE)
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  txt <- gsub("format(score, digits = 5)", "format(score, digits = 17)", txt, fixed = TRUE)
  writeLines(txt, file.path(work, basename(f)))
}

suppressMessages({library(ape); library(phangorn); library(R.utils); library(FNN)})
for (f in c("rrf_rates.R", "rrf_times.R", "rrf_rates_times.R", "corrtest.R", "ddbd.R", "tree2table.R")) {
  source(file.path(work, f))
}
owd <- getwd(); setwd(work)
invisible(suppressWarnings(corrtest(tree.name = tree.name, type = "NEWICK", outgroup = og,
                                    sister.resample = 0, filename = "fp_sr0")))
set.seed(42)
invisible(suppressWarnings(corrtest(tree.name = tree.name, type = "NEWICK", outgroup = og,
                                    sister.resample = 50, filename = "fp_sr50")))
setwd(owd)
cat(readLines(file.path(work, "fp_sr0_corrtest.txt"))[1], "\n")
cat(readLines(file.path(work, "fp_sr50_corrtest.txt"))[1], "\n")
