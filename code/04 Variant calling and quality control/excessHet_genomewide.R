#!/usr/bin/env Rscript

suppressMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
infile  <- args[1]
outfile <- args[2]

cat("Reading input...\n")
dat <- fread(infile, showProgress = TRUE)

setnames(dat,
         c("V1","V2","V3","V4","V5","V6","V7","V8"),
         c("CHROM","POS","ID","called","ref","het","hom","miss"))

cat("Filtering invalid sites...\n")
dat <- dat[called > 0]

dat[, p_alt := (2*hom + het) / (2*called)]
dat[, Hobs  := het]
dat[, q_exp := 2 * p_alt * (1 - p_alt)]
dat[, Hexp  := q_exp * called]

cat("Computing binomial p-values...\n")
dat[, p_binom := 1.0]

idx <- which(dat$q_exp > 0 & dat$q_exp < 1)

dat$p_binom[idx] <- 1 - pbinom(dat$het[idx] - 1,
                               size = dat$called[idx],
                               prob = dat$q_exp[idx])

cat("Multiple testing correction...\n")
dat[, p_bonf_binom := p.adjust(p_binom, method="bonferroni")]
dat[, p_bh_binom   := p.adjust(p_binom, method="BH")]

dat[, signif_bonf_binom := p_bonf_binom < 0.05]
dat[, signif_bh_binom   := p_bh_binom   < 0.05]

cat("Selecting significant sites...\n")
res <- dat[signif_bonf_binom | signif_bh_binom]

cat("Writing result...\n")
fwrite(
  res[, .(CHROM, POS, ID,
          called, ref, het, hom, miss,
          p_alt, Hobs, Hexp,
          p_binom, p_bonf_binom, p_bh_binom,
          signif_bonf_binom, signif_bh_binom)],
  file = outfile,
  sep = "\t"
)

cat("Done. Significant sites:", nrow(res), "\n")
