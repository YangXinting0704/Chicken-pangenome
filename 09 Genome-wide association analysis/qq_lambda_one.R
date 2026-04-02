args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]

dat <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)

p <- dat$p_wald
p <- p[!is.na(p) & p > 0 & p < 1]

if (length(p) < 1000) {
  quit(save = "no")
}

# lambda GC
chisq <- qchisq(1 - p, df = 1)
lambda_gc <- median(chisq) / qchisq(0.5, 1)

out_jpg <- sub("\\.assoc\\.txt$", ".QQ.jpg", infile)
out_txt <- sub("\\.assoc\\.txt$", ".lambdaGC.txt", infile)

# QQ（jpg）
jpeg(out_jpg, width = 1200, height = 1200, res = 200)

obs <- -log10(sort(p))
exp <- -log10(ppoints(length(p)))

plot(
  exp, obs,
  pch = 16, cex = 0.4,
  xlab = "Expected -log10(P)",
  ylab = "Observed -log10(P)",
  main = paste0(
    basename(infile),
    "\nGC inflation factor = ",
    round(lambda_gc, 3)
  )
)

abline(0, 1, col = "red", lwd = 2)
dev.off()

# Output the lambda result
write.table(
  data.frame(
    trait = basename(infile),
    lambda_GC = round(lambda_gc, 4),
    n_variants = length(p)
  ),
  out_txt,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
