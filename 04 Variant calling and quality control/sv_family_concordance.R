#!/usr/bin/env Rscript
suppressMessages({
  library(parallel)
})

### ============ parameter ============
ped_file  <- "pedigree.txt"
geno_file <- "genotypes.tsv"
out_file  <- "SV_family_concordance.tsv"

ncore      <- 32
block_size <- 2000

### ============ pedigree ============
ped <- read.table(ped_file, header=TRUE, sep="\t",
                  stringsAsFactors=FALSE, na.strings=c("NA","","."))
colnames(ped) <- c("ID","Dam","Sire")

### ============ genotype ============
geno_raw <- read.table(geno_file, header=TRUE, sep="\t",
                       stringsAsFactors=FALSE, check.names=FALSE)

sv_id <- geno_raw[[1]]
geno_raw <- geno_raw[,-1]
sample_id <- colnames(geno_raw)

### ============ pedigree alignment ============
ped <- ped[ped$ID %in% sample_id, ]
ped <- ped[match(sample_id, ped$ID), ]
stopifnot(all(ped$ID == sample_id))

nsample <- length(sample_id)
nsv <- length(sv_id)

### ============ Convert GT to an integer matrix ============
# 0/0 -> 0, 0/1 -> 1, 1/1 -> 2, ./. -> NA
convert_gt <- function(x) {
  x[x %in% c("./.", ".")] <- NA
  x[x == "0/0"] <- 0
  x[x == "0/1"] <- 1
  x[x == "1/1"] <- 2
  as.integer(x)
}

geno_mat <- matrix(NA_integer_, nrow=nsv, ncol=nsample)
for (i in seq_len(nsample)) {
  geno_mat[,i] <- convert_gt(geno_raw[[i]])
}
rm(geno_raw)
storage.mode(geno_mat) <- "integer"

### ============ Build the sample pair ============
idx <- seq_len(nsample)
all_pairs <- t(combn(idx, 2))

dam <- ped$Dam
sire <- ped$Sire

same_dam  <- dam[all_pairs[,1]] == dam[all_pairs[,2]] & !is.na(dam[all_pairs[,1]])
same_sire <- sire[all_pairs[,1]] == sire[all_pairs[,2]] & !is.na(sire[all_pairs[,1]])

within_idx <- which(same_dam | same_sire)
bg_idx     <- which(!(same_dam | same_sire))

within_A <- all_pairs[within_idx,1]
within_B <- all_pairs[within_idx,2]
bg_A     <- all_pairs[bg_idx,1]
bg_B     <- all_pairs[bg_idx,2]

rm(all_pairs)

### ============ block processing function ============
process_block <- function(sv_idx, geno_mat,
                          within_A, within_B,
                          bg_A, bg_B) {

  res <- matrix(NA_real_, nrow=length(sv_idx), ncol=5)

  for (i in seq_along(sv_idx)) {
    g <- geno_mat[sv_idx[i], ]

    # within-family
    ga <- g[within_A]; gb <- g[within_B]
    ok <- !(is.na(ga) | is.na(gb))
    w_pairs <- sum(ok)
    w_conc  <- sum(ga[ok] == gb[ok])

    # background
    ga <- g[bg_A]; gb <- g[bg_B]
    ok <- !(is.na(ga) | is.na(gb))
    b_pairs <- sum(ok)
    b_conc  <- sum(ga[ok] == gb[ok])

    pval <- NA
    if (w_pairs > 0 && b_pairs > 0) {
      mat <- matrix(c(w_conc,
                      w_pairs - w_conc,
                      b_conc,
                      b_pairs - b_conc),
                    nrow=2, byrow=TRUE)
      pval <- fisher.test(mat)$p.value
    }

    res[i,] <- c(w_pairs, w_conc, b_pairs, b_conc, pval)
  }
  res
}

### ============ parallel running ============
blocks <- split(seq_len(nsv),
                ceiling(seq_len(nsv) / block_size))

cl <- makeCluster(ncore)
clusterExport(cl,
  c("geno_mat",
    "within_A","within_B",
    "bg_A","bg_B",
    "process_block"),
  envir=environment()
)

res_list <- parLapply(cl, blocks, process_block,
                      geno_mat=geno_mat,
                      within_A=within_A,
                      within_B=within_B,
                      bg_A=bg_A,
                      bg_B=bg_B)
stopCluster(cl)

res <- do.call(rbind, res_list)
colnames(res) <- c("within_pairs","within_concord",
                   "bg_pairs","bg_concord","fisher_p")

### ============ Summary output ============
out <- data.frame(
  SV_ID = sv_id,
  within_pairs   = res[,1],
  within_concord = res[,2],
  within_rate    = res[,2] / res[,1],
  bg_pairs       = res[,3],
  bg_concord     = res[,4],
  bg_rate        = res[,4] / res[,3],
  fisher_p       = res[,5],
  FDR            = p.adjust(res[,5], "fdr"),
  stringsAsFactors = FALSE
)

write.table(out, out_file, sep="\t",
            quote=FALSE, row.names=FALSE)

message("Finished! Output: ", out_file)
