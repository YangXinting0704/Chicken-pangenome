#!/usr/bin/env Rscript
#transcriptome splicing map cluster to genes
suppressMessages(library(argparser, quietly=TRUE))
suppressMessages(library(dplyr, quietly=TRUE))
suppressMessages(library(foreach, quietly=TRUE))
# leafcutter functions:
get_intron_meta <- function(introns) {
  intron_meta <- do.call(rbind, strsplit(introns,":"))
  colnames(intron_meta) <- c("chr","start","end","clu")
  intron_meta <- as.data.frame(intron_meta, stringsAsFactors=FALSE)
  intron_meta$start <- as.numeric(intron_meta$start)
  intron_meta$end <- as.numeric(intron_meta$end)
  intron_meta
}
#' Work out which gene each cluster belongs to. Note the chromosome names used in the two inputs must match.
map_clusters_to_genes <- function(intron_meta, exons_table) {  
  gene_df <- foreach(chr = sort(unique(intron_meta$chr)), .combine = rbind) %dopar% {  
    intron_chr <- intron_meta[intron_meta$chr == chr, ]  
    exons_chr <- exons_table[exons_table$chr == chr, ]  
    exons_chr$temp <- exons_chr$start  
    intron_chr$temp <- intron_chr$end  
    three_prime_matches <- inner_join(intron_chr, exons_chr, by="temp")  
    exons_chr$temp <- exons_chr$end 
    intron_chr$temp <- intron_chr$start 
    five_prime_matches <- inner_join(intron_chr, exons_chr, by="temp")  
    all_matches <- rbind(three_prime_matches, five_prime_matches)[, c("clu", "gene_name")]  
    all_matches <- all_matches[!duplicated(all_matches), ]  
    if (nrow(all_matches) == 0) return(NULL)  
    all_matches$clu <- paste(chr, all_matches$clu, sep=':')  
    all_matches}  
  clu_df <- gene_df %>% group_by(clu) %>% summarize(genes = paste(gene_name, collapse = ","))  
  class(clu_df) <- "data.frame"  
  clu_df}
p <- arg_parser("LeafCutter: map clusters to genes")
p <- add_argument(p, "intron_counts_file", help="Intron counts file from LeafCutter, typically <prefix>_perind.counts.gz",default = "perind.counts.gz")
p <- add_argument(p, "exon_file", help="File listing all unique exons in annotation. Must have columns: chr, start, end, strand, gene_id[, gene_name].",default = "exons.txt")
p <- add_argument(p, "output_name", help="tissue_perind.counts.gz")
p <- add_argument(p, "--output_dir", short="-o", help="Output directory", default = "/rawdata/")
argv <- list(  
  intron_counts_file = "perind.counts.gz",  
  exon_file = "exons.txt",  
  output_name = "tissue_perind.counts.gz",  
  output_dir = "/rawdata/"  
)
cat("LeafCutter: mapping clusters to genes\n")
library(data.table)  
library(doParallel)
"%&%" = function(a,b){paste0(a,b)}
cl <- makeCluster(2)  
clusterEvalQ(cl, library(dplyr))  
registerDoParallel(cl)  
intron_counts <- read.table(argv$intron_counts_file, header=TRUE, check.names=FALSE, row.names=1)
intron_meta <- get_intron_meta(rownames(intron_counts))
exon_table <- read.table(argv$exon_file, header=TRUE, stringsAsFactors=FALSE)
stopifnot(is.element('gene_id', colnames(exon_table)))
exon_table[, 'gene_name'] <- exon_table[, 'gene_id']
m <- map_clusters_to_genes(intron_meta, exon_table)
write.table(m, file.path(argv$output_dir, argv$output_name), sep = "\t", quote=FALSE, row.names=FALSE)
