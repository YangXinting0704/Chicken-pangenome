#!/usr/bin/env Rscript
options(warn = -1)
suppressMessages(library(stringr))
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(edgeR))
suppressMessages(library(argparser)) 
suppressMessages(library(purrr))
p <- arg_parser("Merge featureCounts(linux version) and calculate FPKM/TPM")
p <- add_argument(p, "--input_path", help="input: a directory containing the counts matrix named with '<sample>.tsv'", type="character",default = "/transcriptom/rawdata/tissue/1.gene/data")
p <- add_argument(p, "--pattern", help="limit pattern of input files using regular expression in R language",type="character",default = ".*\\.tsv$")
p <- add_argument(p, "--output_path", help="output: an existent directory", type="character",default = "/home/bailu/transcriptom/rawdata/tissue/1.gene")
p <- add_argument(p, "--prefix", help="give the file of output matrix a prefix like '<output_prefix>_genes.*'", type="character",default = "tissue",short = "-f")
argv <- parse_args(p)
path <- argv$input_path
pattern <- argv$pattern
output_path <- argv$output_path
output_prefix <- argv$prefix
file_name <- dir(path = path,pattern = pattern)
if(is_empty(file_name)){stop("\nCould not find files with pattern of '",pattern,"'!")}
file <- paste0(path,"/",file_name)
# merge all count matrix
df <- read.table(file[1], header = T,comment.char = "#") %>% select(c(1,6,7))
colnames(df)[3] <- basename(file[1]) %>% str_remove("\\.\\w+$")
cat("1 count matrix has merged!\n")
for (i in 2:length(file)) {
  df_tmp <- read.table(file[i], header = T,comment.char = "#") %>% select(c(1,7))
  colnames(df_tmp)[2] <- basename(file[i]) %>% str_remove("\\.\\w+$")
  df <- df %>% full_join(df_tmp,by="Geneid")
  cat(i," count matrixs have merged!\n")
}
count_mat <- df %>% select(-2) %>% tibble::column_to_rownames(var = "Geneid")
gene_length <- df %>% select(1,2) %>% tibble::column_to_rownames(var = "Geneid")
filter_count_mat <- count_mat[rowSums(count_mat)>0,]
fpkm <- rpkm(count_mat, gene.length = gene_length$Length) %>% as.data.frame()
fpkm2tpm <- function(fpkm) {
  if ((is.matrix(fpkm) | is.data.frame(fpkm)) & all(fpkm >= 0)) {
    tpm <- t(t(fpkm) / colSums(fpkm)) * 10^6
  }
  return(tpm)
}
tpm <- fpkm2tpm(fpkm = fpkm) %>% as.data.frame()

out_count_mat <- count_mat %>% rownames_to_column(var = "gene_id")
count_file <- paste0(output_path, "/", output_prefix, "genes_expression.counts")
write.table(out_count_mat, file = count_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
out_tpm <- tpm %>% rownames_to_column(var = "gene_id")
tpm_file <- paste0(output_path, "/", output_prefix, "genes_expression.tpm")
write.table(out_tpm, file = tpm_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)


#quality control and normalization
library(preprocessCore)
library(limma)
TPM_tmp0<-read.table("gene_expression.tpm",header=T)
row.names(TPM_tmp0)<-TPM_tmp0$ID
TPM_tmp0<-TPM_tmp0[,-1]
count_0.1<-rowSums(TPM_tmp0>0.1)  
nsamples<- n
expr_matrix00<-TPM_tmp0[count_0.1>=(0.2*nsamples),]   
expr_matrix00_qn <- normalizeQuantiles(as.matrix(expr_matrix00))   
INT<-function(x){ qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
expr_matrix00_qn_ind<-t(apply(expr_matrix00_qn,MARGIN=1,FUN=INT))  
colnames(expr_matrix00_qn_ind)<-colnames(expr_matrix00)
rownames(expr_matrix00_qn_ind)<-rownames(expr_matrix00)
write.table(expr_matrix00_qn_ind ,"gene.expression.qn_ind.txt",sep="\t",row.names=T,quote =FALSE)
