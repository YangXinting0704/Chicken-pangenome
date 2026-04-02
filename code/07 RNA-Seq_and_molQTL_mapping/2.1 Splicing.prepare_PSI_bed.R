#!/usr/bin/env Rscript
library(dplyr)  
library(readr)
library(tibble)
library(purrr)
library(data.table)
setwd("/transcriptom/rawdata/2.splicing/")
file_list <- list.files(pattern = "tissue_perind.counts.gz.qqnorm_")
merged_data <- purrr::map_dfr(  
  file_list,  
  ~ fread(.x, header = TRUE) %>%   
     mutate(`#Chr` = as.character(`#Chr`), ID = as.character(ID))  
)
fwrite(merged_data, "tissue.splicing_qqnorm.PSI.txt", sep = "\t", row.names = FALSE, col.names = TRUE)  

#Finally, replace the splicing position information with that of the corresponding gene.