library(readr)
library(dplyr)

bedFormat <- function(file, sample, outpath) {
  df <- read_tsv(file, col_names = T)
  df <- df %>%
    select(1:3) %>%
    arrange(chrom, start) %>%
    mutate(X4 = seq(1, nrow(.))) %>%
    mutate(X4 = paste(sample, X4, sep = "_"))
  write_tsv(df, stringr::str_glue("{outpath}/{sample}.ecc.bed"), col_names = F)
}


args <- commandArgs(trailingOnly = TRUE)
##
## args[1] is the input file
## args[2] is the sample name
## args[3] is the output path

bedFormat(args[1], args[2], args[3])
