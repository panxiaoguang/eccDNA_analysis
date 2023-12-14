library(readr)
library(dplyr)

get_filter_svs <- function(x) {
  df <- read_tsv(stringr::str_glue("merged/{x}.merged.allsv.txt"), col_types = cols(
    "c", "d", "d", "c", "d",
    "d", "c", "d", "c", "c",
    "c", "c", "c"
  ))
  noTRA <- df %>% filter(svclass != "TRA")
  multiCaller <- noTRA %>%
    filter(stringr::str_detect(svmethod, ","))
  TRA <- df %>%
    filter(svclass == "TRA")
  filterSV <- bind_rows(multiCaller, TRA)
  write_tsv(filterSV, stringr::str_glue("merged/{x}.merged.allsv.filter.txt"))
}

samples <- read_tsv("samples", col_names = F, col_types = cols("c")) %>% pull(X1)

lapply(samples, get_filter_svs)
