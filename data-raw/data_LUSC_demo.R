## code to prepare `data_LUSC_demo` dataset goes here
data_LUSC_demo <- readr::read_csv("./demo/data/LUSC.csv")
usethis::use_data(data_LUSC_demo, overwrite = TRUE)
