## code to prepare `data_LUSC_demo` dataset goes here
## Original preprocessing code available upon request
data_LUSC_demo <- readr::read_csv("./demo/data/LUSC.csv")
usethis::use_data(data_LUSC_demo, overwrite = TRUE)
