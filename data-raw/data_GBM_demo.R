## code to prepare `data_GBM_demo` dataset goes here
## Original preprocessing code available upon request
data_GBM_demo <- readr::read_csv("./demo/data/GBM.csv")
usethis::use_data(data_GBM_demo, overwrite = TRUE)
