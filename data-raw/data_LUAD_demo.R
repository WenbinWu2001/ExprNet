## code to prepare `data_LUAD_demo` dataset goes here
## Original preprocessing code available upon request
data_LUAD_demo <- readr::read_csv("./demo/data/LUAD.csv")
usethis::use_data(data_LUAD_demo, overwrite = TRUE)
