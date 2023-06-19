## code to prepare `data_LGG_demo` dataset goes here
## Original preprocessing code available upon request
data_LGG_demo <- readr::read_csv("./demo/data/LGG.csv")
usethis::use_data(data_LGG_demo, overwrite = TRUE)
