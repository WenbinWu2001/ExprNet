## code to prepare `network_demo` dataset goes here
## Original preprocessing code available upon request
network_demo <- igraph::read_graph("./demo/network_info/network", format = "edgelist")
usethis::use_data(network_demo, overwrite = TRUE)
