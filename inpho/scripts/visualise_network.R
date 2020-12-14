# Visualise interactively the obtained network 
# PAUVERT Charlie
# 2020-12-14
library(igraph)
library(visNetwork)

g <- graph_from_data_frame(read.delim(snakemake@input[[1]], header=F))
data <- toVisNetworkData(g, idToLabel=F)
data$nodes$title <- paste0("<p>", data$nodes$id, "</p>")

visNetwork(nodes = data$nodes, edges = data$edges, height = "500px") %>%
  visEdges(arrows = "to") %>% 
  visOptions(highlightNearest = list(enabled = T, hover = T)) %>% 
  visSave(file = snakemake@output[[1]])

