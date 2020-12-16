# Visualise interactively the obtained network 
# PAUVERT Charlie
# 2020-12-14
library(igraph)
library(visNetwork)

# Read the gathered influences
net <- as.matrix(read.delim(snakemake@input[[1]], header=F))
# Trim whitespaces in name and decode percent-encoding characters
net <- apply(trimws(net), 1:2, URLdecode)
# Built network
g <- graph_from_data_frame(net)
data <- toVisNetworkData(g, idToLabel=F)
data$nodes$title <- paste0("<p>", data$nodes$id, "</p>")

visNetwork(nodes = data$nodes, edges = data$edges, height = "500px",
           main = paste0("Philosophers influences built from InPho.\n",
                         vcount(g), " thinkers and ",
                         ecount(g), " influences.")) %>%
  visEdges(arrows = "to") %>% 
  visOptions(highlightNearest = list(enabled = T, hover = T)) %>% 
  visSave(file = snakemake@output[[1]])

