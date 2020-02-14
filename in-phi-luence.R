# Automatic generation of network of influences of philosophers
# based on Wikipedia pages
# PAUVERT Charlie
# 2020-02-11


library(rvest)

# PetScan engine to search for articles in wikipedia
# using the query https://petscan.wmflabs.org/?psid=15433482
# that searched articles using the template Infobox philosopher and 
# belonging to the category Philosophy of science
phi_list<-read.csv("phi-science-list.csv", stringsAsFactors = F)$title

get_influences<-function(foo, filename){
    # Extract the name and construct the link 
  #   with the suffix https://en.wikipedia.org/wiki/
  url<-paste0("https://en.wikipedia.org/wiki/",foo)
  # Extract the name and keep as label
  # label<-gsub("_"," ",foo)# later
  # Go to the page (if it exists)
  html<-read_html(url)
  # Find the infobox and influences (if it exists)
  boxes<-html %>%
    html_node("table.infobox") %>% # Get the wikipedia infobox
    html_nodes("tr.note") # Find all boxes
  id<-boxes %>%
    html_node("div.NavHead") %>% xml_text() == "Influences"
  # True False NA vector indicating the influences
  if(!any(na.exclude(id))){
    message("No influences documented, skipping: ",foo)
    return(NA)
  }
  message(foo)
  # Extract the influences
  influences<-boxes %>% html_nodes("ul.NavContent") %>% # Get the boxes that contains a list  
    `[[`(which(na.exclude(id))) %>% # Extract the Influences
    html_nodes("a") %>% # Find the links
    html_attr("href") %>% # Extract the links instead of the text!
    .[!grepl("#|CITESHORT|redlink",.)] %>%  # Remove potential references links, missing refs or inexistant pages
    gsub(".*/wiki/","",.) %>% # Remove the url prefix totally
    sapply(.,URLdecode,USE.NAMES = F) # Decode special characters %C3%A9 -> é
  # Append the Influences into a file
  write.table(
    cbind(influences,foo),file = filename,row.names = F,col.names = F,append = T
    )
  system("sleep .7")
  return(1)
}

out<-lapply(phi_list,get_influences, filename="influences.txt")
# Network generated on Tue Feb 11 16:18:27 2020

# Generate a network from the edges list
library(igraph)

# Read the fetch edge table
inf_list<-as.matrix(read.table("influences.txt", stringsAsFactors = F))

# Duplicates entry are still present, e.g. 
# Georg Wilhelm Friedrich Hegel
# Georg_Hegel
# G._W._F._Hegel
# and that needs to be homogeneized.
# Perform a request for each article, save the article title 
#   and create a named list to convert back synonyms to proper nomenclature
synonyms<-sapply( inf_list %>% c() %>% unique(),# From the unique names
                  function(foo){
                    tryCatch(
                      paste0("https://en.wikipedia.org/wiki/",foo) %>%
                      read_html() %>%
                      html_node("h1") %>% html_text() %>% gsub(" ","_", .),
                      error = function(e){foo})
                  })
# Convert names to erase duplicates
inf_list<-apply(inf_list, 1:2,function(foo) synonyms[foo])
# Graph creation
(g<-graph_from_edgelist(inf_list, directed = T))


# Add the "inlist" attribute
V(g)$inlist<-V(g)$name %in% phi_list

# Export the graph to GraphML for future visualisation
write_graph(g, "influences.graphml", format = "graphml")

# Get the 5 best influencers
degree(g,mode = "out") %>% sort(decreasing = T) %>% .[1:5]
# Immanuel_Kant           Aristotle Ludwig_Wittgenstein               Plato 
# 24                  20                  20                  16 
# David_Hume 
# 12 

# Visualisation with ggraph
library(ggraph)
ggraph(g, layout = "igraph",algorithm="fr") +
  geom_edge_link(arrow = arrow(length = unit(2, 'mm')), 
                 end_cap = circle(1, 'mm')) + 
  geom_node_point() +
  geom_node_point(data = function(x){x[x$name=="Ludwig_Wittgenstein",]},color="orange") +
  coord_fixed() + theme_graph()

# Interactive visualisation
library(visNetwork)
library(dplyr)

vis_data<-toVisNetworkData(g)
vis_data$nodes <- vis_data$nodes %>% mutate(
  title = paste0("<a href='https://en.wikipedia.org/wiki/",label,
                 "'><b>",gsub("_"," ",label),"</b></a>"))

visNetwork(vis_data$nodes, vis_data$edges) %>%
  visEdges(arrows = "to") %>% # Add directionality
  visIgraphLayout() %>%
  visOptions(highlightNearest = T)

# Visualise only the philosophers listed
g %>% delete.vertices(V(g)$name %>% .[!. %in% phi_list]) %>% visIgraph()