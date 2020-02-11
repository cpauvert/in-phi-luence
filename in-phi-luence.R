# Automatic generation of network of influences of philosophers
# based on Wikipedia pages
# PAUVERT Charlie
# 2020-02-11


library(rvest)
library(dplyr)

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
  if(!any(id)){
    message("No influences documented, skipping: ",foo)
    return(NA)
  }
  message(foo)
  # Extract the influences
  influences<-boxes %>% html_nodes("ul.NavContent") %>% # Get the boxes that contains a list  
    `[[`(which(na.exclude(id))) %>% # Extract the Influences
    html_nodes("a") %>% # Find the links
    html_text()
  # Remove potential references links
  influences<-influences[!grepl("\\[",influences)]
  # Homogeneize with foo by adding underscore
  influences<-gsub(" ","_",influences)
  # Append the Influences into a file
  write.table(
    cbind(foo,influences),file = filename,row.names = F,col.names = F,append = T
    )
  system("sleep 1")
  return(1)
}

out<-lapply(phi_list[1:10],get_influences, filename="influences.txt")

