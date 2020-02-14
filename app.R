library(shiny)
library(igraph)
library(dplyr)
library(visNetwork)

# Define server side required to visualise the network
server <- function(input, output,session) {
   
  data<-reactive({
    # Read the network object
    g<-read_graph("influences.graphml", format = "graphml")
    vis_data<-toVisNetworkData(g)
    vis_data$nodes <- vis_data$nodes %>% mutate(
      title = paste0("<a target='_blank' href='https://en.wikipedia.org/wiki/",label,
                     "'><b>",gsub("_"," ",label),"</b></a>"))
    vis_data
  })
   output$network <- renderVisNetwork({
     visNetwork(data()$nodes, data()$edges) %>%
       visEdges(arrows = "to") %>% # Add directionality
       visIgraphLayout() %>%
       visOptions(highlightNearest = T)
   })
   
   observe({
     updateSelectInput(session, "selnode",
                       choices = sort(data()$nodes$id))
   })
   
   observeEvent(input$focus,{
     visNetworkProxy("network") %>%
       visFocus(id = input$selnode) %>%
       visSelectNodes(id = input$selnode)
   })
   
   observeEvent(input$reset,{
     visNetworkProxy("network") %>%
       visFit() %>%
       visUnselectAll()
   })
}

# Define UI for the visualisation
ui <- fluidPage(
  
  # Application title
  titlePanel("Influences of Philosophers of science"),
  
  # Sidebar with node search enabled
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "selnode",
                  label = "Search Philosophers",
                  choices = "",
                  multiple = FALSE),
    actionButton("focus", "Focus"),
    actionButton("reset", "Reset view")),
    # Show the network
    mainPanel(
      visNetworkOutput("network",height = "500px", width = "auto")
    )
  )
)


# Run the application 
shinyApp(ui = ui, server = server)

