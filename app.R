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
     visNetwork(data()$nodes, data()$edges,
                main = list(
                  text=paste(nrow(data()$nodes),"philosophers and",nrow(data()$edges),"links"),
                  style="font-family:Helvetica,Arial,sans-serif;font-size=14px;text-align:center"),
                footer = list(
                  text="Network generated on 2020-02-14.",
                  style="font-family:Helvetica,Arial,sans-serif;font-size=8px;text-align:right")
                ) %>%
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
   inph_list<-isolate({
     g<-read_graph("influences.graphml", format = "graphml")
     foo<-list("from"=as_adj_list(g, mode = "in"),
          "to"=as_adj_list(g, mode = "out"))
     lapply(foo, function(dir) lapply(dir, function(x) x %>% as_ids() %>% sort()))
   })
   output$influencers<-renderTable(
     data.frame("Influencers"=inph_list[["from"]][[input$selnode]]),
     striped = TRUE, hover = TRUE, align = "c")
   output$influencees<-renderTable(
     data.frame("Influencees"=inph_list[["to"]][[input$selnode]]),
     striped = TRUE, hover = TRUE, align = "c")
   output$philosopher<-renderText(input$selnode)
}

# Define UI for the visualisation
ui <- fluidPage(
  # Application title
  titlePanel("In-phi-luence: a network view of philosophers"),
  # Sidebar with node search enabled
  sidebarLayout(
    sidebarPanel(
      h4("What is it?"),
      p("Philosophers of science are here represented as",em("nodes"),"of the network.",
        "The",em("links"), "between nodes indicates the", strong("influences"),"one philosopher had",
        "on another one. The direction of the link follows knowledge flow."),
      helpText(strong("Example:"),"Aristotle was influenced by Plato"),
      br(),hr(),
      selectInput(inputId = "selnode",
                  label = "Find Philosophers in the network",
                  choices = "",
                  multiple = FALSE),
    actionButton("focus", "Focus"),
    actionButton("reset", "Reset view"),
    br(),
    helpText("Note: Select a philosopher and click on focus",
             "to highlight its influences. ",
             "Use 'Reset' to get back to the full network view.")),
    # Show the network
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Network",
                           visNetworkOutput("network",height = "500px", width = "auto")
                           ),
                  tabPanel("Details", fluidRow(
                    h3(textOutput("philosopher")),
                    column(3, tableOutput("influencers")),
                    column(3, tableOutput("influencees"))
                  ))),
    )),
  br(),
    fluidRow(
      column(4, h4("How?"),
             p("Wikipedia articles were first fetched if they had a Philosopher Infobox",
               "and if they belong to the category of Philosophy of science.",
               "Mentions of influences in the infobox were mined, collected and gathered into a network."),
             p("The code to fetch and clean data is available at",
               a(href="https://github.com/cpauvert/in-phi-luence","https://github.com/cpauvert/in-phi-luence"))),
    column(4,
           h4("Why?"),
           p("This app allows to visualise knowledge flow as a whole.",
             "It is also a way to promote the content of Wikipedia articles which are accessible",
             "once hovering on the nodes."),
           p("Scholars and students in philosophy could explore the network",
             "and even spot missing links that should be indicated in Wikipedia articles.")),
    column(4 ,h4("What's next?"),p("New features could include:"),
           tags$ul(
             tags$li("Expand the search to all philosophers"),
             tags$li("Intersect the network with one built from articles in different languages"),
             tags$li("Restrict the visualisation to philosophers only, not influences like novelists")
           ),p("Suggestions and remarks are welcomed (en/fr)",
           a(href="https://github.com/cpauvert/in-phi-luence/issues","here.")))
  )
)


# Run the application 
shinyApp(ui = ui, server = server)

