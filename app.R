library(shiny)
library(igraph)
library(dplyr)
library(visNetwork)
library(DT)

# Define server side required to visualise the network
server <- function(input, output,session) {
  net<-reactive({
    # Data sources dependent input file
    g_src<-switch (input$data_src,
            "wk" = "influences.graphml",
            "inpho" = "inpho/inpho.graphml")
    # Read the network object
    read_graph(g_src, format = "graphml")
  })
  data<-reactive({
    link_src<-switch (input$data_src,
                      "wk" = "https://en.wikipedia.org/wiki/",
                      "inpho" = "https://www.inphoproject.org/entity?redirect=true&q=")
    vis_data<-toVisNetworkData(net())
    vis_data$nodes <- vis_data$nodes %>% mutate(
      title = paste0("<a target='_blank' href='",link_src, label,
                     "'><b>",gsub("_"," ",label),"</b></a>"))
    vis_data
  })
   output$network <- renderVisNetwork({
     visNetwork(data()$nodes, data()$edges) %>%
       visNodes(shape = "diamond", color = list(background = "darkgray", border = "white", highlight = "black"), borderWidth = 2) %>%
       visEdges(arrows = list(to = list(enabled = TRUE, scaleFactor = 0.5)), color = "darkgray", width = 2) %>%
       visIgraphLayout() %>%
       visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T))
   })
   
   metrics <- reactive({
     data.frame(
       Philosophers=V(net()) %>% as_ids() %>% gsub("_", " ", .),
       In=degree(net(), mode = "in"),
       Out=degree(net(), mode = "out")
       )
   })
   output$table <- DT::renderDataTable(metrics(),
                                       rownames = FALSE,
                                       caption = paste(
                                         "Network metrics at the node (philosopher) level.",
                                         "The In degree indicates the number of incoming links to the network node,",
                                         "or the number of influences for the philosopher.",
                                         "Conversely, the Out degree indicates the number of influencees stemming from the philosopher.")
   )
   observe({
     updateSelectInput(session, "selnode",
                       choices = data()$nodes$id %>%
                         sort() %>%
                         gsub("_", " ",.))
   })
   
   observeEvent(input$show,{
     visNetworkProxy("network") %>%
       visSelectNodes(id = gsub(" ", "_", input$selnode))
   })
   
   observeEvent(input$unselect,{
     visNetworkProxy("network") %>%
       visUnselectAll()
   })
   inph_list<-reactive({
     adj_lists<-list(
       "influencers" = E(net())[to(gsub(" ", "_", input$selnode))] %>% # Get the influencers
         tail_of(net(), .), # hence the tail of the arrow
       "influencees" = E(net())[from(gsub(" ", "_", input$selnode))] %>% # Get the influencees
         head_of(net(), .) # hence the head of the arrow
     )
     lapply(adj_lists, function(x){
       x %>% as_ids() %>%
         sort() %>% gsub("_", " ",.)}
       )
   })
   # Lists solution from https://stackoverflow.com/a/50414101
   output$influencers<-renderUI(
     lapply(inph_list()[["influencers"]], function(x) tags$li(x))
   )
   output$influencees<-renderUI(
     lapply(inph_list()[["influencees"]], function(x) tags$li(x))
   )
   output$philosopher<-renderText(input$selnode)
   output$title<-renderText(
     paste(nrow(data()$nodes), "philosophers and",
           nrow(data()$edges), "influences links")
   )
}

# Define UI for the visualisation
ui <- navbarPage(
  # Application title
  "In-phi-luence",
  tabPanel("Visualise influences",
  # Sidebar with node search enabled
  sidebarLayout(
    sidebarPanel(
      h4("A network view of philosophers"),
      p("Philosophers are here represented as",em("nodes"),"of the network.",
        "The",em("links"), "between nodes indicates the", strong("influence"),"one philosopher had",
        "on another one. The direction of the link follows knowledge flow."),
      helpText(strong("Example:"),"Aristotle was influenced by Plato"),
      hr(),
      radioButtons("data_src", "Select a source for the influences data",
                   choices = c("The Free Encyclopedia Wikipedia (en)" = "wk",
                               "The Internet Philosophy Ontology Project" = "inpho")),
      strong(""),
      selectInput(inputId = "selnode",
                  label = "Explore freely the interactive network or find a philosopher in the list:",
                  choices = "",
                  multiple = FALSE),
    actionButton("show", "Show"),
    actionButton("unselect", "Unselect"),
    br(),
    helpText("Note: Search for a philosopher and click on `Show`",
             "to highlight its position in the network.")),
    # Show the network
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Network",
                           h4(textOutput("title")),
                           visNetworkOutput("network",height = "500px", width = "auto")
                           ),
                  tabPanel("Table",
                           br(),
                           DT::dataTableOutput("table")),
                  tabPanel("Influences list",
                           h4(tags$em(textOutput("philosopher"))),
                           strong("was influenced by:"), tags$ul(uiOutput("influencers")),
                           strong("has influenced:"), tags$ul(uiOutput("influencees"))
                  ))
    ))),
  tabPanel("Materials and methods",
           fluidPage(
             fluidRow(
               h3("Rationale"),
               p("This app allows to visualise knowledge flow as a whole.",
                 "It is also a way to promote the content of knowledge sources that are accessible",
                 "once hovering on the nodes in the network.",
                 "Scholars and students in philosophy could explore the network",
                 "and even spot missing links that should be indicated in Wikipedia articles.",
                 "The original idea stems from my partner wondering whether it was possible to visualise connections",
                 "between philosophers and their schools of thoughts."),
               h3("Influences data sources")),
             fluidRow(
               column(6, h4("The Free Encyclopedia Wikipedia (en)"),
                      p(a(href="https://en.wikipedia.org/", "Wikipedia"), "articles were first fetched if they had a Philosopher Infobox",
                        "and if they belong to the category of Philosophy of science.",
                        "Mentions of influences in the infobox were mined, collected and gathered into a network.",
                        "Additional philosophers, novelist and school of thoughts were also added after following",
                        "the influences linked in the listed Wikipedia pages."),
                      p("The currently displayed network was generated on 2020-02-14."),
                      p("Wikipedia articles were gathered with the powerful search tool",
                        a(href="https://petscan.wmflabs.org/","PetScan"),".",
                        "R code was then used to fetch and clean data and is available at",
                        a(href="https://github.com/cpauvert/in-phi-luence","https://github.com/cpauvert/in-phi-luence"))),
               column(6, h4("The Internet Philosophy Ontology Project (InPho)"),
                      p("This", a(href="https://www.inphoproject.org/", "amazing scholarly resource"), "compile ontologies",
                        "on philosophers and which are then made accessible through API or files.",
                        "Monthly archives of the InPho ontologies were fetched. Automatic mining of the ontologies extracted",
                        "relevant properties (such as", code("has_influenced"), "or",code("was_influenced") ,
                        ") and all results were gather into a network."),
                      p("The currently displayed network was generated on 2020-12-16."),
                      p("R and Python code were used to fetch and clean data and are availaible at",
                        a(href="https://github.com/cpauvert/in-phi-luence/inpho","https://github.com/cpauvert/in-phi-luence/inpho"))
               )),
             fluidRow(h3("What's next?"),
                      p("Future avenues include an update of the sources used, a comparison between the two generated networks.",
                        "The latter could help precise missing influences in Wikipedia articles for instance.",
                        "Planned milestones are listed",a(href="https://github.com/cpauvert/in-phi-luence/milestones","here"),".",
                        "Suggestions and remarks are very welcomed (en/fr)",
                        a(href="https://github.com/cpauvert/in-phi-luence/issues","here."))
             )))
)


# Run the application 
shinyApp(ui = ui, server = server)

