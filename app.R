# Load packages
lapply(c("shiny"), require, character.only = T)

##################################################################################
### Setup page with a navigation bar
ui <- navbarPage( title = "NDOW Upland Harvest",
                  
                  ###INSTRUCTIONS
                  tabPanel(
                    title = "Instructions",
                    
                    p("This is where instructions will go!"),
                    
                    downloadButton("downloadallData", "Download Model Estimates")
                  ),
                  
                  ### Harvest Forecast
                  tabPanel(
                    title = "Harvest Forecast",
                    
                    selectInput("alltoplot", "Which dataset would you like to see forcast?",
                                choices = c("Hunter Effort", "Total Harvest")),
                    
                    conditionalPanel(
                      condition = "input$alltoplot == 'Hunter Effort'",
                      
                      plotOutput("allHfc.plot"),
                      
                      plotOutput("allHcor.plot")
                      
                    ),
                    
                    conditionalPanel(
                      condition = "input$alltoplot == 'Total Harvest'",
                      
                      plotOutput("allNfc.plot"),
                      
                      plotOutput("allNcor.plot")
                      
                    ),
                    
                  ),
                  
                  ### Species - Specific Information
                  tabPanel(
                    title = "Species-Specific",
                    
                    
                    selectInput("speciestoplot", "Which species would you like to see data for?",
                                choices = c("Hunter Effort", "Total Harvest", "Chukar Site Abundance", "Sage Grouse Wing-Bee")),
                    
                    plotOutput("oneNfc.plot"),
                    
                    plotOutput("oneHfc.plot"),
                    
                    plotOutput("oneNcor.plot"),
                    
                    conditionalPanel(
                      condition = "input$speciestoplot == 'Sage Grouse'",
                      plotOutput("sg.plot")
                      
                      ),
                    
                    conditionalPanel(
                      input$speciestoplot == "input$speciestoplot == 'Chukar'",
                      plotOutput("chuk.plot")
                      
                    ),
                    
                  )
)


##################################################################################
### Server Instructions
server <-  function(input,output){
  
  output$allNfc.plot <- renderPlot(
    ggplot(data = read.csv("./www/out_N_all.csv"), aes(x = Year, y = Estimate, group = Species)) +
      geom_line(aes(color = Species)) +
      facet_wrap(vars(Region))
  )
  
  output$allHfc.plot <- renderPlot(
    ggplot(data = read.csv("./www/out_H_all.csv"), aes(x = Year, y = Estimate, group = Species)) +
      geom_line(aes(color = Species)) +
      facet_wrap(vars(Region))
  )
  
  ### "Download/Input Data"
  # Reactive value for selected dataset ----
  species.plot <- reactive({
    
    switch(input$speciestoplot,
           "Blue Grouse" = "blgr",
           "California Quail" = "caqu",
           "Chukar" = "chuk",
           
           )
  })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, ".csv", sep = "")
    },
    content = function(file) {
      write.csv(read.csv("FullDataDownload.csv"), file, row.names = F)
    }
  )
}


##################################################################################
### Construct Shiny App
shinyApp(ui = ui, server = server)

# ################################################################################
#https://shiny.rstudio.com/images/shiny-cheatsheet.pdf
#https://vimeo.com/131218530?embedded=true&source=video_title&owner=22717988
#https://ecoforecast.org/workshops/r-shiny-seminar-series/
#https://mastering-shiny.org/
# ### Example code from the RStudio Tutorial Video
# 
# 
# ## CSS Options
# ## Can stylize page using premade or adjusted css files which can be found online or made
# 
# 
# ### User Interface
# #can alternatively use fixedPage(), dashboardPage(), or navbarPage()/navbarMenu()
# ui <- fluidPage(
#   sidebarLayout(
#     sidebarPanel(),
#     mainPanel()
#   ),
#   
#   fluidRow(column(3), #number is width, out of 12 units
#            column(5)),
#   fluidRow(4, offset = 8), #same for offset
#   
#   wellPanel(), #visually combines objects together in output
#   #Many types of panels
#   tabPanel("Tab 1", ...), #Acts as its own ui, needs a title, then list any elements using...
#   tabsetPanel(
#     tabPanel("tab 1", "contents"),
#     tabPanel("tab 2", "contents"),
#     tabPanel("tab 3", "contents")
#   ), #combines panels together
#   
#   navlistPanel(), #alternative config for tabs using a sidebar
#   
#   fluidh1("Upland Harvest Trends"),
#   p("Here are the", tags$strong("Instructions")),
#   
#   hr(), #line across to seperate sections
#   
#   img(height = 1000,
#       width = 1000,
#       src = "https://www.URL.com/r.jpg"), #can also link to file just by name.ext in app directory folder when placed in "www" folder
#   
#   #Input() functions
#   sliderInput(inputId = "num", 
#               label = "Choose a number",
#               value = 25,
#               min = 1,
#               max = 100),
#   
#   #Output() fuynctions
#   plotOutput(outoutId = "hist"),
#   
#   actionButton(inputID = "clicks",
#                label = "Click Me")
# )
# 
# ### Server Instructions
# server <-  function(input,output){
#   
#   data <- reactive({
#     rnorm(input$num)
#   })
#   
#   output$hist <- renderPlot({
#     hist(rnorm(data()))
#   })
#   
#   output$stats <- renderPrint({
#     summary(data())
#   })
#   
#   observeEvent(input$go, {
#     print(input$clicks)
#   })
#   
#   reactiveEvent()
#   
#   reactiveValues(data = rnorm(100))
# }
# 
# ### Construct Shiny App
# shinyApp(ui = ui, server = server)