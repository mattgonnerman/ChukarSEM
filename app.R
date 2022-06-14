#https://daniel-gibson.shinyapps.io/RLP_MODEL/?_ga=2.73890838.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/JSM_MODEL/?_ga=2.76979217.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/CANV_REDH_REL/?_ga=2.111113217.54103719.1653590001-1179135978.1653590001


# Load packages
lapply(c("shiny","ggplot2", "dplyr"), require, character.only = T)

#Fix for unaligned checkboxes
tweaks <- tags$head(tags$style(HTML(
  ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
          }
         .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
          }
        "
)))
  

##################################################################################
### Setup page with a navigation bar
ui <- navbarPage( title = "NDOW Upland Harvest",
                  tweaks,
                  
                  ###INSTRUCTIONS
                  tabPanel(
                    title = "Instructions",
                    
                    p("This is where instructions will go!\n
                      \n
                      [Forestcasting Summary] \n
                      \n
                      [Species Specific Summary]"),
                    
                    downloadButton("downloadallData", "Download Model Estimates")
                  ),
                  
                  ### Harvest Forecast
                  tabPanel(
                    title = "Forecast",
                    
                    selectInput("forecastset", "Which dataset would you like to see forcast?",
                                choices = c("Birds per Hunter" = "BPH", 
                                            "Total Harvest" = "N", 
                                            "Hunter Participation" = "H")),
                    
                    selectInput("spphighlight", "Species:",
                                choices = c("All",
                                            "BLGR",
                                            "CAQU",
                                            "CHUK",
                                            "HUPA",
                                            "PHEA",
                                            "RABB",
                                            "SAGR"))
                  ),
                  
                  ### Past Years
                  tabPanel(
                    title = "Past Years",
                    
                    sidebarLayout(
                      sidebarPanel(#Specify Plot Settings
                        selectInput("pastset", "Which dataset would you like to see past estimates for?",
                                    choices = c("Birds per Hunter" = "BPH", 
                                                "Total Harvest" = "N", 
                                                "Hunter Participation" = "H")),
                        
                        checkboxGroupInput("pastspp", "Species:",
                                           choices = c("BLGR","CAQU","CHUK","HUPA","PHEA","RABB","SAGR"),
                                           selected = c("CHUK", "SAGR"), inline = TRUE, width = "350px"),
                        
                        sliderInput("pastyears", "Years", min=1976, max=2017, 
                                    value=c(1976, 2017), sep=""),),
                      mainPanel(#Outplut Plot
                        plotOutput("pastplot"),
                        
                        
                        #Survey specific graphs for Grouse and Chukar
                        conditionalPanel(
                          condition = "input$speciestoplot == 'SAGR'",
                          plotOutput("sg.plot")
                          
                        ),
                        
                        conditionalPanel(
                          condition = "input$speciestoplot == 'CHUK'",
                          plotOutput("chuk.plot")
                          
                        )
                      )
                    )
                  )
)


##################################################################################
### Server Instructions
server <-  function(input,output){
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","PHEA","RABB","SAGR")
  colors.constant <- c("#0033a0","#015b6e","#016f55","#01833b","#80931e","#c09b0f","#ffa300")
  
  
  past.data <- reactive(read.csv(paste("./www/out_", input$pastset, "_all.csv", sep = "")) %>%
                          mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
                                 Year = Year + 1975,
                                 Species = species.constant[as.numeric(Species)]) %>%
                          filter(Year >= input$pastyears[1], Year <= input$pastyears[2],
                                 Species %in% input$pastspp))
  
  colors <- reactive(setNames(colors.constant[seq(1,7, length.out = floor(length(unique(past.data()$Species))))],
                              species.constant[which(species.constant %in% unique(past.data()$Species))]))
  

  
  output$pastplot <- renderPlot(
    ggplot(data = past.data(), aes(x = Year, y = Estimate, group = Species)) +
      geom_line(aes(color = Species), size = 1.3) +
      facet_wrap(vars(Region)) +
      theme_classic(base_size = 18) +
      scale_color_manual(values = colors()) +
      theme(axis.title = element_blank(), legend.position = "bottom")
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