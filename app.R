#https://daniel-gibson.shinyapps.io/RLP_MODEL/?_ga=2.73890838.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/JSM_MODEL/?_ga=2.76979217.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/CANV_REDH_REL/?_ga=2.111113217.54103719.1653590001-1179135978.1653590001


# Load packages
lapply(c("shiny", "shinyWidgets", "ggplot2", "dplyr"), require, character.only = T)

# Load Data from Model
load("./www/NDOW_SEM_app_objects.R")

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
  ), setBackgroundColor(color = "blue") ) )

  

##################################################################################
### Setup page with a navigation bar
ui <- navbarPage( tweaks,
                  
                  title = "NDOW Upland Harvest", selected = "Instructions",
                  ###INSTRUCTIONS
                  tabPanel(
                    title = "Instructions",
                    
                    # p("This is where instructions will go!", <br>,"[Forestcasting Summary]", <br>,
                    #   "[Species Specific Summary]"),
                    
                  ),
                  
                  ### Harvest Forecast
                  tabPanel(
                    title = "Forecast",
                    
                    sidebarLayout(
                      sidebarPanel(selectInput("forecastset", "Which dataset would you like to see forcast?",
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
                                                           "SAGR")),
                                   
                                   
                                   
                                   
                                   sliderInput("econ", "Economic Predictor", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("bbs", "Predator Predictor", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("awssi", "Winter Severity", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("pdsi", "Drought Conditions", min=-2, max=2, 
                                               value = 0, sep=""),
                      ),
                      
                      mainPanel(#Output Plot
                        plotOutput("forecastplot"),
                      )
                      
                    )
                    
                    
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
                                    value=c(1976, 2017), sep=""),
                        
                        #Do we want a switch to show correlation among species? May need to calculate rho for BPH
                        # checkboxInput("rhoswitch", "Plot Species Correlations", F),
                        
                        downloadButton("downloadallData", "Download Model Estimates")
                      ),
                      
                      mainPanel(#Output Plot
                        plotOutput("pastplot"),
                        )
                      
                      # conditionalPanel(
                      #   condition = "input.rhoswitch == T",
                      #   plotOutput("pastplot")
                      # )
                      
                    )
                  )
)



##################################################################################
### Server Instructions
server <-  function(input,output){
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","PHEA","RABB","SAGR")
  colors.constant <- c("#0033a0","#015b6e","#016f55","#01833b","#80931e","#c09b0f","#ffa300")
  
  
  past.data <- reactive(
    read.csv(paste("./www/out_", input$pastset, "_all.csv", sep = "")) %>%
      mutate(Region = ifelse(as.numeric(Region) == 1 ,"Eastern", "Western"),
             Year = Year + 1975,
             Species = species.constant[as.numeric(Species)]) %>%
      filter(Year >= input$pastyears[1], Year <= input$pastyears[2],
             Species %in% input$pastspp))
  
  colors <- reactive(setNames(colors.constant[seq(1,7, length.out = floor(length(unique(past.data()$Species))))],
                              species.constant[which(species.constant %in% unique(past.data()$Species))]))
  
  
  forecast.plot <- reactive({
    econ.var <- input$econ
    winter.var <-  input$awssi
    predator.var <- input$bbs
    drought.var <- input$pdsi
    
    forecast.calculate <- data.frame(Year = 1:15, Est = econ.var + winter.var + predator.var + drought.var + 1:15)
    forecast.calculate
  })
  
  output$forecastplot <- renderPlot(
    ggplot(data = forecast.plot(), aes(x = Year, y = Est)) +
      geom_point(size = 3)
  )
  
  
  output$pastplot <- renderPlot(
    ggplot(data = past.data(), aes(x = Year, y = Estimate, group = Species)) +
      geom_line(aes(color = Species), size = 1.3) +
      facet_wrap(vars(Region)) +
      theme_classic(base_size = 18) +
      scale_color_manual(values = colors()) +
      theme(axis.title = element_blank(), legend.position = "bottom")
    )
  
  output$corplot <- renderPlot(
    ggplot(data = past.data(), aes(x = Year, y = Estimate, group = Species)) +
      geom_line(aes(color = Species), size = 1.3) +
      facet_wrap(vars(Region)) +
      theme_classic(base_size = 18) +
      scale_color_manual(values = colors()) +
      theme(axis.title = element_blank(), legend.position = "bottom")
  )
  
  
  
  
  
  ### "Download/Input Data"
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