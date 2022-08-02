#https://daniel-gibson.shinyapps.io/RLP_MODEL/?_ga=2.73890838.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/JSM_MODEL/?_ga=2.76979217.54103719.1653590001-1179135978.1653590001
#https://daniel-gibson.shinyapps.io/CANV_REDH_REL/?_ga=2.111113217.54103719.1653590001-1179135978.1653590001


# Load packages
lapply(c("shiny", "shinyWidgets", "ggplot2", "dplyr", "mvtnorm", "R.utils", "tidyr", "lubridate"), require, character.only = T)

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
                                   
                                   # sliderInput("years.f.plot", "Years to Display", min =1976, max = (1975 + ncol(appobject$N)), 
                                   #             value = c(1976:(1975 + ncol(appobject$N))), sep=""),

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
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","SAGR")
  colors.constant <- c("#0033a0","#015b6e","#01833b","#c09b0f","#ffa300")
  colors <- reactive(setNames(colors.constant[seq(1,length(species.constant), length.out = floor(length(unique(past.data()$Species))))],
                              species.constant[which(species.constant %in% unique(past.data()$Species))]))
  
  d.plot.input <- reactive({
    if(input$forecastset == "BPH"){
      d.input.array <- appobject$N.data/(1+appobject$H.data)
    }else{
      if(input$forecastset == "H"){
        d.input.array <- 1000*appobject$H.data
      }else{
        d.input.array <- 1000*appobject$N.data
      }
    }
    
    as.data.frame.table(d.input.array) %>% 
      dplyr::rename(Species = Var1, Year = Var2, Region = Var3, Est = Freq) %>%
      filter(!is.na(Est)) %>%
      mutate(Year = as.numeric(Year),
             Region = as.numeric(Region),
             Species = as.numeric(Species)) %>%
      mutate(Year = Year + 1975)  %>%
      mutate(Species = recode(Species, 
                              '5' = "SAGR",
                              '1' = "BLGR",
                              '3' = "CHUK",
                              '4' = 'HUPA',
                              '2' = 'CAQU'))
    # %>%
    #   filter(Year %in% input$years.f.plot[1]:input$years.f.plot[2])
  })
  
  forecast.plot.input <- reactive({
    #Prep Inputs
    H.reg <- appobject$H.reg
    N.reg <- appobject$N.reg
    Latent <- appobject$latent.trend
    n.species <- nrow(H.reg[[1]])
    econ.var <- input$econ
    winter.var <-  input$awssi
    predator.var <- input$bbs
    drought.var <- input$pdsi
    
    d.start <- 1976
    f.start <- min(which(is.na(appobject$N.data[1,,1]))) + 1975
    d.end <- f.start - 1
    f.end <- ncol(appobject$N.data) + 1975
    y.b <- f.start-1976 #years of data before forecast
    n.year.f <- length(f.end:d.start) #If you want to show more than next year, use this as another dimension in array
    n.samples <- 1000
    
    #Function for formatting forecast graph inputs
    sem_out_comb <- function(f.arr){
      as.data.frame(wrap(f.arr, map=list(NA, 2))) %>%
        mutate(Species = rep(1:n.species, n.year.f),
               Year = rep(1:n.year.f, each = n.species)) %>%
        pivot_longer(cols = 1:2, names_to = "Region", names_prefix = "V") %>%
        as.data.frame()
    }
    #Empty arrays to fill
    mu.N <- mu.H <- array(NA, dim = c(n.species, 2, n.year.f, n.samples))
    BPH.Est <- BPH.LCL <- BPH.UCL <- N.Est <- N.LCL <- N.UCL <- H.Est <- H.LCL <- H.UCL <- array(NA, dim = c(n.species, 2, n.year.f))
    
    for(r in 1:2){
      for(s in 1:n.species){
        for(t in 1:n.year.f){
          #Hunter Effort
          mu.H[s,r,t,] <- rnorm(n.samples, H.reg[[r]]$A.hunt[s], H.reg[[r]]$A.hunt.sd[s]) +
            econ.var*rnorm(n.samples, H.reg[[r]]$B.econ.hunt[s], H.reg[[r]]$B.econ.hunt.sd[s]) +
            rnorm(n.samples, as.numeric(Latent$mean[[r]][s,t]), as.numeric(Latent$sd[[r]][s,t]))
          
          H.Est[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.5)))
          H.LCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.075)))
          H.UCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(.935)))
          
          # Total Harvest
          mu.N[s,r,t,] <- rnorm(n.samples, N.reg[[r]]$A.harv[s], N.reg[[r]]$A.harv.sd[s]) +
            rnorm(n.samples, N.reg[[r]]$B.hunter[s], N.reg[[r]]$B.hunter.sd[s])*rnorm(n.samples, as.numeric(Latent$mean[[r]][s,t]), as.numeric(Latent$sd[[r]][s,t]))+
            winter.var * rnorm(n.samples, N.reg[[r]]$B.ws.harv[s], N.reg[[r]]$B.ws.harv.sd[s]) +
            drought.var * rnorm(n.samples, N.reg[[r]]$B.pdsi.harv[s], N.reg[[r]]$B.pdsi.harv.sd[s]) +
            predator.var * rnorm(n.samples, N.reg[[r]]$B.bbs.harv[s], N.reg[[r]]$B.bbs.harv.sd[s])
          
          N.Est[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.5)))
          N.LCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.075)))
          N.UCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(.935)))
        }
      }
    }
    
    #Birds Per Hunter
    mu.BPH <- mu.N/mu.H
    for(r in 1:2){
      for(s in 1:n.species){
        for(t in 1:n.year.f){
          BPH.Est[s,r,t] <- 1000*exp(quantile(mu.BPH[s,r,t,], c(.5)))
          BPH.LCL[s,r,t] <- 1000*exp(quantile(mu.BPH[s,r,t,], c(.075)))
          BPH.UCL[s,r,t] <- 1000*exp(quantile(mu.BPH[s,r,t,], c(.935)))
        }
      }
    }
    
    N.df1 <- sem_out_comb(N.Est) %>% dplyr::rename(Est = value)
    N.df2 <- sem_out_comb(N.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(N.df1, ., by = c("Species", "Region", "Year"))
    N.df <- sem_out_comb(N.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(N.df2, ., by = c("Species", "Region", "Year")) %>% 
      mutate(Year = Year + 1975) %>%
      filter(Year > d.end) %>%
      mutate(Species = recode(Species, 
                              '5' = "SAGR",
                              '1' = "BLGR",
                              '3' = "CHUK",
                              '4' = 'HUPA',
                              '2' = 'CAQU'))
    
    H.df1 <- sem_out_comb(H.Est) %>% dplyr::rename(Est = value)
    H.df2 <- sem_out_comb(H.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(H.df1, ., by = c("Species", "Region", "Year"))
    H.df <- sem_out_comb(H.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(H.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975)  %>%
      filter( Year > d.end)%>%
      mutate(Species = recode(Species, 
                              '5' = "SAGR",
                              '1' = "BLGR",
                              '3' = "CHUK",
                              '4' = 'HUPA',
                              '2' = 'CAQU'))
    
    BPH.df1 <- sem_out_comb(BPH.Est) %>% dplyr::rename(Est = value)
    BPH.df2 <- sem_out_comb(BPH.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(H.df1, ., by = c("Species", "Region", "Year"))
    BPH.df <- sem_out_comb(BPH.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(H.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975) %>%
      filter( Year > d.end)%>%
      mutate(Species = recode(Species, 
                              '5' = "SAGR",
                              '1' = "BLGR",
                              '3' = "CHUK",
                              '4' = 'HUPA',
                              '2' = 'CAQU'))
    
    if(input$forecastset == "BPH"){
      BPH.df 
    }else{
      if(input$forecastset == "H"){
        H.df
      }else{
        N.df
      }
    }
  })
  
  output$forecastplot <- renderPlot(
    ggplot(data = forecast.plot.input(), aes(x = Year, color = Species)) +
      geom_errorbar(aes(ymin = UCL, ymax = LCL)) +
      geom_point(aes(y = Est, color = Species)) +
      geom_line(data = d.plot.input(), aes(y = Est)) +
      facet_wrap(vars(Region), scales = "free_y", nrow = 2)
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
 
  # ### "Download/Input Data"
  # # Downloadable csv of selected dataset ----
  # output$downloadData <- downloadHandler(
  #   filename = function() {
  #     paste(input$dataset, ".csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(read.csv("FullDataDownload.csv"), file, row.names = F)
  #   }
  # )
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