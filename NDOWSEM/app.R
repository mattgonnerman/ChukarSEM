library(shinyWidgets)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(R.utils)
library(tidyr)
library(lubridate)
library(scales)
library(plotly)
library(DT)
library(viridis)
library(leaflet)
library(sf)

# Load Data from Model
load("./www/NDOW_SEM_app_objects.R")
load("./www/NDOW_ChukOnly_app_objects.R")
lastyear <- ncol(appobject$N)
lastyear.ch <- ncol(appobject.CH$N)

#Fix for unaligned checkboxes
tweaks <- tags$head(tags$style(HTML(
  ".checkbox-inline { 
                    margin-left: 0px;
                    margin-right: 10px;
          }
         .checkbox-inline+.checkbox-inline {
                    margin-left: 0px;
                    margin-right: 10px;
          }"),
  
  setBackgroundColor(color = "#12263f")),
  
  setSliderColor(rep("#12263f",10),c(1:10)), # Slider Colors, numbers correspond to number of sliders. Just put more than you need
  
  chooseSliderSkin("Flat"),
  
  tags$head(tags$style(HTML(".nav-tabs {background-color: white}"))), #Set background color of tabs to improve readability

  tags$head(tags$style("#forecasttable table {background-color: white;
                       border: 1px solid '#12263f';
                       border-collapse: collapse;
                       }",
                       media="screen", type="text/css")),
  
  tags$head(tags$style("#chuk.forecasttable table {background-color: white;
                       border: 1px solid '#12263f';
                       border-collapse: collapse;
                       }",
                       media="screen", type="text/css")),
  
  tags$style(type = 'text/css', ".tab-panel{ background-color: white}")
  
  
  
  )

tabletweaks <- "
                    .dataTables_wrapper .dataTables_length, .dataTables_wrapper .dataTables_filter, .dataTables_wrapper .dataTables_info, .dataTables_wrapper .dataTables_processing, .dataTables_wrapper .dataTables_paginate {
                    color: #ffffff;
                    }

                    thead {
                    color: black;
                    background-color: white;
                    }
                    
                    .dataTables_filter input {
                            color: black;
                    }
                           
                    .dataTables_length select {
                           color: black !important; 
                    }

                   "



##################################################################################
### Setup page with a navigation bar
ui <- navbarPage( tweaks,
                  
                  title = "NDOW Upland Harvest", selected = "Instructions",
                  ###INSTRUCTIONS
                  tabPanel(
                    title = "Instructions",
                    
                    sidebarPanel(
                      
                      uiOutput("Report"),
                      
                      br(),
                      
                      textOutput("appfor"),
                      tags$head(tags$style("#appfor{color: black;
                                 font-size: 16px;
                                 font-style: bold;
                                 }")),
                      textOutput("appforbody"),
                      tags$head(tags$style("#appforbody{color: black;
                                 font-size: 12px;
                                 }")),
                      br(),
                      textOutput("modelwork"),
                      tags$head(tags$style("#modelwork{color: black;
                                 font-size: 16px;
                                 font-style: bold;
                                 ")),
                      textOutput("modelworkbody"),
                      tags$head(tags$style("#modelworkbody{color: black;
                                 font-size: 12px;
                                 ")),
                      br(),
                      
                      textOutput("appwork"),
                      tags$div(tags$ol(
                        tags$li("Select whether to display estimates for birds per hunter, hunter participation, or total harvest."),
                        tags$li("Select which species (counties for CHUK model) you wish to display estimates for."),
                        tags$li("Select which range of years to display estimates for."),
                        tags$li("Select what confidence level you wish error to be presented in."),
                        tags$li("Set expected values for predictor covariates. Values are initially set as the most recent observed for each covariate"),
                        tags$li("Download raw estimates.")),
                        style = "font-size: 12px;
                      color: black;"),
                      tags$head(tags$style("#appwork{color: black;
                                 font-size: 16px;
                                 font-style: bold;
                                 "))),
                    
                    mainPanel(img(src='exampleimage.jpg', height= 700))
                  ),
                  
                  
                  ### Harvest Forecast
                  tabPanel(
                    title = "Forecast",
                    
                    sidebarLayout(
                      sidebarPanel(selectInput("forecastset", "Which dataset would you like to see forcast?",
                                               choices = c("Birds per Hunter" = "BPH", 
                                                           "Total Harvest" = "N", 
                                                           "Hunter Participation" = "H")),
                                   
                                   checkboxGroupInput("spphighlight", "Species:",
                                                      choices = c("Chukar" = "CHUK",
                                                                  "California Quail" = "CAQU",
                                                                  "Blue Grouse" = "BLGR",
                                                                  "Hungarian Partridge" = "HUPA",
                                                                  "Pheasant" = "PHEA",
                                                                  "Sage Grouse" = "SAGR"), 
                                                      selected = c("CHUK"), inline = T),
                                   
                                   sliderInput("years.f.plot", "Years to Display", min =1976, max = (1975 + lastyear),
                                               value = c(2011,(1975 + lastyear)), sep=""),
                                   selectInput("conflevel", "Confidence Interval",
                                               choices = c("95%" = .95, 
                                                           "85%" = .85, 
                                                           "50%" = .5, 
                                                           "Hide Confidence Interval" = 0)),
                                   
                                   conditionalPanel(
                                     condition = "input.forecastset == 'H'|input.forecastset == 'BPH'",
                                     
                                     sliderInput("econ", HTML("Economic Strength (Low \u2192 High)"), min=-3, max=3, 
                                                 value = 0, sep="", step = 0.1)
                                   ),
                                   conditionalPanel(
                                     condition = "input.forecastset == 'N'|input.forecastset == 'BPH'",
                                     
                                     sliderInput("bbs", HTML("Predator Abundance (Low \u2192 High)"), min=-3, max=3, 
                                                 value = 0, sep="", step = 0.1)
                                   ),
                                   conditionalPanel(
                                     condition = "input.forecastset == 'N'|input.forecastset == 'BPH'",
                                     
                                     sliderInput("awssi", HTML("Winter Severity (Low \u2192 High)"), min=-3, max=3, 
                                                 value = 0, sep="", step = 0.1)
                                   ),
                                   
                                   downloadButton("downloadData", "Download Displayed Estimates")
                      ),
                      
                      mainPanel(#Output Plot
                        tabBox(
                          tabPanel(id = "plot.full", title = "Plot", plotOutput("forecastplot")),
                          tabPanel(id = "table.full", title = "Table", DT::dataTableOutput('forecasttable'),
                                   tags$style(HTML(tabletweaks))),
                        )
                        , width = "800px"
                      )
                    ) 
                  ),
                  
                  ### Chukar-Only Model
                  tabPanel(
                    title = "Chukar County-Level Model",

                    sidebarLayout(
                      sidebarPanel(

                        selectInput("chuk.forecastset", "Which dataset would you like to see forcast?",
                                    choices = c("Birds per Hunter" = "BPH",
                                                "Total Harvest" = "N",
                                                "Hunter Participation" = "H")),
                        
                        checkboxGroupInput("countyhighlight", "Counties:",
                                           choices = appobject.CH$county_order,
                                           selected = "Carson City", inline = T),
                        
                        conditionalPanel(condition="input.chuktabs == 'Plot'|input.chuktabs == 'Table'",
                                         sliderInput("years.chuk.plot", "Years to Display", min =1990, max = (1989 + lastyear.ch),
                                                     value = c(2011,(1989 + lastyear.ch)), sep=""),
                                         selectInput("c.conflevel", "Confidence Interval",
                                                     choices = c("95%" = .95,
                                                                 "85%" = .85,
                                                                 "50%" = .5,
                                                                 "Hide Confidence Interval" = 0))
                        ),
                        
                        conditionalPanel(condition="input.chuktabs == 'Map'",
                                         sliderInput("years.chuk.map", "Year to Display", min =1990, max = (1989 + lastyear.ch),
                                                     value = (1989 + lastyear.ch), sep="")
                        ),
                        
                        conditionalPanel(
                          condition = "input.forecastset == 'H'|input.forecastset == 'BPH'",
                          
                          sliderInput("c.econ", HTML("Economic Strength (Low \u2192 High)"), min=-3, max=3,
                                      value = 0, sep="", step = 0.1)
                        ),
                        conditionalPanel(
                          condition = "input.forecastset == 'N'|input.forecastset == 'BPH'",
                          
                          sliderInput("c.bbs", HTML("Predator Abundance (Low \u2192 High)"), min=-3, max=3,
                                      value = 0, sep="", step = 0.1)
                        ),
                        conditionalPanel(
                          condition = "input.forecastset == 'N'|input.forecastset == 'BPH'",
                          
                          sliderInput("c.awssi", HTML("Winter Severity (Low \u2192 High)"), min=-3, max=3,
                                      value = 0, sep="", step = 0.1)
                        ),
                        
                        downloadButton("downloadData.c", "Download Displayed Estimates")
                      ),

                      mainPanel(
                        tabBox(id = "chuktabs",
                          tabPanel(id = "plot.chuk", title = "Plot", plotOutput("chukplot", height = "710px", width = "800px")),
                          tabPanel(id = "map.chuk", title = "Map", leafletOutput("chukmap", height = "600px", width = "600px")),
                          tabPanel(id = "table.chuk", title = "Table", DT::dataTableOutput('chuk.forecasttable'),
                                   tags$style(HTML(tabletweaks)))
                        )
                        , width = "800px")
                    )
                  )
)



##################################################################################
### Server Instructions
server <-  function(input,output,session){
  
  output$Report<- renderUI({
    tags$a(href="report.pdf", target="_blank", 
           "To read a detailed description of the model and results, please follow this link to view our full report.",
           style = "font-size: 16px; color: black;")
  })
  
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","PHEA","SAGR")
  colors.constant <- c("#E8A323","#23E8A3","#A323E8","#074A36","#36074A", "#4A3607")
  colors <- reactive(setNames(colors.constant[which(species.constant %in% forecast.plot.input()$Species)],
                              species.constant[which(species.constant %in% forecast.plot.input()$Species)]))
  
  colors.codes.chuk <- viridis(length(appobject.CH$county_order))
  colors.chuk <- reactive(setNames(colors.codes.chuk[which(appobject.CH$county_order %in% CHUK.forecast.plot.input()$County)],
                                   appobject.CH$county_order[which(appobject.CH$county_order %in% CHUK.forecast.plot.input()$County)]))
  
  output$tmptxt <- renderText({"Under Construction"})
  
  output$appfor <- renderText({"What is this app for?"})
  output$appforbody <- renderText({"Upland game birds are popular quarry among hunters but traditionally have received less monitoring attention than other game such as waterfowl or big game.  The dynamics of many upland game populations are poorly understood as a result, which both complicates management and provides a dearth of information to stakeholders, including hunters, about population status.  In Nevada, annual surveys of hunters provide information on hunter effort and total harvest for each of the state’s upland game species.  For a few species, specifically chukar and greater sage-grouse, detailed historic population survey data also exist.  If upland game bird populations respond similarly to background environmental factors (e.g., harsh winter weather) among species, it is likely this suite of data can be used to infer historic population dynamics and predict future patterns based on environment-covariate relationships and shared variance among data sources. "})
  output$modelwork <- renderText({"How does the model work?"})
  output$modelworkbody <- renderText({"We employed a Bayesian SEM framework to estimate the correlated variation in total harvest (N) and hunter participation (H) for upland game birds. We accounted for correlation among species in two ways. 1) We used a multivariate normal distribution to describe N and H, as it allowed for estimating the correlation structure, and 2) used common priors among species to define beta coefficients used to describe variation in N and H. We linked variation between H and N by allowing N to vary according to the latent trend in H, absent variation associated with identified covariates (e.g., economic conditions). The model is comprised of sub-models describing variation in hunter participation, total harvest, chukar site abundance, and variation in predictor variables. We applied this model in two ways; 1) as a regional multi-species approach applied to six species of upland game birds in the northeast and northwest regions of Nevada and 2) a county-specific single-species approach applied to chukar exclusively. "})
  output$appwork <- renderText({"How to use this app? (reference image on right)"})

  observe({
    find.na <- data.frame(A = appobject$data.all$une,
                          B = appobject$data.all$gas,
                          C = appobject$data.all$pdi,
                          D = appobject$data.all$res)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(appobject$pred.econ[t,2])
    updateSliderInput(session, "econ", value = val)
    
    find.na <- data.frame(A = appobject$data.all$ravens,
                          B = appobject$data.all$rthawk,
                          C = appobject$data.all$nharr,
                          D = appobject$data.all$pfal)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(appobject$pred.bbs[t,2])
    updateSliderInput(session, "bbs", value = val)
    
    find.na <- as.data.frame(appobject$data.all$awssi) %>% mutate(Mean = (Eastern + Western )/ 2)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(find.na[t,3])
    updateSliderInput(session, "awssi", value = val)
    

    find.na <- data.frame(A = appobject.CH$data.all$une,
                          B = appobject.CH$data.all$gas,
                          C = appobject.CH$data.all$pdi,
                          D = appobject.CH$data.all$res)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(appobject.CH$pred.econ[t,2])
    updateSliderInput(session, "c.econ", value = val)
    
    find.na <- data.frame(A = appobject.CH$data.all$ravens,
                          B = appobject.CH$data.all$rthawk,
                          C = appobject.CH$data.all$nharr,
                          D = appobject.CH$data.all$pfal)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(appobject.CH$pred.bbs[t,2])
    updateSliderInput(session, "c.bbs", value = val)
    
    find.na <- as.data.frame(appobject.CH$data.all$awssi) %>% mutate(Mean = (Eastern + Western )/ 2)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(find.na[t,3])
    updateSliderInput(session, "c.awssi", value = val)
  })

                              
  forecast.plot.input.nodelay <- reactive({
    #Prep Inputs
    H.reg <- appobject$H.reg
    N.reg <- appobject$N.reg
    Latent <- appobject$latent.trend
    n.species <- nrow(H.reg[[1]])
    econ.var <- input$econ
    winter.var <-  input$awssi
    predator.var <- input$bbs
    
    d.start <- 1976
    f.start <- min(which(is.na(appobject$N.data[1,,1]))) + 1975
    d.end <- f.start - 1
    f.end <- ncol(appobject$N.data) + 1975
    y.b <- f.start-1976 #years of data before forecast
    n.year.f <- length(f.end:d.start) #If you want to show more than next year, use this as another dimension in array
    n.samples <- 10000
    
    
    #Function for formatting forecast graph inputs
    sem_out_comb <- function(f.arr){
      as.data.frame(wrap(f.arr, map=list(NA, 2))) %>%
        mutate(Species = rep(1:n.species, n.year.f),
               Year = rep(1:n.year.f, each = n.species)) %>%
        pivot_longer(cols = 1:2, names_to = "Region", names_prefix = "V") %>%
        as.data.frame()
    }
    
    ### Prepare Real Data
    if(input$forecastset == "BPH"){
      d.input.array <- appobject$N.data/(1+appobject$H.data)
    }else{
      if(input$forecastset == "H"){
        d.input.array <- 1000*appobject$H.data
      }else{
        d.input.array <- 1000*appobject$N.data
      }
    }
    
    d.input <- as.data.frame.table(d.input.array) %>% 
      dplyr::rename(Species = Var1, Year = Var2, Region = Var3, Real = Freq) %>%
      filter(!is.na(Real)) %>%
      mutate(Year = as.numeric(Year),
             Region = as.numeric(Region),
             Species = as.numeric(Species)) %>%
      mutate(Year = Year + 1975)  %>%
      mutate(Species = recode(Species, 
                              '1' = 'BLGR',
                              '2' = 'CAQU',
                              '3' = 'CHUK',
                              '4' = 'HUPA',
                              '5' = 'PHEA',
                              '6' = 'SAGR'))
    
    #Empty arrays to fill
    latent <- mu.N <- mu.H <- array(NA, dim = c(n.species, 2, n.year.f, n.samples))
    BPH.Est <- BPH.LCL <- BPH.UCL <- N.Est <- N.LCL <- N.UCL <- H.Est <- H.LCL <- H.UCL <- array(NA, dim = c(n.species, 2, n.year.f))
    
    for(r in 1:2){
      for(s in 1:n.species){
        for(t in 1:n.year.f){
          #Latent time trend
          latent[s,r,t,] <- rnorm(n.samples, as.numeric(Latent$mean[[r]][s,t]), as.numeric(Latent$sd[[r]][s,t]))
          
          #Hunter Effort
          mu.H[s,r,t,] <- rnorm(n.samples, H.reg[[r]]$A.hunt[s], H.reg[[r]]$A.hunt.sd[s]) +
            ifelse(is.na(appobject$pred.econ[t,2]), econ.var, appobject$pred.econ[t,2])*rnorm(n.samples, H.reg[[r]]$B.econ.hunt[s], H.reg[[r]]$B.econ.hunt.sd[s]) +
            latent[s,r,t,]
          
          CI <- (1-as.numeric(input$conflevel))/2
          
          H.Est[s,r,t] <- round(1000*exp(quantile(mu.H[s,r,t,], c(.5))))
          H.LCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(CI)))
          H.UCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(1-CI)))
          
          # Total Harvest
          mu.N[s,r,t,] <- rnorm(n.samples, N.reg[[r]]$A.harv[s], N.reg[[r]]$A.harv.sd[s]) +
            rnorm(n.samples, N.reg[[r]]$B.hunter[s], N.reg[[r]]$B.hunter.sd[s])*latent[s,r,t,]+
            ifelse(is.na(appobject$data.all$awssi[t,r]), winter.var, appobject$data.all$awssi[t,r]) * rnorm(n.samples, N.reg[[r]]$B.ws.harv[s], N.reg[[r]]$B.ws.harv.sd[s]) +
            ifelse(is.na(appobject$pred.bbs[t,2]), predator.var, appobject$pred.bbs[t,2]) * rnorm(n.samples, N.reg[[r]]$B.bbs.harv[s], N.reg[[r]]$B.bbs.harv.sd[s])
          
          N.Est[s,r,t] <- round(1000*exp(quantile(mu.N[s,r,t,], c(.5))))
          N.LCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(CI)))
          N.UCL[s,r,t] <- 1000*exp(quantile(mu.N[s,r,t,], c(1-CI)))
        }
      }
    }
    
    #Birds Per Hunter
    mu.BPH <- exp(mu.N - mu.H)
    for(r in 1:2){
      for(s in 1:n.species){
        for(t in 1:n.year.f){
          BPH.Est[s,r,t] <- quantile(mu.BPH[s,r,t,], c(.5))
          BPH.LCL[s,r,t] <- quantile(mu.BPH[s,r,t,], c(CI))
          BPH.UCL[s,r,t] <- quantile(mu.BPH[s,r,t,], c(1-CI))
        }
      }
    }
    
    N.df1 <- sem_out_comb(N.Est) %>% dplyr::rename(Est = value)
    N.df2 <- sem_out_comb(N.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(N.df1, ., by = c("Species", "Region", "Year"))
    N.df <- sem_out_comb(N.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(N.df2, ., by = c("Species", "Region", "Year")) %>% 
      mutate(Year = Year + 1975) %>%
      mutate(Species = recode(Species, 
                              '1' = 'BLGR',
                              '2' = 'CAQU',
                              '3' = 'CHUK',
                              '4' = 'HUPA',
                              '5' = 'PHEA',
                              '6' = 'SAGR'))
    
    H.df1 <- sem_out_comb(H.Est) %>% dplyr::rename(Est = value)
    H.df2 <- sem_out_comb(H.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(H.df1, ., by = c("Species", "Region", "Year"))
    H.df <- sem_out_comb(H.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(H.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975)  %>%
      mutate(Species = recode(Species, 
                              '1' = 'BLGR',
                              '2' = 'CAQU',
                              '3' = 'CHUK',
                              '4' = 'HUPA',
                              '5' = 'PHEA',
                              '6' = 'SAGR'))
    
    BPH.df1 <- sem_out_comb(BPH.Est) %>% dplyr::rename(Est = value)
    BPH.df2 <- sem_out_comb(BPH.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(BPH.df1, ., by = c("Species", "Region", "Year"))
    BPH.df <- sem_out_comb(BPH.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(BPH.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975) %>%
      mutate(Species = recode(Species, 
                              '1' = 'BLGR',
                              '2' = 'CAQU',
                              '3' = 'CHUK',
                              '4' = 'HUPA',
                              '5' = 'PHEA',
                              '6' = 'SAGR'))
    
    if(input$forecastset == "BPH"){
      f.input <- BPH.df 
    }else{
      if(input$forecastset == "H"){
        f.input <- H.df
      }else{
        f.input <- N.df
      }
    }
    
    merge(f.input, d.input, by = c("Species", "Region", "Year"), all = T) %>%
      mutate(LCL.use = ifelse(is.na(Real) & Year > 2000, LCL, NA),
             UCL.use = ifelse(is.na(Real) & Year > 2000, UCL, NA),
             Graph.Y = ifelse(is.na(Real), Est, Real),
             Region = as.factor(ifelse(Region == 1, "East", "West")),
             LCL.use = ifelse(Year == vline.intercept(), Real, LCL.use),
             UCL.use = ifelse(Year == vline.intercept(), Real, UCL.use)) %>%
      filter(Year %in% input$years.f.plot[1]:input$years.f.plot[2]) %>% 
      filter(Species %in% input$spphighlight)
  })
  
  forecast.plot.input <- forecast.plot.input.nodelay %>% debounce(100)
  
  CHUK.forecast.plot.input.nodelay1 <- reactive({
    #Prep Inputs
    c.H.reg <- appobject.CH$H.reg
    c.N.reg <- appobject.CH$N.reg
    c.Latent <- appobject.CH$latent.trend
    county_order <- appobject.CH$county_order
    county_reg <- appobject.CH$county_reg
    n.counties <- nrow(c.H.reg)
    n.year.chuk <- ncol(appobject.CH$H.data)
    c.econ.var <- input$c.econ
    c.winter.var <-  input$c.awssi
    c.predator.var <- input$c.bbs
    c.n.samples <- 10000

    #Empty arrays to fill
    c.latent <- c.mu.N <- c.mu.H <- array(NA, dim = c(n.counties, n.year.chuk, c.n.samples))
    c.BPH.Est <- c.BPH.LCL <- c.BPH.UCL <- c.N.Est <- c.N.LCL <- c.N.UCL <- c.H.Est <- c.H.LCL <- c.H.UCL <- matrix(NA, nrow = n.counties, ncol = n.year.chuk)

    c.Latent.mean <- c.Latent %>% dplyr::select(County, Year, Latent) %>%
      pivot_wider(names_from = Year, values_from = Latent) %>% dplyr::select(-County)
    c.Latent.sd <- c.Latent %>% dplyr::select(County, Year, Latent.sd) %>%
      pivot_wider(names_from = Year, values_from = Latent.sd) %>% dplyr::select(-County)

    county_reg <- c(1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 2)
    
    for(s in 1:n.counties){
      for(t in 1:n.year.chuk){
        #Latent time trend
        c.latent[s,t,] <- rnorm(c.n.samples, as.numeric(c.Latent.mean[s,t]), as.numeric(c.Latent.sd[s,t]))

        #Hunter Effort
        c.mu.H[s,t,] <- rnorm(c.n.samples, c.H.reg$A.hunt[s], c.H.reg$A.hunt.sd[s]) +
          ifelse(is.na(appobject.CH$pred.econ[t,2]), c.econ.var, appobject.CH$pred.econ[t,2])*rnorm(c.n.samples, c.H.reg$B.econ.hunt[s], c.H.reg$B.econ.hunt.sd[s]) +
          c.latent[s,t,]

        c.CI <- (1-as.numeric(input$c.conflevel))/2

        c.H.Est[s,t] <- round(1000*exp(quantile(c.mu.H[s,t,], c(.5))))
        c.H.LCL[s,t] <- 1000*exp(quantile(c.mu.H[s,t,], c(c.CI)))
        c.H.UCL[s,t] <- 1000*exp(quantile(c.mu.H[s,t,], c(1-c.CI)))

        # Total Harvest
        c.mu.N[s,t,] <- rnorm(c.n.samples, c.N.reg$A.harv[s], c.N.reg$A.harv.sd[s]) +
          rnorm(c.n.samples, c.N.reg$B.hunter[s], c.N.reg$B.hunter.sd[s])*c.latent[s,t,]+
          ifelse(is.na(appobject.CH$data.all$awssi[t,county_reg[s]]), c.winter.var, appobject.CH$data.all$awssi[t,county_reg[s]]) * rnorm(c.n.samples, c.N.reg$B.ws.harv[s], c.N.reg$B.ws.harv.sd[s]) +
          ifelse(is.na(appobject.CH$pred.bbs[t,2]), c.predator.var, appobject.CH$pred.bbs[t,2]) * rnorm(c.n.samples, c.N.reg$B.bbs.harv[s], c.N.reg$B.bbs.harv.sd[s])

        c.N.Est[s,t] <- round(1000*exp(quantile(c.mu.N[s,t,], c(.5))))
        c.N.LCL[s,t] <- 1000*exp(quantile(c.mu.N[s,t,], c(c.CI)))
        c.N.UCL[s,t] <- 1000*exp(quantile(c.mu.N[s,t,], c(1-c.CI)))
      }
    }

    c.N.df1 <- as.data.frame(c.N.Est) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(Est = value)
    c.N.df2 <- as.data.frame(c.N.LCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(LCL = value) %>%
      merge(c.N.df1, ., by = c("County", "Year"))
    c.N.df <- as.data.frame(c.N.UCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(UCL = value) %>%
      merge(c.N.df2, ., by = c("County", "Year")) %>%
      mutate(County = factor(County, levels = 1:n.counties, labels = county_order)) %>%
      mutate(Year = 1989 + as.numeric(as.character(Year))) %>%
      arrange(County, Year)

    c.H.df1 <- as.data.frame(c.H.Est) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(Est = value)
    c.H.df2 <- as.data.frame(c.H.LCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(LCL = value) %>%
      merge(c.H.df1, ., by = c("County", "Year"))
    c.H.df <- as.data.frame(c.H.UCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(UCL = value) %>%
      merge(c.H.df2, ., by = c("County", "Year")) %>%
      mutate(County = factor(County, levels = 1:n.counties, labels = county_order)) %>%
      mutate(Year = 1989 + as.numeric(as.character(Year))) %>%
      arrange(County, Year)


    #Birds Per Hunter
    c.mu.BPH <- exp(c.mu.N - c.mu.H)
    for(s in 1:n.counties){
      for(t in 1:n.year.chuk){
        c.BPH.Est[s,t] <- quantile(c.mu.BPH[s,t,], c(.5))
        c.BPH.LCL[s,t] <- quantile(c.mu.BPH[s,t,], c(c.CI))
        c.BPH.UCL[s,t] <- quantile(c.mu.BPH[s,t,], c(1-c.CI))
      }
    }

    c.BPH.df1 <- as.data.frame(c.BPH.Est) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(Est = value)
    c.BPH.df2 <- as.data.frame(c.BPH.LCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(LCL = value) %>%
      merge(c.BPH.df1, ., by = c("County", "Year"))
    c.BPH.df <- as.data.frame(c.BPH.UCL) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(cols = 1:n.year.chuk, names_to = "Year", names_prefix = "V") %>%
      dplyr::rename(UCL = value) %>%
      merge(c.BPH.df2, ., by = c("County", "Year")) %>%
      mutate(County = factor(County, levels = 1:n.counties, labels = county_order)) %>%
      mutate(Year = 1989 + as.numeric(as.character(Year))) %>%
      arrange(County, Year)


    if(input$chuk.forecastset == "BPH"){
      chuk.f.input <- c.BPH.df
    }else{
      if(input$chuk.forecastset == "H"){
        chuk.f.input <- c.H.df
      }else{
        chuk.f.input <- c.N.df
      }
    }

    ### Prepare Real Data
    if(input$chuk.forecastset == "BPH"){
      chuk.input.df <- (1000*appobject.CH$N.data)/(1+(1000*appobject.CH$H.data))
    }else{
      if(input$chuk.forecastset == "H"){
        chuk.input.df <- 1000*appobject.CH$H.data
      }else{
        chuk.input.df <- 1000*appobject.CH$N.data
      }
    }

    chuk.d.input <- as.data.frame(chuk.input.df) %>%
      mutate(County = 1:n.counties) %>%
      pivot_longer(names_to = "Year", cols = 1:ncol(chuk.input.df), values_to = "Real") %>%
      filter(!is.na(Real)) %>%
      mutate(Year = as.numeric(Year),
             County = as.numeric(County)) %>%
      mutate(County = factor(County, levels = 1:n.counties, labels = county_order))


    ### Merge Real Values and Estimates
    merge(chuk.f.input, chuk.d.input, by = c("County","Year"), all = T)  %>%
      filter(County %in% input$countyhighlight)
  })
  
  CHUK.forecast.plot.input.nodelay <- reactive({CHUK.forecast.plot.input.nodelay1() %>%
    filter(Year %in% input$years.chuk.plot[1]:input$years.chuk.plot[2])})
  
  CHUK.forecast.plot.input <- CHUK.forecast.plot.input.nodelay %>% debounce(100)
  
  vline.intercept <- reactive({
    1975 + min(which(is.na(appobject$N.data[1,,1]))) - 1
  })
  
  linetext.pos <- reactive({
    reg1 <- forecast.plot.input() %>% filter(Region == 1)
    reg2 <- forecast.plot.input() %>% filter(Region == 2)
    
    c(max(as.numeric(reg1$UCL.use), na.rm = T),max(as.numeric(reg2$UCL.use), na.rm = T))
  })
  
  ylabel.rec <- reactive({
    ifelse(input$forecastset == "BPH", "Birds Per Hunter",
           ifelse(input$forecastset == "N", "Total Harvest", "Hunter Participation"))
  })
  
 estrealtext <- reactive({
   reg1 <- forecast.plot.input() %>% filter(Region == "East")
   reg2 <- forecast.plot.input() %>% filter(Region == "West")

   data.frame(Species = NA,
              Region = c("East", "West", "East", "West"),
              Text = c("\n\nEstimated", "\n\nEstimated", "Observed\n\n", "Observed\n\n"),
              X =  1975 + min(which(is.na(appobject$N.data[1,,1]))) - 1,
              Y = c(.8*max(as.numeric(reg1$UCL.use), na.rm = T),.8*max(as.numeric(reg2$UCL.use), na.rm = T),.8*max(as.numeric(reg1$UCL.use), na.rm = T),.8*max(as.numeric(reg2$UCL.use), na.rm = T)))
 })
  
  output$forecastplot <- renderPlot(
    ggplot(data = forecast.plot.input(), aes(x = Year, fill = Species, color = Species)) +
      geom_ribbon(aes(ymin = UCL.use, ymax = LCL.use), alpha = .4, linetype = "dashed") +
      geom_line(aes(y = Graph.Y, color = Species), size = 1.4) +
      geom_vline(xintercept = vline.intercept(), size = 1.5, linetype = "dashed") +
      geom_text(data = estrealtext(), aes(x = X, y = Y, label = Text),
                color = "black", angle=90, alpha = .5) +
      facet_wrap(vars(Region), scales = "free_y", nrow = 2) +
      theme_classic(base_size = 18) + 
      labs(y = ylabel.rec()) +
      scale_x_continuous(breaks= pretty_breaks(n = 4)) +
      scale_color_manual(values = colors()) +
      scale_fill_manual(values = colors()) +
      theme(legend.position = "right",
            plot.background = element_rect(colour = "#a4a9ac", fill=NA, size=2))+
      guides(color=guide_legend(nrow=6, byrow=TRUE),
             fill=guide_legend(nrow=6, byrow=TRUE)),
    
    height = 700,
    width = 800
  )
  
  
  c.ylabel.rec <- reactive({
    ifelse(input$chuk.forecastset == "BPH", "Birds Per Hunter",
           ifelse(input$chuk.forecastset == "N", "Total Harvest", "Hunter Participation"))
  })

  output$chukplot <- renderPlot(
    ggplot(data = CHUK.forecast.plot.input(), aes(x = Year, fill = County, color = County)) +
      geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = .4, linetype = "dashed") +
      geom_line(aes(y = Est, color = County), size = 1.4) +
      geom_point(aes(y = Real)) +
      theme_classic(base_size = 18) +
      labs(y = c.ylabel.rec()) +
      scale_x_continuous(breaks= pretty_breaks(n = 4)) +
      scale_color_manual(values = colors.chuk()) +
      scale_fill_manual(values = colors.chuk()) +
      theme(legend.position = "right",
            plot.background = element_rect(colour = "#a4a9ac", fill=NA, size=2))+
      guides(color=guide_legend(nrow=7, byrow=TRUE),
             fill=guide_legend(nrow=7, byrow=TRUE)),

    height = 700,
    width = 800
  )
  
  
  forecasttable <- reactive({
    if(input$forecastset == "BPH"){
      forecast.plot.input() %>%
        dplyr::select(Year, Value = Graph.Y, Species, Region) %>%
        group_by(Region) %>%
        mutate(Value = round(Value, digits = 2)) %>%
        pivot_wider(names_from = "Year", values_from = "Value") %>%
        arrange(Region, Species)
    }else{
      forecast.plot.input() %>%
        dplyr::select(Year, Value = Graph.Y, Species, Region) %>%
        mutate(Value = as.integer(Value)) %>%
        group_by(Region) %>%
        mutate(Value = round(Value, digits = 2)) %>%
        pivot_wider(names_from = "Year", values_from = "Value") %>%
        arrange(Region, Species)
    }
    
  })
  
  chuk.forecasttable <- reactive({
    if(input$chuk.forecastset == "BPH"){
      CHUK.forecast.plot.input() %>%
        dplyr::select(Year, Value = Est, County) %>%
        mutate(Value = round(Value, digits = 2)) %>%
        pivot_wider(names_from = "Year", values_from = "Value") %>%
        arrange(County)
    }else{
      CHUK.forecast.plot.input() %>%
        dplyr::select(Year, Value = Est, County) %>%
        mutate(Value = as.integer(Value)) %>%
        mutate(Value = round(Value, digits = 2)) %>%
        pivot_wider(names_from = "Year", values_from = "Value") %>%
        arrange(County)
    }

  })
  
  output$forecasttable <- DT::renderDataTable(forecasttable(),
                                              options = list(
                                                scrollY = T,
                                                scrollX = T,
                                                columnDefs = list(list( targets = "_all", width = '15px'))
                                              )
                                      )
  
  output$chuk.forecasttable <- DT::renderDataTable(chuk.forecasttable(),
                                              options = list(
                                                scrollY = T,
                                                scrollX = T,
                                                columnDefs = list(list( targets = "_all", width = '15px'))
                                              )
  )
  
  ### CHUKAR MAP
  chuk.county.sf <- reactive({
    st_read("./www/nv_county.shp") %>% 
      dplyr::select(County = NAME, geometry) %>%
      filter(County %in% appobject.CH$county_order) %>%
      arrange(County) %>%
      merge(., CHUK.forecast.plot.input.nodelay1() %>% filter(Year == input$years.chuk.map), by = "County") %>%
      st_transform(4326) %>%
      mutate(Est = round(Est, 3))
  })
  
  output$chukmap <- renderLeaflet({
    pal <- colorNumeric(
      palette = "magma",
      domain = CHUK.forecast.plot.input()$Est)
    
    
    
   leaflet(chuk.county.sf()) %>%
      setView(lng=-117.224121, lat=38.376019 , zoom=6)%>%
      addTiles() %>%
      addPolygons(stroke = F, smoothFactor = .8, fillOpacity = .8,
                  color = ~pal(Est), label = ~Est) %>%
     addLegend("bottomright", pal = pal, values = ~Est,
               title = input$chuk.forecastset,
               opacity = 1)
  })

  ### "Download/Input Data"
  saveestimates <- reactive({forecast.plot.input() %>% 
      dplyr::select(-LCL.use, -UCL.use, -Graph.Y)
    })
  
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("NDOW_SEM_", input$forecastset ,".csv", sep = "")
    },
    content = function(file) {
      write.csv(saveestimates(),
                file, row.names = F)
    }
  )
}


##################################################################################
### Construct Shiny App
shinyApp(ui = ui, server = server)