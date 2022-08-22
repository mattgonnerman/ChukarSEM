library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(R.utils)
library(tidyr)
library(lubridate)
library(scales)
library(plotly)
library(DT)

# Load Data from Model
load("./www/NDOW_SEM_app_objects.R")
load("./www/NDOW_ChukOnly_app_objects.R")
lastyear <- ncol(appobject$N)

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
  
  setSliderColor(rep("#12263f",5),c(1:5)),
  
  chooseSliderSkin("Flat"),
  tags$head(tags$style("#forecasttable table {background-color: white;
                       border: 1px solid '#12263f';
                       border-collapse: collapse;
                       }",
                       media="screen", type="text/css"))
  )



##################################################################################
### Setup page with a navigation bar
ui <- navbarPage( tweaks,
                  
                  title = "NDOW Upland Harvest", selected = "Instructions",
                  ###INSTRUCTIONS
                  tabPanel(
                    title = "Instructions",
                    
                    textOutput("appinstrhead"),
                    tags$div(tags$ul(
                      tags$li("ENSURE MODEL IS UP TO DATE [SEE BELOW]"),
                      tags$li("FORECAST ESTIMATES"),
                      tags$li("DOWNLOAD ESTIMATES"),
                      tags$li("CHUKAR ONLY MODEL")),
                      style = "font-size: 15px;
                      color: white;"),
                    tags$head(tags$style("#appinstrhead{color: white;
                                 font-size: 20px;
                                 font-style: bold;
                                 }")),
                    
                    textOutput("datainsthead"),
                    tags$div(tags$ul(
                      tags$li("UPDATE DATA SHEETS"),
                      tags$li("RERUN MODEL"),
                      tags$li("UPDATE SHINY")),
                      style = "font-size: 15px;
                      color: white;"),
                    tags$head(tags$style("#datainsthead{color: white;
                                 font-size: 20px;
                                 font-style: bold;
                                 "))
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
                                   conditionalPanel(
                                     condition = "input.forecastset == 'N'|input.forecastset == 'BPH'",
                                     
                                     sliderInput("pdsi", HTML("Drought Severity (High \u2192 Low)"), min=-3, max=3, 
                                                 value = 0, sep="", step = 0.1)
                                   ),
                                   
                                   downloadButton("downloadData", "Download Displayed Estimates")
                      ),
                      
                      mainPanel(#Output Plot
                        plotOutput("forecastplot", width = "800px", height = "710px"),
                        DT::dataTableOutput('forecasttable')
                      )
                      
                    )
                    
                    
                  ),
                  
                  ### Chukar-Only Model
                  tabPanel(
                    title = "Chukar-Only Model",
                    mainPanel(#Output Plot
                      textOutput("tmptxt")
                    ),
                    tags$head(tags$style("#tmptxt{color: white;
                                 font-size: 20px;
                                 font-style: italic;
                                 }"
                                         )
                              )
                  )
)



##################################################################################
### Server Instructions
server <-  function(input,output,session){
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","PHEA","SAGR")
  colors.constant <- c("#E8A323","#23E8A3","#A323E8","#074A36","#36074A", "#4A3607")
  colors <- reactive(setNames(colors.constant[which(species.constant %in% forecast.plot.input()$Species)],
                              species.constant[which(species.constant %in% forecast.plot.input()$Species)]))
  
  output$tmptxt <- renderText({"Under Construction"})
  
  output$appinstrhead <- renderText({"Using this app:"})
  output$datainsthead <- renderText({"Updating model with new data:"})

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
    
    find.na <- as.data.frame(appobject$data.all$pdsi) %>% mutate(Mean = (Eastern + Western )/ 2)
    t <- max(sapply(find.na, function(x) max(which(!is.na(x)))))
    val <- as.numeric(find.na[t,3])
    updateSliderInput(session, "pdsi", value = val)
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
            econ.var*rnorm(n.samples, H.reg[[r]]$B.econ.hunt[s], H.reg[[r]]$B.econ.hunt.sd[s]) +
            latent[s,r,t,]
          
          CI <- (1-as.numeric(input$conflevel))/2
          
          H.Est[s,r,t] <- round(1000*exp(quantile(mu.H[s,r,t,], c(.5))))
          H.LCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(CI)))
          H.UCL[s,r,t] <- 1000*exp(quantile(mu.H[s,r,t,], c(1-CI)))
          
          # Total Harvest
          mu.N[s,r,t,] <- rnorm(n.samples, N.reg[[r]]$A.harv[s], N.reg[[r]]$A.harv.sd[s]) +
            rnorm(n.samples, N.reg[[r]]$B.hunter[s], N.reg[[r]]$B.hunter.sd[s])*latent[s,r,t,]+
            winter.var * rnorm(n.samples, N.reg[[r]]$B.ws.harv[s], N.reg[[r]]$B.ws.harv.sd[s]) +
            drought.var * rnorm(n.samples, N.reg[[r]]$B.pdsi.harv[s], N.reg[[r]]$B.pdsi.harv.sd[s]) +
            predator.var * rnorm(n.samples, N.reg[[r]]$B.bbs.harv[s], N.reg[[r]]$B.bbs.harv.sd[s])
          
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
  
  CHUK.forecast.plot.input <- reactive({
    #Prep Inputs
    c.H.reg <- appobject.CH$H.reg
    c.N.reg <- appobject.CH$N.reg
    c.Latent <- appobject.CH$latent.trend
    county_order <- appobject.CH$county_order
    county_reg <- appobject.CH$county_reg
    n.counties <- nrow(c.H.reg)
    n.year.chuk <- ncol(appobject.CH$H.data)
    # c.econ.var <- input$c.econ
    # c.winter.var <-  input$c.awssi
    # c.predator.var <- input$c.bbs
    c.econ.var <- 1
    c.winter.var <-  1
    c.predator.var <- 1
    c.n.samples <- 100
    
    #Empty arrays to fill
    c.latent <- c.mu.N <- c.mu.H <- array(NA, dim = c(n.counties, n.year.chuk, c.n.samples))
    c.BPH.Est <- c.BPH.LCL <- c.BPH.UCL <- c.N.Est <- c.N.LCL <- c.N.UCL <- c.H.Est <- c.H.LCL <- c.H.UCL <- matrix(NA, nrow = n.counties, ncol = n.year.chuk)
    
    c.Latent.mean <- c.Latent %>% dplyr::select(County, Year, Latent) %>%
      pivot_wider(names_from = Year, values_from = Latent) %>% dplyr::select(-County)
    c.Latent.sd <- c.Latent %>% dplyr::select(County, Year, Latent.sd) %>%
      pivot_wider(names_from = Year, values_from = Latent.sd) %>% dplyr::select(-County)
    
    for(s in 1:n.counties){
      for(t in 1:n.year.chuk){
        #Latent time trend
        c.latent[s,t,] <- rnorm(c.n.samples, as.numeric(c.Latent.mean[s,t]), as.numeric(c.Latent.sd[s,t]))
        
        #Hunter Effort
        c.mu.H[s,t,] <- rnorm(c.n.samples, c.H.reg$A.hunt[s], c.H.reg$A.hunt.sd[s]) +
          c.econ.var*rnorm(c.n.samples, c.H.reg$B.econ.hunt[s], c.H.reg$B.econ.hunt.sd[s]) +
          c.latent[s,t,]
        
        # c.CI <- (1-as.numeric(input$conflevel))/2
        c.CI <- (1-as.numeric(.95))/2
        
        c.H.Est[s,t] <- round(100*exp(quantile(c.mu.H[s,t,], c(.5))))
        c.H.LCL[s,t] <- 100*exp(quantile(c.mu.H[s,t,], c(c.CI)))
        c.H.UCL[s,t] <- 100*exp(quantile(c.mu.H[s,t,], c(1-c.CI)))
        
        # Total Harvest
        c.mu.N[s,t,] <- rnorm(c.n.samples, c.N.reg$A.harv[s], c.N.reg$A.harv.sd[s]) +
          rnorm(c.n.samples, c.N.reg$B.hunter[s], c.N.reg$B.hunter.sd[s])*c.latent[s,t,]+
          c.winter.var * rnorm(c.n.samples, c.N.reg$B.ws.harv[s], c.N.reg$B.ws.harv.sd[s]) +
          c.predator.var * rnorm(c.n.samples, c.N.reg$B.bbs.harv[s], c.N.reg$B.bbs.harv.sd[s])
        
        c.N.Est[s,t] <- round(100*exp(quantile(c.mu.N[s,t,], c(.5))))
        c.N.LCL[s,t] <- 100*exp(quantile(c.mu.N[s,t,], c(c.CI)))
        c.N.UCL[s,t] <- 100*exp(quantile(c.mu.N[s,t,], c(1-c.CI)))
      }
    }
    
    chuk_out_comb <- function(f.df){
      as.data.frame(f.df) %>%
        mutate(County = 1:n.counties) %>%
        pivot_longer(cols = "County", names_to = "Year", names_prefix = "V") %>%
        as.data.frame()
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
    
    
    if(input$forecastset == "BPH"){
      chuk.f.input <- c.BPH.df 
    }else{
      if(input$forecastset == "H"){
        chuk.f.input <- c.H.df
      }else{
        chuk.f.input <- c.N.df
      }
    }
    
    ### Prepare Real Data
    if(input$forecastset == "BPH"){
      chuk.input.df <- (100*appobject.CH$N.data)/(1+100*appobject.CH$H.data)
    }else{
      if(input$forecastset == "H"){
        chuk.input.df <- 100*appobject.CH$H.data
      }else{
        chuk.input.df <- 100*appobject.CH$N.data
      }
    }
    
    chuk.d.input <- as.data.frame(chuk.input.df) %>% 
      mutate(County = 1:n.counties, ) %>%
      pivot_longer(names_to = "Year", cols = 1:ncol(chuk.input.df), values_to = "Real") %>%
      filter(!is.na(Real)) %>%
      mutate(Year = as.numeric(Year),
             County = as.numeric(County)) %>% 
      mutate(County = factor(County, levels = 1:n.counties, labels = county_order))
    
    
    ### Merge Real Values and Estimates
    merge(chuk.f.input, chuk.d.input, by = c("County","Year"), all = T) %>%
      filter(Year %in% input$years.chuk.plot[1]:input$years.chuk.plot[2]) %>% 
      filter(County %in% input$countyhighlight)
  })
  
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
      # geom_text(aes(x=vline.intercept(), label="\n\nEstimated", y=linetext.pos()),
      #           color = "black", angle=90, alpha = .5) +
      # geom_text(aes(x=vline.intercept(), label="Observered\n\n", y=linetext.pos()),
      #           color = "black", angle=90, alpha = .5) +
      geom_text(data = estrealtext(), aes(x = X, y = Y, label = Text),
                color = "black", angle=90, alpha = .5) +
      facet_wrap(vars(Region), scales = "free_y", nrow = 2) +
      theme_classic(base_size = 18) + 
      labs(y = ylabel.rec()) +
      scale_x_continuous(breaks= pretty_breaks(n = 4)) +
      scale_color_manual(values = colors()) +
      scale_fill_manual(values = colors()) +
      theme(legend.position = c(.15,.9),
            plot.background = element_rect(colour = "#a4a9ac", fill=NA, size=2))+
      guides(color=guide_legend(nrow=3, byrow=TRUE),
             fill=guide_legend(nrow=3, byrow=TRUE)),
    
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
  
  output$forecasttable <- DT::renderDataTable(forecasttable(),
                                              options = list(
                                                scrollY = T,
                                                scrollX = T,
                                                columnDefs = list(list( targets = "_all", width = '15px'))
                                              )
                                      )
  
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