# Install and Load packages
packages <- c("shinyWidgets", "ggplot2", "dplyr", "mvtnorm", "R.utils", "tidyr", "lubridate", "scales")
# '%notin%' <- Negate('%in%')
# for(i in 1:length(packages)){
#   if(packages[i] %notin% installed.packages()){
#     install.packages(packages[i])
#   }
# }

lapply(packages, library, character.only = T)

library(shinyWidgets)
library(ggplot2)
library(dplyr)
library(mvtnorm)
library(R.utils)
library(tidyr)
library(lubridate)
library(scales)
library(plotly)

# Load Data from Model
load("./www/NDOW_SEM_app_objects.R")
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
                                                                  "Blue Grouse" = "BLGR",
                                                                  "Hungarian Partridge" = "HUPA",
                                                                  "Pheasant" = "PHEA",
                                                                  "Sage Grouse" = "SAGR"), 
                                                      selected = c("CHUK"), inline = T),
                                   
                                   sliderInput("years.f.plot", "Years to Display", min =1976, max = (1975 + lastyear),
                                               value = c(2011,(1975 + lastyear)), sep=""),
                                   
                                   sliderInput("econ", "Economic Predictor", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("bbs", "Predator Predictor", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("awssi", "Winter Severity", min=-2, max=2, 
                                               value = 0, sep=""),
                                   sliderInput("pdsi", "Drought Conditions", min=-2, max=2, 
                                               value = 0, sep=""),
                                   
                                   downloadButton("downloadData", "Download Displayed Estimates")
                      ),
                      
                      mainPanel(#Output Plot
                        plotOutput("forecastplot"),
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
server <-  function(input,output){
  
  species.constant <- c("BLGR","CAQU","CHUK","HUPA","SAGR")
  colors.constant <- c("#0033a0","#015b6e","#01833b","#c09b0f","#ffa300")
  colors <- reactive(setNames(colors.constant[seq(1,length(species.constant), length.out = floor(length(unique(past.data()$Species))))],
                              species.constant[which(species.constant %in% unique(past.data()$Species))]))
  
  output$tmptxt <- renderText({"Under Construction"})
  
  output$appinstrhead <- renderText({"Using this app:"})
  output$datainsthead <- renderText({"Updating model with new data:"})
                              
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
                              '2' = 'CHUK',
                              '3' = 'HUPA',
                              '4' = 'PHEA',
                              '5' = 'SAGR'))
    
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
    mu.BPH <- exp(mu.N - mu.H)
    for(r in 1:2){
      for(s in 1:n.species){
        for(t in 1:n.year.f){
          BPH.Est[s,r,t] <- quantile(mu.BPH[s,r,t,], c(.5))
          BPH.LCL[s,r,t] <- quantile(mu.BPH[s,r,t,], c(.075))
          BPH.UCL[s,r,t] <- quantile(mu.BPH[s,r,t,], c(.935))
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
                              '2' = 'CHUK',
                              '3' = 'HUPA',
                              '4' = 'PHEA',
                              '5' = 'SAGR'))
    
    H.df1 <- sem_out_comb(H.Est) %>% dplyr::rename(Est = value)
    H.df2 <- sem_out_comb(H.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(H.df1, ., by = c("Species", "Region", "Year"))
    H.df <- sem_out_comb(H.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(H.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975)  %>%
      mutate(Species = recode(Species, 
                              '1' = 'BLGR',
                              '2' = 'CHUK',
                              '3' = 'HUPA',
                              '4' = 'PHEA',
                              '5' = 'SAGR'))
    
    BPH.df1 <- sem_out_comb(BPH.Est) %>% dplyr::rename(Est = value)
    BPH.df2 <- sem_out_comb(BPH.LCL) %>% dplyr::rename(LCL = value) %>%
      merge(BPH.df1, ., by = c("Species", "Region", "Year"))
    BPH.df <- sem_out_comb(BPH.UCL) %>% dplyr::rename(UCL = value) %>%
      merge(BPH.df2, ., by = c("Species", "Region", "Year"))%>% 
      mutate(Year = Year + 1975) %>%
      mutate(Species = recode(Species,  
                              '1' = 'BLGR',
                              '2' = 'CHUK',
                              '3' = 'HUPA',
                              '4' = 'PHEA',
                              '5' = 'SAGR'))
    
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
  
  vline.intercept <- reactive({
    1975 + min(which(is.na(appobject$N.data[1,,1]))) - 1
  })
  
  linetext.pos <- reactive({
    max(as.numeric(.6 * forecast.plot.input()$UCL.use), na.rm = T)
  })
  
  output$forecastplot <- renderPlot(
    ggplot(data = forecast.plot.input(), aes(x = Year, fill = Species, color = Species)) +
      geom_ribbon(aes(ymin = UCL.use, ymax = LCL.use), alpha = .4, linetype = "dashed") +
      geom_line(aes(y = Graph.Y, color = Species), size = 1.4) +
      geom_vline(xintercept = vline.intercept(), size = 1.5, linetype = "dashed") +
      geom_text(aes(x=vline.intercept(), label="\n\nEstimated", y=linetext.pos()), 
                color = "black", angle=90, alpha = .5) +
      geom_text(aes(x=vline.intercept(), label="Observered\n\n", y=linetext.pos()),
                color = "black", angle=90, alpha = .5) +
      facet_wrap(vars(Region), scales = "free_y", nrow = 2) +
      theme_classic(base_size = 18) + 
      scale_x_continuous(breaks= pretty_breaks(n = 4)) +
      theme(legend.position = "bottom",
            axis.title.y = element_blank()),
    
    height = 700
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