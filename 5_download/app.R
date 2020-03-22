# 5_download/app.R
library(quantmod) # contains the findValleys function, which maybe we should just extract and put verbatim in a source file instead of loading this whole thing...?
library(minpack.lm) # contains the nlsLM function, which we use for our fitting
library(modelr) # used in both the data modeling and the analysis model fitting 
library(SciViews) # contains the ln function used in the data modeling
library(signal) # contains the savistky golay filter (savgolfilt), used to generate the first derivative data in both data modeling and analysis model fitting  


library(shinyBS) # drop-down panels
library(tidyverse) #  handling data structures and plotting

# source("support_scripts/upload_formatters.R")
# source("support_scripts/layout_handling.R")
# source("support_scripts/plotting.R")
# source("support_scripts/analysis.R")

library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(rhandsontable) # user-interactive tables 
library(shiny) # for shiny web-apps

# User interface -------------------------------------------------------------------------------
ui <- fluidPage(
                                                 tabPanel("Download results", value = "downloads_tab",
                                                          # Sidebar layout with input and output definitions
                                                          sidebarLayout(
                                                              # p("Uploading instructions", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                              # p("Upload raw RFU data by either of the following two methods:", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"),
                                                              # Sidebar panel for inputs 
                                                              sidebarPanel(p("Select and preview downloads", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center") %>% strong(), #HTML("<h2><strong>Select and preview downloads</h2></strong>"),
                                                                           p("Plots can be downloaded directly from the analysis window.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"), #HTML("<i><h3>To download plots, right-click on them directly in the analysis window and select 'save'. </i></h3>"),
                                                                           # Input: Choose dataset 
                                                                           selectInput("dataset1", p("Quick results, averaged by user-defined replicates", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left") %>% strong(),
                                                                                       choices = c("Tma by dRFU",
                                                                                                   "Tma by best fit",
                                                                                                   "Replicate-averaged raw data",
                                                                                                   "RFU data with fits")

                                                                                       
                                                                                       
                                                                                       ),
                                                                           textInput("dataset1_download_name", p("Set file name for Quick Result", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left") , value = Sys.Date()),
                                                                           downloadButton("downloadData1", "Download quick result"),
                                                                           p("----", style = "font-family: 'Avenir Next'; font-size: 20px; color: white",align = "center"),
                                                                           
                                                                           selectInput("dataset2", p("Supplemental files", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left") %>% strong(),
                                                                                       choices = c("Tma by dRFU, no replicate averaging",
                                                                                                   "Tma by best fit, no replicate averaging",
                                                                                                   "Reformatted raw data",
                                                                                                   "Normalized raw data",
                                                                                                   "First derivative of raw data") 
                                                                                       ), 
                                                                           
                                                                           textInput("dataset2_download_name", p("Set file name for Supplemental File", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"), value = Sys.Date()),
                                                                           downloadButton("downloadData2", "Download supplemental file"),
                                                                           p("----", style = "font-family: 'Avenir Next'; font-size: 20px; color: white",align = "center"),
                                                                           
                                                                           selectInput("dataset3", p("R outputs (.rmd)", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left") %>% strong(),
                                                                                       choices = c("Labeled, nested data" ,
                                                                                       "All fitted models, full results" ,
                                                                                       "All fitted models, BIC rankings",
                                                                                       "All fitted models, Tmas" ,
                                                                                       "All fitted models, all outputs")), 
                                                                           
                                                                           textInput("dataset3_download_name", p("Set file name for rmd", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"), value = Sys.Date()),
                                                                           downloadButton("downloadData3", "Download supplemental file"),
                                                                           p("----", style = "font-family: 'Avenir Next'; font-size: 20px; color: white",align = "center")
                                                              ), # end sidebarPanel 
                                                              
                                                              # Main panel for displaying outputs 
                                                              mainPanel(
                                                                        tabsetPanel(
                                                                            tabPanel("Preview: quick downloads", dataTableOutput("table_set1"), style = "overflow-x: scroll;"),
                                                                            tabPanel("Preview: supplemental files", dataTableOutput("table_set2"), style = "overflow-x: scroll;")
                                                                        )
                                                              ) # end main panel
                                                          ) # end sidebarLayout
                                                 )
                                                 
) # end fluidPage, end ui

# ================================================================================================================================================================================================
# ================================================================================================================================================================================================
# ================================================================================================================================================================================================

server <- function(input, output, session) {
    # Reactive value for selected dataset ----
    datasetInput1 <- reactive({
        switch(input$dataset1,
               "Tma by dRFU" = head(mtcars),
               "Tma by best fit" = head(mtcars),
               "Replicate-averaged raw data" = head(mtcars),
               "RFU data with fits" = head(mtcars)
               )
    })

    datasetInput2 <- reactive({
        switch(input$dataset2,
               "Tma by dRFU, no replicate averaging"= head(mtcars) ,
                 "Tma by best fit, no replicate averaging"= head(mtcars) ,
                 "Reformatted raw data"= head(mtcars) ,
                 "Normalized raw data"= head(mtcars) ,
                 "First derivative of raw data" = head(mtcars) 
             
        )
    })
    
    datasetInput3 <- reactive({ # raw outputs from DSFworld
        switch(input$dataset3,
               "Labeled, nested data" = values$df,
               "All fitted models, full results" = values$df_models,
               "All fitted models, BIC rankings" = values$df_BIC_models,
               "All fitted models, Tmas" = values$df_tm_models,
               "All fitted models, all outputs" = values$model_list
        )
    })
    
    # make some download formats
    
    df_mean_sd <- reactive({values_df %>% # the raw data
                            unnest(data) %>%
                            select(c(condition, Temperature, mean, sd)) %>%
                            distinct(.keep_all = TRUE) %>%
                            pivot_wider(names_from = condition, values_from = c(mean, sd)) 
                            })
    
    df_norm_mean_sd <- reactive({ values_df %>% # the normalized data
                                  unnest(data) %>%
                                  select(c(condition, Temperature, mean_norm, sd_norm)) %>%
                                  distinct(.keep_all = TRUE) %>%
                                  pivot_wider(names_from = condition, values_from = c(mean_norm, sd_norm))
                                 })
    
    # # Table of selected dataset ----
    output$table_set1 <- renderDataTable(datasetInput1(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))

    output$table_set2 <- renderDataTable(datasetInput2(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))


    # Downloadable csv of selected dataset ----
    output$downloadData1 <- downloadHandler(
        filename = function() {
            paste(input$dataset1_download_name, "-", input$dataset1, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(datasetInput1(), file, row.names = FALSE)
        }
    )

    output$downloadData2 <- downloadHandler(
        filename = function() {
            paste(input$dataset2_download_name, "-", input$dataset2, ".csv", sep = "")
        },
        content = function(file) {
            write.csv(datasetInput2(), file, row.names = FALSE)
        }
    )
    
    output$downloadData3 <- downloadHandler(
        filename = function() {
            paste(input$dataset3_download_name, "-", input$dataset3, ".rds", sep = "")
        },
        content = function(file) {
            write_rds(datasetInput3(), file)
        }
    )
 
}

shinyApp(ui, server)