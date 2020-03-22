# 2_add_layouts/app.R
library(shinyBS) # drop-down panels
library(tidyverse) #  handling data structures and plotting

source("support_scripts/upload_formatters.R")
source("support_scripts/layout_handling.R")

library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(rhandsontable) # user-interactive tables 
library(shiny) # for shiny web-apps 

df_sample <- read.csv("sample_data_file.csv")

ui <- navbarPage( useShinyalert(),
                 # Data Analysis --------------------------------------------------------------------------------
                 tabPanel(p("to data analysis", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "data_analysis_mother", # end tab panel (tabset, div, main still remaining)
                          tabsetPanel(id = "inTabset_analysis", # tabset for all analysis sub-tabs
                                      # analyze and visualize --------------------------- 
                                      tabPanel(p("2 | analyze and visualize", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "analysis_tab",
                                               tags$head( # det the slider aesthetic
                                                   tags$style(
                                                       ".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: grey; border-color: transparent;}",
                                                       ".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: grey; border-color: transparent;}"
                                                   )
                                               ),
                                               sidebarLayout(
                                                   sidebarPanel(
                                                       bsCollapse(id = "plot_aes", open = "Panel 2",
                                                                  bsCollapsePanel(p("Set plate layout and replicates", style = "font-family: 'Avenir Next'; font-size: 16px; color: black",align = "center"), 
                                                                                  bsCollapsePanel(p("Method 1 - upload layout", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), 
                                                                                                  p("Our favorite approach to DSF data analysis.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                  p("Use the template below to create a layout file for your experiment. Each plate in the layout file defines a new experimental variable (e.g. compound, pH, concentration), with the varible name provided in the first column of the layout file. You can define any number of variables by adding additional plates to the layout file. Using this method, data can be visualized by user-defined variables (e.g. color by concentration).", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                  p("Layouts are connected to data by well name, so your data must have a 'well' column to use this feature.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                  p("For more information, see the instructions tab.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                  downloadButton("sample_layout_file", p("Download layout template", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center")), #
                                                                                                  p("...", style = "font-family: 'Avenir Next'; font-size: 12px; color: white",align = "center"),
                                                                                                  fileInput("layout_file", p("Upload your csv layout file", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                            accept = c(
                                                                                                                "text/csv",
                                                                                                                "text/comma-separated-values,text/plain",
                                                                                                                ".csv")
                                                                                                  )),
                                                                                  p("...", style = "font-family: 'Avenir Next'; font-size: 15px; color: white",align = "center"),
                                                                                  p("Method 2 - edit manually", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"),
                                                                                  rHandsontableOutput("r_table"),
                                                                                  p("...", style = "font-family: 'Avenir Next'; font-size: 15px; color: white",align = "center"),
                                                                                  uiOutput("handson_update_button"),
                                                                                  #actionButton("submit_handson_names", p("Update names", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center")),
                                                                                  style = "default")
                                                       )),
                                                   mainPanel( 
                                                              p("The layout-combined data is shown in the data table below. The data itself is nested as a single entry per row of this dataframe; that row is not rendered in the table."),
                                                              DT::dataTableOutput('df_layout_table'),
                                                              verbatimTextOutput("df_struc_print")
                                                              # df_struc_print <- renderPrint({ str(values$df) })
                                                       )
                                                   
                                               )
                                      )
                                      # end tabpanel
                          )) # end tabset Panel (contains all "analysis sub-panels)
                 
) # end 


# Define server logic required to draw a histogram
server <- function(session, input, output) {
    values <- reactiveValues() 
    ## just in the layout standalone
    values$data_raw <- df_sample
    output$df_struc_print <- renderPrint({ str(values$df) }) # this is used only the applet, to print out the structures
 
    
####### begin 2_layouts applet server   ##### 
####### basic data formatting post-upload
    observeEvent(values$data_raw, { # ultimately, observe the transfer to the analysis page
        print("observe2")
        req(values$data_raw) # but leave this requirement as is
        
        tryCatch({
            values$df <- nest_raw(values$data_raw) %>% # active dataframe, used for plotting and calculations
                         add_standardized_wells()
            
            values$df_1 <- nest_raw(values$data_raw) %>% # original dataframe, used when a "clean slate" is needed, e.g. if a new layout file is uploaded
                add_standardized_wells()
            
        }, error = function(e) {
            shinyalert("Please ensure that your data is formatted correctly", "In the 'upload data' tab, you data should be displayed with Temperature in the first column, and RFU data in the columns to the right.")
        }
        )
    })
    
########## End render GUI elements for the analysis page ######
    
####### begin layout handling and updating
    values$r_table_update <- 1 # initialize the r_table update counter
    output$df_layout_table <- DT::renderDataTable({ values$df %>% select(-data) })

    # GUI elements
    output$handson_update_button <- renderUI({ # for layout
        req(values$df)
        actionButton("submit_handson_names", p("Update names from manual table (above)", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),  width = '100%')
    })
    
    # data handling
    layout <- reactive({
        req(input$layout_file)
        
       tryCatch({
            make_layout(input$layout_file$datapath) %>% # all columns are characters
                add_standardized_wells()
            
        }, error = function(e) {
            shinyalert("Please ensure that your layout is formatted correctly", "Layout file should be uploaded asa  UTF-8 csv. A layout template can be downloaded below.")
        })
    })

    observeEvent( layout(), { 
      req(layout()) # don't do this if the layout doesn't exist
        tryCatch({
            values$df <- join_layout_nest(values$df_1, layout() ) 
            
        }, error = function(e) {
          shinyalert("Please ensure that your layout is formatted correctly", "Layout file should be uploaded asa  UTF-8 csv. A layout template can be downloaded below.")
        })
        
    })
    
    observeEvent(input$submit_handson_names, { # when r_table is updated
        req(input$r_table)
      
        values$df<- hot_to_r(input$r_table) %>% # update the layout
                            as_tibble() %>%
                            ensure_standardized_wells() %>%
                            join_layout_nest( values$df_1, . )
        
        write_rds(values$df, "values_df_with_layout.rds")
        
       values$r_table_update <- values$r_table_update + 1 # trigger the re-rendering of the r_table
        
    }, ignoreInit = TRUE, ignoreNULL = TRUE)
 
  output$r_table <- renderRHandsontable({
              
            tryCatch({
                req(values$df)
                trigger <- values$r_table_update
                
                if (is.null(input$layout_file) == FALSE) { # if there is a layout file
                    layout_vars <- names(layout())[!c(names(layout()) %in% c("well_", "well_f_", "row_", "col_", "row", "column"))]
                    
                    handson_df <- values$df %>%
                        select( one_of( layout_vars )) # this will always include "condition"
                    
                } else { # if no layout file
                    handson_df <- values$df %>%
                        select(well, condition)
                }
                rhandsontable( handson_df, height = 200, useTypes = TRUE, stretch = "all") %>% hot_col(c("well"), readOnly = TRUE)
                
            }, error = function(e) {
              shinyalert("Please ensure that your layout is formatted correctly", "Layout file should be uploaded asa  UTF-8 csv. A layout template can be downloaded below.")
            })
            
           
        })
  ### end 2_layouts applet server
    
} # end server
# Run the application 
shinyApp(ui = ui, server = server)
