# 3_plotting/app.R
library(shinyBS) # drop-down panels
library(tidyverse) #  handling data structures and plotting

source("support_scripts/upload_formatters.R")
source("support_scripts/layout_handling.R")
source("support_scripts/plotting.R")

library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(rhandsontable) # user-interactive tables 
library(shiny) # for shiny web-apps 

ui <- navbarPage("",
                 #tabPanel(p("to data analysis", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "right")),
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
                                                       ),
                                                       bsCollapse(id = "plot_aes", open = "Panel 2",
                                                                  bsCollapsePanel(p("Make plots", style = "font-family: 'Avenir Next'; font-size: 16px; color: black", align = "center"),
                                                                                  # update_plot
                                                                                  uiOutput("trigger_df_1"),
                                                                                  uiOutput("update_plot"),
                                                                                  radioButtons("facet", "Make sub-plots",
                                                                                               c("Single plot" = "none",
                                                                                                 "Subset by one variable" = "wrap",
                                                                                                 "Subset by two variables" = "grid")),
                                                                                  uiOutput("wrap_by"), # reactive selector for the graph?
                                                                                  uiOutput("grid_rows"),
                                                                                  uiOutput("grid_cols"),
                                                                                  uiOutput("color_by"),
                                                                                  checkboxInput("use_linetypes", "Vary line types", FALSE),
                                                                                  uiOutput("linetype_by"),
                                                                                  uiOutput("mean_or_each"),
                                                                                  radioButtons("mean_or_each", "Show mean?",
                                                                                               c("Each replicate" = "each",
                                                                                                 "Mean" = "mean")),
                                                                                  
                                                                                  radioButtons("fix_free", "Equalize y axes?",
                                                                                               c("Equal" = "fixed",
                                                                                                 "Independent" = "free")),
                                                                                  bsCollapsePanel(h5("Edit plot labels"),
                                                                                                  textInput("plot_title", "Plot title", "Raw RFU Data"),
                                                                                                  textInput("legend_title", "Legend title", "Condition"),
                                                                                                  uiOutput("linetype_title"),
                                                                                                  textInput("y_title", "y-axis title", "RFU"),
                                                                                                  textInput("x_title", "x-axis title", "Temperature (ºC)"),
                                                                                                  numericInput("text_size", "Plot text size", 10, min = 4, max = 20),
                                                                                                  checkboxInput("hide_legend", "Hide legend", FALSE)),
                                                                                  style = "default")
                                                       )
                                                       
                                                   ),
                                                   
                                                   mainPanel(
                                                       wellPanel(
                                                           tags$script("$(document).on('shiny:connected', function(event) {
                                                                                var myWidth = $(window).width();
                                                                                Shiny.onInputChange('shiny_width',myWidth)
                                                                                
                                                                                });"),
                                                           
                                                           tags$script("$(document).on('shiny:connected', function(event) {
                                                                                var myHeight = $(window).height();
                                                                                Shiny.onInputChange('shiny_height',myHeight)
                                                                                
                                                                                });"),
                                                           
                                                           id = "facet_plot_re_examined",
                                                           
                                                           tags$style(HTML('#q1 {margin-top: 30px}')),
                                                           splitLayout(cellWidths = c("20%", "80%"), 
                                                                       plotDownloadUI("plot1"), 
                                                                       textInput("plot_download_name", "Downloaded plot name", value = "dsfworld_plot")
                                                           ),
                                                           plotOutput("data", height = "auto") %>% withSpinner(color="#525252"), style = ("overflow-y:scroll; max-height: 600px") 
                                                       ))
                                               )
                                      )
                                      # end tabpanel
                          )) # end tabset Panel (contains all "analysis sub-panels)
                 
) # end 

# Define server logic required to draw a histogram
server <- function(session, input, output) {
    values <- reactiveValues() 
    values$data_raw <- readRDS("values_df_with_layout.rds")
    values$df <- readRDS("values_df_with_layout.rds")
    values$df_1 <- readRDS("values_df_with_layout.rds")
    values$plot_chosen <- "initial"

    ######## begin data layout and handling #########   
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
    ######## end data layout and handling #########   
    
    ########## Render GUI elements for the analysis page ######
    output$update_plot <- renderUI({
        req(values$df)
        actionButton("update_plot", p("Update plot", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center") %>% strong(),  width = '100%')
    })
    ########## End render GUI elements for the analysis page ######
    
    ######## eval selections, to pass to the plotting function #########
    output$wrap_by <- renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        if (input$facet == "none" | input$facet == "grid") return(NULL)
        varSelectInput("wrap_by", label = "Sub-plot by",  # as select input
                       data = values$df %>% get_layout_vars()
        ) # this is able to take a reactive value
    })
    
    output$grid_rows <- renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        if (input$facet == "none" | input$facet == "wrap") return(NULL)
        varSelectInput("grid_rows", label = "Sub-plot grid, rows",  # as select input
                       data = values$df %>% get_layout_vars()
                       ) # this is able to take a reactive value
    })
    
    output$grid_cols <- renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        req(input$grid_rows)
        if (input$facet == "none" | input$facet == "wrap") return(NULL)
        
        tryCatch({
            col_data <- values$df %>% 
                get_layout_vars() %>%
                select(-!!input$grid_rows)
        }, error = function(e) {
            col_data <- NULL
        })
        
        varSelectInput("grid_cols", label = "Sub-plot grid, columns",  # as select input
                        data = col_data,
                       selected = "- ")
    }) 
    
    output$color_by <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        varSelectInput("color_by", label = "Color",  # as select input
                       data = values$df %>% get_layout_vars(),
                       selected = "-"
                       ) # this is able to take a reactive value
    })  
    
    output$linetype_by <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        if (input$use_linetypes == FALSE) return(NULL)
        varSelectInput("linetype_by", label = "Line types",  # as select input
                       data = values$df %>% get_layout_vars(),
                       selected = "-"
                       ) # this is able to take a reactive value
    })
    
    output$linetype_title <-renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        if (input$use_linetypes == FALSE) return(NULL)
        textInput("linetype_title", "Line type legend title", "Condition2")
    })
    
    plot_title_d <- reactive({input$plot_title})
    plot_legend_d <- reactive({input$legend_title})
    plot_legend_linetype <- reactive({ input$linetype_title })
    
    legend_position <- reactive({ 
        if (input$hide_legend == TRUE ) { legend_pos <- "none"
        } else { legend_pos <- "right"}
        legend_pos})
    
    ######## end eval selections, to pass to the plotting function #########  
    
    # update and re-render the plot only when the update-plot button is clicked!
    #init_plot <- eventReactive( input$trigger_df_1, { # when new data is uploaded
    plot_initial <- eventReactive( values$data_raw,  { # only when the "update plot" button is clicked, update the plot 
        print("plot changed")
        req(values$df) # only render the plot if there is data
        
        unnest(values$df) %>%
                ggplot(aes(x = Temperature, y = value, group = well)) +
                    geom_line(size = 0.5, alpha = 0.7) +
                    labs(title = "Raw RFU Data", x = "Temperature (ºC)", y = "RFU") +
                    theme_bw() +
                    dsfworld_default +
                    theme(  text = element_text(size = 10*1.25),
                            axis.text = element_text(size = 10),
                            plot.title = element_text(lineheight=.8, face="bold", size = 10*1.5)) 
    })
    
    plot_updated <- eventReactive( input$update_plot,  { # only when the "update plot" button is clicked, update the plot 
        print("plot changed")
            req(values$df) # only render the plot if there is data

            df_RFU_plot <- unnest(values$df) %>%
                plyr::mutate("-" = rep("", nrow(.))) %>%
                plyr::mutate("- " = rep("", nrow(.)))
 
            facet_func(df = df_RFU_plot,# reacts to the appearance and changes to the dataframe, to the uploading of format files
                       mean_or_each = input$mean_or_each,
                       color_by = !!input$color_by,
                       linetype_by = !!input$linetype_by,
                       use_linetypes = input$use_linetypes,
                       facet = input$facet,
                       facet_by = !!input$wrap_by,
                       facet_rows = !!input$grid_rows,
                       facet_cols = !!input$grid_cols,
                       set_title = plot_title_d(),
                       legend_title = plot_legend_d(),
                       legend_linetype_title = plot_legend_linetype(),
                       fix_free = input$fix_free,
                       text_size = input$text_size,
                       legend_position = legend_position(),
                       x_title = input$x_title,
                       y_title = input$y_title)
        })
    
    observeEvent(values$data_raw, { values$plot_chosen <- "initial"  })    
    observeEvent(input$update_plot, { values$plot_chosen <- "updated"  })
    observeEvent(values$show_model_plot, { values$show_model_plot <- "model"  })  
  
chosen_plot <- reactive({
        if (values$plot_chosen == "updated") { # TRUE when the "update plot" button was clicked more recently than new data uploads or "show model plot"
            plot_updated()
        } else if (values$plot_chosen == "model") { # TRUE when the "show model plot" button was clicked more recently than new data uploads or "update model plot"
            plot2()
        } else { # if its the initial plot
            plot_initial()  # this is the default
        }
    })

plot_height <- eventReactive(input$update_plot, {
                if (input$facet == "none") {
                    height <- 400 
                } else {
                    # adapted from https://github.com/rstudio/shiny/issues/650
                    h_dyn <- gg_facet_nrow_ng(plot_updated()) * ((session$clientData$output_data_width-100)/(gg_facet_ncol_ng(plot_updated())))*(1/1.618)
                    
                    if ( h_dyn < session$clientData$output_data_width * (1/1.618) ) { 
                        height <- session$clientData$output_data_width * (1/1.618)
                    } else { height <- h_dyn
                    }
                } 
     height
    })

output$data <- renderPlot({ # is there a way to implement renderCachedPlot that would be worthwhile here?
         chosen_plot()

     }, height = function() 
         if (values$plot_chosen == "initial") { 400 
         } else {
             plot_height() 
             }  
     ) 

} # end server
# Run the application 
shinyApp(ui = ui, server = server)
