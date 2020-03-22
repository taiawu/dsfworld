# perhaps we could put the loading of packages into an observe event, to keep them all froam loading right at the start. would this make the first page faster?
# 1_upload/app.R
library(shinyBS) # drop-down panels
library(shinyalert) # pop-up windows
library(tidyverse) #  handling data structures and plotting
library(shiny) # for shiny web-apps 
source("support_scripts/upload_formatters.R") # scripts written to analyze and visualize DSF data

# generate vectors of possible well names, for matching to layouts
WELLS1 <- make_well_names("ROWS", "1")
wells_any <- c(WELLS1, # e.g. A1 .... P24
               make_well_names("ROWS", "01"), # e.g. A01 .... P24
               make_well_names("rows", "1"), # e.g. a1 .... p24
               make_well_names("rows", "01") # e.g. a01 .... p24
            )

ui <- navbarPage(
                useShinyalert(),
                                      tabPanel(
                                          p("1 | upload data", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "uploads_tab",
                                               sidebarLayout( # Sidebar layout with input and output definitions
                                                   sidebarPanel(# Sidebar panel for inputs
                                                       #uploading---------------------------
                                                       fileInput("uploaded_file", p("Browse or drag-and-drop raw (RFU) DSF data", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), # Input: Select a file
                                                                 multiple = FALSE,
                                                                 accept = c("text/csv",
                                                                            "text/comma-separated-values,text/plain",
                                                                            ".csv",
                                                                            ".xls",
                                                                            ".xlsx")),
                                              
                                                       bsTooltip("help1", HTML("To analyze data, please upload it as a .tsv, .csv, .xls, or .xlsx, formatted with Temperature in the first column and raw fluorescence measurements in the remaining columns. A correctly-formatted example file can be downloaded at left. Minor reformatting can be done after uploading using the Reformatting assistance options. DSFworld can accept and reformat data files exactly as they are exported from the instruments listed in under Supported Reformatting (at left). See the Instructions tab for more information. Incompatible with Explorer versions 9 and earlier."),
                                                                 "right", options = list(container = "body"), trigger = "hover"),
                                                       bsCollapse(id = "upload_help", open = "Panel 1",
                                                                  bsCollapsePanel(p("Uploading instructions", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                                  p("Upload raw RFU data by either of the following two methods:", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"),
                                                                                  p("(i) exactly as exported from the instruments listed under 'reformat raw from instrument', below", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"),
                                                                                  p("(ii) as formatted in the example file, (download it below): a UTF-8 csv file, with Temperature in the first column and RFU data in the columns to the right. ", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"),
                                                                                  p(" ", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                  p("After uploading, your data will appear in its current format in a table at right. Minor adjustments can then be made, such as replacing cycle number values with Temperature.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                  p("To use the plate-layout capabilities (e.g. setting replicates and making custom plots) in the analysis window, each data column must be named by well. Most instruments automatically export data with wells as column names, but if necessary, you can artifically write well names onto your data under 'alter delminiters, headers', below", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                  downloadButton("download_sample_input", "Download example file", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff")
                                                                  )),
                                                       bsCollapse(id = "file_parse_types", open = "Panel 1",
                                                                  bsCollapsePanel(p("Alter delimiters, headers", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                                  checkboxInput("header", "Header", TRUE), # Input: Checkbox if file has header
                                                                                  checkboxInput("name_to_well", "Overwrite column names with wells", FALSE), # Input: Checkbox if file has header
                                                                                  radioButtons("sep", "Separator",  # Input: Select separator
                                                                                               choices = c(Comma = ",",
                                                                                                           Semicolon = ";",
                                                                                                           Tab = "\t"),
                                                                                               selected = ","),
                                                                                  radioButtons("quote", "Quote",  # Input: Select quotes
                                                                                               choices = c(None = "",
                                                                                                           "Double Quote" = '"',
                                                                                                           "Single Quote" = "'"),
                                                                                               selected = '"'))),
                                                       bsCollapse(id = "instrument_reformatting", open = "Panel 1",
                                                                  bsCollapsePanel(p("Reformat raw from instrument", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                                  p("Select your instrument from the list below, and upload data exactly as exported.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "left"),
                                                                                  radioButtons("reformat", "", # Input: Select type of reformatting necessary
                                                                                               choices = c(None = "none",
                                                                                                           Biorad = "biorad",
                                                                                                           Stratagene = "stratagene",
                                                                                                           quantStudio = "quantStudio",
                                                                                                           qTower = "qTower" # will have to figure out how to deal with multiple reader errors
                                                                                               ),
                                                                                               selected = "none"))),
                                                       bsCollapse(id = "cycle_to_T_panel", open = "Panel 1",
                                                                  bsCollapsePanel(p("Convert cycle number to temperature", style = "font-family: 'Avenir Next'; font-size: 14px; color: black", align = "center"),
                                                                                  checkboxInput("cycle_to_T", "Convert cycle number to temperature? (if yes, specify below)", FALSE),
                                                                                  textInput(inputId="start_T", label="Starting Temp (C)", value = 25),
                                                                                  textInput(inputId="increment_T", label="Increase per cycle (C)", value = 1)
                                                                                  )
                                                                  ),
                                                       actionButton('jumpToAnalysis', p("Analyze", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                    icon("chart-area"), width = '100%', style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff")


                                                   ),  # end sidebar panel

                                                   # Main panel for displaying outputs
                                                   mainPanel(
                                                       tags$style(type='text/css', "#instructions {font-size: 18px; line-height: +2;} "),
                                                       #HTML("instructions"),
                                                       dataTableOutput("input_file"), style = "overflow-x: scroll;"
                                                   ) # end main panel
                                               ) # end sidebarLayout

                                      )
                 )

server <- function(input, output) {
    values <- reactiveValues() # initalize the reactive values container (is this effectively a class? is shiny OOP-like...?) wish i'd realized that earlier)
    
    ########### data uploading  ----------
    # download a sample file
    output$download_sample_input <- downloadHandler(
        filename = function() {
            paste('dsfworld_upload_format.csv', sep='')
        },
        content = function(file) {
            read_csv("sample_file.csv")
            write.csv(read_csv("sample_file.csv"), file, row.names = FALSE)
        }
    )
    
    data_raw <- reactive({
        print("uploading raw file")
        req(input$uploaded_file)
        file <- input$uploaded_file$datapath
        
        tryCatch({
            
            if (input$reformat == "none" ) {          df <- read.table(input$uploaded_file$datapath, # file
                                                                                 header = input$header, # colnames
                                                                                 sep = input$sep,
                                                                                 quote = input$quote,
                                                                                 stringsAsFactors =  FALSE)
            
            } else if (input$reformat == "biorad") { df <- read.table(input$uploaded_file$datapath, # file
                                                                      header = input$header, # colnames
                                                                      sep = input$sep,
                                                                      quote = input$quote,
                                                                      stringsAsFactors =  FALSE) %>% 
                                                                      format_biorad()
            
            } else if (input$reformat == "stratagene") { df <- read.table(input$uploaded_file$datapath, # file
                                                                          header = input$header, # colnames
                                                                          sep = input$sep,
                                                                          quote = input$quote,
                                                                          stringsAsFactors =  FALSE) %>% 
                                                                          format_stratagene()
            
            } else if (input$reformat == "quantStudio") { df <- read_quantStudio(input$uploaded_file$datapath)
            } else if (input$reformat == "qTower") { df <- read_qTower(input$uploaded_file$datapath) # take path bc read.table will error on this file
            }
            
            df <- df %>% # in case someone has a file that reads in as characters
                mutate_if(is.factor, as.character) %>% # make any factors characters
                mutate_all(as.numeric) # make all numeric
            
            if (input$name_to_well == TRUE) {
                df <- df %>%
                    set_names(c("Temperature", WELLS1[c(1:(ncol(.)-1))]))
            }
            # cycle number to temperature
            if (input$cycle_to_T == TRUE) {
                Temps_calc <- cycle_to_T_func(as.numeric(input$start_T), as.numeric(input$increment_T), df)
                df_cycle <- dplyr::bind_cols(Temps_calc, df)[-2] 
                names(df_cycle)[1] <- "Temperature"
                return(df_cycle)
                    
            } else {
                names(df)[1] <- "Temperature"
                df
            }
        },   
        error = function(e){
            shinyalert("This file needs pre-formatting", "Please select your instrument from 'Reformat raw from instrument'. If you don't see your instrument there, please format your data as shown in the downloadable template and upload again.")
            values$data_raw <<- NULL
        }
        )
    }) # read the input file
    
    observeEvent(data_raw(), {
        values$data_raw <- data_raw()

        # set the following values based on the data
        tryCatch({
            low_T <- isolate( data_raw()$Temperature %>% min() )
            high_T <- isolate( data_raw()$Temperature %>% max() )
            n_meas <- isolate( data_raw() %>% nrow() )
            
            n2r <<- make_temp_n2r(range(low_T:high_T)) #make_temp_n2r(range(values$data$Temperature)) # an example of how this could be used
            win3d <<- floor(3/((n2r(1) - n2r(0))/n_meas))
                      if ( win3d < 5 ) { win3d <<- 5 }
            
            sgfilt_nest <<- sgfilt_set_n(n_ = find_sgolay_width( win3d ))
        },   
        error = function(e){
            print("win3 errored! setting win3d to 7")
            win3d <<- 7
            sgfilt_nest <<- sgfilt_set_n( n_ = find_sgolay_width( 7 ) )
        })
    }) # write to values
    
    observeEvent(data_raw(), { # ultimately, observe the transfer to the analysis page
        req(values$data_raw) # but leave this requirement as is
        tryCatch({
            values$df <- values$data_raw %>% # this is the active dataframe, used for plotting and calculations
                gather(well, value, -Temperature) %>%
                group_by(well) %>%
                mutate(value_norm = BBmisc::normalize(value, method = "range", range = c(0,1)), ###### if we do this as a mutate, it ignores the groups!!!!!!!
                       Temperature_norm = BBmisc::normalize(Temperature, method = "range", range = c(0,1)))  %>%
                nest %>%
                plyr::mutate(new_names = well)
            
            values$df_1 <- values$data_raw %>% # this is the original dataframe, used when a "clean slate" is needed, e.g. if a new layout file is uploaded
                gather(well, value, -Temperature) %>%
                group_by(well) %>%
                mutate(value_norm = BBmisc::normalize(value, method = "range", range = c(0,1)), ###### if we do this as a mutate, it ignores the groups!!!!!!!
                       Temperature_norm = BBmisc::normalize(Temperature, method = "range", range = c(0,1)))  %>%
                nest %>%
                plyr::mutate(new_names = well)
            }, error = function(e){
                shinyalert("Please ensure that your data is formatted correctly", "In the 'upload data' tab, you data should be displayed with Temperature in the first column, and RFU data in the columns to the right.")
            }
        )
    })
    
    
    ### the first steps of the analysis app, to make sure they stitch together ok
        output$input_file <- renderDataTable({
            req(input$uploaded_file)
            tryCatch(
                values$data_raw,
                error = function(e){
                    shinyalert("File needs pre-formatting!", "Please select your instrument from 'Supported Reformatting'. Or, if you don't see your instrument there, please format your data as shown in the downloadable template and upload again.")
                }
            )
        }, options = list(scrollX = TRUE, scrollY = 500, scrollCollapse = TRUE, paging = FALSE, dom = 't')
        )
        
        
        ##### end data uploads applet
}

shinyApp(ui, server)
