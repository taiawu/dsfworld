# 0_integrated_modules/app.R
library(quantmod) # contains the findValleys function, which maybe we should just extract and put verbatim in a source file instead of loading this whole thing...?
library(minpack.lm) # contains the nlsLM function, which we use for our fitting
library(modelr) # used in both the data modeling and the analysis model fitting 
library(SciViews) # contains the ln function used in the data modeling
library(signal) # contains the savistky golay filter (savgolfilt), used to generate the first derivative data in both data modeling and analysis model fitting  

library(shinyBS) # drop-down panels
library(tidyverse) #  handling data structures and plotting

source("scripts/data_modeling.R") # computational models and all associate plots

source("scripts/upload_formatting.R") # functions to assist with raw data uploading
source("scripts/layout_handling.R") # functions to apply layouts to data
source("scripts/plotting.R") # functions to make all plots displayable in the analysis window
source("scripts/analysis.R") # perform Tma analyses by dRFU or curve fitting

library(shinyalert) # pop-up error messages
library(shinycssloaders) # spinning plot loading icon
library(rhandsontable) # user-interactive tables 
library(shiny) # for shiny web-apps


named_mods <- c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred") %>% # determine the model names
    set_names("Fit 1", "Fit 2", "Fit 3", "Fit 4")

# generate vectors of possible well names, for matching to layouts
WELLS1 <- make_well_names("ROWS", "1") # create well names, used in the uploading page 
wells_any <- c(WELLS1, # e.g. A1 .... P24
               make_well_names("ROWS", "01"), # e.g. A01 .... P24
               make_well_names("rows", "1"), # e.g. a1 .... p24
               make_well_names("rows", "01") # e.g. a01 .... p24
)

# Define server logic required to draw a histogram
server <- function(session, input, output) {
    ########### inter-website navigations 
    observeEvent(input$jumpToAnalysis, {
        updateTabsetPanel(session, "inTabset_analysis",
                          selected = "analysis_tab")
    })
    
    ########### interactive modeling ----------
    df_model <-  reactive({
        #make_model_df(start_T = (273+25), end_T = (273+95), dHu_ = 270, T_half_ = (55+273), dCp_ = 8, Ea_ = 100, T_star_ = 85+273, v_ = 1, nat_dye = 0, unf_dye = 1, fin_dye = 1, decay_rate = 0.8) -> df
        make_model_df(start_T = (273+25), end_T = (273+95), 
                      dHu_ = input$dHu_, 
                      T_half_ = input$T_half_ + 273, #55+273, 
                      dCp_ = input$dCp_, 
                      Ea_ = input$Ea_, 
                      T_star_ = input$T_star_ + 273, 
                      v_ = input$v_, 
                      nat_dye =input$nat_dye_, 
                      unf_dye = input$unf_dye_, 
                      fin_dye = input$fin_dye_, 
                      decay_rate = input$decay_rate_
        )
    })
    
    p_model <- reactive(make_model_plot(df_model()))
    output$plot_model <- renderPlot(p_model())
    ####### interactive modeling 
    
    
    ########## data analysis #########
    values <- reactiveValues() # initalize the reactive values container
    
    ###### data uploading #####
                # inputs: 
                        # input$uploaded_file # the file uploaded by the user
                # outputs: 
                        ## used beginning of the data handling structures
                        # data_raw()
                        # values$data_raw # raw data, retained without labels, layouts, or any modifications past reformatting
                        # values$df  # the active, nested dataframe to which layouts are added. Passed to layouts, plots, and analysis
                        # values$df_1 # the nested dataframe created directly from data_raw() and used to reset values$df_1
                        
                        ## used for model fitting
                        # low_T 
                        # high_T 
                        # n_meas 
                        # n2r 
                        # win3d 
                        # sgfilt_nest 
    
    ####### data upoading functions ####
    output$download_sample_input <- downloadHandler(
        filename = function() {
            paste('dsfworld_sample_raw_data_format.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_sample_raw_data_format.csv")
            write.csv(read_csv("dsfworld_sample_raw_data_format.csv"), file, row.names = FALSE)
        }
    )
    
    output$sample_layout_file <- downloadHandler(
        filename = function() {
            paste('dsfworld_example_layout.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_example_layout.csv")
            write.csv(read_csv("dsfworld_example_layout.csv"), file, row.names = FALSE, fileEncoding = "UTF-8")
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
    
    move_to_analysis <- FALSE # set this as false to start
    # observeEvent({ id$inTabset_analysis
    #                input$jumpToAnalysis }, { # ultimately, observe the transfer to the analysis page
    # open_analysis_tab <- FALSE
    # if ( input$inTabset_analysis == "analysis_tab" ) {
    #     open_analysis_tab <- TRUE   
    # }
    
   # observeEvent(input$inTabset_analysis , {
      panel <- reactive({ input$inTabset_analysis  }) 
      observeEvent(input$inTabset_analysis, {
          print("panel changed")
          print(panel())
      })
     # print("panel change")
     # print(panel())
        # print(input$inTabset_analysis)
        # panel_change <- panel_change + 1
   # })
     # list(input$test1,input$test2)
   #data_raw_counter <- reactive({ 0})  #reactive({0
   counter <- reactiveValues(data_upload = 0,
                             data_process = 0)
   
   toListen <- reactive({
       list(input$inTabset_analysis, counter$data_upload)
   })
   
    observeEvent( toListen()        , {
         print("tirggered win3d calc")
         print(input$inTabset_analysis)
         req(input$inTabset_analysis == "analysis_tab") # hwo to require that something be true, rather than existing? # only proceed if the user is on the analysis tab
         req(data_raw())
         print("passed the req tabset")
        #values$data_raw <- data_raw()
        
        # set the following values based on the data
        tryCatch({
            low_T <- isolate( data_raw()$Temperature %>% min() )
            high_T <- isolate( data_raw()$Temperature %>% max() )
            n_meas <- isolate( data_raw() %>% nrow() )
            
            n2r <<- make_temp_n2r(range(low_T:high_T)) #make_temp_n2r(range(values$data$Temperature)) # an example of how this could be used
            win3d <<- floor(3/((n2r(1) - n2r(0))/n_meas))
            if ( win3d < 5 ) { win3d <<- 5 }
            sgfilt_nest <<- sgfilt_set_n(n_ = find_sgolay_width( win3d ))
            counter$data_process <- counter$data_process + 1  #data_raw_counter() <- data_raw_counter() + 1 ##### THIS SEEMS TO FAIL FIX THIS FIRST
            move_to_analysis <<- TRUE
        },   
        error = function(e) {
            print("win3 errored! setting win3d to 7")
            win3d <<- 7
            sgfilt_nest <<- sgfilt_set_n( n_ = find_sgolay_width( 7 ) )
            move_to_analysis <<- FALSE
        }, warning = function(w) {
            print("win3 warning! setting win3d to 7")
            win3d <<- 7
            sgfilt_nest <<- sgfilt_set_n( n_ = find_sgolay_width( 7 ) ) 
            move_to_analysis <<- FALSE
        })
    }) # write to values

    
    #observeEvent({ data_raw() }, { # ultimately, observe the transfer to the analysis page
    observeEvent(counter$data_process, {
        req(data_raw()) # but leave this requirement as is
        if (move_to_analysis == FALSE) {
            shinyalert("Please ensure that your data is formatted correctly", "In the 'upload data' tab, you data should be displayed with Temperature in the first column, and RFU data in the columns to the right.")  
        }
       # print(input$inTabset_analysis)
        #req(input$inTabset_analysis == "analysis_tab")
        req(move_to_analysis == TRUE)
        tryCatch({
            print("assigning values$df")
                    values$df <- nest_raw(data_raw()) %>%  add_standardized_wells() # REVISIT1 # active dataframe, used for plotting and calculations
                    values$df_1 <- nest_raw(data_raw()) %>%  add_standardized_wells()  # REVISIT1 # needs an option for when there are non-well names in the data
                    min(unnest(values$df)$Temperature) # this fails if there is something wrong with data_raw()
                    }, error = function(e) {
            shinyalert("Please ensure that your data is formatted correctly", "In the 'upload data' tab, you data should be displayed with Temperature in the first column, and RFU data in the columns to the right.")
        })
    })
    
    ### the first steps of the analysis app, to make sure they stitch together ok
    output$input_file <- renderDataTable({
        req(input$uploaded_file)
        tryCatch(
            data_raw(),
            error = function(e) {
                shinyalert("File needs pre-formatting!", "Please select your instrument from 'Supported Reformatting'. Or, if you don't see your instrument there, please format your data as shown in the downloadable template and upload again.")
            }
        )
    }, options = list(scrollX = TRUE, scrollY = 500, scrollCollapse = TRUE, paging = FALSE, dom = 't')
    )
    
    ##### end data uploads applet
    
    
    ###### layouts #####
                # inputs: 
                    # input$layout_file # the user-uploaded layout file
                    # values$df_1 # the nested raw data to which the layout information is applied

                # outputs: 
                    # values$df, modified # adds the layout information to values$df_1 to generate values$df
    
    
    ####### layout functions #####
    values$r_table_update <- 1 # initialize the r_table update counter

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
       # write_rds(values$df, "values_df_with_layout.rds")

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
    #### end 2_layouts applet server
    
    
    
    ###### plotting #####
    # inputs: 
    # values$df #
    
    # inputs from the GUI for evals for the plot
    # input$facet                   
    # input$mean_or_each,
    # input$color_by,
    # input$linetype_by,
    # input$use_linetypes,
    # input$facet,
    # input$wrap_by,
    # input$grid_rows,
    # input$grid_cols,
    # plot_title_d(),
    # plot_legend_d(),
    # plot_legend_linetype(),
    # input$fix_free,
    # nput$text_size,
    # legend_position(),
    # input$x_title,
    # input$y_title
    
    # outputs: 
    # plot_initial # a simple, quick-rendering plot
    # plot_updated # the user-made plot from the eval statements
    # fit_plots # the BIC plots of all visualized fits
    # fit_plots_best # the best fits only 
    
    # values$plot_chosen # a string indicating which plot the user has chosen to display 
    # chosen_plot # the plot which will actually be displayed in the window
    
    # output$final_plot
    # output$download_plot # a GUI element
   
    ############### plotting functions #######
    output$update_plot <- renderUI({
        req(values$df)
        actionButton("update_plot", p("Update plot", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center") %>% strong(),  width = '100%')
    })
    
    # eval sections
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
    
    ##end eval selections
    
    # update and re-render the plot only when the update-plot button is clicked!
    #init_plot <- eventReactive( input$trigger_df_1, { # when new data is uploaded
    plot_initial <- eventReactive( values$df,  { # only when the "update plot" button is clicked, update the plot 
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
    
    fit_plots <- eventReactive( input$show_BIC_plot, {
        print("plotting all fits")
        plot_all_fits_shiny(values$df_models_p, values$df_BIC_models_p ) 
    })
    
    fit_plots_best <- eventReactive( input$show_best_fits, {
        print("plotting best fits")
        plot_best_fits_shiny(values$df_models_p, values$df_BIC_best)
    })
    
    values$plot_chosen <- "initial" # the starting value
    observeEvent(values$data_raw, { values$plot_chosen <- "initial"  })    
    observeEvent(input$update_plot, { values$plot_chosen <- "updated"  })
    observeEvent(input$show_BIC_plot, { values$plot_chosen <- "all_model"  })  
    observeEvent(input$show_best_fits, { values$plot_chosen <- "best_model"  })
    
    chosen_plot <- reactive({
        values$plot_chosen
        if (values$plot_chosen == "updated") { # TRUE when the "update plot" button was clicked more recently than new data uploads or "show model plot"
            plot_updated()
        } else if (values$plot_chosen == "all_model") { # TRUE when the "show model plot" button was clicked more recently than new data uploads or "update model plot"
            fit_plots()
        } else if (values$plot_chosen == "best_model") {
            fit_plots_best()
        } else { # if its the initial plot
            plot_initial()  # this is the default
        }
    })
    
    plot_height <- reactive({
        values$plot_chosen
        if (values$plot_chosen == "updated") { # TRUE when the "update plot" button was clicked more recently than new data uploads or "show model plot"
                            #plot_updated()
                            if (input$facet == "none") {
                                height <- 400  # the initial plot has a 400 px height
                            } else {
                                # adapted from https://github.com/rstudio/shiny/issues/650
                                h_dyn <- gg_facet_nrow_ng(plot_updated()) * ((session$clientData$output_plot_width-100)/(gg_facet_ncol_ng(plot_updated())))*(1/1.618)
                                if ( h_dyn < session$clientData$output_plot_width * (1/1.618) ) { # if the calulcated height fits within the screen
                                    height <- session$clientData$output_plot_width * (1/1.618) # if the calulcated height fits within the screen
                                } else { height <- h_dyn }
                            }
                            height # return this
            
        } else if (values$plot_chosen == "all_model") { # TRUE when the "show model plot" button was clicked more recently than new data uploads or "update model plot"
                    # fit_plots()
                    h_dyn <- gg_facet_nrow_ng(fit_plots()) * ((session$clientData$output_plot_width-100)/(gg_facet_ncol_ng(fit_plots())))*(1/1.618)
                    if ( h_dyn < session$clientData$output_plot_width * (1/1.618) ) { # if the calulcated height fits within the screen
                        height <- session$clientData$output_plot_width * (1/1.618) # if the calulcated height fits within the screen
                    } else { height <- h_dyn }
                    height
        } else if (values$plot_chosen == "best_model") {
            # fit_plots_best()
            h_dyn <- gg_facet_nrow_ng(fit_plots_best()) * ((session$clientData$output_plot_width-100)/(gg_facet_ncol_ng(fit_plots_best())))*(1/1.618)
            if ( h_dyn < session$clientData$output_plot_width * (1/1.618) ) { # if the calulcated height fits within the screen
                height <- session$clientData$output_plot_width * (1/1.618) # if the calulcated height fits within the screen
            } else { height <- h_dyn }
            height
        } else { # if its the initial plot
                    400 # return this #plot_initial()  # this is the default
        }
    })
    
    output$plot <- renderPlot({ # is there a way to implement renderCachedPlot that would be worthwhile here?
        chosen_plot()
    }, height = function() 
        plot_height()
    )
    
    output$final_plot <- renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        sliderInput("trim_ends", 
                    p("Restrict fits to a temperature range", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"), 
                    min = min(unnest(values$df)$Temperature), max = max(unnest(values$df)$Temperature),
                    value = c(min(unnest(values$df)$Temperature),
                              max(unnest(values$df)$Temperature)),
                    step = 1)
    }) # trim the ends off of the data to improve fitting
    
    ## download the plot
    output$download_plot <- downloadHandler(
        filename = function() { paste(Sys.Date(), '-dsfworld_plot.pdf', sep='') },
        content = function(file) {
            ggsave(file, plot = chosen_plot(), device = "pdf")
        }
    )

    
    ###### analysis #####
    # inputs: 
    # values$df # converted to values$df_fit, which is 
    
    
    # outputs: 
    ################# analysis functions #####
    
    observeEvent(values$df, {values$df_fit <- values$df} ) # if values$df changes, re-do the fitting. When layouts are updated this will trigger unnecessary re-fitting, but the fitting is fast enough that the extra computations are worth the added simplicty of doing it this way
    
    ## by dRFU
    observeEvent( values$df_fit, {
        # starting from values$df gives this calculation a fresh slate should the user re-format their data multiple times
        values$df_tms <- values$df_fit %>% #df_int %>% # add the first derivative Tms
            plyr::mutate(sgd1 = purrr::map(data, sgfilt_nest, m_ = 1)) %>% # add the first derivative data
            plyr::mutate(dRFU_tma = as_vector(purrr::map2(data, sgd1, Tm_by_dRFU)))
        
        if ("condition" %in% names(values$df_tms)) {
            # make the averaged table
            values$tm_table_dRFU <- values$df_tms %>%
                select(condition, dRFU_tma)  %>%
                group_by(condition)  %>%
                summarise( mean_tm = mean(dRFU_tma) ,
                           sd_tm = sd(dRFU_tma)) %>%
                mutate_if(is.numeric, round, 2)
        } else {
            values$tm_table_dRFU <- values$df_tms %>%
                select(well, dRFU_tma)  %>%
                group_by(well)  %>%
                summarise( mean_tm = mean(dRFU_tma) ,
                           sd_tm = sd(dRFU_tma)) %>%
                mutate_if(is.numeric, round, 2)
        }
        
    })
    
    tma_by_dRFU <- reactive({
        req(values$df_tms)
        if( "dRFU_tma" %in% names(values$df_tms)) { # if there are tms to report, render a table
            df <- values$tm_table_dRFU  %>%
                set_names( c("Condition", "Tma", "SD"))
            
        } else { df <- NULL }
        
        df   
    })
    
    output$tm_table_render <- DT::renderDataTable({
        req(values$df_tms)
        tma_by_dRFU()
    },
    options = list(scrollX = TRUE, scrollY = 200, scrollCollapse = TRUE, paging = FALSE, dom = 'tr')) #scroller = TRUE, dom = 'tr'
    
    # download these values
    output$download_dRFU_tma <- downloadHandler(
        filename = function() {
            paste0(Sys.Date(), '-dsfworld_tma_by_dRFU.csv', sep='')
        },
        content = function(file) {
            write.csv(tma_by_dRFU(), file, row.names = FALSE)
        }
    )
    
    # ## by model fitting
    # set initial values for smoothing and normalizing
    output$trim_ends <- renderUI({ # this is reactive by nature of being a render call? it can accept, therefore, rt(), which is a reactive expression. Can we
        req(values$df)
        sliderInput("trim_ends", 
                    p("Restrict fits to a temperature range", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"), 
                    min = min(unnest(values$df)$Temperature), max = max(unnest(values$df)$Temperature),
                    value = c(min(unnest(values$df)$Temperature),
                              max(unnest(values$df)$Temperature)),
                    step = 1)
    }) # trim the ends off of the data to improve fitting
    
    observeEvent(input$trim_ends, {print("input$trim_ends")
        print(input$trim_ends)
        print(str(input$trim_ends))
        
        if ( (input$trim_ends[2] - input$trim_ends[1]) < 25 ) {
            updateSliderInput(session, "trim_ends", 
                              
                              value = c(min(unnest(values$df)$Temperature),
                                        max(unnest(values$df)$Temperature)),
                              step =  1)
            shinyalert("Not enough data points remaining", "Please increase range to at least 25 measurements")
        } else {
            values$df_fit <- values$df %>%
                mutate(data = map(data, ~ filter(., Temperature %>% between(input$trim_ends[1], input$trim_ends[2])))) 
        }
    })
    
    observeEvent(values$df_fit, { # data_raw(), { ### CHANGE THIS BACK TO data_raw() for integration!!!!!!!!!
        tryCatch({ # REVISIT-- should this be called on data_raw() for any reason?
            low_T <-  values$df_fit %>% unnest() %>% .$Temperature %>% min() 
            high_T <- values$df_fit %>% unnest() %>% .$Temperature %>% max() 
            n_meas <- values$df_fit %>% unnest() %>% .$Temperature %>% unique() %>% length() 
            
            n2r <<- make_temp_n2r(range(low_T:high_T)) #
            win3d <<- floor(3/((n2r(1) - n2r(0))/n_meas))
            if ( win3d < 5 ) { win3d <<- 5 }
            peak_finder_nest <<- make_peak_finder_nest( win3d )
            sgfilt_nest <<- sgfilt_set_n(n_ = find_sgolay_width( win3d ))
            
            outlist1 <- list(
                "n2r" = n2r,
                "win3d" = win3d,
                "peak_finder_nest" = peak_finder_nest,
                "sgfilt_nest" = sgfilt_nest
            )
            #write_rds(outlist1, "outlist1.rds")
        },   
        error = function(e) {
            print("win3 errored! setting win3d to 7")
            n2r <<- make_temp_n2r(range(25:95)) ### REVISIT
            win3d <<- 7
            sgfilt_nest <<- sgfilt_set_n( n_ = find_sgolay_width( 7 ) )
            peak_finder_nest <<- make_peak_finder_nest( win3d )
        })
    }) # write to values
    
    # fit the requested models
    observeEvent(values$df_fit, {
        values$start_pars <- get_start_pars(values$df_fit)
        # write_rds(values$df_fit, "values$df_fit.rds")
        # write_rds(values$start_pars, "values_start_pars.rds")
        
        # first, fit the s1 model
        # this will over-write the summary dataframes, re-setting the models to follow
        values$s1_list <- model_all(s1_model, "s1_pred", values$start_pars, win3d)
        
        values$model_list <- list("s1_list" = values$s1_list)
        values$df_models <- values$s1_list$df_models
        values$df_BIC_models <- values$s1_list$df_BIC
        values$df_tm_models <- values$s1_list$tm_table_models
        values$df_BIC_models_p <- values$df_BIC_models %>%
                                        cond_df_BIC_for_plot (  ) 
       #  write_rds(values$df_tm_models, "values_df_tm_models_s1.rds")
       # "tm_models_all" = tm_table_models$df_tma,

        values$df_tm_models_table <- values$df_tm_models %>%
            dplyr::filter( which_model == "s1_pred"  ) %>%
            plyr::mutate( which_model = grep_and_gsub(.$which_model, c("s1_pred", "s1_d_pred", "s2_pred","s2_d_pred"), c("Fit 1", "Fit 2", "Fit 3", "Fit 4"), c("Other")))  %>% # move this to later, for the for-display table only!
            set_names(c("Condition", "Model", "Tma 1", "Tma 1 SD", "Tma 2", "Tma 2 SD")) %>%
            discard(~all(is.na(.x))) #%>% # get rid of the columns which are all NA (true if the model is not selected)
        # write_rds(values$df_tm_models_table, "values_df_tm_models_table.rds")
        # write_rds(values$s1_list, "values_s1_list.rds")
        model_name_all <- c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred") # doesn't need to be in the server or the observer but is fast enough to justify, since it makes the next step clearer
        model_name_true <- reactive({model_name_all[c(input$s1, input$s1_d, input$s2, input$s2_d)]})
        
        values$df_models_filt <- values$df_models %>% dplyr::filter(which_model %in% model_name_true())
        values$df_BIC_models_filt <- values$df_BIC_models %>% dplyr::filter(which_model %in% model_name_true())
        
        values$df_models_p <- cond_df_model_for_plot( values$df_models_filt, values$df_BIC_models_filt  )
        
        values$df_BIC_models_p <- values$df_BIC_models_filt %>%
            cond_df_BIC_for_plot (  ) # adds is_min #readRDS("../4_analyze/values_df__BIC_models.rds")
        
        values$df_BIC_best <-  cond_df_BIC_for_plot ( values$df_BIC_models_filt   ) %>%
            filter(is_min == TRUE) %>%
            select(c(well, condition, which_model))
        
        # if new data is uploaded, reset all of the buttons as well. perhaps we should set these to watch values$data (unnamed), so it doesn't get over-written by renaming, but i'd need to think more carefully about how to incorporate the names downstream....
        updateButton(session, "s1",  value = TRUE)
        updateButton(session, "s1_d",  value = FALSE)
        updateButton(session, "s2",  value = FALSE)
        updateButton(session, "s2_d",  value = FALSE)
    })
    
    # ahandle model selections
    observeEvent( { input$s1
        input$s1_d
        input$s2
        input$s2_d }, {

            req(values$df_models)

            # fit any newly requested models
            if (input$s1_d == TRUE) { # if the button for a model is clicked
                if ("s1_d_pred" %in% values$df_models$which_model == FALSE) { # if it hasn't already been fit, then fit it and append the values to the summary tibbles
                    values$s1_d_list <- model_all(s1_d_model, "s1_d_pred", values$start_pars, win3d)
                    values$model_list <- c(values$model_list, "s1_d_list" = values$s1_d_list)
                    values$df_models <- values$df_models %>% bind_rows(values$s1_d_list$df_models)
                    values$df_BIC_models <- values$df_BIC_models %>% bind_rows(values$s1_d_list$df_BIC)
                    values$df_tm_models <- values$df_tm_models %>% bind_rows(values$s1_d_list$tm_table_models)
                }}

            if (input$s2 == TRUE) {
                if ("s2_pred" %in% values$df_models$which_model == FALSE) {
                    values$s2_list <- model_all(s2_model, "s2_pred", values$start_pars, win3d)
                    values$model_list <- c(values$model_list, "s2_list" = values$s2_list)
                    values$df_models <- values$df_models %>% bind_rows(values$s2_list$df_models)
                    values$df_BIC_models <- values$df_BIC_models %>% bind_rows(values$s2_list$df_BIC)
                    values$df_tm_models <- values$df_tm_models %>% bind_rows(values$s2_list$tm_table_models)

                }}

            if (input$s2_d == TRUE) {
                if ("s2_d_pred" %in% values$df_models$which_model == FALSE) {
                    values$s2_d_list <- model_all(s2_d_model, "s2_d_pred", values$start_pars, win3d)
                    values$model_list <- c(values$model_list, "s2_d_list" = values$s2_d_list)
                    values$df_models <- values$df_models %>% bind_rows(values$s2_d_list$df_models)
                    values$df_BIC_models <- values$df_BIC_models %>% bind_rows(values$s2_d_list$df_BIC)
                    values$df_tm_models <- values$df_tm_models %>% bind_rows(values$s2_d_list$tm_table_models)
                }}
            # write_rds(values$model_list , "values_model_list.rds")
            # write_rds(values$df_tm_models, "values_df_tm_models.rds")
            # # write_rds(values$df_models, "values_df_models.rds")
            # # write_rds(values$df_BIC_models, "values_df__BIC_models.rds")
            # # write_rds(values$df_tm_models, "values_df_tm_models.rds")
            # # write_rds(values$df_models, "values_df_models.rds")
            # write_rds(values$s1_list, "values_s1_list.rds")
            # write_rds(values$s1_d_list, "values_s1_d_list.rds")
            # write_rds(values$s2_list, "values_s2_list.rds")
            # write_rds(values$s2_d_list, "values_s2_d_list.rds")

            #update the tm table for display df_tm_models_table <- df_tm_models %>%
            model_name_all <- c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred") # doesn't need to be in the server or the observer but is fast enough to justify, since it makes the next step clearer
            model_name_true <- reactive({model_name_all[c(input$s1, input$s1_d, input$s2, input$s2_d)]})

            values$df_tm_models_table <- values$df_tm_models %>%
                dplyr::filter( which_model %in% model_name_true()  ) %>%
                plyr::mutate( which_model = grep_and_gsub(.$which_model, c("s1_pred", "s1_d_pred", "s2_pred","s2_d_pred"), c("Fit 1", "Fit 2", "Fit 3", "Fit 4"), c("Other")))  %>% # move this to later, for the for-display table only!
                set_names(c("Condition", "Model", "Tma 1", "Tma 1 SD", "Tma 2", "Tma 2 SD")) %>%
                discard(~all(is.na(.x)))

            values$df_models_filt <- values$df_models %>% dplyr::filter(which_model %in% model_name_true())
            values$df_BIC_models_filt <- values$df_BIC_models %>% dplyr::filter(which_model %in% model_name_true())

            values$df_models_p <- cond_df_model_for_plot( values$df_models_filt, values$df_BIC_models_filt  )

            values$df_BIC_models_p <- values$df_BIC_models_filt %>%
                cond_df_BIC_for_plot (  ) # adds is_min #readRDS("../4_analyze/values_df__BIC_models.rds")

            values$df_BIC_best <-  cond_df_BIC_for_plot ( values$df_BIC_models_filt   ) %>%
                filter(is_min == TRUE) %>%
                select(c(well, condition, which_model))

            # update which models are available for plotting
            mods_available <- named_mods[c(input$s1, input$s1_d, input$s2, input$s2_d)] # the original named_mods is created outside the server
            updateRadioButtons(session, "choose_model_tm",
                               choices = mods_available,
                               selected = mods_available[1]
            )
        })

    observeEvent( {input$show_BIC_plot
        input$show_best_fits
        input$trim_ends}, {
            model_name_all <- c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred") # doesn't need to be in the server or this observer but is fast enough to justify, since it makes the next step clearer
            model_name_true <- reactive({model_name_all[c(input$s1, input$s1_d, input$s2, input$s2_d)]})

            values$df_models_filt <- values$df_models %>% dplyr::filter(which_model %in% model_name_true())
            values$df_BIC_models_filt <- values$df_BIC_models %>% dplyr::filter(which_model %in% model_name_true())

            values$df_models_p <- cond_df_model_for_plot( values$df_models_filt, values$df_BIC_models_filt  )

            values$df_BIC_models_p <- values$df_BIC_models_filt %>%
                cond_df_BIC_for_plot(  ) # adds is_min #readRDS("../4_analyze/values_df__BIC_models.rds")

            # values$df_BIC_best <-  cond_df_BIC_for_plot ( values$df_BIC_models_filt   ) %>%
            #     filter(is_min == TRUE) %>%
            #     select(c(well, condition, which_model))
        })

    # render the model table
    output$tm_table_render_models <- DT::renderDataTable({ ### new for models
        req(values$df_tm_models_table)
        values$df_tm_models_table
    },
    options = list(scrollX = TRUE, scrollY = 200, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))

    # download the model tms
    output$download_fit_tma <- downloadHandler(
        filename = function() {
            paste0(Sys.Date(), '-dsfworld_tma_by_model_fitting.csv', sep='')
        },
        content = function(file) {
            write.csv(values$df_tm_models_table, file, row.names = FALSE)
        }
    )

    # display the model plot, with all  comonents.
    ## choose the best model
    output$show_BIC_plot_button <- renderUI({
        req(values$df)
        actionButton("show_BIC_plot",
                     p("Display/update fit plot.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),  width = '100%')
    })

    output$show_best_fits_button <- renderUI({
        req(values$df)
        actionButton("show_best_fits",
                     p("Plot selected fits", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),  width = '100%')
    })

    output$best_model_table <- DT::renderDataTable( {
        tryCatch({
            values$fit_sel <-  subset(
                nearPoints(values$df_models_p, #values$df_BIC_models_p,
                           input$plot_dblclick,
                           threshold = 1000, # set large,so anywhere in the plot area will work
                           allRows = TRUE),
                selected_ == TRUE) %>%
                select(c(well, condition, which_model)) %>%
                distinct(which_model, .keep_all = TRUE) %>%
                arrange(condition, well)

            values$df_BIC_best <<- values$df_BIC_best %>%
                filter(! well %in% values$fit_sel$well ) %>% # remove the wells to be overwritten
                full_join(values$fit_sel) %>%
                arrange(condition, well)

            values$df_BIC_display <<- values$df_BIC_best %>%
                mutate(`selected fit` = recode(which_model,
                                               s1_pred = "Fit 1",
                                               s1_d_pred = "Fit 2",
                                               s2_pred = "Fit 3",
                                               s2_d_pred = "Fit 4")) %>%
                select(c(well, condition, `selected fit`))

            values$df_BIC_display # this will render in the table
        }, error = function(e) {

            values$df_BIC_display <<- values$df_BIC_best %>%
                mutate(`selected fit` = recode(which_model,
                                               s1_pred = "Fit 1",
                                               s1_d_pred = "Fit 2",
                                               s2_pred = "Fit 3",
                                               s2_d_pred = "Fit 4")) %>%
                select(c(well, condition, `selected fit`))

            write_rds(values$df_BIC_display, "values_df_BIC_display.rds")
            values$df_BIC_display # this will render in the table
        })


    }, options = list(scrollX = TRUE, scrollY = 200, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))

    ##### end 4_analyze appler server 
    

    ##### downloads
    # Reactive value for selected dataset ----
    # make some download formats
    
    df_mean_sd <- reactive({ 
    df <- values$df %>% # the raw data
          unnest(data) # 
    if (all(c("condition", "Temperature", "mean", "sd") %in% names(df))  ) {
        df %>% 
            select(c(condition, Temperature, mean, sd)) %>%
            distinct(.keep_all = TRUE) %>%
            pivot_wider(names_from = condition, values_from = c(mean, sd))
        } else {
            df %>% 
                select(c(condition, Temperature, value)) %>%
                distinct(.keep_all = TRUE) %>%
                pivot_wider(names_from = condition, values_from = c(value))
    }
    })
    
    df_norm_mean_sd <- reactive({ values$df %>% # the normalized data
            unnest(data) %>%
            select(c(condition, Temperature, mean_norm, sd_norm)) %>%
            distinct(.keep_all = TRUE) %>%
            pivot_wider(names_from = condition, values_from = c(mean_norm, sd_norm))
    })
    
    df_value <- reactive({ values$df %>% # the normalized data
        unnest(data)  %>%
        select(c(condition,well,  Temperature, value))%>%
        distinct(.keep_all = TRUE) %>%
        pivot_wider(names_from = c(well, condition), values_from = c(value))})
    
    df_value_norm <- reactive({ values$df %>% # the normalized data
        unnest(data)  %>%
        select(c(condition, well,  Temperature, value_norm))%>%
        distinct(.keep_all = TRUE) %>%
        pivot_wider(names_from = c(well, condition), values_from = c(value_norm)) })
    
    pred_models <- reactive({ make_wide_preds( values$df_models ) })
    
    df_BIC_models_for_download <- reactive({
        values$df_BIC_models_p  %>%
                                rename(`Best fit` = "is_min") %>%
                                mutate(Fit = recode(which_model,
                                                s1_pred = "Fit 1",
                                                s1_d_pred = "Fit 2",
                                                s2_pred = "Fit 3",
                                                s2_d_pred = "Fit 4")) %>%
                                select(c(well, condition, BIC, Fit, `Best fit`))
            
        
    })
    
    tma_by_dRFU_no_rep_avg <- reactive({
        if ("condition" %in% names(values$df_tms)) {
            # make the averaged table
            values$df_tms %>%
                select(well, condition, dRFU_tma)  
        
        } else {
            values$df_tms %>%
                select(well, dRFU_tma)  
        }
    })
    
    datasetInput1 <- reactive({
        switch(input$dataset1,
               "Tma by dRFU" = tma_by_dRFU(),
               "Tma by best fit" = values$df_tm_models_table, 
               "Replicate-averaged raw data" = df_mean_sd(),
               "Replicate-averaged normalized data" = df_norm_mean_sd(),
               "RFU data with fits" = pred_models() 
        )
    })
    
    datasetInput2 <- reactive({
        switch(input$dataset2,
               "Tma by dRFU, no replicate averaging"= tma_by_dRFU_no_rep_avg(), #head(mtcars) ,
               "Tma by best fit, no replicate averaging"= df_BIC_models_for_download(), # ,
               "Reformatted raw data"= data_raw(), #head(mtcars) ,
               "Normalized raw data"= df_value_norm(),
               "First derivative of raw data" = df_value_norm()
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

    
    # # Table of selected dataset ----
    output$table_set1 <- renderTable(datasetInput1(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))#renderDataTable(datasetInput1(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))
    
    output$table_set2 <- renderTable(datasetInput2(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr')) #renderDataTable(datasetInput2(), options = list(scrollX = TRUE, scrollY = 400, scrollCollapse = TRUE, paging = FALSE, dom = 'tr'))
    
    
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
    
    ##### end downloads
    
    #### instructions panel download layouts
    output$download_ex_layout_96_well <- downloadHandler(
        filename = function() {
            paste('dsfworld_one_variable_layout_96_well.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_one_variable_layout_96_well.csv")
            write.csv(read_csv("dsfworld_one_variable_layout_96_well.csv"), file, row.names = FALSE)
        }
    )
    
    output$download_ex_layout_1 <- downloadHandler(
        filename = function() {
            paste('dsfworld_one_variable_layout.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_one_variable_layout.csv")
            write.csv(read_csv("dsfworld_one_variable_layout.csv"), file, row.names = FALSE)
        }
    )
    
    output$download_ex_layout_2 <- downloadHandler(
        filename = function() {
            paste('dsfworld_two_variable_layout.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_two_variable_layout.csv")
            write.csv(read_csv("dsfworld_two_variable_layout.csv"), file, row.names = FALSE)
        }
    )
    
    output$download_ex_layout_3 <- downloadHandler(
        filename = function() {
            paste('dsfworld_four_variable_layout.csv', sep='')
        },
        content = function(file) {
            read_csv("dsfworld_four_variable_layout.csv")
            write.csv(read_csv("dsfworld_four_variable_layout.csv"), file, row.names = FALSE)
        }
    )
    
    ## download the paper
    output$download_paper <- downloadHandler(
        filename = "20200318_dsfworld_preprint.pdf",
        content = function(file) {
            file.copy("20200318_dsfworld_preprint.pdf", file)
        }
    )
    
    output$download_SI <- downloadHandler(
        filename = "20200318_dsfworld_preprint_SI.pdf",
        content = function(file) {
            file.copy("20200318_dsfworld_preprint_SI.pdf", file)
        }
    )
    
    
    } # end server

###### GUI #####
ui <- navbarPage(useShinyalert(),
                 tabPanel(p("-", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"),
                          column(3),
                          column(6,
                                 p("Welcome to DSF world", style = "font-family: 'Avenir Next'; font-size: 30px; color: white",align = "center"),
                                 tags$hr(style="border-color: black;"),
                                 p("Welcome to DSF world", style = "font-family: 'Avenir Next'; font-size: 50px",align = "center"),
                                 tags$hr(style="border-color: black;"),
                                 shiny::div(tags$img(src = "dye_request_image_v0_small.png", width = "100%"), style = "text-align: center;"),
                                 p("Welcome to DSF world", style = "font-family: 'Avenir Next'; font-size: 15px; color: white",align = "center"),
                                 p("This website was created and is maintained by the Gestwicki lab at UCSF.", style = "font-family: 'Avenir Next'; font-size: 12px",align = "center"),
                                 p("email - dsfworlducsf@gmail.com", style = "font-family: 'Avenir Next'; font-size: 12px",align = "center"),
                                 p("twitter - @GestwickiLab", style = "font-family: 'Avenir Next'; font-size: 12px",align = "center")
                          ),
                          column(3)
                 ),
                 tabPanel( p("to interactive modeling", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "interactive_modeling_mother",
                           tabsetPanel(
                               tabPanel(p("...", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "modeling_welcome",
                                        column(3),
                                        column(6,
                                               p("   ", style = "font-family: 'Avenir Next'; font-size: 30px",align = "center"),
                                               p("Welcome to interactive modeling!", style = "font-family: 'Avenir Next'; font-size: 30px",align = "center"),
                                               p("DSF is a simple readout for a complex phenomenon. These interactive models showcase how the processes underlying DSF--thermodynamics and kinetics of unfolding, and relative dye detection of various protein states--can ultimately effect the data obtained.", style = "font-family: 'Avenir Next'; font-size: 20px",align = "left")
                                        ),
                                        column(3)
                               ), # end tabPanel
                               # DSF data modeler  -----------------------------------------------------------------------
                               tabPanel(p("interactive modeling", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), value = "model_plot",
                                        tags$head( # set the slider aesthetic
                                            tags$style(
                                                ".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-2 .irs-single, .js-irs-2 .irs-bar-edge, .js-irs-2 .irs-bar {background: grey; border-color: transparent;}",
                                                
                                                ".js-irs-3 .irs-single, .js-irs-3 .irs-bar-edge, .js-irs-3 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-4 .irs-single, .js-irs-4 .irs-bar-edge, .js-irs-4 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-5 .irs-single, .js-irs-5 .irs-bar-edge, .js-irs-5 .irs-bar {background: grey; border-color: transparent;}",
                                                
                                                ".js-irs-6 .irs-single, .js-irs-6 .irs-bar-edge, .js-irs-6 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-7 .irs-single, .js-irs-7 .irs-bar-edge, .js-irs-7 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-8 .irs-single, .js-irs-8 .irs-bar-edge, .js-irs-8 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-9 .irs-single, .js-irs-9 .irs-bar-edge, .js-irs-9 .irs-bar {background: grey; border-color: transparent;}",
                                                ".js-irs-10 .irs-single, .js-irs-9 .irs-bar-edge, .js-irs-10 .irs-bar {background: grey; border-color: transparent;}"
                                                
                                            )
                                        ),
                                        sidebarLayout(
                                            sidebarPanel(
                                                p("Tune model parameters below", style = "font-family: 'Avenir Next'; font-size: 20px; color: black",align = "center"),
                                                # # "Tune model parameters",
                                                bsCollapse(id = "thermo_pars", open = "Panel 1",
                                                           bsCollapsePanel(p("Thermodynamic parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                           sliderInput("T_half_", "Thermodynamic melting temperature (C)", min = 25, max = 95, value = 55, step = 1),
                                                                           bsTooltip("T_half_", "At this temperature, the equilibrium ratio of folded to reversibly unfolded states is 1:1.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("dHu_", "Enthalpy of unfolding (kJ/mol)", min = 1, max = 600, value = 250, step = 10),
                                                                           bsTooltip("dHu_", "Change in enthalpy between the folded and reversibly unfolded state. This model assumes the enthalpy of reversibly and irreversibly unfolded states are equal, a common simplification.",
                                                                                     "right", options = list(container = "body", color = "white")),
                                                                           
                                                                           sliderInput("dCp_", "Change in heat capacity with unfolding (kJ/mol)", min = 0.1, max = 50, value = 8, step = 0.5),
                                                                           bsTooltip("dCp_", "When a protein unfolds, the heat capacity of the solution changes. The magnitude of this change influences the Tm of the protein.",
                                                                                     "right", options = list(container = "body", color = "white"))
                                                                           , style = "default")
                                                           
                                                           
                                                ),
                                                
                                                bsCollapse(id = "kinetic_pars", open = "Panel 2",
                                                           bsCollapsePanel(p("Kinetic parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                           sliderInput("Ea_", "Activation energy of unfolding (kJ/mol)", min = 10, max = 1000, value = 200, step = 10),
                                                                           bsTooltip("Ea_", "The energy barrier between the reversibly and irreversibly unfolded states. In this model, this value changes with temperature according to the Arrhenius Equation.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("T_star_", "Temperature at which irreversible unfolding becomes significant, T* (C)", min = 25, max = 95, value = 55, step = 1),
                                                                           bsTooltip("T_star_", "Specifically, the temperature at which the rate constant of irreversible unfolding is 1/min.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("v_", "DSF experiment heating rate (C/min)", min = 0.1, max = 15, value = 1, step = 0.1),
                                                                           bsTooltip("v_", "The thermocycling ramp rate used in the simulated DSF experiment.",
                                                                                     "right", options = list(container = "body")),
                                                                           style = "default")
                                                ),
                                                
                                                bsCollapse(id = "dye_pars", open = "Panel 3",
                                                           bsCollapsePanel(p("Dye parameters", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), 
                                                                           # the heading above causes this warning: Warning in if (getAttribs(panels[[i]])$value %in% open) { : the condition has length > 1 and only the first element will be used
                                                                           sliderInput("nat_dye_", "Folded state", min = 0, max = 1, value = 0, step = 0.1),
                                                                           bsTooltip("nat_dye_", "The degree of dye binding and activation observed to the folded state of the protein. Detection of the folded state underlies many high background issues with DSF for hydrophobic proteins.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("unf_dye_", "Reversibly unfolded state", min = 0, max = 1, value = 1, step = 0.1),
                                                                           bsTooltip("unf_dye_", "The degree of dye binding and activation observed to the reverisibly unfolded state of the protein. Low detection of this state may underlie the invisibility of some proteins in DSF.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("fin_dye_", "Irreversibly unfolded state", min = 0, max = 1, value = 1, step = 0.1),
                                                                           bsTooltip("fin_dye_", "The degree of dye binding and activation observed to the irreversibly unfolded state of the protein. This parameter is non-zero for most proteins.",
                                                                                     "right", options = list(container = "body")),
                                                                           
                                                                           sliderInput("decay_rate_", "Temperature sensitivity of dye activation", min = 0, max = 1, value = 0.85, step = 0.1),
                                                                           bsTooltip("decay_rate_", "If all dye binding sites remained constant over the experiment, this is the rate at which fluorescence would decrease with temperature. Empirically this is close to 0.8, and is related to the general temperature-sensitivity of fluorophore quantum yields and the strength of hydrophobic effect.",
                                                                                     "right", options = list(container = "body"))
                                                                           
                                                                           ,style = "default")
                                                )
                                                
                                            ),
                                            mainPanel(
                                                p("Welcome to DSF world", style = "font-family: 'Avenir Next'; font-size: 8px; color: white",align = "center"),
                                                shiny::div(tags$img(src = "20191108_dsfworld_modeling_scheme_v0.png", width = "100%"), style = "text-align: center;"),
                                                p("Welcome to DSF world", style = "font-family: 'Avenir Next'; font-size: 8px; color: white",align = "center"),
                                                plotOutput("plot_model", width = "100%", height = "400px")
                                            ))
                                        
                               ), # end fluidRow
                               
                               tabPanel(p("about the model", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"), value = "learn_about_model",
                                        
                                        p("about the model", style = "font-family: 'Avenir Next'; font-size: 25px; color: white",align = "center"),
                                        shiny::div(tags$iframe(src = "dsfworld_about_the_model.pdf", width = "100%", height = "500px"), style = "text-align: center;"),
                                        
                               ))),
                 
                 # Data Analysis --------------------------------------------------------------------------------
                 tabPanel(p("to data analysis", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "data_analysis_mother", # end tab panel (tabset, div, main still remaining)
                          tabsetPanel(id = "inTabset_analysis", # tabset for all analysis sub-tabs, used by the "jumpToAnalysis" button 
                                      tabPanel(p("1 | upload data", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "uploads_tab",
                                               
                                               ###### begin UI from uploads applet 
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
                                                       dataTableOutput("input_file") %>% withSpinner(color="#525252"), style = "overflow-x: scroll;"
                                                   ) # end main panel
                                               ) # end sidebarLayout
                                               
                                               ###### end UI from uploads 
                                               
                                      ), # end tabpanel
                                      # analyze and visualize --------------------------- 
                                      tabPanel(p("2 | analyze and visualize", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "analysis_tab",
                                               
                                               ####### 2_layouts applet GUI
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
                                                                                                  downloadButton("sample_layout_file", p("Download example layout", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center")), #
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
                                                       ),
                                                       # tm determination ui --------------------------- 
                                                       bsCollapse(id = "tm_table", open = "Panel 2",
                                                                  bsCollapsePanel(p("Find apparent Tms", style = "font-family: 'Avenir Next'; font-size: 16px; color: black",align = "center"),
                                                                                  bsCollapsePanel(p("By dRFU", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                                                  DT::dataTableOutput("tm_table_render"), #style = "height:400px;"
                                                                                                  p("-----", style = "font-family: 'Avenir Next'; font-size: 10px; color: white",align = "center"),
                                                                                                  downloadButton('download_dRFU_tma', "Download Tmas by dRFU"),
                                                                                                  style = "default"
                                                                                  ),
                                                                                  
                                                                                  bsCollapsePanel(p("By sigmoid fitting", style = "font-family: 'Avenir Next'; font-size: 14px; color: black",align = "center"),
                                                                                                  #wellPanel(
                                                                                                  p("Select the models you would like to fit to your data below.", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center") %>% strong(),
                                                                                                  splitLayout(cellWidths = c("25%", "25%", "25%", "25%"), 
                                                                                                              p("Fit 1", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "center"),
                                                                                                              p("Fit 2", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "center"),
                                                                                                              p("Fit 3", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "center"),
                                                                                                              p("Fit 4", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "center")
                                                                                                  ),
                                                                                                  
                                                                                                  splitLayout(cellWidths = c("25%", "25%", "25%", "25%"), 
                                                                                                              bsButton("s1", label = tags$img(src = "s1_v1.png",
                                                                                                                                              width = "100%"),
                                                                                                                       block = TRUE, type = "toggle", value = TRUE), 
                                                                                                              bsButton("s1_d", label = tags$img(src = "s1_id_v1.png",
                                                                                                                                                width = "100%"),
                                                                                                                       block = TRUE, type = "toggle", value = FALSE),
                                                                                                              bsButton("s2", label = tags$img(src = "s2_v1.png",
                                                                                                                                              width = "100%"),
                                                                                                                       block = TRUE, type = "toggle", value = FALSE), 
                                                                                                              bsButton("s2_d", label = tags$img(src = "s2_id_v1.png",
                                                                                                                                                width = "100%"),
                                                                                                                       block = TRUE, type = "toggle", value = FALSE)
                                                                                                  ),
                                                                                                  
                                                                                                  bsTooltip("s1", "Fit 1: One sigmoid with decay",
                                                                                                            "right", options = list(container = "body")),
                                                                                                  bsTooltip("s1_d", "Fit 2: One sigmoid with decay and starting-temperature fluorescence",
                                                                                                            "right", options = list(container = "body")),
                                                                                                  bsTooltip("s2", "Fit 3: Two sigmoids with decays",
                                                                                                            "right", options = list(container = "body")),
                                                                                                  bsTooltip("s2_d", "Fit 4: Two sigmoids with decays and starting-temperature fluorescence",
                                                                                                            "right", options = list(container = "body")),
                                                                                                  p(" ", style = "font-family: 'Avenir Next'; font-size: 8px; color: black",align = "center"),
                                                                                                  
                                                                                                  
                                                                                                  
                                                                                                  DT::dataTableOutput("tm_table_render_models"), #style = "height:400px;"
                                                                                                  p("-----", style = "font-family: 'Avenir Next'; font-size: 10px; color: white",align = "center"),
                                                                                                  downloadButton('download_fit_tma', "Download Tmas from fits"),
                                                                                                  p("-----", style = "font-family: 'Avenir Next'; font-size: 30px; color: white",align = "center"),
                                                                                                  uiOutput("trim_ends"),  
                                                                                                  p(" ", style = "font-family: 'Avenir Next'; font-size: 8px; color: black",align = "center"),
                                                                                                  bsCollapsePanel(p("Select the best fit for each dataset", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center"),
                                                                                                                  uiOutput("show_BIC_plot_button"),
                                                                                                                  DT::dataTableOutput("best_model_table"),
                                                                                                                  p("", style = "font-family: 'Avenir Next'; font-size: 8px; color: black",align = "center"),
                                                                                                                  p("Selecting fits", style = "font-family: 'Avenir Next'; font-size: 12px; color: black",align = "center") %>% strong(),
                                                                                                                  
                                                                                                                  p("For each condition, the fit option (1-4) with the lowest Bayesian Information Criterion (BIC) is selected by default. This is meant to maximize model quality without over-fitting (see 'About the analysis').", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "left"),
                                                                                                                  
                                                                                                                  p("", style = "font-family: 'Avenir Next'; font-size: 8px; color: black",align = "center"),
                                                                                                                  
                                                                                                                  p("However, you can select a fit manually by double-clicking on the desired fit in the 'All fits' plot.", style = "font-family: 'Avenir Next'; font-size: 10px; color: black",align = "center"),
                                                                                                                  uiOutput("show_best_fits_button")
                                                                                                  )
                                                                                  ), style = "default"
                                                                  )
                                                       ), width = 5
                                                       
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
                                                           downloadButton('download_plot', "Download plot"),
                                                           
                                                           plotOutput("plot", 
                                                                      height = "auto", 
                                                                      dblclick = dblclickOpts(
                                                                          id = "plot_dblclick")) %>% 
                                                               withSpinner(color="#525252"), style = ("overflow-y:scroll; max-height: 600px"),
                                                           verbatimTextOutput("dblclick_info")
                                                           
                                                       ), width = 7)  
                                               )
                                               
                                               ###### end 2_layouts applet GUI
                                      ), # end tabPanel
                                      # Analysis visualzation panel ----------              
                                      tabPanel(p("3 | download results", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "downloads_tab", 
                                               
                                               ##### begin downloads applet GUI
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
                                                           tabPanel("Preview: quick downloads", 
                                                                    tableOutput("table_set1"),  
                                                                    #dataTableOutput("table_set1"), 
                                                                    style = "overflow-x: scroll;overflow-y: scroll;"),
                                                           tabPanel("Preview: supplemental files", 
                                                                    tableOutput("table_set2"),  
                                                                    #dataTableOutput("table_set2"), 
                                                                    style = "overflow-x: scroll;")
                                                       )
                                                   ) # end main panel
                                               ) # end sidebarLayout
                                               ##### end downloads applet GUI
                                      ), # end tabpanel
                                      # About the analysis ------
                                      tabPanel(p("about the analysis", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "about_analysis_tab",
                                               shiny::div(tags$iframe(src = "dsfworld_about_the_analysis.pdf", width = "100%", height = "500px"), style = "text-align: center;")
                                      ), # end tabPanel
                                      
                                      
                                      tabPanel(p("instructions", style = "font-family: 'Avenir Next'; font-size: 15px; color: black",align = "center"), value = "instructions_tab",
                                               
                                               column(3),
                                               column(6,
                                                      tags$div(
                                                          tags$h1("Instructions for DSFworld data analysis"), 
                                                          tags$p("These instructions explain how to use DSFworld to visualize and analyze raw DSF data. For more information on the analyses performed, please see the “about the analysis” tab."),
                                                          tags$h3("1. Upload raw data"), 
                                                          tags$p("To use DSFworld, start by uploading your raw fluorescence versus temperature data. Before moving on to analysis, your uploaded data must be formatted with Temperature in the first column, and fluorescence readings in the columns to the right. DSFworld offers a few tools to help you do this. The uploaded data is displayed in its current format in a table to the right of the grey uploads panel to guide you."),
                                                          # uploading image
                                                          shiny::div(tags$img(src = "well_formatted_data.png", width = "100%"), style = "text-align: center;"),
                                                          tags$h5("DSFworld offers some tools to help you achieve this. "),
                                                          tags$ul(
                                                              tags$li(tags$strong("Reformat raw from instrument."),
                                                                      "If the data was collected on any of the instruments listed under “Reformat raw from instrument” (currently, Biorad, Stratagene, quantStudio, and qTower), raw data can be uploaded exactly as exported from the instrument and DSFworld will automatically format the data for you."), 
                                                              tags$br(),
                                                              tags$li(tags$strong("Pre-format by template."),
                                                                      "If your instrument is not listed under the supported reformatting, format your data prior to uploading with Temperature in the first column, and fluorescence values for each well in the columns to the right. A properly-formatted example file can be downloaded under “Uploading instructions” in the “upload data” tab. Be sure to save your data as a comma-separated value (csv) prior to uploading."),
                                                              tags$br(),
                                                              tags$li(tags$strong("Convert cycle number to temperature."),
                                                                      "Some qPCR instruments export cycle numbers (e.g. 1, 2, 3 ...) for each measurement, instead of the temperature corresponding to that measurement. You can convert these cycle numbers to temperatures under the 'Convert cycle number to temperature' panel by setting a starting temperature, and a number of degrees to increase per cycle."),
                                                              tags$br()
                                                          ),
                                                          
                                                          tags$h5("A few tips on uploading data:"),
                                                          tags$ul(
                                                              tags$li(tags$strong("There is no hard limit to the number of wells you can analyze in a single upload at DSFworld."),
                                                                      "We have tested up to full 384-well plates without trouble. However, DSFworld is not intended for use for high throughput screening. Be aware that for larger datasets (e.g. > 50 wells), displaying plots and tables will take a bit more time."), 
                                                              tags$br(),
                                                              tags$li(tags$strong("We strongly recommend using well names (e.g. A1, A2 ...) for your data columns."),
                                                                      "This is the default for most instruments. While it isn't necessary for analysis, having well names for columns enables DSFworld to assign any experimental variables you define to your data in the 'layouts' portion of analysis. This is useful because it lets you make descriptive visualizations of your data (e.g. color lines based on compound concentration). If your data doesn't have well names but you want to use these features, you can overwrite your column names with well names in the uploads panel--see below for more information."), 
                                                              tags$br(),
                                                              tags$li(tags$strong("If data isn't upload successfully, try re-saving it as a UTF-8 csv."),
                                                                      "(even if it is already a csv!) and re-uploading. You can do this on your computer with Save As > CSV, UTF-8 encoding."), 
                                                          ),
                                                          
                                                          tags$h3("2. Analyze."), 
                                                          tags$p("Once data is properly formatted, it can be visualized and analyzed in the analysis tab. DSFworld has features for both visualizing and analyzing the uploaded data. Apparent melting temperatures are displayed in a drop-down table at left of the analysis tab, from where they can be directly downloaded. These melting temperatures available immediately upon proceeding to the analysis tab, and do not depend on plotting or plate layouts. However, because it is always good to visualize your data before exporting processed results, the following instructions are presented in the order in which we recommend proceeding through analysis."),
                                                          tags$ol(
                                                              tags$li( tags$strong("Define your experimental layout."), 
                                                                       "If you would like, you can define an experiment layout for the uploaded data--that is, what conditions were tested in each well. This is useful because it allows DSFworld to average experimental replicates for you, and more importantly, allows you to to use the plotting features of DSFworld to make plots that visualize the effects of your experimental variables."), 
                                                              tags$p("You can define a single experimental variable by entering the variable in the editable table under 'Set plate layout and replicates' tab, and pressing 'Update names from manual table'. For replicates, enter identical names."),
                                                              tags$p(" For more complicated experiments with more than one variable, you can define any number of experimental variables (e.g. protein, compound name, and compound concentration) by uploading a plate layout as a csv. The point of an experimental layout is to tell DSFworld what conditions are present in each well."),
                                                              tags$p("It's often easiest to get the hang of plate layouts by looking at examples."),
                                                              tags$ul(
                                                                  tags$li(tags$strong("Here's an example of a simple layout with only one variable."),
                                                                          "In this example, DSF was performed on three different mutants of a protein. The name of the experimental variable is entered to the left of the block--in this example, its 'Mutant'."), 
                                                                  downloadButton("download_ex_layout_1", "Download one variable layout example", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                                                  downloadButton("download_ex_layout_96_well", "Download one variable layout example, 96-well", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                                                  tags$p(),
                                                                  tags$br(),
                                                                  tags$li(tags$strong("Here's an example of a layout with two variables."),
                                                                          "SYPRO Orange concentration, and protein concentration. The second variable is defined in a new block, directly below the first. This is a layout for an experiment in which we compared DSF data collected for lysozyme (10 µM SYPRO Orange, 1 µM lysozyme) to three different negative controls: buffer alone, 10 µM SYPRO Orange without protein, and 1 µM lysozyme without dye. Fun fact: this is the experimental layout that we use to test whether a new batch of microtiter plates is compatible with DSF. It's the same one used to generate the data presented in Figure 3a and Supplemental Figure 4 of the paper associated with this website."), 
                                                                  downloadButton("download_ex_layout_2", "Download two variable layout example", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                                                  tags$p(),
                                                                  tags$br(),
                                                                  tags$li(tags$strong("Here's an example of a layout for a more complicated experiment."),
                                                                          "This time there are four variables: SYPRO Orange concentration, compound, compound concentration, and protein concentration. Like before, each experimental variable is defined in its own block, and each new block is pasted directly below the block above it, and the name of that variable defined in the left-most column. This is a layout for an experiment in which we tested the impact of aggregation of four compounds (vemurafenib, miconazole, clotrimazole, and ritonavir) on DSF signal in the absence and presence of lysozyme. Fun fact: this is the experimental layout used to generate some of the results presented in Figure 4 and Supplementary Figures 5-7 of the paper associated with this website."),
                                                                  downloadButton("download_ex_layout_3", "Download four variable layout example", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                                                  tags$p(),
                                                                  tags$br()
                                                              ),
                                                              tags$li(tags$strong("Make custom plots."), 
                                                                      tags$p("The uploaded experimental variables are automatically made available for the creation of custom plots. These variables can be used to create sub-plots, vary color, or vary line-type. Any wells which are identical in all experimental variables are considered replicates, and data can be plotted either as individual replicates, or as means with standard deviations. For sub-plots, the y-axes can be fixed as equal, or allowed to vary between the sub-plots."),
                                                                      tags$p("The following plot aesthetics can also be customized: titles for the plot, legend, and x and y axes can be set; text size can be altered; and the legend (for colors and/or linetypes) can be shown or hidden."),
                                                                      tags$p("To apply changes to a plot, press the 'Update plot' button. The currently displayed plot can be downloaded from the analysis window using the 'Download plot' button.")), 
                                                              
                                                              tags$li(tags$strong("Determine apparent melting temperatures."), "Apparent melting temperatures can be calculated by either the maximum of the first derivative, or by sigmoid fitting. For more information on the analysis methods, see the 'about the analysis' tab."),
                                                              tags$ul(
                                                                  tags$li( tags$strong("Maximum of the first derivative (dRFU)"), 
                                                                           "Apparent melting temperatures are calculated by maximum of the first derivative automatically. Experimental replicates are averaged automatically. Results without replicate averaging can be downloaded from the Downloads tab."), 
                                                                  tags$li( tags$strong("Sigmoid fitting"), 
                                                                           tags$p("DSFworld offers four different models for sigmoid fitting, termed 'Fit 1', 'Fit 2', 'Fit 3', and 'Fit 4'. Upon uploading, apparent melting temperatures are automatically calculated by Fit 1. Additional models can be fit to the data by clicking the buttons for that fit in the 'By sigmoid fitting' drop-down menu."),
                                                                           tags$p("You can trim the temperature range used for the fits using the slider-bar."),
                                                                           tags$p("If multiple fits are applied to the data, DSFworld will automatically display the Fit with the lowest Bayesian Information Criterion (BIC). However, the best fit can also be selected manually. To do this, open the 'Select the best fit for each dataset'  dropdown tab. A table is displayed in this drop-down window which shows the currently-selected fit for each data set. Click the 'Display/update fit plot.' button. A new plot will be shown in the plot panel. This plot is specifically designed to facilitate manual fit selection. Each individual sub-plot shows a single data set, with one Fit option. Each row of sub-plots contains an individual data set, and each column shows a fit option. To select a particular fit option for a given data set, double click on that sub-plot."),
                                                                           tags$p("To plot only the selected fits, press the 'Plot selected fits' button.")
                                                                  )
                                                              )
                                                          ),
                                                          tags$h3("3. Download."),
                                                          tags$p("All plots and replicate-averaged apparent melting temperatures can be downloaded directly from the analysis window. However, if you would like to download additional forms of the data, this can be done in the 'download results' panel. If you would like to download results for a particular fit, you must first select that fit in the data analysis window.")),
                                                      tags$br(),
                                                      tags$br()),
                                               column(3)
                                               
                                               
                                               ) # end tabPanel
                          )), # end tabset Panel (contains all "analysis sub-panels)
                 tabPanel( p(" . . .", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "closing_remarks_tab",
                           tabPanel( p(" . . .", style = "font-family: 'Avenir Next'; font-size: 20px; color: grey",align = "center"), value = "closing_remarks_tab",
                                     column(4,
                                            shiny::div(tags$img(src = "dsfworld_logo_grey.png", width = "400px"), style = "text-align: center;")
                                     ),
                                     column(6,
                                            tags$br(),
                                            tags$br(),
                                            tags$p("DSFworld was created to help users complete more successful DSF experiments."), 
                                            #tags$p("We hope that the interactive models and data analyses offered here can help users develop a strong working relationship with DSF results--both the underlying concepts and real data."),
                                            tags$p("The code for DSFworld is available on GitHub from https://github.com/gestwicki-lab/dsfworld. This includes the full code for this website and all associated scripts, as well as modular mini-web applications for each of the tasks tackled here: interactive modeling, data uploading, data layouts, plotting, analyses, and downloads."),
                                            tags$p("DSFworld is presented in a publication, alongside practical tips for DSF and a deeper discussion of Model 2. You can download a version of that paper, and its supplementary information below."),
                                            downloadButton("download_paper", "Download the full paper", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                            downloadButton("download_SI", "Download the supplementary information", width = '50%',style="font-size: 14px; color: #00000; background-color: #fffff; border-color: #ffff"),
                                            tags$p(),
                                            tags$br(),
                                            tags$p("Thank you to the many beta-testers of DSFworld, especially Ziyang Zhang, Douglas Wassarman, and Sarah Williams."),
                                            tags$p("DSFworld is written in R. we're indebted to the many creators of R, R Studio, and the R packages who made this project possible. Particularly, thank you to Joe Cheng, Hadley Wickham, and their teams for creating R Shiny and the tidyverse."),
                                            
                                            tags$ul(
                                                tags$li("The user interface was created using R Shiny, as well as the packages shinyBS for drop-down panels, shinyalert for pop-up messages, shinycssloaders for busy spinners, and rhandsontable for manual entry and editing of condition names in the analysis window. "), 
                                                tags$br(),
                                                tags$li("The model fitting uses modelr and broom for data and model handling, minpack.lm for fitting, signal for Savistky-Golay filtering, quantmod to assist in the finding of local maximima and minima for starting parameter estimates, and SciViews to perform natural logarithms."),
                                                tags$br(),
                                                tags$br()
                                                )),
                                     column(2))# end tabpanel
                           )
                 )
###### END GUI #####

shinyApp(ui = ui, server = server)
