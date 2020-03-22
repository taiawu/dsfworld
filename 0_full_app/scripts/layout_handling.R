#
# new daughter layout function
df_to_layout <- function(df, layout_type) {
  df_m <-   set_names( df ,  c("type","row",as.numeric( df [1,-c(1,2)]))) %>%
    . [ -1 , -1] %>%
    reshape2::melt( . ,id.vars = "row") %>%
    mutate( . , well = as_vector(map2( . $row,  . $variable, paste0)) ) %>%
    set_names( . , c("row", "column", layout_type, "well"))
  df_m
}

make_layout <- function( filename ) { # from path to raw layout to a final fomatted layout file 
  # read the layout file, and split each layout into an individual 
  layout_list <- data.table::fread( filename, header = TRUE) %>%
    as_tibble() %>%
    split( . ,  . $Type)
  
  # put into a merge-able form
  layout <- df_to_layout(layout_list[[1]], names(layout_list)[[1]])[c(1,2,4)] # initialize the list
  for (i in c(1:length(layout_list))) {
    layout <- layout %>%
      mutate("var" =  as_vector(df_to_layout(layout_list[[i]], layout_type = names(layout_list)[[i]])[3] )) %>% # append the column of interest
      set_names(c(names(layout), names(layout_list)[[i]])) # rename based on the column of interest
  }
  layout <- layout %>% 
            unite("condition", c(4:ncol(.)), remove = FALSE) %>% # create a unique column, used to define groups after averaging
            mutate_if(is.factor, as.character)
            
  layout
}

nest_raw <- function( data_raw ) {
  df <- data_raw %>% # this is the active dataframe, used for plotting and calculations
    gather(well, value, -Temperature) %>% # call whatever column names "well"
    group_by(well) %>%
    mutate(value_norm = BBmisc::normalize(value, method = "range", range = c(0,1)), ###### if we do this as a mutate, it ignores the groups!!!!!!!
           Temperature_norm = BBmisc::normalize(Temperature, method = "range", range = c(0,1)))  %>%
    nest_legacy() %>%
    #plyr::mutate(new_names = well)
    plyr::mutate(condition = well)
}

# parse_well_vec <- function( well_vec ){
#   l <- list(
#     col =  parse_number(well_vec),
#     row = str_extract_all(well_vec, "[A-Z; a-z]", simplify = TRUE) 
#     %>% str_to_upper(locale = "en") 
#     %>% as_vector()
#   )
#   l
# }
parse_well_vec <- function( well_vec ){
  col <- parse_number(well_vec) %>% as.character() # pull out all numbers. if no numbers, 
  row <- str_remove(well_vec, col %>% replace_na(replace = "1"))
  
  l <- list(
    col =  col,
    row = row %>% str_to_upper(locale = "en") %>% as_vector()
  )
  l
}

add_standardized_wells <- function( df, make_factor ) {
  df %>%
    mutate(row_ = parse_well_vec(.$well)$row,
           col_ = parse_well_vec(.$well)$col) %>%
    mutate(well_ = map2(.$row_, .$col_, paste0) %>% as_vector()) %>%
    mutate(well_f_ = factor(.$well_, levels = make_well_names("ROWS", "1"))) # as a factor, so things will order correctly
}

ensure_standardized_wells <- function( df ) {
  if ( all(c("well_", "well_f_", "row_", "col_") %in% names(df)) == FALSE ) { # if standardized well columns are missing
    df_out <- df %>%
      add_standardized_wells()
  } else {
    df_out <- df
  }
  
  df_out 
}

join_layout_nest <- function(by_well, layout) { 
  by_well_ <- ensure_standardized_wells(by_well) # this will always be fresh, un-layout-joined dataframe
  layout_ <- ensure_standardized_wells(layout)
  l_names <- names(layout_)
  dup_cols <- l_names[!l_names %in% c("well_", "well_f_", "row_", "col_")]# columns which are already in the layout
  #dup_cols <- names(layout_)[!c(names(layout_) %in% c("well_", "well_f_", "row_", "col_"))]
  
  common_cols <- c("well","well_", "well_f_", "row_", "col_", "row", "column") %>%
                 .[. %in% names(layout_) ]

  if (!"condition" %in% names(layout_)) { # this should never happen....
    print("no_cond see join_layout_nest in layout_handling.R")
    layout_ <- layout_ %>%
                unite("condition", -one_of(c("well","well_", "well_f_", "row_", "col_", "row", "column")), remove = FALSE)
  } else { # if "condition" is present in the layout
    if (  all(l_names[!l_names == "condition"] %in% common_cols) == TRUE ) { # if "condition" is the only column unique to the layout
      common_cols <- c("well_", "well_f_", "row_", "col_", "row", "column") # retain the well column
      layout_ <- layout_ %>%
        # select(-condition) %>%
        # unite("condition", -one_of(common_cols), remove = FALSE) %>%
        mutate_if(is.factor, as.character)
      
      print("all names in common cols ZZ")
    } else { # if a layout has been added 
    layout_ <- layout_ %>%
               select(-condition) %>%
               unite("condition", -one_of(common_cols), remove = FALSE) %>%
               mutate_if(is.factor, as.character)
    }
  }

    # layout_ <- layout_ %>%
    #            select(-condition) %>%
    #            unite("condition", -one_of(common_cols), remove = FALSE) %>%
    #            mutate_if(is.factor, as.character)
    # %>% # clear the existing condition column
    #            unite("condition", -one_of(common_cols), remove = FALSE) 
  #   %>% # create a unique column, used to define groups after averaging
  #              mutate_if(is.factor, as.character)
   
  # join and re-average the data according to the new conditions
  by_well_ %>%
    unnest_legacy() %>%
    dplyr::select(-one_of(dup_cols )) %>% # if it's already in layout, drop it
    dplyr::left_join(., layout_, by = c("well_", "well_f_", "row_", "col_")) %>%
    group_by(condition, Temperature) %>%
    dplyr::mutate(mean = mean(value),
                  sd = sd(value),
                  mean_norm = mean(value_norm),
                  sd_norm = sd(value_norm)) %>%
    ungroup() %>%
    nest(data = c(Temperature, value, mean, sd, Temperature_norm, value_norm, mean_norm, sd_norm)) # THIS IS NOT nest_legacy
}


dsfworld_default <- theme( # adapted from free amino acids hit call
  text = element_text(size = 10),
  # axis.title.x = element_blank(),
  # axis.title.y = element_blank(),
  axis.text = element_text(size = 8),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 12),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(), 
  strip.background = element_blank(),
  aspect.ratio = (1/1.618)
  # axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  # axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')
)

# https://ggplot2.tidyverse.org/reference/aes.html
# https://ggplot2.tidyverse.org/reference/vars.html

facet_func <- function(df, mean_or_each, color_by, linetype_by, facet, facet_by, facet_rows, facet_cols, set_title, legend_title, legend_linetype_title, use_linetypes, fix_free, text_size, legend_position, x_title, y_title) {
  #df <- dplyr::filter(df, Allele != !!enquo(drop_cols))
  #if (color_by == "Uncolored" ) { color_by <- NULL } else { color_by <- enquo(color_by) }
  color_by <- enquo(color_by) # tell how to color
  linetype_by <- enquo(linetype_by)
  #print(length(unique(df$color_by)))
  
  if (mean_or_each == "mean") {
    p <- ggplot(df, aes(x = Temperature, y = mean, group = well, color = !!color_by)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.5)
  } else if (mean_or_each == "each") {
    p <- ggplot(df, aes(x = Temperature, y = value, group = well, color = !!color_by))
    
  }
  
  if (facet == "wrap") {
    facet_by = enquo(facet_by)
    p <- p + 
      facet_wrap(vars(!!facet_by), scales = fix_free)
  } 
  else if (facet == "grid") {
    facet_rows = enquo(facet_rows)
    facet_cols = enquo(facet_cols)
    p <- p +
      facet_grid(rows = vars(!!facet_rows), cols = vars(!!facet_cols), scales = fix_free)
  }
  
  if (use_linetypes == TRUE ) {
    p <- p +
      geom_line(  aes( linetype = !!linetype_by ))  + ###delete the aes to revert
      theme_bw() + 
      scale_color_viridis_d() +
      labs(title = set_title, color = legend_title, x = x_title, y = y_title, linetype = legend_linetype_title) +
      dsfworld_default +
      theme(  text = element_text(size = text_size*1.25),
              axis.text = element_text(size = text_size),
              plot.title = element_text(lineheight=.8, face="bold", size = text_size*1.5),
              legend.position = legend_position)
  } else {
    p <- p +
      geom_line()  + ###delete the aes to revert
      theme_bw() + 
      scale_color_viridis_d() +
      labs(title = set_title, color = legend_title, x = x_title, y = y_title) +
      dsfworld_default +
      theme(  text = element_text(size = text_size*1.25),
              axis.text = element_text(size = text_size),
              plot.title = element_text(lineheight=.8, face="bold", size = text_size*1.5),
              legend.position = legend_position) 
    
  }
  
  p
}

gg_facet_nrow_ng <- function(p){ # determine the number of rows in a ggplot
  assertive.types::assert_is_any_of(p, 'ggplot')
  p %>%
    ggplot2::ggplot_build() %>%
    magrittr::extract2('layout') %>%
    magrittr::extract2('layout') %>%
    magrittr::extract2('ROW') %>%
    unique() %>%
    length()
}

gg_facet_ncol_ng <- function(p){ # determine the number of cols in a ggplot
  assertive.types::assert_is_any_of(p, 'ggplot')
  p %>%
    ggplot2::ggplot_build() %>%
    magrittr::extract2('layout') %>% 
    magrittr::extract2('layout') %>%
    magrittr::extract2('COL') %>%
    unique() %>%
    length()
}

plotDownloadUI <- function(id, height = 400) {
  ns <- NS(id)
  downloadButton(ns("download_plot"), "Download figure")
}


plotDownload <- function(input, output, session, plotFun) { #https://github.com/czeildi/shiny-modules-examples
  output$plot <- renderPlot({
    plotFun()
  })
  
  output$download_plot <- downloadHandler(
    filename = function() {
      #filename <- paste0("test", ".pdf")#paste0(input$plot_download_name, ".pdf")
      paste0(as.character(input$plot_download_name), ".pdf")
    },
    content = function(file) {
      ggsave(file, width = session$clientData$output_data_width/100, height = height/100, limitsize = FALSE)
    }
  )
}