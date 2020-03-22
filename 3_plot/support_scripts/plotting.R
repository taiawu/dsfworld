get_layout_vars <- function( df ) {
  layout_vars <- names(df) %>%
    .[!. %in% c("well_", "well_f_", "row_", "col_", "row", "column", "data")]
  df %>%
    select( one_of(layout_vars) ) 
  # %>%
  #   plyr::mutate("-" = rep("", nrow(.)))
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