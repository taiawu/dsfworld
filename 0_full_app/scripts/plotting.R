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

## for model fitting plots
# this code to be added to plotting.R
cond_df_model_for_plot <- function( model_df, df_BIC ) {
  model_df %>%
    group_by(well) %>%
    nest() %>%
    ungroup() %>%
    full_join(df_BIC) %>%
    unnest()  %>%
    select(well, condition, Temperature, Temperature_norm, pred, value_norm, which_model, component, BIC) %>%
    pivot_longer(-c(well, Temperature, Temperature_norm, which_model, component, BIC, condition), names_to = "which_value", values_to = "value") %>%
    mutate(which_model_f = factor(.$which_model, levels = c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred"))) 
}

cond_df_BIC_for_plot <- function( df_BIC ) {
  df_BIC %>%
    group_by(well) %>%
    nest() %>%
    ungroup() %>%
    unnest() %>%
    group_by(well) %>%
    mutate(is_min = (BIC == min(BIC))) %>%
    ungroup()%>%
    mutate(which_model_f = factor(.$which_model, levels = c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred"))) 
}

make_model_names <- function( model_names_present ) {
  mods <- tibble(
    model_names = c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred"),
    model_names_f = factor(c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred"), levels = c("s1_pred", "s1_d_pred", "s2_pred", "s2_d_pred")),
    human_read =  c("Fit 1", "Fit 2", "Fit 3", "Fit 4"),
    human_read_f = factor(c("Fit 1", "Fit 2", "Fit 3", "Fit 4"), levels = c("Fit 1", "Fit 2", "Fit 3", "Fit 4"))
  )
  
  mods_present <- mods %>%
    filter(model_names %in% model_names_present) 
  name_vec_f <- mods_present$human_read_f
  names(name_vec_f) <- mods_present$model_names_f
  
  name_vec <- mods_present$human_read
  names(name_vec) <- mods_present$model_names
  
  list("factor" = name_vec_f, "character" = name_vec)
}

match_well_to_cond <- function( df_models ) {
  df <- df_models %>% distinct(well, .keep_all = TRUE)
  wells <- map2(df$well, df$condition, paste, sep = "\n") %>% as_vector()
  names(wells) <- df$well
  wells
}

dsfworld_default_model <- theme( # adapted from free amino acids hit call
  text = element_text(size = 10*1.25),
  axis.text = element_text(size = 8*1.25),
  axis.text.x = element_text(angle = 0, vjust = 0.5),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 12*1.25),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(), 
  strip.background = element_blank(),
  strip.text.y = element_text(angle = 0),
  aspect.ratio = (1/1.618)
)


plot_all_fits_shiny <- function(df_models, df_BIC) { # takes pre-conditioned dataframes 
  
  # create vectors used to label and position data within the plots
  model_names <- make_model_names( df_models$which_model %>% unique() )
  well_names <- match_well_to_cond( df_models )
  mid_temp <- (max(df_models$Temperature) - min(df_models$Temperature))/2 + min(df_models$Temperature)
  
  
  # create the plot
  df_models %>%
    ggplot() +
    geom_line (aes(x = Temperature, # RFU data, fitted and "real"
                   y = value, 
                   linetype = which_value, 
                   color = component, 
                   group = interaction(which_model, component, which_value)
    ), 
    
    alpha = 0.5) +
    geom_text (data = df_BIC, aes(label = paste0("BIC ", round(BIC, 0)), # BIC values for fits
                                  x = mid_temp,
                                  y = 1.3,
                                  group = well,
                                  alpha = is_min),
               size = 3) +
    facet_grid (well~which_model_f,
                labeller = labeller(well = well_names,
                                    which_model_f = model_names$character)
    ) +
    
    scale_alpha_manual (values = c(0.7, 1), guide = FALSE) +
    scale_linetype_manual (values = c("pred" = "solid", "value_norm" = "dashed"),
                           name="",
                           breaks=c( "value_norm", "pred"),
                           labels=c("Data", "Fits")
    ) +
    scale_color_manual (values = c("full_pred" = "#081d58", "initial_decay" = "#edf8b1", "sigmoid_1" = "#253494", "sigmoid_2" = "#41b6c4"),
                        name="Fit component",
                        breaks=c("full_pred", "sigmoid_1", "sigmoid_2", "initial_decay"),
                        labels=c("Final fit", "Sigmoid 1", "Sigmoid 2", "Initial RFU")
    ) +
    guides  (linetype = guide_legend(order = 1),
             color = guide_legend(order = 2))+
    
    ylim (c(0, 1.4))+
    theme_bw() +
    dsfworld_default_model +
    labs (x = "Temperature (ºC)", y = "Normalized RFU")
}

plot_best_fits_shiny <- function(df_models, df_best_in) {
  # create vectors used to label and position data within the plots
  well_names <- match_well_to_cond( df_models )
  mid_temp <- (max(df_models$Temperature) - min(df_models$Temperature))/2 + min(df_models$Temperature)
  
  df_best <- df_best_in %>%
    mutate(best_model = which_model) %>%
    mutate(best_model_human = recode(best_model,
                                     s1_pred = "Fit 1",
                                     s1_d_pred = "Fit 2",
                                     s2_pred = "Fit 3",
                                     s2_d_pred = "Fit 4")) 
  
  df_models_best <- df_models %>% 
    full_join(df_best)  %>%
    dplyr::filter(which_model == best_model) %>%
    mutate(best_model = recode(best_model,
                               s1_pred = "Fit 1",
                               s1_d_pred = "Fit 2",
                               s2_pred = "Fit 3",
                               s2_d_pred = "Fit 4")) 
  df_models_best %>%
    ggplot() +
    geom_line (aes(x = Temperature, # RFU data, fitted and "real"
                   y = value, 
                   linetype = which_value, 
                   color = component, 
                   group = interaction(which_model, component, which_value)
    ), 
    alpha = 0.5) +
    geom_text (data = df_best, 
               aes(label = paste0("Selected: ", best_model_human), # BIC values for fits
                   x = mid_temp,
                   y = 1.3,
                   group = well),
               size = 3) +
    facet_wrap(~well,  #, 
               labeller = labeller(well = well_names)
    ) +
    
    # scale_alpha_manual (values = c(0.7, 1), guide = FALSE) +
    scale_linetype_manual (values = c("pred" = "solid", "value_norm" = "dashed"),
                           name="",
                           breaks=c( "value_norm", "pred"),
                           labels=c("Data", "Fits")
    ) +
    scale_color_manual (values = c("full_pred" = "#081d58", "initial_decay" = "#edf8b1", "sigmoid_1" = "#253494", "sigmoid_2" = "#41b6c4"),
                        name="Fit component",
                        breaks=c("full_pred", "sigmoid_1", "sigmoid_2", "initial_decay"),
                        labels=c("Final fit", "Sigmoid 1", "Sigmoid 2", "Initial RFU")
    ) +
    guides  (linetype = guide_legend(order = 1),
             color = guide_legend(order = 2))+
    
    ylim (c(0, 1.4))+
    theme_bw() +
    dsfworld_default_model +
    labs (x = "Temperature (ºC)", y = "Normalized RFU")
}