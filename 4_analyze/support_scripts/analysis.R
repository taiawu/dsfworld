qT_parse_2 <- function(datafile_with_path, start_T, end_T) {
  df_qt <- read.csv(datafile_with_path, row.names = NULL, skip = 18, header = FALSE, stringsAsFactors = FALSE)  %>%
    set_names( . ,  c("well", c(start_T:end_T)))
  
  channel_rows <- df_qt[[1]] %in% c("FAM", "JOE", "TAMRA", "ROX", "Cy5", "Cy5.5", "Sypro") # determine the rows which contain the channels
  
  channels <- df_qt[channel_rows,] %>% # determine the actual measured channels
    .[[1]]
  
  channel_rep <- c()
  for (i in channels) {
    a <- rep(i, 384)
    channel_rep <- c(channel_rep,a)
  }
  
  df_qt2 <- df_qt[!channel_rows,] %>%
    mutate( . , well_channel = as_vector(map2( . $well,  channel_rep, paste, sep = ":"))) %>% # remove the rows which contain only channels from the dataframe
    .[-1]
  
  df_qt_t <- as_tibble(cbind(Temperature = names( df_qt2 ), t( df_qt2 ))) %>%
    set_names( . , c("Temperature", df_qt2$well_channel))
  
  # list_out <- list(df_qt_2, df_qt_t)
  
  df_qt_t
  df_m <- df_qt_t %>%
    reshape2::melt(id.vars = "Temperature") %>%
    mutate( . , well = as_vector(lapply((strsplit(as.character( . $variable), ':')), function(x) {x[1]} ))) %>%
    mutate( . , channel = as_vector(lapply((strsplit(as.character( . $variable), ':')), function(x) {x[2]} )))
  
  out_list <- list(long = df_m, wide = df_qt_t)
}

# new daughter layout function
df_to_layout <- function(df, layout_type) {
  df_m <-   set_names( df ,  c("type","row",as.numeric( df [1,-c(1,2)]))) %>%
    . [ -1 , -1] %>%
    reshape2::melt( . ,id.vars = "row") %>%
    mutate( . , well = as_vector(map2( . $row,  . $variable, paste0)) ) %>%
    set_names( . , c("row", "column", layout_type, "well"))
  df_m
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

gg_facet_ncol_ng <- function(p){
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



######## determination of Tm by dRFU
##the fitting functions
sgfilt_set_n <- function(n_) { # a closure to set the sg filter length, which will be based on the length of the input data
  function(data, m_) {
    x <- data$value
    out <- sgolayfilt(x, p = 3, n = n_, m = m_)
  }
}

find_sgolay_width <- function(df) { # this works with the true temperature window
  range <- max(df$Temperature) - min(df$Temperature)
  nmeas <- length(unique(df$Temperature))
  out <- floor(3/(range/nmeas)) - floor(3/(range/nmeas))%%2 # ensure it's odd
  
  if (out < 5 ) {out <- 5}
  
  out
}

#sgfilt_nest <- sgfilt_set_n(n_ = find_sgolay_width(df)) # this will ultimately be set by the uploaded data
sgfilt_nest <- sgfilt_set_n(n_ = 13) # this will ultimately be set by the uploaded data

Tm_by_dRFU <- function( data,  sgd1) {
  df <- tibble( x = data$Temperature, 
                y = sgd1)
  
  grid <- tibble( x = seq(min(df$x), max(df$x), by = 0.1) ) %>% modelr::add_predictions(loess(y ~ x, data = df, span = 0.1))
  
  tma <- grid$x[which(grid$pred == max(grid$pred))]
  
  tma # the apparent Tm from the first derivative
}

#######

# model fitting

# load packages
# library(quantmod) # contains the findValleys function, which maybe we should just extract and put verbatim in a source file instead of loading this whole thing...?
# library(minpack.lm) # contains the nlsLM function, which we use for our fitting
# #library(assertive.types) # to figure out how many rows are in a plot
# #library(readr) # to read files into dataframes
# library(modelr)
# library(signal)
# library(tidyverse)
#source("~/Box Sync/papers/paper1_dsfworld/20190929_dsfworld_modeling_update/dsfworld5_uploads_page_alone/support_scripts/dsfworld5_data_analysis.R") # scripts written to analyze and visualize DSF data)

grep_and_gsub <- function(vec_in, g_pattern_vec, rep_TRUE_vec, else_pattern) { # general function, which is likeley present in a different support script, so look at that accordingly
  a <- vec_in # establish initial vector
  for (i in c(1:length( g_pattern_vec ))) {
    b <- grepl( g_pattern_vec[[i]], vec_in )
    a[b] <- rep_TRUE_vec[[i]]
  }
  
  a[ a == vec_in ] <- else_pattern
  
  a
}

#### the fitting functions
sgfilt_set_n <- function(n_) { # a closure to set the sg filter length, which will be based on the length of the input data
  function(data, m_) {
    x <- data$value_norm
    out <- sgolayfilt(x, p = 3, n = n_, m = m_)
  }
}


find_sgolay_width <- function( win3d ) { # this works with the true temperature window
  
  out <- 3 * win3d
  
  if (out < 5 ) {out <- 5}
  
  out
}

make_peak_bool <- function(df, peaks) {
  bool <- rep(FALSE, times  = length(df$Temperature)  ) # This 
  bool[peaks] <- TRUE
  bool
}

make_peak_finder_nest <- function( win3d ) {
  
  function( data, use_loess ) { # assumes a single variable, and the named column of interest  
    df <- data %>% 
      mutate(sgd1 = sgfilt_nest( . , m = 1))  %>%
      mutate(sgd2 = sgfilt_nest( . , m = 2)) 
    
    if (use_loess == TRUE ) {
      sg1 <- loess( df$sgd1 ~ df$Temperature_norm, span = 0.1) %>%
        predict( . )
      sg2 <- loess( df$sgd2 ~ df$Temperature_norm, span = 0.1) %>%
        predict( . )
    } else {
      sg1 <- df$sgd1
      sg2 <- df$sgd2      
    }
    
    peak_loc <- findPeaks(sg1) %>% # find the peaks in the first derivative
      . [sg1[ . ] > 0.0002] %>% # don't take any super tiny peaks
      .[ . > 5] %>% # don't take anything within the first five points
      .[ . < 95] # dont take anything in the last five points 
    
    
    
    valley_loc <- quantmod::findValleys(abs(sg2), thresh = 0.000001) %>%# 0.000001
      . [sg1[ . ] > 0]  %>% # keep only those which do not occur on a [significantly] negative slope
      .[ . > 5] %>% # don't take anything within the first five points
      .[ . < 95] # dont take anything in the last five points
    
    points <- c(peak_loc, valley_loc) # locations of the peaks
    
    
    if (length(points) > 1) { # if there is more than one peak, filter for overly-close peaks
      diffs <- points %>% diff
      points <- points[c(TRUE,(diffs > win3d))] # don't take any points that are within three measurements of one another--maybe fix this to take an n2r argument???
    }
    
    points_final <- points
    
    if (length(points < 2)) {
      points_final <- c(points_final, points_final[[1]] + 1) # if xmid1 and xmid2 start equal, this often causes singular gradient at starting parameters. This makes intuitive sense?
    }
    points_final
  }
  
}

# this is a close which is run once at the beginning of the 
make_temp_n2r <- function( Temps_real) {
  function( Temp_norm ) {
    Temp_norm*(max(Temps_real) - min(Temps_real)) + min(Temps_real)
  }
  
}


s1_d_model <- function(df, par_start) {# runs fine
  xmid1 <- df[par_start,]$Temperature_norm[1]
  if (is.na(xmid1) == TRUE ) {xmid1 <- 0.6} # when wil this be NA, anyway? catch this
  xmid2 <- df[par_start,]$Temperature_norm[2]
  if (is.na(xmid2) == TRUE ) {xmid2 <- 0.62} # when wil this be NA, anyway? catch this
  
  id_d_start <- df$value_norm[1] # the initial value
  
  tryCatch(
    nlsLM(value_norm ~ Asym/(1 + exp((xmid - Temperature_norm)/scal))*exp(d * (Temperature_norm - xmid)) +
            id_d * exp(Temperature_norm* id_b),
          
          data = df,
          start = list(Asym = 1, xmid = xmid1  , scal = 0.03, d = -1,
                       id_d = id_d_start, id_b = -5
          ),
          lower = c(Asym = .1, xmid = 0.1, scal = 0.01, d = -10,
                    id_d = 0.01, id_b = -20
          ),
          control = list(maxiter = 500)
    ),
    warning = function(w) return(NA), error = function(e) return(NA)
  )
}

s1_model <- function(df, par_start) { # runs fine
  xmid1 <- df[par_start,]$Temperature_norm[1]
  if (is.na(xmid1) == TRUE ) {xmid1 <- 0.6} # when wil this be NA, anyway? catch this
  xmid2 <- df[par_start,]$Temperature_norm[2]
  if (is.na(xmid2) == TRUE ) {xmid2 <- 0.62} # when wil this be NA, anyway? catch this
  
  tryCatch(
    nlsLM(value_norm ~ Asym/(1 + exp((xmid - Temperature_norm)/scal))*exp(d * (Temperature_norm - xmid)),
          
          data = df,
          start = list(Asym = 1, xmid = xmid1  , scal = 0.03, d = -1
          ),
          lower = c(Asym = .1, xmid = 0.1, scal = 0.01, d = -10
          ),
          control = list(maxiter = 500)
    ),
    warning = function(w) return(NA), error = function(e) return(NA)
  )
}

s2_d_model <- function(df, par_start) {# runs fine
  xmid1 <- df[par_start,]$Temperature_norm[1]
  if (is.na(xmid1) == TRUE ) {xmid1 <- 0.6} # when wil this be NA, anyway? catch this
  xmid2 <- df[par_start,]$Temperature_norm[2]
  if (is.na(xmid2) == TRUE ) {xmid2 <- 0.62} # when wil this be NA, anyway? catch this
  
  id_d_start <- df$value_norm[1] # the initial value
  
  tryCatch(
    nlsLM(value_norm ~ Asym/(1 + exp((xmid - Temperature_norm)/scal))*exp(d * (Temperature_norm - xmid)) +
            Asym2/(1 + exp((xmid2 - Temperature_norm)/scal2))*exp(d2 * (Temperature_norm - xmid2)) +
            id_d * exp(Temperature_norm* id_b),
          
          data = df,
          start = list(Asym = 1,    xmid = xmid1,  scal = 0.03,  d = -1,
                       Asym2 = 0.5, xmid2 = xmid2, scal2 = 0.03, d2 = -2,
                       id_d = id_d_start, id_b = -5
          ),
          lower = c(Asym = .01, xmid = 0.1, scal = 0.01, d = -10,
                    Asym2 = .01, xmid2 = 0.1, scal2 = 0.01, d2 = -10,
                    id_d = 0.01, id_b = -20
          ),
          
          control = list(maxiter = 500)
    ),
    warning = function(w) return(NA), error = function(e) return(NA)
  )
  
}

s2_model <- function(df, par_start) { # runs fine
  xmid1 <- df[par_start,]$Temperature_norm[1]
  if (is.na(xmid1) == TRUE ) {xmid1 <- 0.6} # when wil this be NA, anyway? catch this
  xmid2 <- df[par_start,]$Temperature_norm[2]
  if (is.na(xmid2) == TRUE ) {xmid2 <- 0.62} # when wil this be NA, anyway? catch this
  
  tryCatch(
    nlsLM(value_norm ~ Asym/(1 + exp((xmid - Temperature_norm)/scal))*exp(d * (Temperature_norm - xmid)) +
            Asym2/(1 + exp((xmid2 - Temperature_norm)/scal2))*exp(d2 * (Temperature_norm - xmid2)),
          
          data = df,
          start = list(Asym = 1,    xmid = xmid1,  scal = 0.03,  d = -1,
                       Asym2 = 0.5, xmid2 = xmid2, scal2 = 0.03, d2 = -2
          ),
          lower = c(Asym = .01, xmid = 0.1, scal = 0.01, d = -10,
                    Asym2 = .01, xmid2 = 0.1, scal2 = 0.01, d2 = -10 
          ),
          
          control = list(maxiter = 500)
    ),
    warning = function(w) return(NA), error = function(e) return(NA)
  )
}

fit_dsf <- function(by_variable, which_model, which_vars) {
  by_variable %>%
    mutate(sgd1 = purrr::map(data, sgfilt_nest, m_ = 1),# fit the model, adding it as a column to the nested dataframe
           peaks = purrr::map(data, peak_finder_nest, use_loess = TRUE), # find the initial peaks
           model = purrr::map2(data, peaks, which_model)
    )
}

handle_fits <- function( df_fit) {
  df_fit %>%# fit the model, adding it as a column to the nested dataframe
    mutate(
      model_pars = purrr::map(model, broom::tidy), # extract model parameters
      xmid_start = map2(data, peaks, make_peak_bool),
      broom_aug = map2(model, data, broom::augment),
      resids = purrr::map2(data, model, add_residuals), # add the residuals to the model, model is the column we created earlier, not the actual model
      predictions = purrr::map2(data, model, add_predictions), # model is the column we created earlier, not the actual model itself
      glance = purrr::map(model, broom::glance)
    ) 
}

# implement models
fit_any_model <- function(by_variable, which_model) {
  # implementation of the fitting functions
  fit_dsf(by_variable, 
          which_model =  which_model, 
          #which_vars = unique(by_variable$variable)) %>% ## changed
          which_vars = unique(by_variable$well)) %>% ## changed
    dplyr::filter(is.na(model) == FALSE) %>%
    
    handle_fits()
  
}

# # s3_BIC <- s3$glance
# extract_model_element <- function(model_df, info, sub_info, pred_name) {
#   model_df %>%
#     #select(variable, info)  %>%
#     select(well, info)  %>%
#     unnest_legacy() %>%
#     #select(variable, sub_info) %>%
#     select(well, sub_info) %>%
#     plyr::mutate(which_model = rep(pred_name, nrow(.))) 
# }

extract_model_element <- function(model_df, info, sub_info, pred_name) {
  model_df %>%
    select(well, condition, info)  %>%
    unnest_legacy() %>%
    select(well, condition, sub_info) %>%
    plyr::mutate(which_model = rep(pred_name, nrow(.))) 
}

## restructure this to break up the determination of starting parameters, which happens only once, and the fitting of the model
# first, get the starting parameter for all models
get_start_pars <- function(by_variable) {
  by_variable %>%
    mutate(
      sgd1 = purrr::map(data, sgfilt_nest, m_ = 1),# fit the model, adding it as a column to the nested dataframe
      peaks = purrr::map(data, peak_finder_nest, use_loess = TRUE)
      # , # find the initial peaks
      # model = purrr::map2(data, peaks, which_model)
    )
}

fit_model_from_pars <- function(by_variable, which_model) { #### THIS NEEDS TO GET by_variable AFTER THE START PARS HAVE BEEN GENERATED.....  this should have a different name; be a different object?
  by_variable <-   by_variable %>%
                    mutate(model = purrr::map2(data, peaks, which_model)) %>%
                    dplyr::filter(is.na(model) == FALSE) %>% # retain only rows which 
                    handle_fits()
  
  by_variable
}

# build individual component predictions from the models
build_sig1 <- function(  data, model_pars, which_model) {
  
  pars <- as.list(model_pars$estimate)
  names(pars) <- model_pars$term
  
  Temperature_norm <- data %>%
    select(Temperature_norm) %>%
    unique() %>%
    arrange()
  Temperature <- data %>%
    select(Temperature) %>%
    unique() %>%
    arrange()
  
  pred_s <- pars$Asym/(1 + exp((pars$xmid - Temperature_norm)/pars$scal))*exp(pars$d * (Temperature_norm - pars$xmid))
  
  tibble( Temperature = as_vector(Temperature),
          Temperature_norm = as_vector(Temperature_norm),
          pred = as_vector(pred_s),
          value_norm = rep(NA, length(pred_s)),
          which_model = rep(which_model, length(pred_s)),
          component = rep("sigmoid_1", length(pred_s))
  )
}

build_sig2 <- function(  data, model_pars, which_model) {
  print("building sig2")
  pars <- as.list(model_pars$estimate)
  names(pars) <- model_pars$term
  
  Temperature_norm <- data %>%
    select(Temperature_norm) %>%
    unique() %>%
    arrange()
  
  Temperature <- data %>%
    select(Temperature) %>%
    unique() %>%
    arrange()
  
  pred_s <- pars$Asym2/(1 + exp((pars$xmid2 - Temperature_norm)/pars$scal2))*exp(pars$d2 * (Temperature_norm - pars$xmid2))
  
  tibble( Temperature = as_vector(Temperature),
          Temperature_norm = as_vector(Temperature_norm),
          pred = as_vector(pred_s),
          value_norm = rep(NA, length(pred_s)),
          which_model = rep(which_model, length(pred_s)),
          component = rep("sigmoid_2", length(pred_s))
  )
}

build_d <- function(  data, model_pars, which_model) {
  pars <- as.list(model_pars$estimate)
  names(pars) <- model_pars$term
  
  Temperature_norm <- data %>%
    select(Temperature_norm) %>%
    unique() %>%
    arrange()
  
  Temperature <- data %>%
    select(Temperature) %>%
    unique() %>%
    arrange()
  
  pred_s <- pars$id_d * exp(Temperature_norm* pars$id_b)
  
  tibble( Temperature = as_vector(Temperature),
          Temperature_norm = as_vector(Temperature_norm),
          pred = as_vector(pred_s),
          value_norm = rep(NA, length(pred_s)),
          which_model = rep(which_model, length(pred_s)),
          component = rep("initial_decay", length(pred_s))
  )
}

build_s1_predictions <- function( df ) {
  which_model <- "s1_pred"
  
  s1 <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_sig1, which_model)) %>%
    unnest_legacy(pred) 
  
  s1
}

build_s1_d_predictions <- function( df ) {
  
  which_model <- "s1_d_pred"
  
  s1 <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_sig1, which_model)) %>%
    unnest_legacy(pred) 
  
  
  sd <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_d, which_model)) %>%
    unnest_legacy(pred)
  
  rbind(s1, sd)
}

build_s2_predictions <- function( df) {
  which_model <- "s2_pred"
  
  s1 <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_sig1, which_model)) %>%
    unnest_legacy(pred) 
  
  s2 <-  df %>%
    mutate(sig1_pred = purrr::map2(data, model_pars, build_sig2, which_model)) %>%
    unnest_legacy(sig1_pred)
  
  rbind(s1, s2)
}

build_s2_d_predictions <- function( df ) {
  which_model <- "s2_d_pred"
  s1 <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_sig1, which_model)) %>%
    unnest_legacy(pred)
  
  s2 <-  df %>%
    mutate(sig1_pred = purrr::map2(data, model_pars, build_sig2, which_model)) %>%
    unnest_legacy(sig1_pred)
  
  
  sd <- df %>%
    mutate(pred = purrr::map2(data, model_pars, build_d, which_model)) %>%
    unnest_legacy(pred)
  
  rbind(s1, s2, sd)
}

# extract_preds <- function( df, pred_name ) {
#   df %>% 
#     unnest_legacy(predictions) %>% 
#     #select( variable, Temperature, Temperature_norm, pred, value_norm ) %>% 
#     select( well, Temperature, Temperature_norm, pred, value_norm ) %>% 
#     #mutate(pred_model = rep(pred_name, nrow(.))) %>% #%>% rename(s1_pred = pred)
#     plyr::mutate(which_model = rep(pred_name, nrow(.))) %>% #%>% rename(s1_pred = pred)
#     plyr::mutate(component = rep("full_pred", nrow(.))) 
# }
# 
# build_predictions <- function(model_name, df) {
#   if (model_name == "s1_pred")   { out <- build_s1_predictions(df)
#   } else if (model_name == "s1_d_pred") { out <- build_s1_d_predictions(df)
#   } else if (model_name == "s2_pred")   { out <- build_s2_predictions(df)
#   } else if (model_name == "s2_d_pred") { out <- build_s2_d_predictions(df) 
#   }
#   out
# }

extract_preds <- function( df, pred_name ) {
  df %>% 
    unnest_legacy(predictions) %>% 
    select( well, Temperature, Temperature_norm, pred, value_norm , condition) %>% 
    plyr::mutate(which_model = rep(pred_name, nrow(.))) %>% 
    plyr::mutate(component = rep("full_pred", nrow(.))) 
}

build_predictions <- function(model_name, df) {
  if (model_name == "s1_pred")   { out <- build_s1_predictions(df)
  } else if (model_name == "s1_d_pred") { out <- build_s1_d_predictions(df)
  } else if (model_name == "s2_pred")   { out <- build_s2_predictions(df)
  } else if (model_name == "s2_d_pred") { out <- build_s2_d_predictions(df) 
  }
  out 
  # %>%
  #   select( well, Temperature, Temperature_norm, pred, value_norm , condition, component) 
}

readable_model_names <- function( model ) {
  out <- model
  for (i in 1:length(out)) {
    if ( out[i] == "s1_d_pred") { out[i] <- "One sigmoid, with initial fluor"
    } else if ( out[i] == "s1_pred") { out[i] <- "One sigmoid, no initial fluor"
    } else if ( out[i] == "s2_d_pred") { out[i] <- "Two sigmoids, with initial fluor"
    } else if ( out[i] == "s2_pred") { out[i] <- "Two sigmoids, no initial fluor"
    }
  }
  out     
}

dsfworld_default_model <- theme( # adapted from free amino acids hit call
  text = element_text(size = 10),
  # axis.title.x = element_blank(),
  # axis.title.y = element_blank(),
  axis.text = element_text(size = 8),
  axis.text.x = element_text(angle = 45, vjust = 0.5),
  legend.position = "right",
  plot.title = element_text(lineheight=.8, face="bold", size = 12),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(), 
  strip.background = element_blank()
  # ,
  # aspect.ratio = (1/1.618)
  # axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
  # axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid')
)

make_model_plot_sel <- function(df_models, df_BIC, fit_vars, show_models, df_min_BIC) {
  # these are the same every time, so it would be most efficient to do them outside of this function, but I think the time added is worth the simplicity in how things are coded in the app itself...
  facet_labels <- c("1 sigmoid \n+with initial fluor.", "1 sigmoid \n no initial fluor.", "2 sigmoids \n with initial fluor.", "2 sigmoids\n  no initial fluor." )
  names(facet_labels) <- c("s1_d_pred", "s1_pred","s2_d_pred", "s2_pred")
  
  square <- with(df_models, {
    tibble( x = c(min(Temperature), max(Temperature), max(Temperature), min(Temperature)),
            y = c(0, 0, 1.5, 1.5))
  })
  
  selected_model <- tibble(well = unique(df_models$well),
                           s1_pred= rep(FALSE, length(unique(df_models$well))),
                           s1_d_pred = rep(FALSE, length(unique(df_models$well))),
                           s2_pred = rep(FALSE, length(unique(df_models$well))),
                           s2_d_pred = rep(FALSE, length(unique(df_models$well)))
  )
  
  for ( i_var in df_min_BIC$well) {
    model_choice <- df_min_BIC[df_min_BIC$well == i_var,]$which_model
    selected_model[selected_model$well == i_var,c(2:ncol(selected_model))] <- c(model_choice == names(selected_model[-1]))
  }
  
  # a boolean dataframe showing which model is selected
  df_min_BIC_for_plot <- selected_model %>%
    gather(which_model, selected, -well) %>%
    dplyr::filter(well %in% fit_vars)  %>%
    dplyr::filter(which_model %in% show_models)
  
  # a dataframe with all of the raw and modeled data, and fit components
  df_for_plot <- df_models %>%
    dplyr::filter(well %in% fit_vars)  %>%
    dplyr::filter(which_model %in% show_models)
  
  #df_for_plot
  
  # a dataframe with the actual BIC values
  BIC_for_plot <- df_BIC %>%
    dplyr::filter(well %in% fit_vars) %>%
    dplyr::filter(which_model %in% show_models)
  
  ggplot(df_for_plot, aes(x = Temperature, y = pred)) +
    geom_polygon(aes(x = x, y = y, fill = selected),  data = merge(df_min_BIC_for_plot, square), alpha = 0.1, show.legend = FALSE) +
    scale_fill_manual( values = c(NA, "black") ) +
    geom_text(data = BIC_for_plot,
              aes(x = ((max(df_for_plot$Temperature) - min(df_for_plot$Temperature))/2)+ min(df_for_plot$Temperature), # center the label in the plot region
                  y = 1.3, label = paste0("BIC ", round(BIC, 0)),
                  group= NULL, color = BIC),
              size = 4, show.legend = FALSE) +
    
    geom_point(aes(Temperature, value_norm, shape = "raw data"), size = 1, alpha = 0.15) +
    geom_point(aes(Temperature, value_norm), size = 1, alpha = 0.15) +
    
    geom_line(aes(linetype = component), color = "black") +
    facet_grid(well~which_model,  labeller = labeller(which_model = facet_labels)) +
    scale_linetype_manual(values=c("full_pred" = "solid", "sigmoid_1" = "dashed", "sigmoid_2" = "dotdash", "initial_decay" =  "dotted"),
                          name="Fits",
                          breaks = c("full_pred", "sigmoid_1", "sigmoid_2", "initial_decay"),
                          labels=c("final fit", "simgoid 1", "sigmoid 2", "initial fluor.")) +
    labs(title = "Comparison of models", y = "Normalized RFU", x = "Temperature (C)") +
    labs( shape = "") +
    theme_bw() +
    dsfworld_default_model +
    ylim( c(0, 1.5) ) +
    scale_color_viridis_c()
}

## Tma from model fits
# function to calculate the Tm by dRFU for the models--different only in the input formats for the function. ideally would be merged with the dRFU Tma function.
Tm_by_dRFU_for_model <- function( data ) {
  df <- tibble( x = data$Temperature, 
                y = data$sgd1_pred)
  
  grid <- tibble( x = seq(min(df$x), max(df$x), by = 0.1) ) %>% modelr::add_predictions(loess(y ~ x, data = df, span = 0.1))
  
  tma <- grid$x[which(grid$pred == max(grid$pred))]
  
  tma # the apparent Tm from the first derivative
}

mean_no_NaN <- function( vec ) {
  a <- mean(vec, na.rm = TRUE)
  
  if( is.nan(a) == TRUE) {
    a <- NA
  }
  a
}

get_condition <- function(data) {
  data$condition %>% unique()
}

add_sig2_col <- function( data ) {
  if ("sigmoid_2" %in% names(data)) {
    data
  } else {
    data %>%
      mutate(sigmoid_2 = map(sigmoid_1, function(x) {NA}) %>% as_vector())
     # mutate(sigmoid_2 = map(.$sigmoid_1, function(x) {NA}) %>% as_vector()) # if the data were not grouped, but it is
  }
}


model_tms_by_dRFU <- function( df_models, win3d ) {
  print("making model tms by drFU, calling sigmoid_2")
  df_pred_sgd1 <- df_models  %>% # the model outcomes, with sgd1 added for all predictions within the nested data
    group_by(well, which_model, component)  %>%
    mutate(sgd1_pred = sgolayfilt(pred, p = 3, n = find_sgolay_width( win3d ), m = 1)) %>%
    nest() %>%
    mutate(tma = map(data, Tm_by_dRFU_for_model) %>% as_vector())
  
  df_tma <- df_pred_sgd1 %>% # all calculated temperatures, with the rest of the data now removed
    mutate(tma = map(data, Tm_by_dRFU_for_model) %>% as_vector()) %>%
    ungroup() %>%
    mutate(condition = map(data, get_condition) %>% as_vector()) %>%
    select(-data) %>%
    dplyr::filter(component %in% c("sigmoid_1", "sigmoid_2")) %>%
    pivot_wider(names_from = component, values_from  = tma)  %>%
    group_by(condition, which_model) %>%
    add_sig2_col() %>% # if the sigmoid_2 column is missing, add in an empty one for it 
    mutate(mean_tma1 = mean_no_NaN(sigmoid_1) %>% round(1),
           mean_tma2 = mean_no_NaN(sigmoid_2) %>% round(1),
           sd_tma1 = sd(sigmoid_1) %>% round(1),
           sd_tma2 = sd(sigmoid_2) %>% round(1)) 
  
  df_tma_mean <- df_tma %>% # the mean tmas only, for downloading and displaying
    distinct(mean_tma1, mean_tma2, .keep_all = TRUE) %>%
    ungroup() %>%
    dplyr::select(c(condition, which_model, mean_tma1, sd_tma1, mean_tma2, sd_tma2))
  
  out <- list("df_pred_sgd1" = df_pred_sgd1,
              "df_tma" = df_tma,
              "df_tma_mean" = df_tma_mean)
}


# get_Tms <- function(model_df, model) {
#   pars <- model_df %>%
#     unnest_legacy(model_pars) %>%
#     dplyr::filter(grepl('xmid', term)) %>%
#     mutate( tma = map(estimate, n2r) %>% 
#               as_vector() %>%
#               round( digits = 1) ) %>%
#     select(well, condition, term, tma) %>% # 20200310
#     #select(well,  term, tma) %>%
#     plyr::mutate(which_model = rep(model, nrow(.)))
#   #pivot_wider(names_from = term, values_from = tma)
# }

# make_model_tm_table <- function( model_tms) {
#   # for the s1 models, add a dummy xmid2 column. clunkly, i know, but a simple fix
#   if (!"xmid2" %in% names(model_tms)) { model_tms <- model_tms %>%  plyr::mutate(xmid2 = rep(NA, nrow(.)))}
#   
#   if ("condition" %in% names(model_tms)) {
#     # make the averaged table
#     tm_table_models <- model_tms %>%
#       select(condition, which_model, xmid, xmid2) %>%
#       group_by(condition, which_model)   %>%
#       summarise( mean_xmid1 = mean(xmid) ,
#                  sd_xmid1 = sd(xmid),
#                  mean_xmid2 = mean(xmid2) ,
#                  sd_xmid2 = sd(xmid2)) %>%
#       mutate_if(is.numeric, round, 2) %>%
#       ungroup()
#     
#   } else {
#     tm_table_models <- model_tms %>%
#       select(well, which_model, xmid, xmid2) %>%
#       group_by(well, which_model)   %>%
#       summarise( mean_xmid1 = mean(xmid) ,
#                  sd_xmid1 = sd(xmid),
#                  mean_xmid2 = mean(xmid2) ,
#                  sd_xmid2 = sd(xmid2)) %>%
#       mutate_if(is.numeric, round, 2) %>%
#       ungroup()
#   }
#   
#   tm_table_models
# }

model_all <- function(which_model, model_name, start_pars, win3d) {
  
  model_fit <- start_pars %>% fit_model_from_pars(which_model = which_model)
  
  df_BIC <- extract_model_element(model_fit, "glance", "BIC", model_name) # these will be created for s1, the default fit, and the fastest one. they will be added to for each additional model

  df_models <- bind_rows( extract_preds(model_fit, model_name), # these will be created for s1, the default fit, and the fastest one. they will be added to for each additional model
                      build_predictions(model_name, model_fit) )

  tm_table_models <- model_tms_by_dRFU( df_models, win3d )

  out_list <- list("model" = model_fit, 
                   "df_BIC" = df_BIC, 
                   "df_models" = df_models, 
                   "tm_table_models" = tm_table_models$df_tma_mean,
                   "tm_models_all" = tm_table_models$df_tma,
                   "df_models_sgd1" = tm_table_models$df_pred_sgd1
                   )
}


########### plots 
facet_func_with_model <- function(df, 
                                  #mean_or_each, # the models are fit to the individual traces, not the means, so this isn't an option here
                                  color_by, 
                                  #linetype_by,  # not an option, since the models will be linetyped automatically
                                  facet, 
                                  facet_by, 
                                  facet_rows, 
                                  facet_cols, 
                                  set_title, 
                                  legend_title, 
                                  legend_linetype_title, 
                                  #use_linetypes, # not an option, since the models will be linetyped automatically
                                  fix_free, 
                                  text_size, 
                                  legend_position, 
                                  x_title, 
                                  y_title,
                                  
                                  show_components, # a boolean, whether to add all model components, or just the overall fit
                                  color_models,
                                  tm_df_gathered,
                                  show_tm
) {
  
  color_by <- enquo(color_by) # tell how to color
  if (show_components == FALSE) { df <- df %>% filter(component == "full_pred")} # show only the full prediction
  
  df <- df %>%
    plyr::mutate(well_component  =  map2(.$well, .$component, paste, sep = "-") %>% as_vector() ) # a grouping variable for the components
  
  
  p <- ggplot(df, aes(x = Temperature,  group = well, color = !!color_by))
  
  if (color_models == TRUE) {
    p <- p + geom_line(aes(y = pred, linetype = component, group = well_component))
  } else if (color_models == FALSE) {
    p <- p + geom_line(aes(y = pred, linetype = component, group = well_component), color = "grey")
  }
  
  p <- p + geom_point(aes(y = value_norm, shape = "raw data"), size = 1, alpha = 0.5) +# add the raw data points
    scale_linetype_manual(values=c("full_pred" = "solid", "sigmoid_1" = "dashed", "sigmoid_2" = "dotdash", "initial_decay" =  "dotted"), # unique to the modeling-overlay plot
                          name="Fits",
                          breaks = c("full_pred", "sigmoid_1", "sigmoid_2", "initial_decay"),
                          labels=c("final fit", "simgoid 1", "sigmoid 2", "initial fluor."))
  
  if (show_tm == TRUE) {
    p <- p + geom_vline(data = tm_df_gathered, aes(xintercept = Tm))
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
  
  p <- p +
    theme_bw() + # here to end, identical to the typical DSFworld plot
    scale_color_viridis_d() +
    labs(title = set_title, color = legend_title, x = x_title, y = y_title, shape = "") +
    dsfworld_default_model +
    theme(  text = element_text(size = text_size*1.25),
            axis.text = element_text(size = text_size),
            plot.title = element_text(lineheight=.8, face="bold", size = text_size*1.5),
            legend.position = legend_position) 
  
  p
}

facet_func2 <- function(df, 
                        mean_or_each, 
                        color_by, 
                        linetype_by, 
                        facet, 
                        facet_by, 
                        facet_rows, 
                        facet_cols, set_title, 
                        legend_title, 
                        legend_linetype_title, 
                        use_linetypes, 
                        fix_free, 
                        text_size, 
                        legend_position, 
                        x_title, 
                        y_title,
                        
                        tm_df_gathered,
                        show_tm
) {
  
  color_by <- enquo(color_by) # tell how to color
  linetype_by <- enquo(linetype_by)
  
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
      dsfworld_default_model +
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
      dsfworld_default_model +
      theme(  text = element_text(size = text_size*1.25),
              axis.text = element_text(size = text_size),
              plot.title = element_text(lineheight=.8, face="bold", size = text_size*1.5),
              legend.position = legend_position) 
    
  }
  
  if (show_tm == TRUE) {
    p <- p + geom_vline(data = tm_df_gathered, aes(xintercept = Tm))
  }
  
  p
}
