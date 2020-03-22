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

parse_well_vec <- function( well_vec ){
  print("parsing well vec")
  l <- list(
    col =  parse_number(well_vec),
    row = str_extract_all(well_vec, "[A-Z; a-z]", simplify = TRUE) 
    %>% str_to_upper(locale = "en") 
    %>% as_vector()
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
  dup_cols <- l_names[!l_names %in% c("well_", "well_f_", "row_", "col_")]
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
      print("all names in common cols ZZ")
    } ### IX THIS--laout needs to have something left to combine
  }
  
    layout_ <- layout_ %>%
               select(-condition) %>%
               unite("condition", -one_of(common_cols), remove = FALSE) %>%
               mutate_if(is.factor, as.character)
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

