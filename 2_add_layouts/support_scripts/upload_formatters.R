library(outliers)
library(readr) # to read files into dataframes
library(signal) # contains the savitsky-golay filter, whichs masks poly() in stats
library(stats)
library(tibble)
library(dplyr)
library(purrr)
library(nnet)
library(readxl) # for read_excel for quantStudio

cycle_to_T_func <- function(start_T, inc_T, df) { # returns a numeric vector
  Temps <- rep(start_T, nrow(df))
  steps <- seq.int(from = 0, to = (nrow(df)-1), by = 1)
  steps2 <- steps*inc_T
  Temps2 <- steps2 + Temps
  df <- data.frame("Temperature" = Temps2)
  df
} 

# qTower and quantStudio formats require both reading and reformatting
# qTower
read_qTower <- function(filename) {
  df <- read_csv(filename) %>%
    . [is.na( . [[11]]) == FALSE,] %>%
    t(.) %>%
    as_tibble(.) %>%
    set_names(., c("Temperature", .[1,][-1])) %>%
    . [-1,] %>%
    mutate_if(is.factor, as.character) %>%
    mutate_all(as.numeric) %>%
    round(3)
  df
}

# quantStudio
read_quantStudio <- function(filename) { # completed
  df_raw <- readxl::read_excel(filename, skip=5, sheet = 4) %>%
    . [is.na(.[3]) == FALSE, c(2,4,5)] %>%
    set_names( . , c("Well", "Temperature", "RFU")) %>% 
    .[-1,] %>%
    pivot_wider(names_from = "Well", values_from = RFU) %>%
    mutate_if(is.factor, as.character) %>%
    mutate_all(as.numeric) %>%
    round(3)
}


format_stratagene <- function(df) {
  # determine if the data frame is a single chunk, or multiple
  chunk_num_raw <- table(df[,1]) # will be null if a single chunk
  d_rows <- which(df[,1] == "Dissociation Curve") # determine start rows of multiple chunks, if present
  df <- df[,-1] # remove the first column, needed only to determine chunk presence and locations
  
  if (length(chunk_num_raw) > 0 ) { # if multiple chunks present
    chunk_num <- (chunk_num_raw[[1]] + 1) # determine how many
    df_final <- df[1:d_rows[1],]
    
    for (i in c(1:(length(d_rows)-1))) { # stitch everything together into a single dataframe 
      df_chunk <- df[(d_rows[i]+1):(d_rows[i+1]),]
      df_final <- dplyr::bind_cols(df_final, df_chunk) # append temperature column  
    }
    df <- df_final # over-write input df 
  }
  
  df <- df[rowSums(is.na(df)) < ncol(df),] # remove rows which are all NAs, once present between chunks
  well_names <- as_vector(df[1,]) # extract well names
  well_names <- c("Temperature" , well_names[c(TRUE, FALSE)]) # remove empty columns from well names
  df <- df[-c(1,2),] #  remove leading empty rows
  temperatures <- data.frame( "Temperatures" = df[1]) # extract temperatures
  
  to_delete <- seq(1, ncol(df), 2) # remove duplicated temperature rows
  df <- df[,-to_delete] # remove duplicated temperature rows
  
  df <- dplyr::bind_cols(temperatures, df) # append temperature column
  df <- set_names(df, nm = well_names) # re-set the names
}

format_none <- function(filepath) { 
}

format_biorad <- function(df) {
  df[,-1]
}


#### making the various well names 
make_well_names <- function(row_style, num_style) {
  if (row_style == "rows") { rows <-  letters[1:16] } else {rows <- LETTERS[1:16] }
  if (num_style == "01") { cols <- c(c("01", "02", "03", "04", "05", "06", "07", "08", "09"), c(10:24)) } else {cols <- c(1:24) }
  
  int <-  list(row = rows,
               col =cols) %>%
    cross() %>% # general all possible combos of rows and columns
    map(lift(paste0)) %>% # pate them together
    
    as_vector() %>%
    as_tibble()  %>%
    mutate(nchar = lapply(.$value, nchar) %>% as_vector(),
           char =  str_extract_all(.$value, "[A-Z; a-z]", simplify = TRUE) 
           %>% str_to_upper(locale = "en") 
           %>% as_vector()) 
  
  
  if (num_style == "01") {
    ordered <- int %>%
      dplyr::arrange(value) %>%
      .[[1]]
  } else {
    ordered <- int %>%
      dplyr::arrange(char) %>%
      .[[1]]
  }
  ordered
}

# e.g. use this way
# well_mappings <- tibble(
#   ROWS1 = make_well_names("ROWS", "1"),
#   ROWS01 = make_well_names("ROWS", "01"),
#   rows1 = make_well_names("rows", "1"),
#   rows01 = make_well_names("rows", "01")
# )