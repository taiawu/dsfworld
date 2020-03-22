# what will ultimately be a sources script
# the purely thermodynamic model
model_K <- function(dHu, T_half, dCp, Temp) {
  dG <- dHu*((T_half - Temp) / T_half) - dCp * (T_half - Temp * (1 - ln(Temp/T_half))) # return this value, the dG at a given temperature
  K <- exp(- (dG / (0.008314 * Temp)))
  #frac_unfold <- (K/(1 + K))
}

model_thermo_unfold <- function(K) {
  (K/(1 + K)) # fraction unfolded at any given time
}

#  Lumry-Eyring models
model_LE <- function(Ea, Temp, T_star, dHu, T_half, v, lowT, dCp, sel_output) {
  R <- 0.008314 # kj/(mol*K)
  
  # calculate the mole fractions of native, unfolded, and final state 
  integrand <- function(Temp) {
    exp(((-Ea/0.008314)*(1/Temp-1/T_star)))*  
      #exp(((-dHu/R)*(1/Temp-1/T_half)))/exp(((-dHu/R)*(1/Temp-1/T_half)+1)) # without dCp
      (exp( dHu*((T_half - Temp) / T_half) - dCp * (T_half - Temp * (1 - ln(Temp/T_half)))) /  # e^(K), recalculated here, which is inefficient, since we have a column for K
         exp(dHu*((T_half - Temp) / T_half) - dCp * (T_half - Temp * (1 - ln(Temp/T_half))) + 1 )) #e^(K + 1), recalculated here, which is inefficient, since we have a column for K
  }
  
  integ <- stats::integrate(integrand, lower = lowT, upper = Temp)
  L <- exp((-1/v)*as.numeric(integ$value))
  
  #K <- calc_K(dHu, Temp, T_half) # thermodynamic
  K <- model_K(dHu, T_half, dCp, Temp)  
  
  xF <- 1 - L # fraction final 
  xU <- (K / (K + 1)) * L # fraction unfolded
  xN <- (1 / (K + 1)) * L # fraction native
  
  output_list <- list(xN, xU, xF)
  
  output_list[[sel_output]]  # return one of the following outputs 
}

# dye detection models, with decay
RFU_detec <- function(detec , decay_rate, nrow_df ) {
  decay <- as.numeric(lapply(c(1:nrow_df), function(x) { 1 - decay_rate*(x/nrow_df) } )) 
  detec_out <- decay * detec
  detec_out
}

# create the full model
make_model_df <- function(start_T, end_T, dHu_, T_half_, dCp_, Ea_, T_star_, v_, nat_dye, unf_dye, fin_dye, decay_rate) {
  
  model_df <- tibble("Temperature" = seq(from = start_T, to = end_T, by = 1)) %>%
    mutate( K = model_K(dHu = dHu_, T_half= T_half_, dCp = dCp_, Temp = .$Temperature)) %>% # calculate K
    mutate( thermo_unfold = model_thermo_unfold(K = .$K)) %>% # thermodynamic model of fraction unfolded
    mutate( thermo_fold = 1 - thermo_unfold) %>% # thermodynamic model for fraction folded
    
    mutate( kin_final = as_vector(Map(model_LE, Temp = .$Temperature, Ea = Ea_, T_star = T_star_, dHu = dHu_, T_half = T_half_, v = v_, lowT = 273, dCp = dCp_, sel_output = 3)))  %>%
    mutate( kin_native = as_vector(Map(model_LE, Temp = .$Temperature, Ea = Ea_, T_star = T_star_, dHu = dHu_, T_half = T_half_, v = v_, lowT = 273, dCp = dCp_, sel_output = 1)))  %>%
    mutate( kin_unfolded = as_vector(Map(model_LE, Temp = .$Temperature, Ea = Ea_, T_star = T_star_, dHu = dHu_, T_half = T_half_, v = v_, lowT = 273, dCp = dCp_, sel_output = 2)))  %>%
    
    mutate( RFU_native = RFU_detec(nat_dye, decay_rate, nrow(.)) * .$kin_native)  %>% # works with either a vector or a double
    mutate( RFU_unfolded = RFU_detec(unf_dye,decay_rate, nrow(.)) * .$kin_unfolded)  %>% # works with either a vector or a double
    mutate( RFU_final = RFU_detec(fin_dye,decay_rate, nrow(.)) * .$kin_final)  %>% # works with either a vector or a double
    mutate( RFU  = RFU_native + RFU_unfolded + RFU_final ) %>%
    mutate( dRFU = sgolayfilt( .$RFU, p = 3, n = 5, m = 1))
  
  model_df
}

calc_Tm <- function(model_df) { # not used
  model_df %>%
    dplyr::filter(dRFU == max(model_df$dRFU)) %>%
    .$Temperature %>%
    as.numeric
}


grep_and_gsub <- function(vec_in, g_pattern_vec, rep_TRUE_vec, else_pattern) { # general
  a <- vec_in # establish initial vector
  for (i in c(1:length( g_pattern_vec ))) {
    b <- grepl( g_pattern_vec[[i]], vec_in )
    a[b] <- rep_TRUE_vec[[i]]
  }
  
  a[ a == vec_in ] <- else_pattern
  
  a
}

raster_theme_1 <- theme( # adapted from free amino acids hit call
  text = element_text(size = 16),
  axis.text = element_text(size = 16),
  axis.text.x = element_text(angle = 0, vjust = -0.5),
  #axis.text.x = element_blank(),
  legend.position = "none",
  plot.title = element_text(lineheight=.8, face="bold", size = 12),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  #panel.border = element_blank(),
  strip.background = element_blank(),
  panel.spacing = unit(0, "mm"),
  strip.text.y = element_text(size = 16)
) 

# make_model_plot <- function(df) {
#   df %>%
#     select(Temperature, thermo_unfold, thermo_fold, kin_final, kin_native, kin_unfolded, RFU) %>%
#     gather(variable, value, -Temperature) %>%
#     mutate(model = grep_and_gsub( vec_in = .$variable, g_pattern_vec = c("thermo", "kin", "RFU"), rep_TRUE_vec = c("thermo", "kin", "kin"), else_pattern = "other")) %>%
#     #mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( expression("Thermodynamic \n model"), expression("Thermo-kinetic \n model"),"DSF")) ) %>%
#     mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( "Thermodynamic" ,"Thermodynamic-kinetic ","DSF")) ) %>%
#     mutate(variable_f = factor(.$variable, levels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_unfolded", "kin_final", "RFU"), labels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_unfolded", "kin_final", "RFU"))) %>%
#     
#     ggplot(aes(x = Temperature - 273, y = value, group = variable_f, fill = variable_f, color = variable_f)) +
#     geom_ribbon(aes(ymin = 0, ymax = value), color = NA, alpha = 0.85) +
#     geom_point(size = 1) +
#     scale_color_manual(values = c(thermo_fold = NA, kin_native =  NA, RFU = "black", thermo_unfold = NA, kin_unfolded = NA, kin_final = NA)) +
#     scale_fill_manual( values = c(thermo_fold = "#d8d8d8", kin_native =  "#d8d8d8", RFU = NA, thermo_unfold = "#b2182b", kin_unfolded = "#b2182b", kin_final = "#2166ac")) +
#     facet_grid(rows = vars(model_f)) +
#     theme_bw()+
#     scale_x_continuous(expand = c(0,0))+
#     scale_y_continuous(expand = c(0,0.1), breaks = c(0, 0.5, 1))+
#     labs(x = expression('Temperature ('*~degree*C*')'), y = "Relative population" ) +
#     raster_theme_1 -> p
#   p # return the plot object
# }

# make_model_plot <- function(df) {
#   df %>%
#     select(Temperature, thermo_unfold, thermo_fold, kin_native, kin_final,kin_unfolded, RFU) %>%
#     gather(variable, value, -Temperature) %>%
#     mutate(model = grep_and_gsub( vec_in = .$variable, g_pattern_vec = c("thermo", "kin", "RFU"), rep_TRUE_vec = c("thermo", "kin", "kin"), else_pattern = "other")) %>%
#     #mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( expression("Thermodynamic \n model"), expression("Thermo-kinetic \n model"),"DSF")) ) %>%
#     mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( "Thermodynamic" ,"Thermodynamic-kinetic ","DSF")) ) %>%
#     mutate(variable_f = factor(.$variable, levels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_unfolded", "kin_final", "RFU"), labels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_unfolded", "kin_final", "RFU"))) %>%
#     
#     ggplot(aes(x = Temperature - 273, y = value, group = variable_f, fill = variable_f, color = variable_f)) +
#     geom_ribbon(aes(ymin = 0, ymax = value), color = NA, alpha = 0.85) +
#     geom_line(size = 1) +
#     scale_color_manual(values = c(thermo_fold = NA, kin_native =  NA, RFU = "orange", thermo_unfold = NA, kin_unfolded = NA, kin_final = NA)) +
#     scale_fill_manual( values = c(thermo_fold = "#EBEBEC", kin_native =  "#EBEBEC", RFU = NA,  # nearly white
#                                   thermo_unfold = "#222222", kin_unfolded = "#222222", # dark grey
#                                   kin_final = "#d8d8d8")) + # light grey
#     facet_grid(rows = vars(model_f)) +
#     theme_bw()+
#     scale_x_continuous(expand = c(0,0))+
#     scale_y_continuous(expand = c(0,0.1), breaks = c(0, 0.5, 1))+
#     labs(x = expression('Temperature ('*~degree*C*')'), y = "Relative population" ) +
#     raster_theme_1 -> p
#   p # return the plot object
# }

make_model_plot <- function(df) {
  df %>%
    select(Temperature, thermo_unfold, thermo_fold, kin_native, RFU, kin_final,kin_unfolded) %>%
    gather(variable, value, -Temperature) %>%
    mutate(model = grep_and_gsub( vec_in = .$variable, g_pattern_vec = c("thermo", "kin", "RFU"), rep_TRUE_vec = c("thermo", "kin", "kin"), else_pattern = "other")) %>%
    #mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( expression("Thermodynamic \n model"), expression("Thermo-kinetic \n model"),"DSF")) ) %>%
    mutate(model_f = factor(.$model, levels = c("thermo", "kin","RFU"), labels = c( "Thermodynamic" ,"Thermodynamic-kinetic ","DSF")) ) %>%
    mutate(variable_f = factor(.$variable, levels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_final", "kin_unfolded", "RFU"), labels = c("thermo_fold", "thermo_unfold", "kin_native", "kin_final","kin_unfolded",  "RFU"))) %>%
    
    ggplot(aes(x = Temperature - 273, y = value, group = variable_f, fill = variable_f, color = variable_f)) +
    geom_ribbon(aes(ymin = 0, ymax = value), color = NA, alpha = 0.85) +
    geom_line(size = 1) +
    scale_color_manual(values = c(thermo_fold = NA, kin_native =  NA, RFU = "orange", thermo_unfold = NA, kin_unfolded = NA, kin_final = NA)) +
    scale_fill_manual( values = c(thermo_fold = "#E6E6E6", kin_native =  "#E6E6E6", RFU = NA,  # nearly white
                                  thermo_unfold = "#000000", kin_unfolded = "#000000", # dark grey
                                  kin_final = "#484848")) + # light grey
    facet_grid(rows = vars(model_f)) +
    theme_bw()+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0.1), breaks = c(0, 0.5, 1))+
    labs(x = expression('Temperature ('*~degree*C*')'), y = "Relative population" ) +
    raster_theme_1 -> p
  p # return the plot object
}