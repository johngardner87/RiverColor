
##########################################################################
### Functions used in making and analyzing river surface reflectance database
###########################################################################

# Function to clean and filter surface reflectance pull
pull_clean <- function(df, dswe = 2, minPix = 5, maxcScore = 900, maxClouds = 50,
                       maxRef = 10000, minRef = 1) {
  
  data <- df %>%
    dplyr::filter_if(is.numeric,all_vars(!is.na(.))) %>%
    dplyr::filter( 
      pixelCount > minPix,
      Cloud_Cover < maxClouds,
      cScore < maxcScore,
      dswe <= dswe,
      hillshadow > 0, 
      blue <= maxRef,
      blue >= minRef,
      green <= maxRef,
      green >= minRef,
      red <= maxRef,
      red >= minRef,
      nir <= maxRef,
      nir >= minRef,
      swir1 <= maxRef,
      swir1 >= minRef,
      swir2 <= maxRef,
      swir2 >= minRef) %>%
    separate(date, into =c("date", "time"), sep="T", remove=T) %>%
    dplyr::mutate(
      sat = factor(sat, levels = c(8,7,5), labels = c( 'LC08','LE07', 'LT05')),
      date_time= ymd_hms(paste0(date, time, sep=" ")),
      date = ymd(date),
      month = as.numeric(month(date)),
      year = year(date),
      hour = hour(date_time)) %>%
    mutate(season = case_when(
      month %in%  9:11 ~ "Fall",
      month %in%  c(12,1,2)  ~ "Winter",
      month %in%  3:5  ~ "Spring",
      TRUE ~ "Summer")) %>%
    mutate(decade= cut(year, breaks = c(1983,1990,1995,2000,2005,2010,2015,2020),
                       labels = c(1990,1995,2000,2005,2010,2015,2020) )) %>%
    rename(LS_ID = 'system:index') %>% 
    dplyr::group_by(ID) %>%
    dplyr::mutate(count =n(),
                  max_year=max(year, na.rm=T),
                  min_year = min(year, na.rm=T),
                  n_years = (max_year - min_year)) %>%
    ungroup() 
  
  return(data)
}

##########################################################################
# Function that calculates band metrics and water color

pull_transform <- function(df, maxRGB=10000, RGB=F) {
  
  if(RGB == T) { 
    
    maxRGB <- maxRGB
    
    data <- df %>%
      filter_at(vars(red, green, blue), all_vars(.< maxRGB))
    
  }else{ 
    data <- df
    
    maxRGB <- df  %>%
      dplyr::select(red, green, blue) %>%
      dplyr::summarise(maxRGB = max(., na.rm=T)) 
    
    maxRGB <- maxRGB$maxRGB
  }
  
  data<- data %>%
    dplyr::mutate(NR = nir/red,
                  BR = blue/red,
                  GR = green/red,
                  SR = swir1/red,
                  BG = blue/green,
                  RG = red/green, 
                  NG = nir/green,
                  SG = swir1/green,
                  BN = blue/nir,
                  GN = green/nir,
                  RN = red/nir,
                  SN = swir1/nir,
                  BS = blue/swir1,
                  GS = green/swir1,
                  RS = red/swir1,
                  NS = nir/swir1,
                  R.GN = red/ (green + nir),
                  R.GB = red/ (green + blue),
                  R.GS = red/ (green + swir1),
                  R.BN = red/ (blue + nir),
                  R.BS = red/ (blue + swir1),
                  R.NS = red/ (nir + swir1),
                  G.BR = green/ (blue + red),
                  G.BN = green / (blue + nir),
                  G.BS = green / (blue + swir1),
                  G.RN = green / ( red + nir),
                  G.RB = green / (red + blue),
                  G.NS = green / (nir + swir1),
                  B.RG = blue / (red + green),
                  B.RN = blue / (red + nir),
                  B.RS = blue / (red + swir1),
                  B.GN = blue / (green + nir),
                  B.GS = blue / (green + swir1),
                  B.NS = blue / (nir + swir1),
                  N.RG = nir / (red + green),
                  N.RB = nir / (red + blue),
                  N.RS = nir / (red + swir1),
                  N.GB = nir / (green + blue),
                  N.GS = nir / (green + swir1),
                  N.BS = nir / (blue  + swir1),
                  
                  GR2 = (green + red) / 2,
                  GN2 = (green + nir) / 2,
                  #blooms
                  BR_G = (blue - red) / green,
                  NS_NR = (nir - swir1) / (red - swir1),
                  fai = nir - (red + (swir1-red)*((830-660)/(1650-660))),
                  # fai = (nir - red) + (red -swir) * (830-660)/(1648-660)
                  N_S= nir - swir1,
                  N_R = nir- red,
                  #
                  ndvi = ((nir-red)/(nir+red)),
                  ndwi = ((green- swir1)/(green + swir1)),
                  ndssi = ((blue - nir)/ (blue + nir)),
                  gn.gn= ((green- nir)/ (green + nir)),
                  hue = rgb2hsv(r=red, g=green, b=blue, maxColorValue = maxRGB)[1,],
                  saturation = rgb2hsv(r=red, g=green,  b=blue, maxColorValue = maxRGB)[2,],
                  bright = rgb2hsv(r=red, g=green,  b=blue, maxColorValue = maxRGB)[3,],
                  bright_tot = (red + green + nir +blue),
                  dw = chroma(R=red, G=green, B=blue),
                  hexcolor = rgb(r=red, g=green, b=blue, maxColorValue = maxRGB)) 

  return(data)
}

#################################################################

# Function used in function above to calculate dominant wavelength in chromaticity 
# colorspace.

chroma <- function(R, G, B) {
  require(colorscience)
  require(tidyverse)
# Converst R,G, and B spectral reflectance to dominant wavelength based
# on CIE chromaticity color space

# see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
# Classification of Inland Water With the Forel-Ule
# Scale: A Case Study of Lake Taihu

# chromaticity.diagram.color.fill()
Xi <- 2.7689*R + 1.7517*G + 1.1302*B
Yi <- 1.0000*R + 4.5907*G + 0.0601*B
Zi <- 0.0565*G + 5.5943*B

x <-  Xi / (Xi + Yi +  Zi)
y <-  Yi / (Xi + Yi +  Zi)
z <-  Zi / (Xi + Yi +  Zi)

# calculate hue angle
alpha <- atan2( (x - (1/3)), (y - (1/3))) * 180/pi

# make look up table for hue angle to wavelength conversion
cie <- cccie31 %>%
  dplyr::mutate(a = atan2( (x - (1/3)), (y - (1/3))) * 180/pi) %>%
  dplyr::filter(wlnm <= 700) %>%
  dplyr::filter(wlnm >=380) 

# find nearest dominant wavelength to hue angle
wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']

return(wl)
}

##############################################################
# Function to add column of chromaticity coordinates and hue angle to dataframe
hueangle <- function(df, R, G, B) {
  
  require(colorscience)

  R = as.data.frame(df)[,R]
  G = as.data.frame(df)[,G]
  B = as.data.frame(df)[,B]
  
  # Add hue angle
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  chroma_x <-  Xi / (Xi + Yi +  Zi)
  chroma_y <-  Yi / (Xi + Yi +  Zi)
  chroma_z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  hue_angle <- atan2( (chroma_x - (1/3)), (chroma_y - (1/3))) * 180/pi
  
  return(cbind(df, chroma_x, chroma_y, chroma_z, hue_angle))
}







# test 3 is chroma x y and z
#attr(test3, "clrsp") <- "CIEXYZ"
#cieplot(test3, cex=0.5)

#test %>%
#  hueangle(., R="red", G="green", B="blue")

#####################################################################
# Function to calculate purity in chromaticity colorspace. Purity is similar
# to saturation is HSV colorspace

purity <- function(R, G, B) {
  require(colorscience)
  # Converst R,G, and B spectral reflectance to dominant wavelength based
  # on CIE chromaticity color space
  
  # see Wang et al 2015. MODIS-Based Radiometric Color Extraction and
  # Classification of Inland Water With the Forel-Ule
  # Scale: A Case Study of Lake Taihu
  
  
  # chromaticity.diagram.color.fill()
  Xi <- 2.7689*R + 1.7517*G + 1.1302*B
  Yi <- 1.0000*R + 4.5907*G + 0.0601*B
  Zi <- 0.0565*G + 5.5943*B
  
  # calculate coordinates on chromaticity diagram
  x <-  Xi / (Xi + Yi +  Zi)
  y <-  Yi / (Xi + Yi +  Zi)
  z <-  Zi / (Xi + Yi +  Zi)
  
  # calculate hue angle
  alpha <- atan2( (x - (1/3)), (y - (1/3))) * 180/pi
  
  # make look up table for hue angle to wavelength conversion
  cie <- cccie31 %>%
    dplyr::mutate(a = atan2( (x - (1/3)), (y - (1/3))) * 180/pi) %>%
    dplyr::filter(wlnm <= 700) %>%
    dplyr::filter(wlnm >=380) #%>%

  # find nearest dominant wavelength to hue angle
  wl <- cie[as.vector(sapply(alpha,function(x) which.min(abs(x - cie$a)))), 'wlnm']

  wl_xy <- wl %>%
    as_data_frame() %>%
    left_join(cie %>%
    dplyr::filter(wlnm %in% wl), by= c("value"= "wlnm"))
  
  xd <- wl_xy %>% pull(x)
    
  yd <- wl_xy %>% pull(y)
  
  P <- sqrt((x - (1/3))^2 + (y - (1/3))^2) / sqrt((xd - (1/3))^2 + (yd - (1/3))^2) 
  
  return(P)
}

####################################################
wlnm2RGB_mod = function(wavelength, Gamma = 0.8, IntensityMax = 1) {
  
  if (IntensityMax < 0) 
    IntensityMax <- 0
  if (IntensityMax > 1) 
    IntensityMax <- 1
  wavelength <- round(wavelength, digits = 0)
  Red <- Green <- Blue <- factorL <- rep(0, length(wavelength))
  w <- which(wavelength > 380 & wavelength <= 440)
  if (length(w) > 0) {
    Red[w] <- -(wavelength[w] - 440)/(440 - 380)
    Green[w] <- 0
    Blue[w] <- 1
  }
  w <- which(wavelength > 440 & wavelength <= 490)
  if (length(w) > 0) {
    Red[w] <- 0
    Green[w] <- (wavelength[w] - 440)/(490 - 440)
    Blue[w] <- 1
  }
  w <- which(wavelength > 490 & wavelength <= 510)
  if (length(w) > 0) {
    Red[w] <- 0
    Green[w] <- 1
    Blue[w] <- -(wavelength[w] - 510)/(510 - 490)
  }
  w <- which(wavelength > 510 & wavelength <= 580)
  if (length(w) > 0) {
    Red[w] <- (wavelength[w] - 510)/(580 - 510)
    Green[w] <- 1
    Blue[w] <- 0
  }
  w <- which(wavelength > 580 & wavelength <= 645)
  if (length(w) > 0) {
    Red[w] <- 1
    Green[w] <- -(wavelength[w] - 645)/(645 - 580)
    Blue[w] <- 0
  }
  w <- which(wavelength > 645 & wavelength <= 780)
  if (length(w) > 0) {
    Red[w] <- 1
    Green[w] <- 0
    Blue[w] <- 0
  }
  w <- which(wavelength > 380 & wavelength <= 420)
  if (length(w) > 0) 
    factorL[w] <- 0.3 + 0.7 * (wavelength[w] - 380)/(420 - 
                                                       380)
  w <- which(wavelength > 420 & wavelength <= 700)
  if (length(w) > 0) 
    factorL[w] <- 1
  w <- which(wavelength > 700 & wavelength <= 780)
  if (length(w) > 0) 
    factorL[w] <- 0.3 + 0.7 * (780 - wavelength[w])/(780 - 
                                                       700)
  R <- Adjust(Red, factorL, Gamma, IntensityMax)
  G <- Adjust(Green, factorL, Gamma, IntensityMax)
  B <- Adjust(Blue, factorL, Gamma, IntensityMax)
  
  output = as_tibble(list(R = R, G = G, B = B)) %>% 
    transmute(color = rgb(R, G, B, maxColorValue = 1))
  
  
  return(as.character(output))
}


###########################################################
# function that calculates the highest mode for hue

hue_mode <- function(data) {

mode <-(data %>%
          pull(hue) %>%
          amps())$Peaks %>%
          as_data_frame() 

colnames(mode) <- c("location", "amp")

out <- mode[which.max(mode$amp), "location"] 

return(out)
}

#########################################################
# function that calculates the mode(s) from a distribution of dominant
# wavelength data, as well as antimodes, and coefficients of bimodality
# 

modesum <-function(data) {
  
  require(modes)
  data_1 <- data %>%
    dplyr::select(-sat)
  
  sum <- data_1 %>%
    dplyr::summarise_all( list(~median(.), ~mean(.), ~sd(.))) 
                         
  mode_df <-(data_1 %>%
    pull(dw) %>%
    amps())$Peaks %>%
    as_tibble() 
  
  colnames(mode_df) <- c("location", "amp")
  
  antimode_df <- (data_1 %>%
                    pull(dw) %>%
                    amps())$Antimode %>%
                    as_tibble()
  colnames(antimode_df) <- c("location", "amp")
#   
  
  if(nrow(mode_df) ==1) {
    antimode <- tibble(location= NA, amp = NA)
    mode2 <- tibble(location= NA, amp = NA)
    mode1 <- mode_df
    
  } else if (nrow(mode_df) > 1 & nrow(antimode_df) == 0) {
    
    mode1 <- mode_df[which.max(mode_df$amp), ] 
    
    mode_df2 <- mode_df %>%
      filter(location != mode1$location)
    
    antimode <- tibble(location= NA, amp = NA)
    mode2 <- mode_df2[which.max(mode_df2$amp), ]
    
  } else if(nrow(mode_df) > 1 & nrow(antimode_df) == 1) {
    
    mode1 <- mode_df[which.max(mode_df$amp), ] 
    
    mode_df2 <- mode_df %>%
      filter(location != mode1$location)
    
    mode2 <- mode_df2[which.max(mode_df2$amp), ]
    
    modes <- bind_rows(mode1, mode2) %>%
      arrange(location)
    
    antimode1 <- antimode_df %>%
      filter(between(location, modes$location[1], modes$location[2])) %>%
      filter(amp == min(amp)) %>%
      filter(row_number()==ceiling(n()/2))

        if(nrow(antimode1)==0) {
        antimode <- tibble(location= NA, amp = NA)
        mode2 <- mode2
        mode1 <- mode_df[which.max(mode_df$amp), ] 
      } else{
        antimode <- antimode1
        mode2 <- mode2
        mode1 <- mode_df[which.max(mode_df$amp), ] 
      }

  } else if(nrow(mode_df) > 1 & nrow(antimode_df) > 1) {

    mode1 <- mode_df[which.max(mode_df$amp), ] 
    
    mode_df2 <- mode_df %>%
      filter(location != mode1$location)
    
    mode2 <- mode_df2[which.max(mode_df2$amp), ]
    
    modes <- bind_rows(mode1, mode2) %>%
      arrange(location)
    
    antimode1 <- antimode_df %>%
      filter(between(location, modes$location[1], modes$location[2])) %>%
      filter(amp == min(amp)) %>%
      filter(row_number()==ceiling(n()/2))
    
    if(nrow(antimode1)==0) {
      antimode <- tibble(location= NA, amp = NA)
      mode2 <- mode2
      mode1 <- mode_df[which.max(mode_df$amp), ] 
    } else{
      antimode <- antimode1
      mode2 <- mode2
      mode1 <- mode_df[which.max(mode_df$amp), ] 
    }
  } else {
      
      antimode <- tibble(location= NA, amp = NA)
      mode2 <- tibble(location= NA, amp = NA)
      mode1 <- tibble(location= NA, amp = NA)
    }
  
  #> 0.55 is bimodal
  bi_coef <- data_1 %>%
    pull(dw) %>%
    bimodality_coefficient()
  
#  bi_amp <- suppressWarnings(data_1 %>%
#    pull(dw) %>%
#    bimodality_amplitude(., fig=F))

  skew <- skewness(data_1 %>% pull(dw))
  
  out <- tibble(dw_mode1=mode1$location ,dw_mode2=mode2$location, 
               dw_mode1_amp=mode1$amp ,dw_mode2_amp=mode2$amp, 
                dw_antimode=antimode$location,  dw_bi_coef = bi_coef,
                 skew=skew) %>%
    mutate(higher_mode = ifelse(mode1$location > mode2$location, 1, 2))
  
  # calculate fraction of samples around each mode
  if( bi_coef > 0.54 & is.finite(out$dw_antimode) & !is.na(out$dw_mode2)) {
    
  x <- data %>%
    dplyr::select(dw) %>%
    mutate(count = n(),
           mode = ifelse(dw > out$dw_antimode, "high", "low")) %>%
    group_by(mode, count) %>%
    summarise(color_count=n()) %>%
    mutate(mode_frac = color_count/count) %>%
    ungroup() %>%
    dplyr::select(mode, mode_frac) %>%
    spread(key = mode, value=mode_frac) %>%
    dplyr::rename_all(paste0, "_frac")
  
  x_sat <-  data %>%
    dplyr::select(dw, sat) %>%
    mutate(count = n(),
           mode = ifelse(dw > out$dw_antimode, "high", "low")) %>%
    group_by(sat) %>%
    mutate(count_sat = n()) %>%
    ungroup() %>%
    group_by(sat, mode, count_sat) %>%
    summarise(color_count=n()) %>%
    mutate(mode_frac = color_count/count_sat) %>%
    mutate(mode_sat = paste(mode, sat, sep="_")) %>%
    ungroup() %>%
    dplyr::select(mode_sat, mode_frac) %>%
    spread(key = mode_sat, value=mode_frac) %>%
    dplyr::rename_all(paste0, "_frac")
  } else { 
    x <- tibble(high_frac =NA, low_frac=NA)
    x_sat <- tibble( high_LC08_frac = NA, high_LE07_frac= NA,
                     high_LT05_frac=NA , low_LC08_frac=NA  ,
                     low_LE07_frac=NA ,low_LT05_frac=NA )
  }
  
  output <- as_tibble(cbind(sum, out, x, x_sat)) 
    
#if(nrow(output) >1) {print("error: too many rows in output")}

# if bimodal..then do fraction in high or low
# output for which mode is higher amp, Main_mode
  # fraction measurements around each mode

  return(output)
}


# nest_ID <- sr_clean_all %>%
#   filter(ID == 10727  ) %>%
#   filter(dw > 400, dw < 699) %>%
#   filter(count > 50 & n_years > 10) %>%
#   filter(Cloud_Cover < 50) %>%
#   dplyr::select(ID,decade, sat, Cloud_Cover, pixelCount,dw,
#                 ndssi, hue, saturation, bright, bright_tot,
#                 fai, BR_G , NS_NR, N_S, N_R,  ndvi) %>%
#   group_by(ID, decade) %>%
#   mutate(count_decade = n())  %>%
#   filter(count_decade > 10) %>%
#   nest() 
# 
# #data1 <- (nest_ID %>% filter(ID==18674))
# 
# data <- nest_ID$data[[1]]
# 
# test <- modesum(data)
#
# ggplot() + geom_density(aes(data$dw))

#######################################################
# Function to summarise color data by year

sum_yr <-function(data) {
  
  require(modes)
  data_1 <- data %>%
    dplyr::select(-sat)
  
  sum <- data_1 %>%
    dplyr::summarise_all( list(~median(.), ~mean(.), ~sd(.))) 
  
  mode_df <-(data_1 %>%
               pull(dw) %>%
              amps())$Peaks %>%
    as_tibble() 
  
  colnames(mode_df) <- c("location", "amp")
   
  mode1 <- mode_df[which.max(mode_df$amp), ] 
  
  out <- tibble(dw_mode1=mode1$location)
  
  output <- as_tibble(cbind(sum, out)) 
   
  return(output)
}

##################################################

colorFracSat <- function(data) {
  
  x <- data %>%
    dplyr::select(dw) %>%
    mutate(count = n(),
           color = case_when(
             dw < 495 ~ "blue",
             dw >= 495 & dw < 560 ~ "green",
             dw >=560 ~ "yellow"  )) %>%
    group_by(color, count) %>%
    summarise(color_count=n()) %>%
    mutate(color_frac = color_count/count) %>%
    ungroup() %>%
    dplyr::select(color, color_frac) %>%
    spread(key = color, value=color_frac) %>%
    dplyr::rename_all(paste0, "_frac")
  
  x_sat <- data %>%
    dplyr::select(dw, sat) %>%
    mutate(color = case_when(
             dw < 495 ~ "blue",
             dw >= 495 & dw < 560 ~ "green",
             dw >=560 ~ "yellow"  )) %>%
    group_by(sat) %>%
    mutate(count_sat = n()) %>%
    ungroup() %>%
    group_by(sat, color, count_sat) %>%
    summarise(color_count=n()) %>%
    mutate(color_frac = color_count/count_sat) %>%
    mutate(color_sat = paste(color, sat, sep="_")) %>%
    ungroup() %>%
    dplyr::select(color_sat, color_frac) %>%
    spread(key = color_sat, value=color_frac) %>%
    dplyr::rename_all(paste0, "_frac")
  
  
out <- cbind(x, x_sat)

return(out)
}

#################################################
# 
# colorFracTest <- function(data) {
#   
#   x <- data %>%
#     dplyr::select(dw) %>%
#     mutate(count = n(),
#            mode = ifelse(dw > out$dw_antimode, "high", "low")) %>%
#     group_by(mode, count) %>%
#     summarise(color_count=n()) %>%
#     mutate(mode_frac = color_count/count) %>%
#     ungroup() %>%
#     dplyr::select(mode, mode_frac) %>%
#     spread(key = mode, value=mode_frac) %>%
#     dplyr::rename_all(paste0, "_frac")
#   
#   x_sat <-  data %>%
#     dplyr::select(dw, sat) %>%
#     mutate(count = n(),
#            mode = ifelse(dw > out$dw_antimode, "high", "low")) %>%
#     group_by(sat) %>%
#     mutate(count_sat = n()) %>%
#     ungroup() %>%
#     group_by(sat, mode, count_sat) %>%
#     summarise(color_count=n()) %>%
#     mutate(mode_frac = color_count/count_sat) %>%
#     mutate(mode_sat = paste(mode, sat, sep="_")) %>%
#     ungroup() %>%
#     dplyr::select(mode_sat, mode_frac) %>%
#     spread(key = mode_sat, value=mode_frac) %>%
#     dplyr::rename_all(paste0, "_frac")
#   
#   out <- cbind(x, x_sat)
#   
#   return(out)
# }

#########################################################
# 
ls_correction <- function(data) {

# bands are quite different but dw is not 
sr_57 <- data %>%
  filter(sat %in% c("LE07", "LT05")) %>%
  filter(between(date, "1999-01-01", "2012-05-01" )) %>%
  # filter to site with enough data
  filter(n_years > 10) %>%
  select(ID, date, sat, count, n_years, blue, red, green, nir, swir1, swir2) %>%
  gather(blue:swir2, green, key='band', value='value') 

# plot distribution compare means
# ggplot(data=sr_57, aes(x=value, color=sat)) +
#   geom_density() +
#   facet_wrap(~band, scales="free")


# do ranking plotting percentiles, joining, and correcting
sr_57_rank  <- sr_57 %>%
  droplevels() %>%
  filter(sat =="LT05") %>%
  group_by(band) %>%
  nest() %>%
  mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
          ret = purrr::invoke_map(tibble, ret)) %>%
  unnest(ret) %>%
  dplyr::select(-data) %>%
  pivot_longer(
    cols= contains("%")
  ) %>%
  mutate(quant = parse_number(name)/100) %>%
  rename(value_5 = value) %>%
  inner_join(sr_57 %>%
               droplevels() %>%
               filter(sat =="LE07") %>%
               group_by(band) %>%
               nest() %>%
               mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
                       ret = purrr::invoke_map(tibble, ret)) %>%
               unnest(ret) %>%
               dplyr::select(-data) %>%
               pivot_longer(
                 cols= contains("%")
               ) %>%
               mutate(quant = parse_number(name)/100) %>%
               rename(value_7 = value) %>%
               dplyr::select(-name),
             by=c("band", "quant")) 

poly_5_trunc <- function(df){
  lm(value_7 ~ poly(value_5, 2, raw=T), data = df %>%
       filter(!quant %in% c(0, 1))  )
}
poly_5_all <- function(df){
  lm(value_7 ~ poly(value_5, 2, raw=T), data = df)
}

## polynomial correction fit
poly_57 <- sr_57_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  nest(-band) %>%
  mutate( model = purrr::map(data, poly_5_trunc)) %>%
  mutate( model_all = purrr::map(data, poly_5_all)) %>%
  mutate( pred = purrr::map2(model, data, predict)) %>%
  mutate( pred_all = purrr::map2(model_all, data, predict)) %>%
  unnest(c(pred, pred_all, data))  %>%
  dplyr::select(-model, -model_all)

png("D:/Dropbox/projects/rivercolor/figs/poly_5_correction.png",
    units="in", width = 6, height=4, res = 300)

ggplot(poly_57) +
   geom_point( aes(x=value_5, y= value_7))+
   geom_point( aes(x=pred_all, y= value_7), color="red", alpha=0.5)+
   geom_point( aes(x=pred, y= value_7), color="blue", alpha=0.5)+
   geom_abline(aes(slope=1, intercept=0)) +
   facet_wrap(~band, scales="free") +
  theme_bw() 

dev.off()
#ggsave("figs/poly_5_correction.png", units="in", width = 6, height=4, dpi = 300)

#
write_csv(poly_57, "D:/Dropbox/projects/rivercolor/out/ls_57_quantile_correction.csv")

##
coef_5 <- sr_57_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  filter(!quant %in% c(0, 1)) %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_5, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LT05") %>%
  mutate(fit = "98_quant")

coef_5_all <- sr_57_rank %>%
  ungroup() %>%
  #  filter(band != "dw") %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_5, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LT05") %>%
  mutate(fit = "all_quant")

########################

sr_78 <- data %>%
  filter(sat %in% c("LE07", "LC08")) %>%
  filter(date > "2013-04-11" ) %>%
  # filter to site with enough data
  filter(n_years > 10) %>%
  select(ID, date, sat, count, n_years, blue, red, green, nir, swir1, swir2) %>%
  gather(blue:swir2, key='band', value='value') 

# plot distribution compare means
# ggplot(data=sr_78, aes(x=value, color=sat)) +
#   geom_density() +
#   facet_wrap(~band, scales="free")

# do ranking plotting percentiles, joining, and correcting
sr_78_rank  <- sr_78 %>%
  droplevels() %>%
  filter(sat =="LC08") %>%
  group_by(band) %>%
  nest() %>%
  mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
          ret = purrr::invoke_map(tibble, ret)) %>%
  unnest(ret) %>%
  dplyr::select(-data) %>%
  pivot_longer(
    cols= contains("%")
  ) %>%
  mutate(quant = parse_number(name)/100) %>%
  rename(value_8 = value) %>%
  inner_join(sr_78 %>%
               droplevels() %>%
               filter(sat =="LE07") %>%
               group_by(band) %>%
               nest() %>%
               mutate( ret = purrr::map(data, ~quantile(.$value, probs = seq(0,1,0.01))),
                       ret = purrr::invoke_map(tibble, ret)) %>%
               unnest(ret) %>%
               dplyr::select(-data) %>%
               pivot_longer(
                 cols= contains("%")
               ) %>%
               mutate(quant = parse_number(name)/100) %>%
               rename(value_7 = value) %>%
               dplyr::select(-name),
             by=c("band", "quant"))  

# ggplot(sr_78_rank, aes(value_8, value_7)) +
#   geom_point() +
#   geom_abline(aes(slope=1, intercept=0), color="red") +
#   geom_smooth(method="lm", color="blue") +
#   facet_wrap(~band, scales = "free")

poly_8_trunc <- function(df){
  lm(value_7 ~ poly(value_8, 2), data = df %>%
       filter(!quant %in% c(0, 1))  )
}
poly_8_all <- function(df){
  lm(value_7 ~ poly(value_8, 2), data = df)
}

poly_78 <- sr_78_rank %>%
  ungroup() %>%
#  filter(band != "dw") %>%
  nest(-band) %>%
  mutate( model = purrr::map(data, poly_8_trunc)) %>%
  mutate( model_all = purrr::map(data, poly_8_all)) %>%
  mutate( pred = purrr::map2(model, data, predict)) %>%
  mutate( pred_all = purrr::map2(model_all, data, predict)) %>%
  unnest(c(pred, pred_all, data)) %>%
  dplyr::select(-model, -model_all)

png("D:/Dropbox/projects/rivercolor/figs/poly_8_correction.png",
    units="in", width = 6, height=4, res = 300)

ggplot(poly_78) +
  geom_point( aes(x=value_8, y= value_7))+
  geom_point( aes(x=pred_all, y= value_7), color="red", alpha=0.5)+
  geom_point( aes(x=pred, y= value_7), color="blue", alpha=0.5)+
  geom_abline(aes(slope=1, intercept=0)) +
  facet_wrap(~band, scales="free") +
  theme_bw()
dev.off()

#ggsave("figs/poly_8_correction.png", units="in", width = 6, height=4, dpi = 300)

#
write_csv(poly_78, "D:/Dropbox/projects/rivercolor/out/ls_78_quantile_correction.csv")

#
coef_8 <- sr_78_rank %>%
  ungroup() %>%
  filter(!quant %in% c(0, 1)) %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_8, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LC08") %>%
  mutate(fit = "98_quant")

#
coef_8_all <- sr_78_rank %>%
  ungroup() %>%
  group_by(band)  %>%
  nest() %>%
  mutate( model = purrr::map(data, ~lm(value_7 ~ poly(value_8, 2, raw=T), data = .) %>%
                               tidy %>%
                               dplyr::select(term, estimate) %>%
                               spread(term, estimate))) %>%
  unnest(model) %>%
  dplyr::select(-data) %>%
  rename(band= 1, intercept=2, coef1=3, coef2=4 )  %>%
  mutate(sat = "LC08") %>%
  mutate(fit = "all_quant")

#
coef_7 <- tibble(band = c("blue", "red", "green", "nir", "swir1", "swir2"), intercept = 0, coef1=1, coef2=0, sat= "LE07")

#
corr_coef <- bind_rows(coef_5, coef_7, coef_8, coef_5_all, coef_8_all) %>%
  ungroup()

return(corr_coef)

}

#####################################################################
# Funntion for Mann-Kendall trend analysis with variance correction
mod_mk <- function(x, col) {
  
  require(modifiedmk)
  x <- x %>%
    pull(col)
  
  y <- mmky(x)
  
  out = tibble(tau = y["Tau"], sen = y["Sen's slope"], new.p = y["new P-value"], old.p = y["old P.value"],
               new.var = y["new.variance"], old.var = y["old.variance"], n=length(x))
  return(out)
}

####################################################################
# Function for classic Mann-Kendall test
mkmk <- function(x, col) {
  
  require(modifiedmk)
  x <- x %>%
    pull(col)
  
  y <- mkttest(x)
  
  out = tibble(tau = y["Tau"], sen = y["Sen's slope"], p = y["P-value"], var = y["Var(S)"],
               s = y["S"], n=length(x))
  return(out)
}

#################################################################

lm_trend <- function(data, col) {
  
  y <- data %>% pull(col) 
  
  x <- data %>% pull(decade) %>%
    as.character() %>%
    as.numeric()
  
  lm <- lm(y ~ x)

  out <- tibble(slope =lm$coefficients[2], intercept=lm$coefficients[1],
               n=length(lm$model$y), slope.SE = summary(lm)$coefficients[2,2],
               r.squared = summary(lm)$r.squared, p =summary(lm)$coefficients[2,4]    )
  return(out)
}

poly_trend <- function(data, col) {
  
  y <- data %>% pull(col) 
  
  x <- data %>% pull(decade) %>%
    as.character() %>%
    as.numeric()
  
  lm <- lm(y ~ poly(x, 2, raw=T), data = data)

  out <- tibble( slope2 =lm$coefficients[3], slope1 =lm$coefficients[2],
                intercept=lm$coefficients[1],
                n=length(lm$model$y), r.squared = summary(lm)$r.squared,
                p = glance(lm)$p.value) %>%
    mutate(dir = ifelse(slope2 > 0, "up", "down"))
return(out)
}
