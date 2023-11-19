library(tidyverse)
library(data.table)
library(ggpubr)
library(purrr)
library(parallel)
rm(list=ls())

# Suppress summarise info from dplyr
options(dplyr.summarise.inform = FALSE)

# Set needed directory paths ----------------------------------------------
data_dir = "./data/"
res_dir = "./results/"
out_dir = data_dir

data_file = "example_sensor_data.txt"
map_file = "XY_to_real_cm_space.txt"

# Set options -------------------------------------------------------------
numcores = 5

# Read in sensor data -----------------------------------------------------
sensors = fread(paste0(data_dir,data_file))

# Read in interpolation model comparisons ---------------------------------
my_files = paste(c("par","tair","hr","vpd"),"_out.txt", sep = "")
my_og = rbindlist(lapply(paste0(res_dir, my_files), fread), use.names = TRUE)
my_units = data.frame(sensor = c("par","tair","hr","vpd"),
                      unit = c("PPFD (µmol m-2 s-1)",
                               "Air temperature (°C)",
                               "Air humidity (%)",
                               "Air VPD (kPa)"))

# Quick sanity check ------------------------------------------------------
c(length(unique(my_og$sensor)), length(unique(my_og$date)))
c(length(unique(sensors$sensor)), length(unique(sensors$date))) #The first element should be 10x (or however many sensors there are)

## Table needed for model fitting
my_best_time = my_og %>% 
  group_by(date, sensor) %>% 
  slice_min(order_by = rmse, n = 1) %>% 
  ungroup()

## What time points have more than one 'best' model ?
my_check = my_best_time %>% 
  group_by(sensor) %>% 
  count(date) %>% 
  ungroup() %>% 
  filter(n>1) %>% 
  mutate(sensor_date = paste0(sensor, "_", date))

## When there is more than one best model, choose the one that was last (prior time point) selected
my_olo = my_best_time %>% 
  mutate(sensor_date = paste0(sensor, "_", date)) %>% 
  filter(!sensor_date%in%my_check$sensor_date)

my_mto = my_best_time %>% 
  mutate(sensor_date = paste0(sensor, "_", date)) %>% 
  filter(sensor_date%in%my_check$sensor_date)

## All time points with more than one best
to_fix = unique(my_mto$sensor_date)

prior_mod_fun = function(x){
  ## The time point with ties for best model 
  current_to_fix = x
  tp = my_mto %>% 
    filter(sensor_date==current_to_fix)
  ## Prior time points
  tp_prior = my_olo %>% 
    ungroup() %>% 
    filter(sensor==unique(tp$sensor),
           date<unique(tp$date))
  ## Pull the last model from prior time points that had no ties. When the models in question are not present in previous time points, set seed and choose one.
  in_prior = tp_prior %>% 
    filter(model%in%tp$model)
  if (nrow(in_prior)==0) {
    set.seed(1234)
    tp_prior = sample(tp$model, 1)
  }else{
    tp_prior = tp_prior %>% 
      mutate(idx = 1:n()) %>% 
      slice_max(order_by = idx, n = 1) %>% 
      pull(model)
  }
  tp_out = tp %>% 
    filter(model==tp_prior)
  return(tp_out)
}

my_mto_fixed = rbindlist(mclapply(to_fix, prior_mod_fun, mc.cores = numcores))
my_best_time = bind_rows(my_olo, my_mto_fixed) %>% 
  select(-c(sensor_date))

# Interpolation -----------------------------------------------------------
my_dat = data.table::fread(paste0(data_dir,data_file)) %>% 
  mutate(date = as.POSIXct(date, tz = "UTC"))
my_map = data.table::fread(paste0(data_dir,map_file))

## Set resolution of grid
my_res = 20 # in cm

my_map = my_map %>% 
  mutate(cx = ifelse(cx==min(cx), cx-my_res, 
                     ifelse(cx==max(cx), cx+my_res, cx))) %>% 
  mutate(cy = ifelse(cy==min(cy), cy-my_res, 
                     ifelse(cy==max(cy), cy+my_res, cy))) %>% 
  mutate(cx = plyr::round_any(cx, my_res),
         cy = plyr::round_any(cy, my_res)) %>% 
  distinct()
my_map[my_map<0] = 0

## Use only center of pots positions
my_map = my_map %>% 
  drop_na() %>% 
  group_by(X,Y) %>% 
  summarise(cx = median(cx),
            cy = median(cy)) %>% 
  ungroup() %>% 
  select(cx:cy,X:Y)

my_sensors = unique(my_og$sensor)
my_dates = unique(my_og$date)
my_mods = unique(my_og$model)

# Define grid boundaries --------------------------------------------------
my_bbox = my_map %>% 
  filter(cx == min(cx) | cx == max(cx) | cy == min(cy) | cy == max(cy)) %>% 
  slice(chull(cx, cy))
my_grid = my_map[,1:2] %>% 
  distinct()

# Create copy of grid as a raster object ----------------------------------
my_grid_raster = my_grid %>% 
  dplyr::mutate(cx = cx, cy = cy, Z = 0) %>% 
  raster::rasterFromXYZ(
    crs = NA)

# Create copy of grid as sp object ----------------------------------------
my_grid_sp = as(sf::st_as_sf(my_grid, coords = c("cx", "cy"), crs = NA), "Spatial")

my_fit_func = function(x){
  
  # Sensor name
  ms = my_sensors[x]
  
  # Data to be used for modeling and interpolation
  md = my_dat %>% 
    filter(grepl(paste0(ms), sensor)) %>% 
    mutate(sensor = ms)
  
  ## Model names and matching string specifying the smoothing basis to use
  bs_gam = tibble(model = c("Thin_plate_spline", "GAM_Penalized","GAM_Duchon","GAM_Gaussian"),
                  my_bs = c("tp","ts","ds","gp"))
  
  mdd = md %>% 
    left_join(., my_best_time, by = c("sensor","date")) %>% 
    left_join(., bs_gam, by =c("model"))
  
  ## Separate based on type of model
  gam_mods = tryCatch({
    mdd %>% 
      mutate(tmp = date) %>% 
      filter(model%in%bs_gam$model) %>% 
      nest(., data = c(date:cy, kdim,idp, my_bs)) %>% 
      # slice(1:30) %>%
      mutate(my_fit = map(.x = data, .f = ~mgcv::gam(value ~ s(cx, cy, k = unique(.x$kdim), bs = unique(.x$my_bs)), data = .x, method = "REML"))) %>% 
      mutate(my_pred = map(.x = my_fit, .f = ~tibble(cx = my_map$cx, 
                                                     cy = my_map$cy,
                                                     value = predict(.x, my_map)))) %>% 
      mutate(data = map(.x = data, .f = ~.x[,c("date","sensor")] %>% distinct())) %>% 
      select(model, data, my_pred) %>% 
      unnest(c(my_pred, data)) %>% 
      left_join(., my_map, by = c("cx","cy")) %>% 
      relocate(value, .after = Y)
  }, error = function(e) {
    gam_mods = NULL
    return(gam_mods)
  })
  
  idp_mods_dat = mdd %>% 
    mutate(tmp = date) %>% 
    filter(model=="Inverse_dist_weighting") %>% 
    nest(., data = c(date:cy, kdim,idp, my_bs))
  
  nn_mods_dat = mdd %>% 
    mutate(tmp = date) %>% 
    filter(model=="Nearest_neighbor") %>% 
    nest(., data = c(date:cy, kdim,idp, my_bs))
  
  my_gstat_fun = function(x, y){
    ## Convert data to sf object then as a SpatialPointsDataFrame
    sf_obj = sf::st_as_sf(x, coords = c("cx", "cy"), crs = NA)
    sp_obj = as(sf_obj, "Spatial")
    
    ## Set the boundaries of the SpatialPointsDataFrame to be that of the full grid
    sp_obj@bbox[1,1] = my_grid_raster@extent@xmin
    sp_obj@bbox[1,2] =  my_grid_raster@extent@xmax
    sp_obj@bbox[2,1] = my_grid_raster@extent@ymin
    sp_obj@bbox[2,2] = my_grid_raster@extent@ymax
    
    ## Info to match gam models
    my_date = unique(x$date)
    my_sensor = unique(x$sensor)
    
    if (y=="IDW") {
      ## Get the idp to use 
      idp_use = unique(as.numeric(x$idp))
      
      ## Fit the model and predict
      fit_IDW = gstat::gstat(formula = value ~ 1, data = sp_obj, set = list(idp = idp_use))
      invisible(capture.output(fit_out <- predict(fit_IDW, my_grid_sp)))
      fit_out = as_tibble(fit_out@coords) %>% 
        mutate(value = fit_out@data$var1.pred) %>%
        setNames(.,c("cx","cy","value")) %>% 
        mutate(model = "Inverse_dist_weighting", 
               sensor = my_sensor, 
               date = my_date) %>% 
        left_join(., my_map, by = c("cx","cy")) %>% 
        select("model", "date", "sensor", "cx", "cy", "X", "Y", "value")  
    }
    if (y=="NN") {
      ## Fit the model and predict
      fit_IDW = gstat::gstat(formula = value ~ 1, data = sp_obj)
      invisible(capture.output(fit_out <- predict(fit_IDW, my_grid_sp)))
      fit_out = as_tibble(fit_out@coords) %>% 
        mutate(value = fit_out@data$var1.pred) %>%
        setNames(.,c("cx","cy","value")) %>% 
        mutate(model = "Nearest_neighbor", 
               sensor = my_sensor, 
               date = my_date) %>% 
        left_join(., my_map, by = c("cx","cy")) %>% 
        select("model", "date", "sensor", "cx", "cy", "X", "Y", "value")  
    }
    return(fit_out)
  }
  
  ## IDP models 
  if (nrow(idp_mods_dat)==0) {
    idp_mods = NULL
  }else{
    idp_mods = rbindlist(mclapply(idp_mods_dat$data, my_gstat_fun, y = "IDW", mc.cores = numcores))
  }
  
  ## Nearest neighbor models
  if (nrow(nn_mods_dat)==0) {
    nn_mods = NULL
  }else{
    nn_mods = rbindlist(mclapply(nn_mods_dat$data, my_gstat_fun, y = "NN", mc.cores = numcores))
  }
  my_check_cc = my_check %>% 
    filter(sensor==ms)
  
  ## What is the RMSE of the model used
  rmse_used = mdd %>% 
    select(date, rmse, model) %>% 
    distinct()
  
  my_out = bind_rows(gam_mods, idp_mods, nn_mods) %>% 
    select(-c(cx,cy)) %>% 
    distinct() %>% 
    left_join(., rmse_used, by = c("date", "model")) %>% 
    arrange(date)
  
  return(my_out)
}

my_write = rbindlist(lapply(seq_along(my_sensors), my_fit_func))
data.table::fwrite(my_write, 
                   paste0(out_dir,"example_sensor_data_interpolation_best_per_time_out.txt"),
                   row.names = FALSE, col.names = TRUE, sep = "\t")

print(paste0("DONE!"))
