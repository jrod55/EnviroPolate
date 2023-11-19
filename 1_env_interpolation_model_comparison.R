library(tidyverse)
library(sf)
library(sp)
library(raster)
library(gstat)
library(mgcv)
library(caret)
library(patchwork)
library(viridis)
rm(list=ls())

# Suppress summarise info from dplyr
options(dplyr.summarise.inform = FALSE)

# Set needed directory paths ----------------------------------------------
data_dir = "./data/"
out_dir = "./results/"
tmp_dir = tmpDir()

data_file = "example_sensor_data.txt"
map_file = "XY_to_real_cm_space.txt"

# Set options -------------------------------------------------------------
n_sens = 10
ncores = 5

# Read in data ------------------------------------------------------------
my_d = data.table::fread(paste0(data_dir,data_file))
my_map = data.table::fread(paste0(data_dir,map_file))

# Define pot boundaries ---------------------------------------------------
my_pots = my_map %>% 
  drop_na() %>% 
  group_by(X,Y) %>% 
  summarise(cx = median(cx),
            cy = median(cy))

# Set resolution of grid --------------------------------------------------
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
my_grid_sp = as(st_as_sf(my_grid, coords = c("cx", "cy"), crs = NA), "Spatial")

# Define variables for data subsetting ------------------------------------
env_vars = c("par","tair","hr","vpd")
env_units = c(" PPFD \n (µmol m-2 s-1)",
              " Air temperature \n (°C)",
              " Air humidity \n (%)",
              " Air VPD \n (kPa)")
env_round = c(10, 5, 5, 0.1)
time_points = unique(my_d$date)
ml = length(time_points)

my_vars = data.frame(v1 = rep(sort(time_points), 4),
                     v2 = c(rep(env_vars[1], ml),
                            rep(env_vars[2], ml),
                            rep(env_vars[3], ml),
                            rep(env_vars[4], ml)),
                     v3 = c(rep(env_units[1], ml),
                            rep(env_units[2], ml),
                            rep(env_units[3], ml),
                            rep(env_units[4], ml)),
                     v4 = c(rep(env_round[1], ml),
                            rep(env_round[2], ml),
                            rep(env_round[3], ml),
                            rep(env_round[4], ml))
)

# GAM modeling step -------------------------------------------------------
mod_fun = function(ms){
  
  ## For testing
  # ms = 1
  
  ## Data subsetting 
  tt = as.POSIXct(my_vars[ms,1])
  cc = as.character(my_vars[ms,2])
  cc_lab = as.character(my_vars[ms,3])
  cc_rnd = as.numeric(my_vars[ms,4])
  
  ## Ranges
  d0 = my_d %>% 
    filter(grepl(paste0(cc), sensor))
  nz_vals = plyr::round_any(d0$value, cc_rnd)
  val_min = min(plyr::round_any(d0$value, cc_rnd, f = floor))
  val_max = max(plyr::round_any(d0$value, cc_rnd, f = ceiling))
  nz_vals = nz_vals[!nz_vals<=val_min]
  val_quants = quantile(nz_vals, probs = c(0.25, 0.50, 0.75, 0.99))
  
  ## Grid limits and breaks
  grd_lims_x = c(0, 1825)
  grd_lims_y = c(0, 1400)
  grd_breaks_x = seq(grd_lims_x[1], grd_lims_x[2], 365)
  grd_breaks_y = seq(grd_lims_y[1], grd_lims_y[2], 350)
  
  ## Color breaks from all data timepoints
  mybreaks = plyr::round_any(ceiling(seq(val_min, val_max, length.out = 5)), cc_rnd)
  mybreaks_alt = c(val_min, as.numeric(val_quants))
  mycols = colorRampPalette(viridis(n = 5, option = "plasma"))
  mybreaks_alt = seq(min(mybreaks_alt), max(mybreaks_alt), length.out = 5)
  
  ## Filtering for a given time point
  d1 = d0 %>% 
    filter(date == tt)
  
  ## Calculate variance among the data points
  d1_variance = var(d1$value)
  
  # Bounds for this given time point
  lower_bound = min(d1$value)
  upper_bound = max(d1$value)
  
  ## Filtering for all points before the given time point
  d2 = d0 %>%
    filter(date <= tt)
  
  ## Plot of the grid
  grid_plot = ggplot() +
    geom_point(data = my_pots, aes(x = cx, y = cy), size = 0.25, shape = 4) +
    geom_point(data = d1,
               mapping = aes(x = cx, y = cy, color = value), size = 10, shape  = 19) +
    scale_color_gradientn(name = paste0(cc_lab),
                          values = scales::rescale(mybreaks_alt),
                          colors = mycols(5),
                          breaks = mybreaks_alt,
                          limits = c(min(mybreaks_alt), max(mybreaks_alt)),
                          na.value = mycols(5)[5])+
    scale_x_continuous(breaks = grd_breaks_x, limits = grd_lims_x)+
    scale_y_continuous(breaks = grd_breaks_y, limits = grd_lims_y)+
    theme_bw()+
    theme(axis.text.x = element_text(size = 10,face = "bold", color = "black"),
          strip.text = element_text(size = 10,face = "bold", color = "black"),
          axis.text.y = element_text(size = 10,face = "bold", color = "black"),
          axis.title =  element_text(size = 10,face = "bold", color = "black"),
          plot.title = element_text(size = 10, face = "bold", color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "PhenoArch X Coord. (cm)",
         y = "PhenoArch Y Coord. (cm)")
  
  ## Line plot cummulative
  my_dates = c(min(d0$date), max(d0$date))
  line_plot = d2 %>%
    group_by(date) %>%
    summarise(value = mean(value)) %>%
    ggplot(., aes(x=date, y = value)) +
    geom_line(color = "black") +
    theme_minimal()+
    scale_y_continuous(breaks = mybreaks_alt,
                       limits = c(min(mybreaks_alt), max(mybreaks_alt))) +
    scale_x_datetime(labels = date_format("%Y-%m-%d %H-%M-%S"),
                     date_breaks = "2 weeks", date_minor_breaks = "1 day",
                     date_labels = "%D", limits = my_dates)+
    theme(axis.text.x = element_text(size = 10,face = "bold", color = "black"),
          strip.text = element_text(size = 10,face = "bold", color = "black"),
          axis.text.y = element_text(size = 10,face = "bold", color = "black"),
          axis.title =  element_text(size = 10,face = "bold", color = "black"),
          plot.title = element_text(size = 10, face = "bold", color = "black")) +
    ggtitle(paste0("PhenoArch ",toupper(cc), " sensor"),
            subtitle = paste0("Date and time: ", tt))+
    labs(x = "Date",
         y = paste0(gsub("\n ", "", cc_lab), "\n Avg. across ", n_sens, " sensors"))
  
  ## Convert data to sf object then as a SpatialPointsDataFrame
  sf_obj = st_as_sf(d1, coords = c("cx", "cy"), crs = NA)
  sp_obj = as(sf_obj, "Spatial")
  
  ## Set the boundaries of the SpatialPointsDataFrame to be that of the full grid
  sp_obj@bbox[1,1] = my_grid_raster@extent@xmin
  sp_obj@bbox[1,2] =  my_grid_raster@extent@xmax
  sp_obj@bbox[2,1] = my_grid_raster@extent@ymin
  sp_obj@bbox[2,2] = my_grid_raster@extent@ymax
  
  ## Model fitting
  rmse = function(x){
    invisible(capture.output(rmse_val <- gstat.cv(x, nfold = n_sens)))
    rmse_val = sqrt(sum(rmse_val$residual^2) / length(sp_obj))
  }
  
  ### Nearest neighbor 
  fit_NN = gstat::gstat(formula = value ~ 1, data = sp_obj)
  fit_NN_out = data.frame(model = "Nearest_neighbor",
                          rmse = rmse(fit_NN))
  
  ### Inverse Distance Weighting...test power functions between range below and select the power that minimized rmse-cv
  my_pwr = seq(1,100,1)
  fit_IDW_test = sapply(seq_along(my_pwr), function(x){rmse(gstat::gstat(formula = value ~ 1, data = sp_obj, set = list(idp = my_pwr[x])))})
  idp_use = my_pwr[which.min(fit_IDW_test)] # First occurrence of the min value. i.e) lowest possible power function in case there are ties
  fit_IDW = gstat::gstat(formula = value ~ 1, data = sp_obj, set = list(idp = idp_use))
  fit_IDW_out = data.frame(model = "Inverse_dist_weighting",
                           rmse = rmse(fit_IDW),
                           idp = idp_use)
  
  ### Generalized Additive Models using different basis functions (bs) and number of basis functions to use in the smoothing (k)
  ## Prepare data for manual cross validation
  d1_cv = d1
  d1_folds = cut(seq(1,nrow(d1_cv)),breaks=n_sens,labels=FALSE)
  
  ## Grid combinations of possible GAM models
  my_ps = expand_grid(my_k = seq(n_sens/2,n_sens-1,1),
                      my_s = c("gp", "ds", "ts", "tp"),
                      my_f = seq(1,n_sens,1))
  
  ## Function for calculating residuals for rmse across all folds of leave-one-out cross validation when using the mgcv package
  rmse_resid = function(x){
    testIndexes = my_ps$my_f[x]
    testData = d1_cv[testIndexes, ]
    trainData = d1_cv[-testIndexes, ]
    my_fit = tryCatch({
      mgcv::gam(value ~ s(cx, cy, k = my_ps$my_k[x], bs = my_ps$my_s[x]), data = trainData, method = "REML")
    }, error = function(e) {
      my_fit = FALSE
      return(my_fit)
    })
    if (is.logical(my_fit)) {
      test_resid = NA
    }else{
      my_fit_pred_test = predict.gam(my_fit, testData[,c("cx","cy")], type = "response")
      test_resid = testData$value - my_fit_pred_test   
    }
    return(test_resid)
  }
  
  ## Interpolating with NN as NULL model
  invisible(capture.output(fit_NN_vals <- predict(fit_NN, my_grid_sp)))
  fit_NN_vals = as_tibble(fit_NN_vals@coords) %>% 
    mutate(value = fit_NN_vals@data$var1.pred) %>%
    setNames(.,c("cx","cy","value")) %>% 
    mutate(model = "Nearest_neighbor")
  
  ## If there is zero variance in the data point set all predictions to NN and choose NN as the best model
  if (d1_variance==0) {
    fit_IDW_vals = fit_NN_vals %>% 
      mutate(model = "Inverse_dist_weighting")
    fit_GAM_gp_vals = fit_NN_vals %>% 
      mutate(model = "GAM_Penalized")
    fit_GAM_ds_vals = fit_NN_vals %>% 
      mutate(model = "GAM_Duchon")
    fit_GAM_ts_vals = fit_NN_vals %>% 
      mutate(model = "GAM_Gaussian")
    fit_TPS_vals = fit_NN_vals %>% 
      mutate(model = "Thin_plate_spline")
    
    plt_dat = bind_rows(fit_NN_vals,
                        fit_IDW_vals,
                        fit_GAM_gp_vals,
                        fit_GAM_ds_vals,
                        fit_GAM_ts_vals,
                        fit_TPS_vals
    ) %>% 
      mutate(value = ifelse(value<lower_bound, lower_bound, value)) %>% 
      mutate(value = ifelse(value>upper_bound, upper_bound, value)) %>% 
      mutate(model = factor(model, levels = c("Nearest_neighbor", "Inverse_dist_weighting", "Thin_plate_spline", "GAM_Gaussian","GAM_Duchon", "GAM_Penalized")))
    
    tab_out = fit_NN_out %>% 
      add_row(model = "Inverse_dist_weighting", rmse = NA) %>% 
      add_row(model = "Thin_plate_spline", rmse = NA) %>% 
      add_row(model = "GAM_Gaussian", rmse = NA) %>% 
      add_row(model = "GAM_Duchon", rmse = NA) %>% 
      add_row(model = "GAM_Penalized", rmse = NA) %>%
      mutate(date = tt,
             sensor = cc)
    
    best_mod = tab_out %>% 
      slice_min(order_by = rmse, n = 1) %>% 
      mutate(cx = grd_lims_x[2], cy = grd_lims_y[2]-55) %>% 
      mutate(model = factor(model, levels = c("Nearest_neighbor", "Inverse_dist_weighting", "Thin_plate_spline", "GAM_Gaussian","GAM_Duchon", "GAM_Penalized")))
    
  }else{
    
    ## Fit IDW models
    invisible(capture.output(fit_IDW_vals <- predict(fit_IDW, my_grid_sp)))
    fit_IDW_vals = as_tibble(fit_IDW_vals@coords) %>% 
      mutate(value = fit_IDW_vals@data$var1.pred) %>% 
      setNames(.,c("cx","cy","value")) %>% 
      mutate(model = "Inverse_dist_weighting")
    
    ## Get CV residuals from GAM models
    fit_GAM_test = sapply(seq_along(my_ps$my_k), rmse_resid)
    
    ## Add 10-fold LOOCV RMSE to the parameters used. Note that when an NA value is returned, which occurs when there is no variance in the LOO-fold, it is excluded from the average calculation
    my_ps1 = my_ps %>% 
      mutate(resid = fit_GAM_test) %>% 
      drop_na(resid) %>% 
      group_by(my_k, my_s) %>% 
      summarise(rmse = sqrt(sum(resid^2)/n())) %>% 
      group_by(my_s) %>% 
      slice_min(order_by = rmse, n = 1)
    
    ## Get the smoothing parameter
    k_use_gp = min(my_ps1 %>% filter(my_s=="gp") %>% pull(my_k))
    k_use_ds = min(my_ps1 %>% filter(my_s=="ds") %>% pull(my_k))
    k_use_ts = min(my_ps1 %>% filter(my_s=="ts") %>% pull(my_k))
    k_use_tp = min(my_ps1 %>% filter(my_s=="tp") %>% pull(my_k))
    
    fit_GAM_gp = mgcv::gam(value ~ s(cx, cy, k = k_use_gp, bs = "gp"), data = d1, method = "REML")
    fit_GAM_ds = mgcv::gam(value ~ s(cx, cy, k = k_use_ds, bs = "ds"), data = d1, method = "REML")
    fit_GAM_ts = mgcv::gam(value ~ s(cx, cy, k = k_use_ts, bs = "ts"), data = d1, method = "REML")
    fit_GAM_tp = mgcv::gam(value ~ s(cx, cy, k = k_use_tp, bs = "tp"), data = d1, method = "REML")
    
    mod_names = tibble(my_s = c("ds","gp","ts", "tp"),
                       model = c("GAM_Duchon","GAM_Gaussian","GAM_Penalized","Thin_plate_spline"))
    
    fit_GAM_out = my_ps1 %>% 
      ungroup() %>% 
      left_join(., mod_names, by = "my_s") %>% 
      mutate(kdim = my_k) %>% 
      dplyr::select(-c(my_s,my_k))
    
    fit_GAM_gp_vals = my_grid %>% 
      mutate(value = predict.gam(fit_GAM_gp, my_grid, type = "response"),
             model = "GAM_Gaussian")
    
    fit_GAM_ds_vals = my_grid %>% 
      mutate(value = predict.gam(fit_GAM_ds, my_grid, type = "response"),
             model = "GAM_Duchon")
    
    fit_GAM_ts_vals = my_grid %>% 
      mutate(value = predict.gam(fit_GAM_ts, my_grid, type = "response"),
             model = "GAM_Penalized")
    
    fit_TPS_vals = my_grid %>% 
      mutate(value = predict.gam(fit_GAM_tp, my_grid, type = "response"),
             model = "Thin_plate_spline")
    
    plt_dat = bind_rows(fit_NN_vals,
                        fit_IDW_vals,
                        fit_GAM_gp_vals,
                        fit_GAM_ds_vals,
                        fit_GAM_ts_vals,
                        fit_TPS_vals
    ) %>% 
      mutate(value = ifelse(value<lower_bound, lower_bound, value)) %>% 
      mutate(value = ifelse(value>upper_bound, upper_bound, value)) %>% 
      mutate(model = factor(model, levels = c("Nearest_neighbor", "Inverse_dist_weighting", "Thin_plate_spline", "GAM_Gaussian","GAM_Duchon", "GAM_Penalized")))
    
    tab_out = bind_rows(fit_GAM_out, fit_IDW_out, fit_NN_out) %>% 
      mutate(date = tt,
             sensor = cc)
    
    best_mod = tab_out %>% 
      slice_min(order_by = rmse, n = 1) %>% 
      mutate(cx = grd_lims_x[2], cy = grd_lims_y[2]-55) %>% 
      mutate(model = factor(model, levels = c("Nearest_neighbor", "Inverse_dist_weighting", "Thin_plate_spline", "GAM_Gaussian","GAM_Duchon", "GAM_Penalized")))
    
  }
  
  mods_plot = ggplot(plt_dat) +
    geom_raster(aes(x = cx, y = cy, fill = value)) +
    geom_text(data = best_mod, aes(x = cx, y = cy, label = "*"), color = "red", size = 5)+
    scale_fill_gradientn(name = paste0(cc_lab), 
                         values = scales::rescale(mybreaks_alt),
                         colors = mycols(5),
                         breaks = mybreaks_alt,
                         limits = c(min(mybreaks_alt), max(mybreaks_alt)),
                         na.value = mycols(5)[5]) +
    scale_x_continuous(breaks = grd_breaks_x, limits = grd_lims_x)+
    scale_y_continuous(breaks = grd_breaks_y, limits = grd_lims_y)+
    theme_bw() +
    theme(axis.text.x = element_text(size = 8,face = "bold", color = "black"),
          strip.text = element_text(size = 10,face = "bold", color = "black"),
          axis.text.y = element_text(size = 8,face = "bold", color = "black"),
          axis.title =  element_text(size = 10,face = "bold", color = "black"),
          plot.title = element_text(size = 10, face = "bold", color = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "PhenoArch X Coord. (cm)",
         y = "PhenoArch Y Coord. (cm)\n ") +
    facet_wrap(~model)+
    theme(legend.position="none")
  
  ## Path where figures for video are saved
  tmp_path = paste0(tmp_dir,cc)
  if (dir.exists(tmp_path)) {
    NULL
  }else{
    dir.create(paste0(tmp_path), recursive = TRUE)
  }
  ggsave(cowplot::plot_grid(line_plot, grid_plot, mods_plot, ncol = 1, rel_heights = c(0.5,0.6, 0.75), align = "hv", axis = "tblr"),
         file = paste0(tmp_path, "/img0",ms, ".png"), width = 8.25, height = 8.5, units = "in", bg = "white", scale = 1, dpi = 100)
  print(paste0("Done with ", cc," time-point ", ms,"/", nrow(my_vars)))
  return(tab_out)
}

## Helper function
my_idx = my_vars %>%
  mutate(idx = 1:n())

my_fun = function(x){
  test_idx = my_idx %>% 
    filter(v2 == paste0(x)) %>% 
    pull(idx)
  # test_idx = test_idx[150:170] ## Uncomment only for testing function
  mv_out = data.table::rbindlist(parallel::mclapply(seq_along(my_vars$v1)[test_idx], mod_fun, mc.cores = ncores), fill = TRUE)
  ## Write out the criteria table
  data.table::fwrite(mv_out,
                     paste0(out_dir,x, "_out.txt"),
                     row.names = FALSE, col.names = TRUE, sep = "\t")
}

## Run the helper function
lapply(env_vars, my_fun)

## Commands to turn individual plots into videos. Requires ffmpeg to be installed.
my_ind_paths = paste0(tmp_dir,env_vars,"/*.png")
my_commands = paste0('cat $(ls -1 ', my_ind_paths, ' | sort -V) | ffmpeg -y -i - -filter_complex "setpts=PTS*2" ', out_dir, env_vars,'.mp4')
for (i in seq_along(my_commands)) {
  system(my_commands[i])  
}

print(paste0("DONE!"))