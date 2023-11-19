library(tidyverse)
library(data.table)
library(ggpubr)
library(lubridate)
library(ggthemes)
rm(list=ls())

# Set needed directory paths ----------------------------------------------
data_dir = "./data/"
res_dir = "./results/"
out_dir = res_dir

data_file = "example_sensor_data.txt"
map_file = "XY_to_real_cm_space.txt"
data_interpolation_file = "example_sensor_data_interpolation_best_per_time_out.txt"

# model_results_dir = "~/projects/expose/results/sensors/"


# Plotting ----------------------------------------------------------------

## Read in predicted values by best model
d1 = fread(paste0(data_dir, data_interpolation_file)) %>% 
  select(model, date, sensor, rmse) %>% 
  distinct()

## Set sensor name and units
my_units = data.frame(sensor = c("par","tair","hr","vpd"),
                      unit = c("PPFD (µmol m-2 s-1)",
                               "Air temperature (°C)",
                               "Air humidity (%)",
                               "Air VPD (kPa)"))

## Keep necessary information and make names nice for plotting
my_df = d1 %>% 
  mutate(model = ifelse(model=="Thin_plate_spline", "Thin-plate_spline", model)) %>% 
  mutate(model = ifelse(model=="Inverse_dist_weighting", "Inverse-dist_weighting", model)) %>% 
  mutate(model = gsub("_","\n", model)) %>% 
  mutate(model = factor(model, levels = c("Nearest\nneighbor", "Inverse-dist\nweighting","Thin-plate\nspline", "GAM\nGaussian","GAM\nDuchon", "GAM\nPenalized"))) %>% 
  left_join(., my_units, by = c("sensor"))
my_dates = c(min(my_df$date), max(my_df$date))

## Make plots
my_theme = theme(
  axis.text.x = element_text(size = 10,face = "bold", color = "black"),
  axis.text.y = element_text(size = 10,face = "bold", color = "black"),
  axis.title.y = element_text(size = 10, color = "black", face = "bold"),
  axis.title.x = element_text(size = 10, color = "black", face = "bold"),
  # axis.ticks.y = element_blank(),
  strip.text.x = element_text(size = 10, color = "black", face = "bold"),
  strip.text.y = element_text(size = 10, color = "black", face = "bold"),
  plot.title = element_text(hjust = 0.5,color = "black", face = "bold"), 
  legend.title = element_text(size = 10, color = "black", face = "bold"),
  legend.text = element_text(size = 10, color = "black", face = "bold"))

### Black and white
my_plt = my_df %>% 
  mutate(bm = "") %>% 
  ggplot(.)+
  geom_tile(aes(x = date, y = model), color = "black")+
  scale_x_datetime(labels = date_format("%Y-%m-%d %H"),
                   date_breaks = "6 hours", 
                   date_labels = "%H", 
                   expand = c(0.01, 0.01))+
  facet_wrap(~unit)+
  theme_classic()+
  my_theme+
  labs(x = "Hour", y = "")+
  theme(legend.position="top",
        legend.title = element_blank())
ggsave(my_plt,
       file = paste0(out_dir,"model_comparison_env_sensors_vs_time.png"), width = 16.25, height = 7.75, units = "in", bg = "white", scale = 1, dpi = 400)

### Color single line
my_plt_c = my_df %>% 
  mutate(bm = "") %>% 
  ggplot(.)+
  geom_tile(aes(x = date, y = bm, fill = model))+ ggsci::scale_fill_lancet()+
  scale_x_datetime(labels = date_format("%Y-%m-%d %H"),
                   date_breaks = "6 hours", 
                   date_labels = "%H", 
                   expand = c(0.01, 0.01))+
  facet_wrap(~unit)+
  theme_classic()+
  my_theme+
  labs(x = "Date", y = "")+
  theme(legend.position="top",
        legend.title = element_blank())
ggsave(my_plt_c,
       file = paste0(out_dir,"model_comparison_env_sensors_vs_time_color.png"), width = 8.25, height = 3.75, units = "in", bg = "white", scale = 1.5)

## Plot showing the proportion of all times that a model is chosen
mod_used = d1 %>% 
  mutate(mod_date_sensor = paste0(model,date,sensor)) %>% 
  pull(mod_date_sensor)

my_files = list.files(path = res_dir, pattern = paste0("_out.txt$"))
my_d = rbindlist(lapply(paste0(res_dir, my_files), fread), use.names = TRUE) %>% 
  mutate(mod_date_sensor = paste0(model,date,sensor)) %>% 
  mutate(is_used = ifelse(mod_date_sensor%in%mod_used, 1, 0))

## For testing
# my_d = my_d %>%
#   filter(date%in%d1$date)

## Summarize model RMSE across all time points and the proportion of times a model is chosen
my_d = my_d %>% 
  left_join(., my_units, by = c("sensor")) %>% 
  mutate(model = ifelse(model=="Thin_plate_spline", "Thin-plate_spline", model)) %>% 
  mutate(model = ifelse(model=="Inverse_dist_weighting", "Inverse-dist_weighting", model)) %>% 
  mutate(model = gsub("_","\n", model)) %>% 
  mutate(model = factor(model, levels = c("Nearest\nneighbor", "Inverse-dist\nweighting", "Thin-plate\nspline", "GAM\nGaussian","GAM\nDuchon", "GAM\nPenalized")))

p1 = my_d %>% 
  ggplot(., aes(model, rmse))+
  geom_violin()+
  facet_wrap(~unit, scales = "free_y") + 
  scale_y_continuous(breaks = scales::pretty_breaks(5), 
                     limits = c(0, NA),
                     expand = expansion(mult = c(0.05, 0.1)))+
  theme_linedraw()+
  my_theme +
  theme(strip.text.x = element_text(size = 10, color = "white", face = "bold"),
        strip.text.y = element_text(size = 10, color = "white", face = "bold"))+
  labs(x = "", y = "LOOCV-RMSE\n ")+
  stat_summary(fun=mean, geom="point", shape=18, size=2, color = "red")

p2 = my_d %>% 
  group_by(model, sensor, unit) %>% 
  summarise(my_f = signif(sum(is_used)/n(),2)) %>% 
  ggplot(., aes(x=model, y=my_f)) +
  geom_bar(stat="identity", fill=NA, color = "black", width = 0.65)+
  geom_text(aes(label=my_f, size=my_f), vjust=1.3, color="black", show.legend = FALSE)+
  facet_wrap(~unit, scales = "free_y") + 
  scale_y_continuous(breaks = scales::pretty_breaks(5), limits = c(0, NA),
                     expand = expansion(mult = c(0.15, 0.025))) +
  theme_classic2()+
  my_theme + 
  theme(strip.text.x = element_text(size = 10, color = "black", face = "bold"),
        strip.text.y = element_text(size = 10, color = "black", face = "bold"))+
  labs(x = "", y = "Proportion of times a model is chosen\n ")

plt_out = cowplot::plot_grid(p1,p2, ncol = 1, labels = "AUTO")
ggsave(plt_out,
       file = paste0(out_dir,"model_comparison_env_sensors.png"),
       width = 11.25, height = 8.5, units = "in", bg = "white", scale = 1.2)


print(paste0("DONE!"))