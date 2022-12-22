## Script for gap analysis of sites with stable water temperatures
## George Valentine
## Oct. 2022

## Steps:
  # 1. Calculate slope at all BKT habitat sites in SE US
  # 2. Filter for 25% lowest slope (most stable temps)
  # 3. Check % overlap with shapefile of conserved lands

# Load packages
library(tidyverse)


## Import data

# Load stream segments that are within protected areas
PA_StreamSegments <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/GIS/PA_flowlines.csv")

# Load PCA coordinates (already calculated for segments in BKT habitat)
pca_coords <- fread("pca_coords_v4.csv")

# Load temp model v4 posterior parameter estimates
StreamTemp_PCA_LM_full_v4_params <- fread("StreamTemp_PCA_LM_full_v4_params.csv")

##
# Create an empty dataframe to store slopes
BKT_hab_betas <- data.frame(COMID = pca_coords$COMID,
                            beta = numeric(length = length(pca_coords$COMID)))

# Save thetas from temperature model
temp_model_thetas <- StreamTemp_PCA_LM_full_v4_params[507:512,2]

# Save sd.beta from temperature model
temp_model_sd.beta <- StreamTemp_PCA_LM_full_v4_params[505,2]

# Save sd from temperature model
temp_model_sd <- StreamTemp_PCA_LM_full_v4_params[514,2]

for (i in 1:nrow(BKT_hab_betas)) {
  # print COMID
  print(BKT_hab_betas[i,1])
  
  # calculate slope at the current site
  mu.beta.i <- temp_model_thetas[1] + temp_model_thetas[2] * pca_coords[i,2] + temp_model_thetas[3] * pca_coords[i,3] + temp_model_thetas[4] * pca_coords[i,4] + temp_model_thetas[5] * pca_coords[i,5] + temp_model_thetas[6] * pca_coords[i,6]
  #print(as.numeric(mu.beta.i))
  
  # save slope
  BKT_hab_betas[i,2] <- rnorm(1, mean = as.numeric(mu.beta.i), sd = as.numeric(temp_model_sd.beta))
}

## WHY ARE THEY >1??

# make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = BKT_hab_betas, 
             aes(x = Long, y = Lat, color = beta)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Predicted Air-Water Temperature Slopes \nat SE Trout Sites",
       color = expression(beta)) +
  scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
  theme_classic()

# Filter for the lowest 25% of slopes
BKT_hab_betas_lowest <- BKT_hab_betas %>% 
  mutate(Percentile = ntile(beta, 100)) %>% 
  filter(Percentile <= 25)

# the best BKT sites in protected areas:
BKT_sites_PAs <- BKT_hab_betas_lowest %>% 
  filter(COMID %in% PA_StreamSegments$COMID)

# % best BKT sites in protected areas:
nrow(BKT_sites_PAs)/nrow(BKT_hab_betas_lowest)
