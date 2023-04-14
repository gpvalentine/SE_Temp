## Stream Temp PCA analysis
# Changes from v4:
  # - get ride of versioning in switch over to git
  # - run analysis on weekly summary data instead of daily summary data
  # - Compare nonlinear model (logistic) to linear model

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")

# Load libraries
library(data.table)
library(devtools)
library(pcaMethods)
install_github("vqv/ggbiplot")
library(ggbiplot)     # for PCA plots
library(factoextra)   # for extracting PCA components
library(jagsUI)
library(tidyverse)
library(MCMCvis)
library(here)
library(AICcmodavg)   # For computing model DICs
library(corrplot)

set.seed(1234)

# Load and format the data
# Load EPA Streamcat data for SE US
StreamCat_Covars <- fread("Data/SE_Temp_StreamCat_Covars_Combined.csv")

# Load national NHDplus data
NHDplus_data <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")

# Combine StreamCat and NHDplus data
# Left join national NHD data to SE StreamCat data so that everything is for the SE US
StreamSegment_covars <- StreamCat_Covars %>% 
  left_join(NHDplus_data)

# Load NS204 temperature monitoring sites
NS204_Sites <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/SE_Temp/Data/NS204_sites.csv")

# Load stream segments that are within BKT habitat patches as defined by the EBTJV
BKT_Habitat_Patch_StreamSegments <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/GIS/NHDplus_Flowlines_BKT_Patches.csv")

# A few segments lie in multiple BKT habitat patches and so have been split into multiple rows
# Filter so that we just have one observation of each segment
BKT_Habitat_Patch_StreamSegments <- BKT_Habitat_Patch_StreamSegments %>% 
  distinct(COMID, .keep_all = T)

# Filter covars to just stream segments in BKT habitat patches and the few temperature monitoring sites not in EBBTJV habitat patches
# We lose a few stream segments in BKT habitat patches because they aren't in the StreamCat dataset. Braided streams and errors in NHDplus that were filtered out of StreamCat
BKT_Habitat_StreamSegment_covars <- StreamSegment_covars %>% 
  filter(COMID %in% c(BKT_Habitat_Patch_StreamSegments$COMID, NS204_Sites$COMID)) %>% 
  left_join(NS204_Sites[,c("COMID", "Lat", "Long")]) %>%  # also join in the lat and long of the sites
  left_join(BKT_Habitat_Patch_StreamSegments[,c("COMID", "Out_Lat", "Out_Long")]) %>% # And join in the coordinates of the outlet of the stream segment
  mutate(Lat = ifelse(!is.na(Lat), Lat, Out_Lat),
         Long = ifelse(!is.na(Long), Long, Out_Long)) %>% 
  dplyr::select(-Out_Long,
                -Out_Lat)

# Replace -9998 (the nodata value) with NA
BKT_Habitat_StreamSegment_covars <- replace(BKT_Habitat_StreamSegment_covars, BKT_Habitat_StreamSegment_covars == -9998, NA)

# Remove columns that have all zeros or all the same value: they cannot be scaled by the PCA and won't be useful as predictors
BKT_Habitat_StreamSegment_covars2 <- BKT_Habitat_StreamSegment_covars %>% 
  filter(StreamOrde < 5) %>% # this filters out "rivers" that were getting super low predicted slopes. See "Script to troubleshoot predicted slope outliers.R"
  #filter(!DamDensCat > 0) %>%  # and filter out sites with dams in their catchments. Dams mimic the effects of groundwater and obscure the true buffering
  select_if(~n_distinct(.) > 2) %>% 
  select_if(is.numeric) %>% # retain only numeric columns
  dplyr::arrange(COMID) %>% 
  dplyr::select(-fid,# and remove several columns with GIS IDs that mess up the PCA
                -OBJECTID,
                -GNIS_ID,
                -REACHCODE,
                -WBAREACOMI,
                -FCODE,
                -Shape_Length,
                -FromNode,
                -ToNode,
                -Hydroseq,
                -LevelPathI,
                -Pathlength,
                -TerminalPa,
                -UpLevelPat,
                -UpHydroseq,
                -DnLevelPat,
                -DnMinorHyd,
                -DnHydroseq,
                -FromMeas,
                -ToMeas)

# find the percent missing data in the dataframe
(sum(is.na(BKT_Habitat_StreamSegment_covars2))/prod(dim(BKT_Habitat_StreamSegment_covars2)))*100

# Center and scale the data
BKT_Habitat_StreamSegment_covars3 <- as.matrix(scale(BKT_Habitat_StreamSegment_covars2[,-1]))


# Bayesian pca can take NA values
pca <- pcaMethods::bpca(BKT_Habitat_StreamSegment_covars3[,-1], # leave out first column, which is COMID
            nPcs = 40)

pca_summary.table <- as.data.frame(summary(pca))[,1:5]

# plot PCA
# biplot(pca,
#          ellipse = T) +
#   theme_classic()

# Save PCA loadings
pca.loadings <- as.data.frame(loadings(pca))

# ten greatest loadings for first five PCs
top.loadings <- data.frame(matrix(NA, nrow = 10, ncol = 0))

for (i in 1:5) {
  top.loadings.i <- pca.loadings %>% 
    dplyr::select(i) %>%
    dplyr::arrange(desc(abs(.))) %>% 
    head(10) %>% 
    rownames_to_column(var = "Covariate")
  
  top.loadings <- cbind(top.loadings, top.loadings.i)
}

# make a table of the top loadings to include in the manuscript
top.loadings.table <- data.frame(rbind(round(pca_summary.table["R2",], 3)*100,
                                       data.frame(PC1 = paste(as.list.data.frame(top.loadings[,1]), collapse = ", "),
                                                  PC2 = paste(as.list.data.frame(top.loadings[,3]), collapse = ", "),
                                                  PC3 = paste(as.list.data.frame(top.loadings[,5]), collapse = ", "),
                                                  PC4 = paste(as.list.data.frame(top.loadings[,7]), collapse = ", "),
                                                  PC5 = paste(as.list.data.frame(top.loadings[,9]), collapse = ", "))))

rownames(top.loadings.table) <- c("$R^2$", "Variables")
#rownames(top.loadings.table) <- NULL

# extract PCA coords for individuals
pca_scores <- scores(pca)

# and subset to just the first five 5
pca_scores <- as.data.frame(pca_scores[,1:5])

# Rename columns for easier interpretation
names(pca_scores) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

pca_scores <- pca_scores %>% 
  cbind(COMID = BKT_Habitat_StreamSegment_covars2$COMID) %>% 
  dplyr::select(COMID, PC1, PC2, PC3, PC4, PC5)
# save PCA scores
#fwrite(pca_coords, "pca_coords_v4.csv")

####################################
# Max stream temperature model

# Load NS204 temperature data
NS204_temps_daily <- fread("Data/NS204_temps_daily.csv")
NS204_temps_weekly <- fread("Data/NS204_temps_weekly.csv")

# add COMIDs to the NS204 site data
# I use the file originally used to get flow data - I just need a file that relates SiteIDs to COMIDs
NS204_SiteID_to_COMID <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/GIS/NS204_Sites_wNHDplusCOMID_EXPORT.csv") %>% 
  dplyr::select(SiteID, COMID)

NS204_temps_daily <- NS204_temps_daily %>% 
  left_join(NS204_SiteID_to_COMID)

NS204_temps_weekly <- NS204_temps_weekly %>% 
  left_join(NS204_SiteID_to_COMID)

NS204_temps_daily_filtered <- NS204_temps_daily %>% 
  filter(COMID %in% pca_scores$COMID, # Filter temp data to the COMIDs for which we have principle components
         !is.na(AirTemp_c_MEAN)) %>% 
  arrange(COMID) %>% 
  mutate(SegmentNo = match(COMID, unique(COMID)))  # Add a SegmentNo column to the filtered data

NS204_temps_weekly_filtered <- NS204_temps_weekly %>% 
  filter(COMID %in% pca_scores$COMID, # Filter temp data to the COMIDs for which we have principle components
         !is.na(AirTemp_c_MEAN)) %>% 
  arrange(COMID) %>% 
  mutate(SegmentNo = match(COMID, unique(COMID)))  # Add a SegmentNo column to the filtered data

# Plot the data from a random segment
NS204_temps_daily_filtered %>% 
  filter(COMID %in% sample(unique(COMID), 1)) %>% 
  ggplot() +
  geom_point(aes(x = AirTemp_c_MEAN,
                y = WaterTemp_c_MEAN)) +
  theme_classic()

NS204_temps_weekly_filtered %>% 
  filter(COMID %in% sample(unique(COMID), 1)) %>% 
  ggplot() +
  geom_point(aes(x = AirTemp_c_MEAN,
                 y = WaterTemp_c_MEAN)) +
  theme_classic()

# for the linear model, we exclude observations where the mean daily air temperature is below 0 degrees
NS204_temps_daily_filtered_LM <- NS204_temps_daily_filtered %>% 
  filter(AirTemp_c_MEAN >= 0)

NS204_temps_weekly_filtered_LM <- NS204_temps_weekly_filtered %>% 
  filter(AirTemp_c_MEAN >= 0)

# Filter pca data to just the NS204 sites
NS204_PCA_scores <- pca_scores %>% 
  filter(COMID %in% NS204_temps_daily_filtered$COMID) %>% 
  arrange(COMID)
  

nCOMIDs <- nrow(NS204_PCA_scores)
nObs.daily <- nrow(NS204_temps_daily_filtered)
nObs.weekly <- nrow(NS204_temps_weekly_filtered)
nObs.daily_LM <- nrow(NS204_temps_daily_filtered_LM)
nObs.weekly_LM <- nrow(NS204_temps_weekly_filtered_LM)

###################################################################

## Daily Mean Temperatures ##
# Write model
sink("Analysis/JAGS_Files/SE_Temp_Daily_Mean_StreamTemp_LM_Full.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:6){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.beta - precision for beta
  # sd parameter for tau.beta
  sd.beta ~ dunif(0, 10)
    
  # tau.beta - precision parameter for all betas
  tau.beta <- 1/sd.beta^2

  for (i in 1:nCOMIDs){

    # alpha.i - intercept for likelihood
    alpha[i] ~ dnorm(0, 0.001)

    # mu.beta.i - mean parameter for dnorm
    mu.beta[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # Beta.i - coefficient for daily air temperature
    beta[i] ~ dnorm(mu.beta[i], tau.beta)
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.daily){
    WaterTemp_Mean_Pred[n] <- alpha[SegmentNo[n]] + beta[SegmentNo[n]] * AirTemp_Mean_Obs[n]
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Posterior predictive checks
  # for (n in 1:nObs.daily){
  #   PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
  # }
  
  # Calculate RMSE
  for (n in 1:nObs.daily){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.daily)
  
    # Calculate Bayesian R^2
  for (n in 1:nObs.daily){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.daily - 1)
  var.pred <- sum(x.pred)/(nObs.daily - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.daily = nObs.daily_LM,
                  NS204_PCA_scores = NS204_PCA_scores,
                  AirTemp_Mean_Obs = NS204_temps_daily_filtered_LM$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_daily_filtered_LM$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_daily_filtered_LM$SegmentNo)

# Set parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "sd.beta", "tau.beta", "theta", "tau", "sd", "RMSE", "R2")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

set.seed(1234)

# Fit model
Daily_Mean_StreamTemp_LM_Full <- jagsUI::jags(data = jags_data,
                             parameters.to.save = jags_params,
                             model.file = "Analysis/JAGS_Files/SE_Temp_Daily_Mean_StreamTemp_LM_Full.jags",
                             n.chains = nc,
                             n.iter = ni,
                             n.burnin = nb, 
                             n.thin = nt,
                             parallel = T)

# Save model summary
Daily_Mean_StreamTemp_LM_Full_params <- MCMCsummary(Daily_Mean_StreamTemp_LM_Full, HPD = T)
# Get DIC
Daily_LM_DIC.val <- DIC(Daily_Mean_StreamTemp_LM_Full)

# ## Posterior predictive check
# Daily_Mean_LM_PPC_data <- data.frame(SegmentNo = NS204_temps_daily_filtered_LM$SegmentNo,
#                                      PPC = MCMCsummary(Daily_Mean_StreamTemp_LM_Full, HPD = T, params = "PPC"))
# 
# Daily_Mean_LM_PPCs <- Daily_Mean_LM_PPC_data %>% 
#   group_by(SegmentNo) %>% 
#   dplyr::summarise(PPC_mean = mean(PPC.mean),
#                    PPC_sd = sd(PPC.mean))
# 
# Daily_Mean_LM_PPC_summ.table <- as.data.frame(apply(Daily_Mean_LM_PPCs[,2:3],2,summary))

# plot to visualize
# plot_data <- Daily_Mean_StreamTemp_LM_Full_params %>% 
#   rownames_to_column("param") %>% 
#   filter(str_detect(param, "WaterTemp_Mean_New")) %>% 
#   select(mean) %>% 
#   cbind(AirTemp_c_Mean = NS204_temps_daily_filtered$AirTemp_c_MEAN,
#         WaterTemp_c_Mean = NS204_temps_daily_filtered$WaterTemp_c_MEAN,
#         COMID = NS204_temps_daily_filtered$COMID)
# 
# plot_data %>% 
#   filter(COMID %in% sample(unique(COMID), 1)) %>% 
#   ggplot() +
#   geom_point(aes(x = AirTemp_c_Mean,
#                  y = WaterTemp_c_Mean)) +
#   geom_line(aes(x = AirTemp_c_Mean,
#                 y = mean),
#             color = "red")

## Weekly Mean Temperatures ##
# Write model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_LM_Full.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:6){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.beta - precision for beta
  # sd parameter for tau.beta
  sd.beta ~ dunif(0, 10)
    
  # tau.beta - precision parameter for all betas
  tau.beta <- 1/sd.beta^2

  for (i in 1:nCOMIDs){

    # alpha.i - intercept for likelihood
    alpha[i] ~ dnorm(0, 0.001)

    # mu.beta.i - mean parameter for dnorm
    mu.beta[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # Beta.i - coefficient for weekly air temperature
    beta[i] ~ dnorm(mu.beta[i], tau.beta)
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Mean_Pred[n] <- alpha[SegmentNo[n]] + beta[SegmentNo[n]] * AirTemp_Mean_Obs[n]
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Posterior predictive checks
  for (n in 1:nObs.weekly){
    PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
  }
  
  # Calculate RMSE
  for (n in 1:nObs.weekly){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.weekly)
  
  # Calculate Bayesian R^2
  for (n in 1:nObs.weekly){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.weekly - 1)
  var.pred <- sum(x.pred)/(nObs.weekly - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly_LM,
                  NS204_PCA_scores = NS204_PCA_scores,
                  AirTemp_Mean_Obs = NS204_temps_weekly_filtered_LM$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_weekly_filtered_LM$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_weekly_filtered_LM$SegmentNo)

# Set parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "sd.beta", "tau.beta", "theta", "tau", "sd", "RMSE", "R2", "PPC")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

set.seed(1234)

# Fit model
Weekly_Mean_StreamTemp_LM_Full <- jagsUI::jags(data = jags_data,
                                             parameters.to.save = jags_params,
                                             model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_LM_Full.jags",
                                             n.chains = nc,
                                             n.iter = ni,
                                             n.burnin = nb,
                                             n.thin = nt,
                                             parallel = T)

# Save model output
Weekly_Mean_StreamTemp_LM_Full_params <- MCMCsummary(Weekly_Mean_StreamTemp_LM_Full, HPD = T, excl = "PPC")
# Get DIC
Weekly_LM_DIC.val <- DIC(Weekly_Mean_StreamTemp_LM_Full)

## Posterior predictive check
Weekly_Mean_LM_PPC_data <- data.frame(SegmentNo = NS204_temps_weekly_filtered_LM$SegmentNo,
                                           PPC = MCMCsummary(Weekly_Mean_StreamTemp_LM_Full, HPD = T, params = "PPC"))

Weekly_Mean_LM_PPCs <- Weekly_Mean_LM_PPC_data %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarise(PPC_mean = mean(PPC.mean),
                   PPC_sd = sd(PPC.mean))

Weekly_Mean_LM_PPC_summ.table <- as.data.frame(apply(Weekly_Mean_LM_PPCs[,2:3],2,summary))

## Plot fit over data to visualize
# Extract intercepts and slopes
Weekly_LM_intercepts <- MCMCsummary(Weekly_Mean_StreamTemp_LM_Full, params = "alpha") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)
  
# Extract slopes
Weekly_LM_slopes <- MCMCsummary(Weekly_Mean_StreamTemp_LM_Full, params = "beta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

## fit
# select a random segment
sample_COMID <- NS204_PCA_scores %>% filter(COMID %in% sample(unique(COMID), 1)) %>% .[,"COMID"]

# Create an empty dataframe to store predicted watertemp data
watertemp_pred_weekly_LM <- data.frame(WaterTemp_Mean_pred = numeric())

# predict water temp at the random segment using the parameters from the model for that segment
for (i in 1:nrow(NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID,])) {
  airtemp_obs <- NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID, "AirTemp_c_MEAN"]
  watertemp_pred_weekly_LM[i,1] <- Weekly_LM_intercepts[Weekly_LM_intercepts$COMID == sample_COMID, "mean"] + Weekly_LM_slopes[Weekly_LM_slopes$COMID == sample_COMID, "mean"] * airtemp_obs[i]
}

# plot
ggplot() +
  geom_point(aes(x = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$AirTemp_c_MEAN,
                 y = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$WaterTemp_c_MEAN)) +
  geom_line(aes(x = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$AirTemp_c_MEAN,
                y = watertemp_pred_weekly_LM$WaterTemp_Mean_pred),
            color = "red") + 
  labs(x = "Air Temp (c)",
       y = "Air Temp(c)") +
  theme_classic()


#########################
### Logistic Model ###
#########################

## following Mohseni et al. (1998), the following nonlinear logistic regression model is used:
# T(s) = epsilon + ((zeta - epsilon)/(1 + e^phi(kappa-T(a))))

# Where: 
#   T(s) = estimated stream temp (C)
#   T(a) = measured air temp (C)
#   epsilon = estimated minimum stream temp (C)
#   zeta = estimated Meanimum stream temp (C)
#   phi = measure of steepest slope of the function (C^-1)
#   beta = steepest slope of the function (C)
#   kappa = Air temp (C) at the function's inflection point

## Daily Mean Temperatures ##
# Create a dataframe of the min and Mean water temps at each site
min_max_waterTemps <- NS204_temps_daily_filtered %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarize(min_waterTemp = min(WaterTemp_c_MEAN, na.rm = T),
            max_waterTemp = max(WaterTemp_c_MEAN, na.rm = T))

# min_max_waterTemps <- c(median(min_max_waterTemps$min_waterTemp),
#                         median(min_max_waterTemps$max_waterTemp))
  

# Write model
sink("Analysis/JAGS_Files/SE_Temp_Daily_Mean_StreamTemp_Logistic_Full.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:6){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.phi - precision for all phis
  sd.phi ~ dunif(0, 10)
  tau.phi <- 1/(sd.phi^2)
  s2.phi <- sd.phi^2
  
  # epsilon_i - estimated median minimum stream temp 
  #epsilon ~ dnorm(min_waterTemp, 1/0.01) # uses a normal distribution informed by the measured minimum temperature
  
  # zeta_i - estimated median maximum stream temp
  #zeta ~ dnorm(max_waterTemp, 1/0.01)

  for (i in 1:nCOMIDs){
  
    # epsilon_i - estimated minimum stream temp at segment i 
    epsilon[i] ~ dnorm(min_waterTemps[i], 1/0.01) # uses a normal distribution informed by the measured minimum temperature
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dnorm(max_waterTemps[i], 1/0.01)
  
    # mu.phi_i - mean parameter for phi_i
    mu.phi[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # phi_i - measure of the steepest slope of the function at segment i
    phi[i] ~ dnorm(mu.phi[i], tau.phi)
    
    # kappa_i - Air temp at the function's inflection point for segment i
    kappa[i] ~ dnorm(20, 0.01)
    
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.daily){
    WaterTemp_Mean_Pred[n] <- epsilon[SegmentNo[n]] + ((zeta[SegmentNo[n]] - epsilon[SegmentNo[n]])/(1 + exp(phi[SegmentNo[n]] * (kappa[SegmentNo[n]] - AirTemp_Mean_Obs[n]))))
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    # New data for PPCs
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Relation to calculate phi from beta (slope), zeta, and epsilon
  for (i in 1:nCOMIDs){
    beta[i] <- (phi[i]*(zeta[i] - epsilon[i]))/4
  }
  
  # Posterior predictive checks
  # for (n in 1:nObs.daily){
  # PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_Pred[n])
  # }
  
  # Calculate RMSE
  for (n in 1:nObs.daily){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.daily)
  
  # Calculate Bayesian R^2
  for (n in 1:nObs.daily){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.daily - 1)
  var.pred <- sum(x.pred)/(nObs.daily - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.daily = nObs.daily,
                  NS204_PCA_scores = NS204_PCA_scores,
                  #min_waterTemp = min_max_waterTemps[1],
                  #Mean_waterTemp = min_max_waterTemps[2],
                  min_waterTemps = min_max_waterTemps$min_waterTemp,
                  max_waterTemps = min_max_waterTemps$max_waterTemp,                  
                  AirTemp_Mean_Obs = NS204_temps_daily_filtered$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_daily_filtered$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_daily_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("theta", "epsilon", "zeta", "mu.phi", "sd.phi", "s2.phi", "phi", "beta", "kappa", "sd", "RMSE", "R2")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

set.seed(1234)

# Fit model
Daily_Mean_StreamTemp_Logistic_Full <- jagsUI::jags(data = jags_data,
                                                    parameters.to.save = jags_params,
                                                    model.file = "Analysis/JAGS_Files/SE_Temp_Daily_Mean_StreamTemp_Logistic_Full.jags",
                                                    n.chains = nc,
                                                    n.iter = ni,
                                                    n.burnin = nb,
                                                    n.thin = nt,
                                                    parallel = T)

# Save model output
Daily_Mean_StreamTemp_Logistic_Full_Params <- MCMCsummary(Daily_Mean_StreamTemp_Logistic_Full, HPD = T)
# Get DIC
Daily_Logistic_DIC.val <- DIC(Daily_Mean_StreamTemp_Logistic_Full)

## Posterior predictive check
# Daily_Mean_Logistic_PPC_data <- data.frame(SegmentNo = NS204_temps_daily_filtered$SegmentNo,
#                                            PPC = MCMCsummary(Daily_Mean_StreamTemp_Logistic_Full, HPD = T, params = "PPC"))
# 
# Daily_Mean_Logistic_PPCs <- Daily_Mean_Logistic_PPC_data %>%
#   group_by(SegmentNo) %>%
#   dplyr::summarise(PPC_mean = mean(PPC.mean),
#                    PPC_sd = sd(PPC.mean))
# 
# Daily_Mean_Logistic_PPC_summ.table <- as.data.frame(apply(Daily_Mean_Logistic_PPCs[,2:3],2,summary))

## Weekly Mean Temperatures ##
# Create a dataframe of the min and max water temps at each site
min_max_waterTemps <- NS204_temps_weekly_filtered %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarize(min_waterTemp = min(WaterTemp_c_MEAN, na.rm = T),
            max_waterTemp = max(WaterTemp_c_MEAN, na.rm = T),
            COMID = first(COMID)) %>% 
  left_join(NS204_Sites)

fwrite(min_max_waterTemps, "Data/weekly_min_max_mean_watertemps.csv", row.names = F)

# plot min and Mean temps in space
US_states <- map_data("state")

# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = min_max_waterTemps,
#              aes(x = Long, y = Lat, color = min_waterTemp)) +
#   coord_map("bonne",
#             lat0 = 40,
#             xlim = c(-84.2, -75.5),
#             ylim = c(34.5, 40.5)) +
#   labs(x = "Long",
#        y = "Lat") +
#   theme_classic()
# 
# ggplot() +
#   geom_polygon(data = US_states,
#                aes(x = long, y = lat, group = group),
#                color = "black", fill = NA) +
#   geom_point(data = min_max_waterTemps,
#              aes(x = Long, y = Lat, color = Mean_waterTemp)) +
#   coord_map("bonne",
#             lat0 = 40,
#             xlim = c(-84.2, -75.5),
#             ylim = c(34.5, 40.5)) +
#   labs(x = "Long",
#        y = "Lat") +
#   theme_classic()

# min_max_waterTemps <- c(median(min_max_waterTemps$min_waterTemp),
#                         median(min_max_waterTemps$max_waterTemp))


# Write model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_Logistic_Full.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:6){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.phi - precision for all phis
  sd.phi ~ dunif(0, 10)
  tau.phi <- 1/(sd.phi^2)
  s2.phi <- sd.phi^2
  
  # epsilon_i - estimated median minimum stream temp 
  #epsilon ~ dnorm(min_waterTemp, 1/0.01) # uses a normal distribution informed by the measured minimum temperature
  
  # zeta_i - estimated median maximum stream temp
  #zeta ~ dnorm(max_waterTemp, 1/0.01)

  for (i in 1:nCOMIDs){

    # epsilon_i - estimated minimum stream temp at segment i 
    epsilon[i] ~ dnorm(min_waterTemps[i], 1/0.01) # uses a normal distribution informed by the measured minimum temperature
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dnorm(max_waterTemps[i], 1/0.01)
    
    # mu.phi_i - mean parameter for phi_i
    mu.phi[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # phi_i - measure of the steepest slope of the function at segment i
    phi[i] ~ dnorm(mu.phi[i], tau.phi)
    
    # kappa_i - Air temp at the function's inflection point for segment i
    kappa[i] ~ dnorm(20, 0.01)
    
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Mean_Pred[n] <- epsilon[SegmentNo[n]] + ((zeta[SegmentNo[n]] - epsilon[SegmentNo[n]])/(1 + exp(phi[SegmentNo[n]] * (kappa[SegmentNo[n]] - AirTemp_Mean_Obs[n]))))
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    # New data for PPCs
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Relation to calculate phi from beta (slope), zeta, and epsilon
  for (i in 1:nCOMIDs){
    beta[i] <- (phi[i]*(zeta[i] - epsilon[i]))/4
  }
  
  # Posterior predictive checks
  for (n in 1:nObs.weekly){
    PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
  }
  
  # Calculate RMSE
  for (n in 1:nObs.weekly){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.weekly)
  
  # Calculate Bayesian R^2
  for (n in 1:nObs.weekly){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.weekly - 1)
  var.pred <- sum(x.pred)/(nObs.weekly - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly,
                  NS204_PCA_scores = NS204_PCA_scores,
                  #min_waterTemp = min_max_waterTemps[1],
                  #max_waterTemp = min_max_waterTemps[2],
                  min_waterTemps = min_max_waterTemps$min_waterTemp,
                  max_waterTemps = min_max_waterTemps$max_waterTemp,
                  AirTemp_Mean_Obs = NS204_temps_weekly_filtered$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_weekly_filtered$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_weekly_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("theta", "epsilon", "zeta", "mu.phi", "sd.phi", "s2.phi", "phi", "beta", "kappa", "sd", "PPC", "RMSE", "R2")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

set.seed(1234)

# Fit model
Weekly_Mean_StreamTemp_Logistic_Full <- jagsUI::jags(data = jags_data,
                                              parameters.to.save = jags_params,
                                              model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_Logistic_Full.jags",
                                              n.chains = nc,
                                              n.iter = ni,
                                              n.burnin = nb,
                                              n.thin = nt,
                                              parallel = T)

# Save model output
Weekly_Mean_StreamTemp_Logistic_Full_Params <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, HPD = T,
                                                          excl = c("PPC"))
# Get DIC
Weekly_Logistic_DIC.val <- DIC(Weekly_Mean_StreamTemp_Logistic_Full)

#MCMCtrace(Weekly_Mean_StreamTemp_Logistic_Full, params = "zeta", pdf = F)

## Posterior predictive check
Weekly_Mean_Logistic_PPC_data <- data.frame(SegmentNo = NS204_temps_weekly_filtered$SegmentNo,
                                       PPC = MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, HPD = T, params = "PPC"))

Weekly_Mean_Logistic_PPCs <- Weekly_Mean_Logistic_PPC_data %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarise(PPC_mean = mean(PPC.mean),
                   PPC_sd = sd(PPC.mean))

Weekly_Mean_Logistic_PPC_summ.table <- as.data.frame(apply(Weekly_Mean_Logistic_PPCs[,2:3],2,summary))

## Plot fit over data to visualize
# Extract parameters
Weekly_Logistic_epsilons <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "epsilon") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_zetas <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "zeta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_phis <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "phi") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_kappas <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "kappa") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_betas <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "beta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

## fit
# select a random segment
sample_COMID <- NS204_PCA_scores %>% filter(COMID %in% sample(unique(COMID), 1)) %>% .[,"COMID"]

# Create an empty dataframe to store predicted watertemp data
watertemp_pred_weekly <- data.frame(WaterTemp_Mean_pred = numeric())

# predict water temp at the random segment using the parameters from the model for that segment
for (i in 1:nrow(NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID,])) {
  airtemp_obs <- NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID, "AirTemp_c_MEAN"]
  watertemp_pred_weekly[i,1] <- Weekly_Logistic_epsilons[Weekly_Logistic_epsilons$COMID == sample_COMID, "mean"] + ((Weekly_Logistic_zetas[Weekly_Logistic_zetas$COMID == sample_COMID, "mean"] - Weekly_Logistic_epsilons[Weekly_Logistic_epsilons$COMID == sample_COMID, "mean"])/(1 + exp(Weekly_Logistic_phis[Weekly_Logistic_phis$COMID == sample_COMID, "mean"] * (Weekly_Logistic_kappas[Weekly_Logistic_kappas$COMID == sample_COMID, "mean"] - airtemp_obs[i]))))
}

# plot
ggplot() +
  geom_point(aes(x = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$AirTemp_c_MEAN,
                 y = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$WaterTemp_c_MEAN)) +
  geom_line(aes(x = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$AirTemp_c_MEAN,
                y = watertemp_pred_weekly$WaterTemp_Mean_pred),
            color = "red") +
  labs(x = "Air Temp (c)",
       y = "Water Temp (c)") +
       #title = sample_COMID) +
  theme_classic()

##### 
# Null linear Model (no PCA)
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_LM_Null.jags")
cat("
model{

  ## Priors
  
  # tau.beta - precision for beta
  # sd parameter for tau.betas
  sd.beta ~ dunif(0, 10)
  s2.beta <- sd.beta^2
    
  # tau.beta - precision parameter for all betas
  tau.beta <- 1/s2.beta
  
  # mu.beta.i - mean parameter for dnorm
  mu.beta ~ dnorm(0, 0.001)

  for (i in 1:nCOMIDs){

    # alpha.i - intercept for likelihood
     alpha[i] ~ dnorm(0, 0.001)

    # Beta.i - coefficient for daily air temperature
    beta[i] ~ dnorm(mu.beta, tau.beta)
  }
  
  ## Observation
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Mean_Pred[n] <- alpha[SegmentNo[n]] + beta[SegmentNo[n]] * AirTemp_Mean_Obs[n]
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Posterior predictive checks
  for (n in 1:nObs.weekly){
    PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
  }
  
  # Calculate RMSE
  for (n in 1:nObs.weekly){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.weekly)
  
  # Calculate Bayesian R^2
  for (n in 1:nObs.weekly){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.weekly - 1)
  var.pred <- sum(x.pred)/(nObs.weekly - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly_LM,
                  AirTemp_Mean_Obs = NS204_temps_weekly_filtered_LM$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_weekly_filtered_LM$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_weekly_filtered_LM$SegmentNo)

# Set parameters to save
jags_params <- c("beta", "mu.beta", "sd.beta", "tau.beta", "s2.beta", "theta", "tau", "sd", "RMSE", "R2", "PPC")

# MCMC settings
ni <- 4000
nc <- 3
nb <- 750
nt <- 1

# Fit model
Weekly_Mean_StreamTemp_LM_Null <- jagsUI::jags(data = jags_data,
                                     parameters.to.save = jags_params,
                                     model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_LM_Null.jags",
                                     n.chains = nc,
                                     n.iter = ni,
                                     n.burnin = nb,
                                     n.thin = nt,
                                     parallel = T)

# Save model output
Weekly_Mean_StreamTemp_LM_Null_params <- MCMCsummary(Weekly_Mean_StreamTemp_LM_Null, HPD = T, excl = "PPC")

## Posterior predictive check
Weekly_Mean_LM_Null_PPC_data <- data.frame(SegmentNo = NS204_temps_weekly_filtered_LM$SegmentNo,
                                     PPC = MCMCsummary(Weekly_Mean_StreamTemp_LM_Null, HPD = T, params = "PPC"))

Weekly_Mean_LM_Null_PPCs <- Weekly_Mean_LM_Null_PPC_data %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarise(PPC_mean = mean(PPC.mean),
                   PPC_sd = sd(PPC.mean))

Weekly_Mean_LM_Null_PPC_summ.table <- as.data.frame(apply(Weekly_Mean_LM_PPCs[,2:3],2,summary))

##### Null logistic model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_Logistic_Null.jags")
cat("
model{

  ## Priors
  
  # tau.phi - precision for all phis
  sd.phi ~ dunif(0, 10)
  s2.phi <- sd.phi^2
  tau.phi <- 1/(s2.phi)

  for (i in 1:nCOMIDs){

    # epsilon_i - estimated minimum stream temp at segment i 
    epsilon[i] ~ dnorm(min_waterTemps[i], 1/0.001) # uses a normal distribution informed by the measured minimum temperature
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dnorm(max_waterTemps[i], 1/0.001)
    
    # mu.phi_i - mean parameter for phi_i
    mu.phi[i] ~ dnorm(0, 0.001)

    # phi_i - measure of the steepest slope of the function at segment i
    phi[i] ~ dnorm(mu.phi[i], tau.phi)
    
    # kappa_i - Air temp at the function's inflection point for segment i
    kappa[i] ~ dnorm(20, 0.01)
    
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Mean_Pred[n] <- epsilon[SegmentNo[n]] + ((zeta[SegmentNo[n]] - epsilon[SegmentNo[n]])/(1 + exp(phi[SegmentNo[n]] * (kappa[SegmentNo[n]] - AirTemp_Mean_Obs[n]))))
    WaterTemp_Mean_Obs[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
    
    # New data for PPCs
    WaterTemp_Mean_New[n] ~ dnorm(WaterTemp_Mean_Pred[n], tau)
  }
  
  ## Derived Quantities
  # Relation to calculate phi from beta (slope), zeta, and epsilon
  for (i in 1:nCOMIDs){
    beta[i] <- (phi[i]*(zeta[i] - epsilon[i]))/4
  }
  
  # Posterior predictive checks
  for (n in 1:nObs.weekly){
    PPC[n] <- step(WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
  }
  
  # Calculate RMSE
  for (n in 1:nObs.weekly){
    SE[n] <- ((WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])^2)
  }
  RMSE <- sqrt(sum(SE)/nObs.weekly)
  
  # Calculate Bayesian R^2
  for (n in 1:nObs.weekly){
    resid[n] <- (WaterTemp_Mean_Obs[n] - WaterTemp_Mean_New[n])
    x.resid[n] <- (resid[n] - mean(resid))^2
    x.pred[n] <- (WaterTemp_Mean_Pred[n] - mean(WaterTemp_Mean_Pred))^2
  }
  
  var.resid <- sum(x.resid)/(nObs.weekly - 1)
  var.pred <- sum(x.pred)/(nObs.weekly - 1)
  R2 <- var.pred/(var.pred + var.resid)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly,
                  NS204_PCA_scores = NS204_PCA_scores,
                  min_waterTemps = min_max_waterTemps$min_waterTemp,
                  max_waterTemps = min_max_waterTemps$max_waterTemp,
                  AirTemp_Mean_Obs = NS204_temps_weekly_filtered$AirTemp_c_MEAN,
                  WaterTemp_Mean_Obs = NS204_temps_weekly_filtered$WaterTemp_c_MEAN,
                  SegmentNo = NS204_temps_weekly_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("epsilon", "zeta", "mu.phi", "sd.phi", "s2.phi", "phi", "beta", "kappa", "sd", "PPC", "RMSE", "R2")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

set.seed(1234)

# Fit model
Weekly_Mean_StreamTemp_Logistic_Null <- jagsUI::jags(data = jags_data,
                                                    parameters.to.save = jags_params,
                                                    model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Mean_StreamTemp_Logistic_Null.jags",
                                                    n.chains = nc,
                                                    n.iter = ni,
                                                    n.burnin = nb,
                                                    n.thin = nt,
                                                    parallel = T)

# Save model output
Weekly_Mean_StreamTemp_Logistic_Null_Params <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Null, HPD = T, excl = "PPC")

## Posterior predictive check
Weekly_Mean_Logistic_Null_PPC_data <- data.frame(SegmentNo = NS204_temps_weekly_filtered$SegmentNo,
                                          PPC = MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Null, HPD = T, params = "PPC"))

Weekly_Mean_Logistic_Null_PPCs <- Weekly_Mean_Logistic_Null_PPC_data %>% 
  group_by(SegmentNo) %>% 
  dplyr::summarise(PPC_mean = mean(PPC.mean),
                   PPC_sd = sd(PPC.mean))

Weekly_Mean_Logistic_Null_PPC_summ.table <- as.data.frame(apply(Weekly_Mean_Logistic_Null_PPCs[,2:3],2,summary))

################################
# Map sites

# filter to the sites at which we did our analysis
NS204_Sites <- NS204_Sites %>% 
  filter(COMID %in% NS204_temps_daily_filtered$COMID) %>% 
  arrange(COMID)

US_states <- map_data("state")

NS204_Sites_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = NS204_Sites, 
             aes(x = Long, y = Lat),
             color = "black") +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat") +
  theme_classic()

# for presentations
ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
geom_point(data = NS204_Sites,
           aes(x = Long, y = Lat),
           color = "black") +
  coord_map("albers",
            parameters = c(29.5, 45.5),
            xlim = c(-83.75, -71),
            ylim = c(33, 41.7)) +
  theme_void()

################################
# Model fits
# RMSE
model_RMSEs.table <- data.frame(Model = c("Daily_LM", "Weekly_LM", "Daily_Logistic", "Weekly_Logistic")) %>%
  cbind(RMSE = rbind(MCMCsummary(Daily_Max_StreamTemp_LM_Full, HPD = T, params = "RMSE")[,"mean"],
              MCMCsummary(Weekly_Max_StreamTemp_LM_Full, HPD = T, params = "RMSE")[,"mean"],
              MCMCsummary(Daily_Max_StreamTemp_Logistic_Full, HPD = T, params = "RMSE")[,"mean"],
              MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, HPD = T, params = "RMSE")[,"mean"])) %>% 
  cbind(R2 = rbind(MCMCsummary(Daily_Max_StreamTemp_LM_Full, HPD = T, params = "R2")[,"mean"],
                     MCMCsummary(Weekly_Max_StreamTemp_LM_Full, HPD = T, params = "R2")[,"mean"],
                     MCMCsummary(Daily_Max_StreamTemp_Logistic_Full, HPD = T, params = "R2")[,"mean"],
                     MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, HPD = T, params = "R2")[,"mean"]))

################################
# is there spatial structure in thermal sensitivity?
library(geoR)

logistic_weekly_betas <- Weekly_Max_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column("param") %>% 
  filter(str_detect(param, "^beta")) %>% 
  cbind(NS204_Sites) # join in site data

logistic_weekly_betas.geo <- as.geodata(logistic_weekly_betas[,c("mean", "Lat", "Long")], coords.col = c(2,3))
logistic_weekly_betas_semivariogram <- geoR::variog(logistic_weekly_betas.geo)
logistic_weekly_betas_semivariogram.plot <- ggplot(data = as.data.frame(logistic_weekly_betas_semivariogram[1:2])) +
  geom_point(aes(x = u,
                 y = v)) +
  labs(x = "Distance Class (DD)",
       y = "Semivariance") +
  theme_classic()

################################
# Visualize thetas
# use 95% highest posterior density credible intervals
thetas.table <- Weekly_Mean_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "theta")) %>% 
  as.data.frame() %>% 
  .[2:6,] %>%  # Remove first theta, which corresponds to intercept
  mutate(PC = c("PC1", "PC2", "PC3", "PC4", "PC5"))

theta_samples.table <- data.frame(sample_val = rbind(as.matrix(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,2]),
                                                     as.matrix(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,3]),
                                                     as.matrix(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,4]),
                                                     as.matrix(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,5]),
                                                     as.matrix(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,6])),
                                  PC = rep(c("PC1", "PC2", "PC3", "PC4", "PC5"),
                                           each = length(Weekly_Mean_StreamTemp_Logistic_Full$sims.list$theta[,1])))


theta_posteriors.plot <- ggplot(theta_samples.table) +
  geom_violin(aes(x = PC,
                  y = sample_val),
              trim = F) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             alpha = 0.5) +
  labs(x = "",
       y = "") +
  theme_classic()

################################
# Derived quantities

# C.phi - variation in slopes explained by landscape characteristics
C.phi <- (1 - ((Weekly_Mean_StreamTemp_Logistic_Full_Params['s2.phi',1])/(Weekly_Mean_StreamTemp_Logistic_Full_Params['sd.phi',1])))

# C - variation in water temp explained by slopes and air temp
#C <- (1 - ((StreamTemp_PCA_LM_full_v4_params['tau',1])/(1/StreamTemp_PCA_LM_v4_null_params['tau',1])))

#################################
# Regression slopes from linear models
# used for comparison to findings from other studies
weekly_linear_slopes.table <- Weekly_Mean_StreamTemp_LM_Full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta"))  # separate out betas

daily_linear_slopes.table <- Daily_Mean_StreamTemp_LM_Full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta"))

# Maximum slopes from nonlinear models
weekly_nonlinear_slopes.table <- Weekly_Mean_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta"))

daily_nonlinear_slopes.table <- Daily_Mean_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^beta"))

# Are linear regression slopes correlated to nonlinear regression slopes?
# Daily
cor(daily_linear_slopes.table$mean, daily_nonlinear_slopes.table$mean, method = "pearson")
# Weekly
cor(weekly_linear_slopes.table$mean, weekly_nonlinear_slopes.table$mean, method = "pearson")

# Intercepts from weekly linear model
weekly_linear_intercepts.table <- Weekly_Mean_StreamTemp_LM_Full_params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "^alpha"))  

#################################
# visualize nonlinear betas (slopes) from weekly model on a map
weekly_nonlinear_slopes.table <- weekly_nonlinear_slopes.table %>% 
  cbind(COMID = NS204_PCA_scores[,1]) # bind in COMIDs

# Join in (uncentered) landscape covariates
weekly_nonlinear_slopes.table <- weekly_nonlinear_slopes.table %>% 
  left_join(BKT_Habitat_StreamSegment_covars2)
  
# what variables are slopes correlated with?
slope_corrs.table <- cor(x = weekly_nonlinear_slopes.table[,9:184],
    y = weekly_nonlinear_slopes.table[,2],
    method = "spearman", 
    use = "pairwise.complete.obs") %>% 
  as.data.frame() %>% 
  rownames_to_column("Variable") %>% 
  arrange(-abs(V1)) %>% 
  head(10)


betas_map.plot <- ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = weekly_nonlinear_slopes.table, 
             aes(x = Long, y = Lat, color = mean)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior Mean Air-Water Temperature Slopes \nat NS204 Sites",
       color = expression(beta)) +
  scale_color_viridis_c() +
  #scale_color_viridis_c(limits=c(min(slopes$mean),max(slopes$mean))) +
  theme_classic()

# for presentations
ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = weekly_nonlinear_slopes.table, 
             aes(x = Long, y = Lat, color = mean)) +
  coord_map("albers",
            parameters = c(29.5, 45.5),
            xlim = c(-83.75, -71),
            ylim = c(33, 41.7)) +
  scale_color_viridis_c() +
  labs(color = "Thermal \nStability") +
  theme_void()

##################
##################
# Make predictions of beta (air-water temp slope) at unsampled sites based on covars
# use all sites of BKT habitat as defined by EBTJV
save.image()
# Filter pca data to just the trout sites
Trout_site_PCA_scores <- pca_scores %>% 
  filter(COMID %in% BKT_Habitat_StreamSegment_covars2$COMID) %>% 
  arrange(COMID)

## Predictions with posterior samples of nonlinear model
# Save parameters from weekly nonlinear temperature model
nonlinear_temp_model_params <- MCMCpstr(Weekly_Max_StreamTemp_Logistic_Full,
         params = c("theta", "sd.phi", "zeta", "epsilon"),
         type = 'chains')

# save the number of samples to try. This is the number of iterations from the model
n.samp <- length(nonlinear_temp_model_params$sd.phi)

# Create an empty dataframe to store predictive samples
nonlinear_beta_pstrs <- as.data.frame(matrix(NA, nrow = length(Trout_site_PCA_scores$COMID), ncol = n.samp))

# add COMID to the beginning of it
nonlinear_beta_pstrs <- data.frame(COMID = Trout_site_PCA_scores$COMID) %>% 
  cbind(nonlinear_beta_pstrs)

# Calculate median values for epsilon and zeta
# Since they represent the min and max water temperatures and we don't have these for unsampled sites, we'll have to make do with estimates from where we do have measurements
# It's probably not too far off, since the measured temperatures come from representative watersheds within the same region
med_epsilon <- median(MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "epsilon")[,"mean"])
med_zeta <- median(MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "zeta")[,"mean"])

for (j in 1:n.samp) {
  for (i in 1:nrow(nonlinear_beta_pstrs)) {
    # calculate mu.phi at the current site
    mu.phi.i <- nonlinear_temp_model_params$theta[1,j] + nonlinear_temp_model_params$theta[2,j] * Trout_site_PCA_scores[i,2] + nonlinear_temp_model_params$theta[3,j] * Trout_site_PCA_scores[i,3] + nonlinear_temp_model_params$theta[4,j] * Trout_site_PCA_scores[i,4] + nonlinear_temp_model_params$theta[5,j] * Trout_site_PCA_scores[i,5] + nonlinear_temp_model_params$theta[6,j] * Trout_site_PCA_scores[i,6]
    
    # calculate phi
    phi.i <- rnorm(1, mean = as.numeric(mu.phi.i), sd = as.numeric(nonlinear_temp_model_params$sd.phi[j]))
    
    # and beta from phi. Save.
    nonlinear_beta_pstrs[i,j+1] <- (phi.i*(med_zeta - med_epsilon))/4
  }
}

# calculate posterior means and left join in site attributes
trout_site_nonlinear_betas <- nonlinear_beta_pstrs %>% 
  mutate(pstr_mean_nonlinear_beta = rowMeans(.[,2:12001])) %>% 
  dplyr::select(COMID,
                pstr_mean_nonlinear_beta) %>% 
  left_join(BKT_Habitat_StreamSegment_covars2)

trout_site_nonlinear_betas %>% 
  filter(pstr_mean_nonlinear_beta > 1) %>% 
  left_join(Trout_site_PCA_scores) %>% 
  left_join(StreamSegment_covars) %>% 
  view()

# what variables are betas correlated with?
nonlinear_beta_corrs <- data.frame(corr = cor(trout_site_nonlinear_betas[,c(2:178)], method = "spearman", use = "pairwise.complete.obs")[1,]) %>% 
  rownames_to_column("param") %>% 
  arrange(-abs(corr)) %>% 
  .[-1,]

# and make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_nonlinear_betas, 
             aes(x = Long, 
                 y = Lat, 
                 color = pstr_mean_nonlinear_beta)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Predicted Air-Water Temperature Slopes \nat SE Trout Sites",
       color = expression(beta)) +
  scale_color_viridis_c() +
  #scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
  theme_classic()

fwrite(trout_site_betas, "C:/Users/georgepv/OneDrive - Colostate/Lab PC Backup/Desktop/trout_site_betas.csv")

## Predictions with point estimates of nonlinear model

# Filter pca data to just the trout sites

# Synch_sites <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/Compiled Trout Data (All sources)/SE_Site_Final.csv")
# 
# Trout_site_PCA_scores <- pca_scores %>% 
#   filter(COMID %in% Synch_sites$COMID) %>% 
#   arrange(COMID)

# Save parameters from weekly nonlinear temperature model
# nonlinear_thetas <- Weekly_Max_StreamTemp_Logistic_Full_Params %>% 
#   .[1:6,1]

nonlinear_thetas <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "theta")[,"mean"]

#nonlinear_thetas <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, func = "median", params = "theta")[,"func"]

#nonlinear_sd.phi <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "sd.phi")[,"mean"]

#nonlinear_sd.phi <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full,  func = "median",params = "sd.phi")[,"func"]

# Calculate median values for epsilon and zeta
# Since they represent the min and max water temperatures and we don't have these for unsampled sites, we'll have to make do with estimates from where we do have measurements
# It's probably not too far off, since the measured temperatures come from representative watersheds within the same region
med_epsilon <- median(MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, func = "median", params = "epsilon")[,"func"])
med_zeta <- median(MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, func = "median", params = "zeta")[,"func"])

epsilon <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "epsilon")[,"mean"]
zeta <- MCMCsummary(Weekly_Mean_StreamTemp_Logistic_Full, params = "zeta")[,"mean"]

trout_site_nonlinear_point_betas <- data.frame(COMID = Trout_site_PCA_scores$COMID,
                                               Beta = NA)

for (i in 1:nrow(trout_site_nonlinear_point_betas)) {
  # calculate phi at the current site
  phi.i <- nonlinear_thetas[1] + nonlinear_thetas[2] * Trout_site_PCA_scores[i,2] + nonlinear_thetas[3] * Trout_site_PCA_scores[i,3] + nonlinear_thetas[4] * Trout_site_PCA_scores[i,4] + nonlinear_thetas[5] * Trout_site_PCA_scores[i,5] + nonlinear_thetas[6] * Trout_site_PCA_scores[i,6]
  
  # calculate phi
  #phi.i <- rnorm(1, mean = as.numeric(mu.phi.i), sd = nonlinear_sd.phi)
  
  # and beta from phi. Save.
  trout_site_nonlinear_point_betas[i,"Beta"] <- (phi.i*(zeta - epsilon))/4
}

# join segment data to slope estimates
trout_site_nonlinear_point_betas <- trout_site_nonlinear_point_betas %>% 
  left_join(BKT_Habitat_StreamSegment_covars2)

# what variables are betas correlated with?
nonlinear_point_beta_corrs <- data.frame(corr = cor(trout_site_nonlinear_point_betas[,c(2:178)], method = "spearman", use = "pairwise.complete.obs")[1,]) %>% 
  rownames_to_column("param") %>% 
  arrange(-abs(corr)) %>% 
  .[-1,]

# are predicted slopes at all correlated with observed slopes?
trout_site_nonlinear_point_betas %>% 
  filter(COMID %in% NS204_temps_weekly_filtered$COMID) %>% 
  left_join(weekly_nonlinear_slopes.table) %>% 
  .[,c("Beta", "mean")] %>% 
  #view()
  cor(., method = "spearman", use = "pairwise.complete.obs")
  
# and make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_nonlinear_point_betas, 
             aes(x = Long, 
                 y = Lat, 
                 color = Beta),
             alpha = 0.5) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       title = "Predicted Slopes using Posterior Mean Parameter Estimates",
       color = expression(beta)) +
  scale_color_viridis_c() +
  #scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
  theme_classic()

# for presentations
ggplot() +
  geom_polygon(data = US_states,
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_nonlinear_point_betas, 
             aes(x = Long, 
                 y = Lat, 
                 color = Beta),
             size = 0.5,
             alpha = 0.5) +
  coord_map("albers",
            parameters = c(29.5, 45.5),
            xlim = c(-83.75, -71),
            ylim = c(33, 41.7)) +
  scale_color_viridis_c() +
  labs(color = "Predicted \nThermal \nStability") +
  theme_void()

## Predictions with linear model

# Save parameters from linear model
linear_model_params <- MCMCpstr(Weekly_Max_StreamTemp_LM_Full,
                              params = c("theta", "sd.beta"),
                              type = 'chains')

# save the number of samples to try. This is the number of iterations from the model
n.samp <- length(linear_model_params$sd.beta)

# Create an empty dataframe to store predictive samples
linear_beta_pstrs <- as.data.frame(matrix(NA, nrow = length(Trout_site_PCA_scores$COMID), ncol = n.samp))

# add COMID to the beginning of it
linear_beta_pstrs <- data.frame(COMID = Trout_site_PCA_scores$COMID) %>% 
  cbind(linear_beta_pstrs)

for (j in 1:n.samp) {
  for (i in 1:nrow(linear_beta_pstrs)) {
    # calculate slope at the current site
    linear_beta_pstrs[i,j+1] <- linear_model_params$theta[1,j] + linear_model_params$theta[2,j] * Trout_site_PCA_scores[i,2] + linear_model_params$theta[3,j] * Trout_site_PCA_scores[i,3] + linear_model_params$theta[4,j] * Trout_site_PCA_scores[i,4] + linear_model_params$theta[5,j] * Trout_site_PCA_scores[i,5] + linear_model_params$theta[6,j] * Trout_site_PCA_scores[i,6]
    
    # save slope
    #linear_beta_pstrs[i,j+1] <- rnorm(1, mean = as.numeric(mu.beta.i), sd = as.numeric(linear_model_params$sd.beta[j]))
  }
}

# calculate posterior means and left join in site attributes
trout_site_linear_betas <- linear_beta_pstrs %>% 
  mutate(pstr_mean_linear_beta = rowMeans(.[,2:12001])) %>% 
  dplyr::select(COMID,
                pstr_mean_linear_beta) %>% 
  left_join(Trout_Sites)

trout_site_linear_betas %>% 
  filter(pstr_mean_linear_beta > 1) %>% 
  left_join(Trout_site_PCA_scores) %>% 
  left_join(StreamSegment_covars) %>% 
  view()

# and make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_linear_betas, 
             aes(x = Long, y = Lat, color = pstr_mean_linear_beta)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Predicted Air-Water Temperature Slopes \nat SE Trout Sites",
       color = expression(beta)) +
  scale_color_viridis_c() +
  #scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
  theme_classic()

## Predictions with point estimates of linear model

# Save parameters from weekly linear temperature model
linear_thetas <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, params = "theta")[,"mean"]

#linear_thetas <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, func = "median", params = "theta")[,"func"]

linear_sd.beta <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, params = "sd.beta")[,"mean"]

#linear_sd.beta <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, func = "median", params = "sd.beta")[,"func"]

trout_site_linear_point_betas <- data.frame(COMID = Trout_site_PCA_scores$COMID,
                                               Beta = NA)

for (i in 1:nrow(trout_site_linear_point_betas)) {
  # calculate slope at the current site
  trout_site_linear_point_betas[i,2] <- linear_thetas[1] + linear_thetas[2] * Trout_site_PCA_scores[i,2] + linear_thetas[3] * Trout_site_PCA_scores[i,3] + linear_thetas[4] * Trout_site_PCA_scores[i,4] + linear_thetas[5] * Trout_site_PCA_scores[i,5] + linear_thetas[6] * Trout_site_PCA_scores[i,6]
  
  # save slope
  #trout_site_linear_point_betas[i,2] <- rnorm(1, mean = as.numeric(mu.beta.i), sd = linear_sd.beta)
}

# join segment data to slope estimates
trout_site_linear_point_betas <- trout_site_linear_point_betas %>% 
  left_join(BKT_Habitat_StreamSegment_covars2)

# what variables are betas correlated with?
linear_point_beta_corrs <- data.frame(corr = cor(trout_site_linear_point_betas[,c(2:178)], method = "spearman", use = "pairwise.complete.obs")[1,]) %>% 
  rownames_to_column("param") %>% 
  arrange(-abs(corr)) %>% 
  .[-1,]

# and make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_linear_point_betas, 
             aes(x = Long, 
                 y = Lat, 
                 color = Beta)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       title = "Predicted Slopes (linear) using Posterior Mean Parameter Estimates",
       color = expression(beta)) +
  scale_color_viridis_c() +
  #scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
  theme_classic()
########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here::here("Results", "v1.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
values <- ls()[str_detect(ls(), ".val")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
save(file = file.path(run_dir, "values.RData"), list = values)
