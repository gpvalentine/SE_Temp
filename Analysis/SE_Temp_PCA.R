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

set.seed(1234)

# Load and format the data
# Load EPA Streamcat data for SE US
StreamCat_Covars <- fread(here("Data", "SE_Temp_StreamCat_Covars_Combined.csv"))

# Load national NHDplus data
NHDplus_data <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/GIS Data/NHDplus/NHDPlusV21_NationalData_Seamless_Geodatabase_Lower48_07/NHDPlusv2.1_National_FlowlineData.csv")

# Combine StreamCat and NHDplus data
# Left join national NHD data to SE StreamCat data so that everything is for the SE US
StreamSegment_covars <- StreamCat_Covars %>% 
  left_join(NHDplus_data)

# Load stream segments that are within BKT habitat patches
BKT_Habitat_Patch_StreamSegments <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/GIS/NHDplus_Flowlines_BKT_Patches.csv")

# Filter covars to just stream segments in BKT habitat patches
# We lose a few stream segments in BKT habitat patches because they aren't in the StreamCat dataset. Braided streams and errors in NHDplus that were filtered out of StreamCat
BKT_Habitat_StreamSegment_covars <- StreamSegment_covars %>% 
  filter(COMID %in% BKT_Habitat_Patch_StreamSegments$COMID)

# Remove columns that have all zeros or all the same value: they cannot be scaled by the PCA and won't be useful as predictors
BKT_Habitat_StreamSegment_covars2 <- BKT_Habitat_StreamSegment_covars %>% 
  na_if(-9998) %>%  # Replace -9998 (the nodata value) with NA
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

# Run the PCA
# pca <- prcomp(BKT_Habitat_StreamSegment_covars2[,-1], # leave out first column, which is COMID
#               center = T,
#               scale. = T)

# Bayesian pca can take NA values
pca <- pcaMethods::bpca(BKT_Habitat_StreamSegment_covars3[,-1],
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

rownames(top.loadings.table) <- c("R^2", "Variables")
  

# extract PCA coords for individuals
#pca_coords <- get_pca_ind(pca)$coord
pca_scores <- scores(pca)

# and subset to just the first five 5
#pca_coords <- as.data.frame(pca_coords[,1:5])
pca_scores <- as.data.frame(pca_scores[,1:5])

# Rename columns for easier interpretation
#names(pca_coords) <- c("PC1", "PC2", "PC3", "PC4", "PC5")
names(pca_scores) <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# add COMIDs and rearrange them to the first column
# pca_coords <- pca_coords %>% 
#   cbind(COMID = BKT_Habitat_StreamSegment_covars2$COMID) %>% 
#   dplyr::select(COMID, PC1, PC2, PC3, PC4, PC5)

pca_scores <- pca_scores %>% 
  cbind(COMID = BKT_Habitat_StreamSegment_covars2$COMID) %>% 
  dplyr::select(COMID, PC1, PC2, PC3, PC4, PC5)
# save PCA coords
#fwrite(pca_coords, "pca_coords_v4.csv")

####################################
# Max stream temperature model

# Load NS204 temperature data
NS204_temps_daily <- fread(here("Data", "NS204_temps_daily.csv"))
NS204_temps_weekly <- fread(here("Data", "NS204_temps_weekly.csv"))

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
         !is.na(AirTemp_c_MAX)) %>% 
  arrange(COMID) %>% 
  mutate(SegmentNo = match(COMID, unique(COMID)))  # Add a SegmentNo column to the filtered data

NS204_temps_weekly_filtered <- NS204_temps_weekly %>% 
  filter(COMID %in% pca_scores$COMID, # Filter temp data to the COMIDs for which we have principle components
         !is.na(AirTemp_c_MAX)) %>% 
  arrange(COMID) %>% 
  mutate(SegmentNo = match(COMID, unique(COMID)))  # Add a SegmentNo column to the filtered data

# Plot the data from a random segment
NS204_temps_daily_filtered %>% 
  filter(COMID %in% sample(unique(COMID), 1)) %>% 
  ggplot() +
  geom_point(aes(x = AirTemp_c_MAX,
                y = WaterTemp_c_MAX)) +
  theme_classic()

NS204_temps_weekly_filtered %>% 
  filter(COMID %in% sample(unique(COMID), 1)) %>% 
  ggplot() +
  geom_point(aes(x = AirTemp_c_MAX,
                 y = WaterTemp_c_MAX)) +
  theme_classic()

# for the linear model, we exclude observations where the max daily air temperature is below 0 degrees
NS204_temps_daily_filtered_LM <- NS204_temps_daily_filtered %>% 
  filter(AirTemp_c_MAX >= 0)

NS204_temps_weekly_filtered_LM <- NS204_temps_weekly_filtered %>% 
  filter(AirTemp_c_MAX >= 0)

# Filter pca data to just the NS204 sites
NS204_PCA_scores <- pca_scores %>% 
  filter(COMID %in% NS204_temps_daily_filtered$COMID) %>% 
  arrange(COMID)

nCOMIDs <- nrow(NS204_PCA_scores)
nObs.daily <- nrow(NS204_temps_daily_filtered)
nObs.weekly <- nrow(NS204_temps_weekly_filtered)
nObs.daily_LM <- nrow(NS204_temps_daily_filtered_LM)
nObs.weekly_LM <- nrow(NS204_temps_weekly_filtered_LM)

## Daily Max Temperatures ##
# Write model
sink("Analysis/JAGS_Files/SE_Temp_Daily_Max_StreamTemp_LM_Full.jags")
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
    WaterTemp_Max_Pred[n] <- alpha[SegmentNo[n]] + beta[SegmentNo[n]] * AirTemp_Max_Obs[n]
    WaterTemp_Max_Obs[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    
    WaterTemp_Max_New[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
  }
  
  ## Posterior Predictive Checks
  # Mean
  mean.WaterTemp_Max_Obs <- mean(WaterTemp_Max_Obs)
  mean.WaterTemp_Max_New <- mean(WaterTemp_Max_New)
  pvalue.mean <- step(mean.WaterTemp_Max_New - mean.WaterTemp_Max_Obs)
  
  # sd
  sd.WaterTemp_Max_Obs <- sd(WaterTemp_Max_Obs)
  sd.WaterTemp_Max_New <- sd(WaterTemp_Max_New)
  pvalue.sd <- step(sd.WaterTemp_Max_New - sd.WaterTemp_Max_Obs)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.daily = nObs.daily_LM,
                  NS204_PCA_scores = NS204_PCA_scores,
                  AirTemp_Max_Obs = NS204_temps_daily_filtered_LM$AirTemp_c_MAX,
                  WaterTemp_Max_Obs = NS204_temps_daily_filtered_LM$WaterTemp_c_MAX,
                  SegmentNo = NS204_temps_daily_filtered_LM$SegmentNo)

# Set parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "sd.beta", "tau.beta", "theta", "tau", "sd", "pvalue.mean", "pvalue.sd")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

# Fit model
Daily_Max_StreamTemp_LM_Full <- jagsUI::jags(data = jags_data,
                             parameters.to.save = jags_params,
                             model.file = "Analysis/JAGS_Files/SE_Temp_Daily_Max_StreamTemp_LM_Full.jags",
                             n.chains = nc,
                             n.iter = ni,
                             n.burnin = nb,
                             n.thin = nt,
                             parallel = T)

# Save model summary
Daily_Max_StreamTemp_LM_Full_params <- MCMCsummary(Daily_Max_StreamTemp_LM_Full, HPD = T)
# Get DIC
DIC(Daily_Max_StreamTemp_LM_Full)

# plot to visualize
plot_data <- Daily_Max_StreamTemp_LM_Full_params %>% 
  rownames_to_column("param") %>% 
  filter(str_detect(param, "WaterTemp_Max_New")) %>% 
  select(mean) %>% 
  cbind(AirTemp_c_MAX = NS204_temps_daily_filtered$AirTemp_c_MAX,
        WaterTemp_c_MAX = NS204_temps_daily_filtered$WaterTemp_c_MAX,
        COMID = NS204_temps_daily_filtered$COMID)

plot_data %>% 
  filter(COMID %in% sample(unique(COMID), 1)) %>% 
  ggplot() +
  geom_point(aes(x = AirTemp_c_MAX,
                 y = WaterTemp_c_MAX)) +
  geom_line(aes(x = AirTemp_c_MAX,
                y = mean),
            color = "red")

## Weekly Max Temperatures ##
# Write model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Max_StreamTemp_LM_Full.jags")
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
    WaterTemp_Max_Pred[n] <- alpha[SegmentNo[n]] + beta[SegmentNo[n]] * AirTemp_Max_Obs[n]
    WaterTemp_Max_Obs[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    
    WaterTemp_Max_New[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
  }
  
  ## Posterior Predictive Checks
  # Mean
  mean.WaterTemp_Max_Obs <- mean(WaterTemp_Max_Obs)
  mean.WaterTemp_Max_New <- mean(WaterTemp_Max_New)
  pvalue.mean <- step(mean.WaterTemp_Max_New - mean.WaterTemp_Max_Obs)

  # sd
  sd.WaterTemp_Max_Obs <- sd(WaterTemp_Max_Obs)
  sd.WaterTemp_Max_New <- sd(WaterTemp_Max_New)
  pvalue.sd <- step(sd.WaterTemp_Max_New - sd.WaterTemp_Max_Obs)
  
  # fit <- sum(res[])
  # fit.new <- sum(res.new[])
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly_LM,
                  NS204_PCA_scores = NS204_PCA_scores,
                  AirTemp_Max_Obs = NS204_temps_weekly_filtered_LM$AirTemp_c_MAX,
                  WaterTemp_Max_Obs = NS204_temps_weekly_filtered_LM$WaterTemp_c_MAX,
                  SegmentNo = NS204_temps_weekly_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("alpha", "beta", "mu.beta", "sd.beta", "tau.beta", "theta", "tau", "sd", "pvalue.mean", "pvalue.sd")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

# Fit model
Weekly_Max_StreamTemp_LM_Full <- jagsUI::jags(data = jags_data,
                                             parameters.to.save = jags_params,
                                             model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Max_StreamTemp_LM_Full.jags",
                                             n.chains = nc,
                                             n.iter = ni,
                                             n.burnin = nb,
                                             n.thin = nt,
                                             parallel = T)

# Save model output
Weekly_Max_StreamTemp_LM_Full_params <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, HPD = T)
# Get DIC
DIC(Weekly_Max_StreamTemp_LM_Full)

## Plot fit over data to visualize
# Extract intercepts and slopes
Weekly_LM_intercepts <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, params = "alpha") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)
  
# Extract slopes
Weekly_LM_slopes <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, params = "beta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

## fit
# select a random segment
sample_COMID <- NS204_PCA_scores %>% filter(COMID %in% sample(unique(COMID), 1)) %>% .[,"COMID"]

# Create an empty dataframe to store predicted watertemp data
watertemp_pred_weekly_LM <- data.frame(WaterTemp_MAX_pred = numeric())

# predict water temp at the random segment using the parameters from the model for that segment
for (i in 1:nrow(NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID,])) {
  airtemp_obs <- NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID, "AirTemp_c_MAX"]
  watertemp_pred_weekly_LM[i,1] <- Weekly_LM_intercepts[Weekly_LM_intercepts$COMID == sample_COMID, "mean"] + Weekly_LM_slopes[Weekly_LM_slopes$COMID == sample_COMID, "mean"] * airtemp_obs[i]
}

# plot
ggplot() +
  geom_point(aes(x = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$AirTemp_c_MAX,
                 y = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$WaterTemp_c_MAX)) +
  geom_line(aes(x = NS204_temps_weekly_filtered_LM[NS204_temps_weekly_filtered_LM$COMID == sample_COMID]$AirTemp_c_MAX,
                y = watertemp_pred_weekly_LM$WaterTemp_MAX_pred),
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
#   zeta = estimated maximum stream temp (C)
#   phi = measure of steepest slope of the function (C^-1)
#   beta = steepest slope of the function (C)
#   kappa = Air temp (C) at the function's inflection point

## Daily Max Temperatures ##
# Create a dataframe of the min and max water temps at each site
min_max_waterTemps <- NS204_temps_daily_filtered %>% 
  group_by(SegmentNo) %>% 
  summarize(min_waterTemp = min(WaterTemp_c_MAX, na.rm = T),
            max_waterTemp = max(WaterTemp_c_MAX, na.rm = T))

# Write model
sink("Analysis/JAGS_Files/SE_Temp_Daily_Max_StreamTemp_Logistic_Full.jags")
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

  for (i in 1:nCOMIDs){

    # epsilon_i - estimated minimum stream temp at segment i 
    epsilon[i] ~ dnorm(min_waterTemps[i], 1/0.001) # uses a normal distribution informed by the measured minimum temperature
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dnorm(max_waterTemps[i], 1/0.001)
    
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
    WaterTemp_Max_Pred[n] <- epsilon[SegmentNo[n]] + ((zeta[SegmentNo[n]] - epsilon[SegmentNo[n]])/(1 + exp(phi[SegmentNo[n]] * (kappa[SegmentNo[n]] - AirTemp_Max_Obs[n]))))
    WaterTemp_Max_Obs[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    
    # New data for PPCs
    WaterTemp_Max_New[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
  }
  
  # Relation to calculate phi from beta (slope), zeta, and epsilon
  for (i in 1:nCOMIDs){
    beta[i] <- (phi[i]*(zeta[i] - epsilon[i]))/4
  }
  
  ## Posterior Predictive Checks
  # Mean
  mean.WaterTemp_Max_Obs <- mean(WaterTemp_Max_Obs)
  mean.WaterTemp_Max_New <- mean(WaterTemp_Max_New)
  pvalue.mean <- step(mean.WaterTemp_Max_New - mean.WaterTemp_Max_Obs)
  
  # sd
  sd.WaterTemp_Max_Obs <- sd(WaterTemp_Max_Obs)
  sd.WaterTemp_Max_New <- sd(WaterTemp_Max_New)
  pvalue.sd <- step(sd.WaterTemp_Max_New - sd.WaterTemp_Max_Obs)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.daily = nObs.daily,
                  NS204_PCA_scores = NS204_PCA_scores,
                  min_waterTemps = min_max_waterTemps$min_waterTemp,
                  max_waterTemps = min_max_waterTemps$max_waterTemp,
                  AirTemp_Max_Obs = NS204_temps_daily_filtered$AirTemp_c_MAX,
                  WaterTemp_Max_Obs = NS204_temps_daily_filtered$WaterTemp_c_MAX,
                  SegmentNo = NS204_temps_daily_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("theta", "epsilon", "zeta", "mu.phi", "sd.phi", "s2.phi", "phi", "beta", "kappa", "sd", "pvalue.mean", "pvalue.sd")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

# Fit model
Daily_Max_StreamTemp_Logistic_Full <- jagsUI::jags(data = jags_data,
                                                    parameters.to.save = jags_params,
                                                    model.file = "Analysis/JAGS_Files/SE_Temp_Daily_Max_StreamTemp_Logistic_Full.jags",
                                                    n.chains = nc,
                                                    n.iter = ni,
                                                    n.burnin = nb,
                                                    n.thin = nt,
                                                    parallel = T)

# Save model output
Daily_Max_StreamTemp_Logistic_Full_Params <- MCMCsummary(Daily_Max_StreamTemp_Logistic_Full, HPD = T)
# Get DIC
DIC(Daily_Max_StreamTemp_Logistic_Full)

## Weekly Max Temperatures ##
# Create a dataframe of the min and max water temps at each site
min_max_waterTemps <- NS204_temps_weekly_filtered %>% 
  group_by(SegmentNo) %>% 
  summarize(min_waterTemp = min(WaterTemp_c_MAX, na.rm = T),
            max_waterTemp = max(WaterTemp_c_MAX, na.rm = T))

# Write model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Max_StreamTemp_Logistic_Full.jags")
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

  for (i in 1:nCOMIDs){

    # epsilon_i - estimated minimum stream temp at segment i 
    epsilon[i] ~ dnorm(min_waterTemps[i], 1/0.001) # uses a normal distribution informed by the measured minimum temperature
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dnorm(max_waterTemps[i], 1/0.001)
    
    # mu.phi_i - mean parameter for phi_i
    mu.phi[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # phi_i - measure of the steepest slope of the function at segment i
    phi[i] ~ dnorm(mu.phi[i], tau.phi)
    #phi[i] ~ dnorm(0, 0.01) # testing to see if phi would change if I removed the model on the mean
    
    # kappa_i - Air temp at the function's inflection point for segment i
    kappa[i] ~ dnorm(20, 0.01)
    
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Max_Pred[n] <- epsilon[SegmentNo[n]] + ((zeta[SegmentNo[n]] - epsilon[SegmentNo[n]])/(1 + exp(phi[SegmentNo[n]] * (kappa[SegmentNo[n]] - AirTemp_Max_Obs[n]))))
    WaterTemp_Max_Obs[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    
    # New data for PPCs
    WaterTemp_Max_New[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
  }
  
  # Relation to calculate phi from beta (slope), zeta, and epsilon
  for (i in 1:nCOMIDs){
    beta[i] <- (phi[i]*(zeta[i] - epsilon[i]))/4
  }
  
  ## Posterior Predictive Checks
  # Mean
  mean.WaterTemp_Max_Obs <- mean(WaterTemp_Max_Obs)
  mean.WaterTemp_Max_New <- mean(WaterTemp_Max_New)
  pvalue.mean <- step(mean.WaterTemp_Max_New - mean.WaterTemp_Max_Obs)
  
  # sd
  sd.WaterTemp_Max_Obs <- sd(WaterTemp_Max_Obs)
  sd.WaterTemp_Max_New <- sd(WaterTemp_Max_New)
  pvalue.sd <- step(sd.WaterTemp_Max_New - sd.WaterTemp_Max_Obs)
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs.weekly = nObs.weekly,
                  NS204_PCA_scores = NS204_PCA_scores,
                  min_waterTemps = min_max_waterTemps$min_waterTemp,
                  max_waterTemps = min_max_waterTemps$max_waterTemp,
                  AirTemp_Max_Obs = NS204_temps_weekly_filtered$AirTemp_c_MAX,
                  WaterTemp_Max_Obs = NS204_temps_weekly_filtered$WaterTemp_c_MAX,
                  SegmentNo = NS204_temps_weekly_filtered$SegmentNo)

# Set parameters to save
jags_params <- c("theta", "epsilon", "zeta", "mu.phi", "sd.phi", "s2.phi", "phi", "beta", "kappa", "sd", "pvalue.mean", "pvalue.sd")

# MCMC settings
ni <- 5000
nc <- 3
nb <- 1000
nt <- 1

# Fit model
Weekly_Max_StreamTemp_Logistic_Full <- jagsUI::jags(data = jags_data,
                                              parameters.to.save = jags_params,
                                              model.file = "Analysis/JAGS_Files/SE_Temp_Weekly_Max_StreamTemp_Logistic_Full.jags",
                                              n.chains = nc,
                                              n.iter = ni,
                                              n.burnin = nb,
                                              n.thin = nt,
                                              parallel = T)

# Save model output
Weekly_Max_StreamTemp_Logistic_Full_Params <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, HPD = T)
# Get DIC
DIC(Weekly_Max_StreamTemp_Logistic_Full)

MCMCtrace(Weekly_Max_StreamTemp_Logistic_Full, params = "zeta", pdf = F)

## Plot fit over data to visualize
# Extract parameters
Weekly_Logistic_epsilons <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "epsilon") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_zetas <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "zeta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_phis <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "phi") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_kappas <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "kappa") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

Weekly_Logistic_betas <- MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "beta") %>% 
  dplyr::select(mean) %>% 
  cbind(COMID = NS204_PCA_scores$COMID)

## fit
# select a random segment
sample_COMID <- NS204_PCA_scores %>% filter(COMID %in% sample(unique(COMID), 1)) %>% .[,"COMID"]

# Create an empty dataframe to store predicted watertemp data
watertemp_pred_weekly <- data.frame(WaterTemp_MAX_pred = numeric())

# predict water temp at the random segment using the parameters from the model for that segment
for (i in 1:nrow(NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID,])) {
  airtemp_obs <- NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID, "AirTemp_c_MAX"]
  watertemp_pred_weekly[i,1] <- Weekly_Logistic_epsilons[Weekly_Logistic_epsilons$COMID == sample_COMID, "mean"] + ((Weekly_Logistic_zetas[Weekly_Logistic_zetas$COMID == sample_COMID, "mean"] - Weekly_Logistic_epsilons[Weekly_Logistic_epsilons$COMID == sample_COMID, "mean"])/(1 + exp(Weekly_Logistic_phis[Weekly_Logistic_phis$COMID == sample_COMID, "mean"] * (Weekly_Logistic_kappas[Weekly_Logistic_kappas$COMID == sample_COMID, "mean"] - airtemp_obs[i]))))
}

# plot
ggplot() +
  geom_point(aes(x = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$AirTemp_c_MAX,
                 y = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$WaterTemp_c_MAX)) +
  geom_line(aes(x = NS204_temps_weekly_filtered[NS204_temps_weekly_filtered$COMID == sample_COMID]$AirTemp_c_MAX,
                y = watertemp_pred_weekly$WaterTemp_MAX_pred),
            color = "red") +
  labs(x = "Air Temp (c)",
       y = "Water Temp (c)",
       title = sample_COMID) +
  theme_classic()

##### 
# Null Model (no PCA)
sink("Analysis/NS204_Max_StreamTemp_PCA_LM_null_v4.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:5){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.beta - precision for beta
  # sd parameter for tau.betas
  sd.beta ~ dunif(0, 10)
    
  # tau.beta - precision parameter for all betas
  tau.beta <- 1/sd.beta^2
  
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
  for (n in 1:nObs){
    NS204_summarized_temps_filtered[n,3] ~ dnorm(WaterTemp_Max[n], tau)
    WaterTemp_Max[n] <- alpha[SiteNo[n]] + beta[SiteNo[n]] * NS204_summarized_temps_filtered[n,2]
  }
}
", fill = TRUE)
sink()

# Bundle data
jags_data <- list(nCOMIDs = nCOMIDs,
                  nObs = nObs,
                  NS204_summarized_temps_filtered = data.matrix(NS204_summarized_temps_filtered),
                  SiteNo = NS204_summarized_temps_filtered$SiteNo)

# Set parameters to save
jags_params <- c("beta", "mu.beta", "sd.beta", "tau.beta", "theta", "tau", "sd")

# MCMC settings
ni <- 4000
nc <- 3
nb <- 750
nt <- 1

# Fit model
StreamTemp_PCA_LM_null_v4 <- jagsUI::jags(data = jags_data,
                                     parameters.to.save = jags_params,
                                     model.file = "NS204_Max_StreamTemp_PCA_LM_null_v4.jags",
                                     n.chains = nc,
                                     n.iter = ni,
                                     n.burnin = nb,
                                     n.thin = nt,
                                     parallel = T)

StreamTemp_PCA_LM_v4_null_params <- MCMCsummary(StreamTemp_PCA_LM_null_v4, HPD = T)

################################
# Map sites
NS204_Sites <- fread("C:/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Temperature Modeling Paper/SE_Temp/Data/NS204_sites.csv")

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
thetas <- Weekly_Max_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "theta")) %>% 
  as.data.frame() %>% 
  .[2:6,] %>%  # Remove first theta, which corresponds to intercept
  mutate(PC = c("PC1", "PC2", "PC3", "PC4", "PC5"))


ggplot(thetas) +
  geom_point(aes(x = PC,
                 y = mean)) +
  geom_linerange(aes(x = PC,
                     ymin = `95%_HPDL`,
                     ymax = `95%_HPDU`)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             alpha = 0.5) +
  labs(x = "",
       y = "") +
  theme_classic()

################################
# Derived quantities

# C.phi - variation in slopes explained by landscape characteristics
C.phi <- (1 - ((1/Weekly_Max_StreamTemp_Logistic_Full_Params['tau.phi',1])/(1/Weekly_Max_StreamTemp_Logistic_Full_Params['tau.phi',1])))

# C - variation in water temp explained by slopes and air temp
C <- (1 - ((1/StreamTemp_PCA_LM_full_v4_params['tau',1])/(1/StreamTemp_PCA_LM_v4_null_params['tau',1])))

#################################
# visualize betas (slopes) from full model on a map
slopes <- Weekly_Max_StreamTemp_Logistic_Full_Params %>% 
  rownames_to_column(., "param") %>% 
  filter(str_detect(param, "beta")) %>% # separate out mu.betas
  cbind(NS204_Sites) %>% # bind in COMIDs
  select(mean, COMID, Lat, Long)

# Join in (uncentered) landscape covariates
slopes <- slopes %>% 
  left_join(BKT_Habitat_StreamSegment_covars2)
  
# what variables are slopes correlated with?
cor(slopes[,c(1,5:7,14:129)], method = "spearman", use = "pairwise.complete.obs")[1,]


ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = slopes, 
             aes(x = Long, y = Lat, color = mean)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior Mean Air-Water Temperature Slopes \nat NS204 Sites",
       color = expression(beta)) +
  #scale_color_viridis_c(limits=c(0,1)) +
  scale_color_viridis_c(limits=c(min(slopes$mean),max(slopes$mean))) +
  theme_classic()

###############################
# Evaluate model with cross validation

# Validation data - inverse of training data
NS204_summarized_temps_filtered_VALIDATION <- NS204_summarized_temps %>% 
  select(COMID, # filter for the columns of interest for this model
         AirTemp_c_MAX_Measured,
         WaterTemp_c_MAX_Measured) %>% 
  filter(COMID %in% pca.ind$COMID, # Filter temp data to the COMIDs for which we have principle components
         !is.na(AirTemp_c_MAX_Measured)) %>% 
  arrange(COMID) %>% 
  mutate(SiteNo = match(COMID, unique(COMID))) %>%  # Add a SiteNo column to the filtered data
  anti_join(NS204_summarized_temps_filtered)

nObs_VAL <- nrow(NS204_summarized_temps_filtered_VALIDATION)

## Data for prediction model
# air temps
AirTemp_Max_Measured <- NS204_summarized_temps_filtered_VALIDATION$AirTemp_c_MAX_Measured

# alpha - 
alpha <- StreamTemp_PCA_LM_full_v4_params %>% 
  rownames_to_column() %>% 
  .[1:168,2]
  

# mu.betas - posterior mean slopes for each site from the full model
beta <- StreamTemp_PCA_LM_full_v4_params %>% 
  rownames_to_column() %>% 
  .[169:336,2]

# sd - posterior mean tau (precision) from full model converted to sd
# tau <- StreamTemp_PCA_LM_full_v4_params %>% 
#   .['tau',2]
#sd <- 1/sqrt(tau)
sd <- StreamTemp_PCA_LM_full_v4_params %>% 
  .['sd',1]

# SiteNo
SiteNo <- NS204_summarized_temps_filtered_VALIDATION$SiteNo

# Empty dataframe of watertemp predictions
WaterTemp_Max <- list(rep_len(NA, length.out = nObs_VAL))
WaterTemp_Max_Predicted <- rep_len(NA, length.out = nObs_VAL)

for (n in 1:nObs_VAL){
  print(n)
  WaterTemp_Max[[n]] <- alpha[SiteNo[n]] + beta[SiteNo[n]] * AirTemp_Max_Measured[n] # intercept/alpha is set to 0 here
  WaterTemp_Max_Predicted[n] <- rnorm(n = 1, mean = WaterTemp_Max[[n]], sd = sd)
}

# Evaluate the accuracy of the model, and save the model RMSE
accuracy <- accuracy(WaterTemp_Max_Predicted, NS204_summarized_temps_filtered_VALIDATION$WaterTemp_c_MAX_Measured)
MAX_watertemp_RMSE <- accuracy[2]

##################
##################
# Make predictions of beta (air-water temp slope) at unsampled sites based on covars
# use SE sites from my trout project because I already have coordinates and COMIDs for them

# Load trout site data
Trout_Sites <- fread("/Users/georgepv/OneDrive - Colostate/SE Eco-Hydrology Project/Data/Trout/Trout Data Working/Compiled Trout Data (All sources)/SE_Site_Final.csv")

# Filter to just SE area (south of Mason-Dixon line)
Trout_Sites <- Trout_Sites %>% 
  filter(Lat <= 39.71666666666667) %>% 
  arrange(COMID)

# Filter pca data to just the trout sites
Trout_site_PCA_scores <- pca_scores %>% 
  filter(COMID %in% Trout_Sites$COMID) %>% 
  arrange(COMID)

# Save parameters from temperature model
temp_model_params <- MCMCpstr(Weekly_Max_StreamTemp_Logistic_Full,
         params = c("theta", "sd.phi", "zeta", "epsilon"),
         type = 'chains')

# save the number of samples to try. This is the number of iterations from the model
n.samp <- length(temp_model_params$sd.phi)

# Create an empty dataframe to store predictive samples
beta_pstrs <- as.data.frame(matrix(NA, nrow = length(Trout_site_PCA_scores$COMID), ncol = n.samp))

# add COMID to the beginning of it
beta_pstrs <- data.frame(COMID = Trout_site_PCA_scores$COMID) %>% 
  cbind(beta_pstrs)

# Calculate median values for epsilon and zeta
# Since they represent the min and max water temperatures and we don't have these for unsampled sites, we'll have to make do with estimates from where we do have measurements
# It's probably not too far off, since the measured temperatures come from representative watersheds within the same region
med_epsilon <- median(MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "epsilon")[,"mean"])
med_zeta <- median(MCMCsummary(Weekly_Max_StreamTemp_Logistic_Full, params = "zeta")[,"mean"])

for (j in 1:n.samp) {
  for (i in 1:nrow(beta_pstrs)) {
    # calculate mu.phi at the current site
    mu.phi.i <- temp_model_params$theta[1,j] + temp_model_params$theta[2,j] * Trout_site_PCA_scores[i,2] + temp_model_params$theta[3,j] * Trout_site_PCA_scores[i,3] + temp_model_params$theta[4,j] * Trout_site_PCA_scores[i,4] + temp_model_params$theta[5,j] * Trout_site_PCA_scores[i,5] + temp_model_params$theta[6,j] * Trout_site_PCA_scores[i,6]
    
    # calculate phi
    phi.i <- rnorm(1, mean = as.numeric(mu.phi.i), sd = as.numeric(temp_model_params$sd.phi[j]))
    
    # and beta from phi. Save.
    beta_pstrs[i,j+1] <- (phi.i*(med_zeta - med_epsilon))/4
  }
}

# calculate posterior means and left join in site attributes
trout_site_betas <- beta_pstrs %>% 
  mutate(c = rowMeans(.[,2:12001])) %>% 
  dplyr::select(COMID,
                pstr_mean_beta) %>% 
  left_join(Trout_Sites)

trout_site_betas %>% 
  filter(pstr_mean_beta > 1) %>% 
  left_join(Trout_site_PCA_scores) %>% 
  left_join(StreamSegment_covars) %>% 
  view()

# and make a map
ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = trout_site_betas, 
             aes(x = Long, y = Lat, color = pstr_mean_beta)) +
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

########################################################
# Export plots to the results folder

# Save the directory to which to save results files
run_dir <- here("results", "v1.0")

plots <- ls()[str_detect(ls(), ".plot")]
tables <- ls()[str_detect(ls(), ".table")]
save(file = file.path(run_dir, "plots.RData"), list = plots)
save(file = file.path(run_dir, "tables.RData"), list = tables)
