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

summary(pca)

# plot PCA
biplot(pca,
         ellipse = T) +
  theme_classic()

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

# Save loadings
fwrite(top.loadings, here("Results/v1.0","PCA_top5_loadings.csv"))

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
  mutate(SiteNo = match(COMID, unique(COMID)))  # Add a SiteNo column to the filtered data

NS204_temps_weekly_filtered <- NS204_temps_weekly %>% 
  filter(COMID %in% pca_scores$COMID, # Filter temp data to the COMIDs for which we have principle components
         !is.na(AirTemp_c_MAX)) %>% 
  arrange(COMID) %>% 
  mutate(SiteNo = match(COMID, unique(COMID)))  # Add a SiteNo column to the filtered data

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
    WaterTemp_Max_Pred[n] <- alpha[SiteNo[n]] + beta[SiteNo[n]] * AirTemp_Max_Obs[n]
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
                  SiteNo = NS204_temps_daily_filtered$SiteNo)

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
    WaterTemp_Max_Pred[n] <- alpha[SiteNo[n]] + beta[SiteNo[n]] * AirTemp_Max_Obs[n]
    WaterTemp_Max_Obs[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    
    WaterTemp_Max.new[n] ~ dnorm(WaterTemp_Max_Pred[n], tau)
    # res[n] <- WaterTemp_Max_Obs[n] - WaterTemp_Max_Pred[n]
    # res.new[n] <- WaterTemp_Max.new[n] - WaterTemp_Max_Pred[n]
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
                  SiteNo = NS204_temps_weekly_filtered$SiteNo)

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
                                             model.file = "Analysis/SE_Temp_Weekly_Max_StreamTemp_LM_Full.jags",
                                             n.chains = nc,
                                             n.iter = ni,
                                             n.burnin = nb,
                                             n.thin = nt,
                                             parallel = T)

# Save model output
Weekly_Max_StreamTemp_LM_Full_params <- MCMCsummary(Weekly_Max_StreamTemp_LM_Full, HPD = T)
# Get DIC
DIC(Weekly_Max_StreamTemp_LM_Full)

#########################
### Logistic Model ###
#########################

## following Mohseni et al. (1998), the following nonlinear logistic regression model is used:
# T(s) = epsilon + ((zeta - epsilon)/1 + e^beta(kappa-T(a)))

# Where: 
#   T(s) = estimated stream temp (C)
#   T(a) = measured air temp (C)
#   epsilon = estimated minimum stream temp (C)
#   zeta = estimated maximum stream temp (C)
#   beta = steepest slope of the function (C^-1)
#   kappa = Air temp (C) at the function's inflection point

## Weekly Max Temperatures ##
# Write model
sink("Analysis/JAGS_Files/SE_Temp_Weekly_Max_StreamTemp_Logistic_Full.jags")
cat("
model{

  ## Priors
  
  # theta.m - coefficients for PCA values
  for (m in 1:6){
   theta[m] ~ dnorm(0, 0.01)
  }
  
  # tau.beta - precision for all betas
  sd.beta ~ dunif(0, 10)
  tau.beta <- 1/sd.beta^2

  for (i in 1:nCOMIDs){

    # epsilon_i - estimated minimum stream temp at segment i
    epsilon[i] ~ dgamma((5^2/5^2), (5/5^2))
    
    # zeta_i - estimated maximum stream temp at segment i
    zeta[i] ~ dgamma((25^2/10^2), (25/10^2))
    
    # mu.beta_i - mean parameter for dnorm
    mu.beta[i] <- theta[1] + theta[2] * NS204_PCA_scores[i,2] + theta[3] * NS204_PCA_scores[i,3] + theta[4] * NS204_PCA_scores[i,4] + theta[5] * NS204_PCA_scores[i,5] + theta[6] * NS204_PCA_scores[i,6]

    # beta_i - steepest slope of the function at segment i
    beta[i] ~ dnorm(mu.beta[i], tau.beta)
    
    # kappa_i - Air temp at the function's inflection point for segment i
    kappa[i] ~ dnorm(20, 0.01)
    
  }
  
  # tau - prescision parameter for water temperature
  sd ~ dunif(0,10)
  tau <- 1/(sd^2)
  
  ## Process
  for (n in 1:nObs.weekly){
    WaterTemp_Max_Pred[n] <- epsilon[[SiteNo]] + ((zeta[SiteNo[n]] - epsilon[[SiteNo]])/1 + exp(beta[SiteNo[n]] * (kappa[SiteNo[n]] - AirTemp_Max_Obs[n])))
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
                  nObs.weekly = nObs.weekly,
                  NS204_PCA_scores = NS204_PCA_scores,
                  AirTemp_Max_Obs = NS204_temps_weekly_filtered$AirTemp_c_MAX,
                  WaterTemp_Max_Obs = NS204_temps_weekly_filtered$WaterTemp_c_MAX,
                  SiteNo = NS204_temps_weekly_filtered$SiteNo)

# Set parameters to save
jags_params <- c("theta", "epsilon", "zeta", "mu.beta", "sd.beta", "beta", "kappa", "sd", "pvalue.mean", "pvalue.sd", "WaterTemp_Max_New")

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

StreamTemp_PCA_LM_v4_null_params <- as.data.frame(StreamTemp_PCA_LM_null_v4$summary)
view(StreamTemp_PCA_LM_v4_null_params) 

################################
# Visualize thetas
# use 95% highest posterior density credible intervals
thetas <- StreamTemp_PCA_LM_full_v4_params %>% 
  rownames_to_column(., "param") %>% 
  .[508:512,2] %>% 
  as.data.frame() %>% 
  mutate(PC = c("PC1", "PC2", "PC3", "PC4", "PC5"),
         Pstr_0.025 = as.data.frame(MCMCpstr(StreamTemp_PCA_LM_full_v4, params = "theta", func = function(x)hdi(x, 0.95)))[2:6,1],
         Pstr_0.975 = as.data.frame(MCMCpstr(StreamTemp_PCA_LM_full_v4, params = "theta", func = function(x)hdi(x, 0.95)))[2:6,2]) %>% 
  rename(., Pstr_Mean = 1) %>% 
  relocate(PC, .before = everything())

ggplot(thetas) +
  geom_point(aes(x = PC,
                 y = Pstr_Mean)) +
  geom_linerange(aes(x = PC,
                     ymin = Pstr_0.025,
                     ymax = Pstr_0.975)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             alpha = 0.5) +
  labs(x = "",
       y = "") +
  theme_classic()

################################
# Derived quantities

# C.beta - variation in slopes explained by landscape characteristics
C.beta.mean <- (1 - ((1/StreamTemp_PCA_LM_full_v4_params['tau.beta',1])/(1/StreamTemp_PCA_LM_v4_null_params['tau.beta',1])))
C.beta.025 <- (1 - ((1/StreamTemp_PCA_LM_full_v4_params['tau.beta',3])/(1/StreamTemp_PCA_LM_v4_null_params['tau.beta',3])))
C.beta.975 <- (1 - ((1/StreamTemp_PCA_LM_full_v4_params['tau.beta',7])/(1/StreamTemp_PCA_LM_v4_null_params['tau.beta',7])))

# C - variation in water temp explained by slopes and air temp
C <- (1 - ((1/StreamTemp_PCA_LM_full_v4_params['tau',1])/(1/StreamTemp_PCA_LM_v4_null_params['tau',1])))

#################################
# visualize betas (slopes) from full model on a map
slopes <- StreamTemp_PCA_LM_full_v4_params %>% 
  rownames_to_column(., "param") %>% 
  .[169:336,2] %>% # separate out mu.betas
  cbind(NS204_PCA_scores[,1]) %>% # bind in COMIDs
  as.data.frame()

names(slopes) <- c("beta", "COMID")

slopes <- slopes %>% 
  left_join(TempModel_Data)
  
# what variables are slopes correlated with?
cor(slopes[,c(1,5:7,14:129)], method = "spearman", use = "pairwise.complete.obs")[1,]
  
US_states <- map_data("state")

ggplot() +
  geom_polygon(data = US_states, 
               aes(x = long, y = lat, group = group),
               color = "black", fill = NA) +
  geom_point(data = slopes, 
             aes(x = Long, y = Lat, color = beta)) +
  coord_map("bonne",
            lat0 = 40,
            xlim = c(-84.2, -75.5),
            ylim = c(34.5, 40.5)) +
  labs(x = "Long",
       y = "Lat",
       #title = "Posterior Mean Air-Water Temperature Slopes \nat NS204 Sites",
       color = expression(beta)) +
  #scale_color_viridis_c(limits=c(0,1)) +
  scale_color_viridis_c(limits=c(min(slopes$beta),max(slopes$beta))) +
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
  filter(COMID %in% Trout_sites$COMID) %>% 
  arrange(COMID)

# Save parameters from temperature model
temp_model_params <- MCMCpstr(StreamTemp_PCA_LM_full_v4,
         params = c("theta", "sd.beta"),
         type = 'chains')

# save the number of samples to try. This is the number of iterations from the model
n.samp <- length(temp_model_params$sd.beta)

# Create an empty dataframe to store predictive samples
beta_pstrs <- as.data.frame(matrix(NA, nrow = length(Trout_site_PCA_scores$COMID), ncol = n.samp))

# add COMID to the beginning of it
beta_pstrs <- data.frame(COMID = Trout_site_PCA_scores$COMID) %>% 
  cbind(beta_pstrs)

for (j in 1:n.samp) {
  for (i in 1:nrow(beta_pstrs)) {
    # calculate slope at the current site
    mu.beta.i <- temp_model_params$theta[1,j] + temp_model_params$theta[2,j] * Trout_site_PCA_scores[i,2] + temp_model_params$theta[3,j] * Trout_site_PCA_scores[i,3] + temp_model_params$theta[4,j] * Trout_site_PCA_scores[i,4] + temp_model_params$theta[5,j] * Trout_site_PCA_scores[i,5] + temp_model_params$theta[6,j] * Trout_site_PCA_scores[i,6]
    
    # save slope
    beta_pstrs[i,j+1] <- rnorm(1, mean = as.numeric(mu.beta.i), sd = as.numeric(temp_model_params$sd.beta[j]))
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
