# Simple Occupancy Model for FWCE 409
# Colton Padilla
# 4/24/2023

# Last Updated: 4/24/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Install packages
install.packages(c("unmarked", "ggplot2", "MuMIn"))

# Load the packages
library(unmarked)
library(ggplot2)
library(MuMIn)

# Load the data
data <- read_csv("./Data/Sample_Data.csv")

# View the head of the data
head(dat)

#-----------------------------------------------------#
# Step 1 - Set up the dataframe in unmarked format ####
#-----------------------------------------------------#

# Detection history
det <- data[, 2:5]

# Site covariates
  # These are covariates that are intrinsic to the site and don't change
  # between surveys. Examples would be:
  # Slope, Aspect, Landcover
site_covs <- data[, 6:8]

# Survey covariates
  # These could be the same as site covariates (i.e. visual obstruction from trees)
  # or they could vary by survey. Examples:
  # Date, Temperature, Rainfall
survey_covs <- list(surv = data[, 9:12])

# Now create the object in unmarked format
umf <- unmarkedFrameOccu(y = det, # Detection history
                         siteCovs = site_covs, # Site covariates
                         obsCovs = survey_covs) # Survey covariates

# View the object
summary(umf)

#--------------------------------------------#
# Step 2 - Fit some models to select from ####
#--------------------------------------------#

# Fit the null model
null <- occu(formula = ~ 1 # Detection is first
                       ~ 1, # Occupancy is second
             data = umf) # The data frame that we created

# Fit a site covariate
site <- occu(formula = ~ 1
                       ~ cov1,
             data = umf)

# Fit a survey covariate
surv <- occu(formula = ~ surv
                       ~ 1,
             data = umf)

# Fit a survey covariate and a site covariate
site_surv <- occu(formula = ~ surv
                            ~ cov2,
                  data = umf)

#---------------------------------------#
# Step 3 - Model Selection using AIC ####
#---------------------------------------#

# Use unmarked for model selection------#
# Create a model list
fit <- fitList("psi(.) p(.)" = null,
               "psi(cov1) p(.)" = site,
               "psi(.) p(surv)" = surv,
               "psi(cov2) p(surv)" = site_surv)

# Now do model selection
aic_unmarked <- modSel(fit)
aic_unmarked # View

# Use MuMIn -------------#
# Create a list of models
modlist <- list(null, site, surv, site_surv)

# Now create an AICtable
aic_MuMIn <- model.sel(modlist)
aic_MuMIn # View

# View the top model summary
summary(site)

#-------------------------------------------------------#
# Step 4 - Get Occupancy and Detection Probabilities ####
#-------------------------------------------------------#

# Get the estimate for normal occupancy
backTransform( # The backtransformation function
  linearComb(site, # The model
             coefficients = c(1,0), # First coefficient
             "state") # Occupancy
)

# Create the estimates and 95% intervals
occ <- data.frame(
  parameter = "Occupancy_1",
  estimate = 0.498, # Mean estimate
  lower = 0.498 - 1.96 * 0.075, # 95% lower
  upper = 0.498 + 1.96 * 0.075 # 95% upper
)

# Get estimate for occupancy at cov1 mean
backTransform(
  linearComb(site, 
             coefficients = c(1,1), # Second coefficient
             "state")
)

# Create the estimates and 95% intervals
occ1 <- data.frame(
  parameter = "Occupancy_cov1",
  estimate = 0.925, # Mean estimate
  lower = 0.925 - 1.96 * 0.0367, # 95% lower
  upper = 0.925 + 1.96 * 0.0367 # 95% upper
)

# Get detection estimate
  # We don't need the Linear combination function because we are not using a 
  # covariate for survey in this model
backTransform(site, type = "det")

# Create estimates and 95% intervals
det <- data.frame(
  parameter = "Detection_1",
  estimate = 0.502,
  lower = 0.502 - 1.96 * 0.024,
  upper = 0.502 + 1.96 * 0.024
)

# Combine all those estimates
param_est <- do.call(rbind, list(occ, occ1, det))
param_est # View

#---------------------------------------------#
# Step 5 - Plot your estimates graphically ####
#---------------------------------------------#

# Using ggplot this time
ggplot(data = param_est, # Data
       mapping = aes(x = parameter, # X-axis is parameter
                     y = estimate)) + # Y-axis is the mean estimate
  geom_point(size = 2) + # Points
  geom_errorbar( # Errorbars
    mapping = aes(ymin = lower, # 95% Lower
                  ymax = upper), # 95% Upper
    width = 0.2, # Width of horizontal bars
    linewidth = 1) + # Width of the lines themselves
  theme_bw() + # Black and white theme
  labs(x = "Parameter", # X-axis label
       y = "Probability", # Y-axis label
       title = "Probability of Detection and Occupancy") + # Title the plot
  theme(axis.title.y = element_text(margin = margin(r = 10)), # Space Y title
        axis.title.x = element_text(margin = margin(t = 10))) # Space X title


#--------------------------------------------------------#
# Step 5 - Congratulations you ran an occupancy model ####
#--------------------------------------------------------#

# These can be way way more complicated than this. This is just a simple model.
# Unmarked can do many different models including the N-mixture to estimate 
# abundance rather than occupancy.