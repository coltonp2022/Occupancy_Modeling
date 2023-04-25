# Simple Occupancy Model for FWCE 409
# Colton Padilla
# 4/24/2023

# Last Updated: 4/24/2023

#-----------------#
# Housekeeping ####
#-----------------#

# Install packages
install.packages(c("unmarked", "ggplot2", "MuMIn", "glmmTMB"))

# Load the packages
library(unmarked)
library(ggplot2)
library(MuMIn)

# Load the data
data <- read.csv("./Data/Sample_Data.csv")

# View the head of the data
head(data)

#-----------------------------------------------------#
# Step 1 - Set up the dataframe in unmarked format ####
#-----------------------------------------------------#

# Detection history
det <- data[, 2:5]

# Site covariates
  # These are covariates that are intrinsic to the site and don't change
  # between surveys. Examples would be:
  # Slope, Aspect, Elevation
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

# NOTE FOR GRADUATE STUDENTS:
# Unmarked models are supported by MuMIn. You can run the dredge on them to run
# all possible combinations of your covariates. Here is code to do it. 
# It takes forever if you have a lot of covariates.

# Create the global model
global <- occu(formula = ~ surv
                         ~ cov1 + cov2 + cov3,
               data = umf)

# Dredge
aic <- dredge(global)
aic # View the aic 

# View the dredge function and its arguments if needed
#?dredge

# CAUTION:
# Depending on your study design, the dredge may violate your design. If you are
# doing a priori model sets, this is a violation of that method. This is used 
# for exploratory analysis or for finding the absolute best model of everything.
# This method is highly criticized in the paper that is the gold standard for
# model selection. Here is the quote from it:

#"“Let the computer find out” is a poor strategy and usually reflects the fact
# that the researcher did not bother to think clearly about the problem of 
# interest and its scientific setting (Burnham and Anderson, 2002)."

# Only use the dredge if you have solid backing and reasoning to do it. 
# If you think about your question in detail and get to the heart of what you
# actually need/want to answer, you likely will not need this approach at all.

#-------------------------------------------------------#
# Step 4 - Get Occupancy and Detection Probabilities ####
#-------------------------------------------------------#

# Estimate occupancy using prediction
occ <- predict(site, type = "state")
occ$cov1 <- data$cov1 # Add the values for covariates
occ$site <- data$site # Add site variable
head(occ) # View

# Get detection estimate
backTransform(site, type = "det")

#-------------------------------------------------------#
# Step 5 - Plot your occupancy estimates graphically ####
#-------------------------------------------------------#

# Using ggplot to plot covariate effect
ggplot(data = occ, # Data
       mapping = aes(x = cov1, # X-axis is parameter
                     y = Predicted)) + # Y-axis is mean estimate
  geom_line(color = "darkgreen") + # Line
  geom_ribbon( # Shading
    mapping = aes(ymin = lower, # 95% Lower
                  ymax = upper), # 95% Upper
    alpha = 0.5, # Width of horizontal bars
    fill = "darkgreen", # Color the ribbon
    linewidth = 1) + # Width of the lines themselves
  theme_bw() + # Black and white theme
  labs(x = "Covariate 1 Value", # X-axis label
       y = "Probability of Occupancy") + # Y-axis title
  theme(axis.title.y = element_text(margin = margin(r = 10)), # Space Y title
        axis.title.x = element_text(margin = margin(t = 10))) # Space X title

# Plot by site
ggplot(data = occ, # Data
       mapping = aes(x = site, # X-axis is parameter
                     y = Predicted)) + # Y-axis is mean estimate
  geom_point(color = "darkgreen") + # Line
  geom_errorbar( # Shading
    mapping = aes(ymin = lower, # 95% Lower
                  ymax = upper), # 95% Upper
    alpha = 0.5, # Width of horizontal bars
    color = "darkgreen", # Color the ribbon
    linewidth = 1) + # Width of the lines themselves
  theme_bw() + # Black and white theme
  geom_hline(aes(yintercept = 1)) + # Put an intercept at 1
  labs(x = "Site Number", # X-axis label
       y = "Probability of Occupancy") + # Y-axis title
  theme(axis.title.y = element_text(margin = margin(r = 10)), # Space Y title
        axis.title.x = element_text(margin = margin(t = 10))) # Space X title


#--------------------------------------------------------#
# Step 5 - Congratulations you ran an occupancy model ####
#--------------------------------------------------------#

# These can be way way more complicated than this. This is just a simple model.
# Unmarked can do many different models including the N-mixture to estimate 
# abundance rather than occupancy. This package is great because it does 
# everything within R itself rather than another program. Another package is 
# RPresence but it uses a software called PRESENCE in the background similar to
# RMark and Mark. I would say unmarked and RPresence each have pros and cons. 
# Whichever you find to be more intuitive is probably the way to go. There was
# a graduate class through USGS and the Co-op for RPresence which was great. 