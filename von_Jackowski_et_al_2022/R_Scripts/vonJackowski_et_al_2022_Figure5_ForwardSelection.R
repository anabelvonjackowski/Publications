
############################################################################################
##  von Jackowski et al., 2022 EM 10.1111/1462-2920.16036 ## #Forward selection of explanatory variables
############################################################################################
# based on the tutorial from David Zeleny Lab http://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel

library(vegan); packageVersion("vegan")
library(beepr)

# model containing only species matrix and intercept
BAC.rda.0 <- rda(BAC.otu ~ 1, data = BAC.env, scale = TRUE) 

#stepwise selection
BAC.rda.sel.os <- ordistep(BAC.rda.0, scope = formula (BAC.rda.all), direction = 'both') 
beep(sound = 1, expr = NULL)

#variance partitioning 
# varpart(BAC.otu, ~ Temperature_C, ~ Salinity_PSU, ~ Arabinose, ~ Rhamnose, 
# data = BAC.env[,c("Temperature_C", "Salinity_PSU", "Arabinose","Rhamnose")])

# Df    AIC      F Pr(>F)   
# - Glycine        1 496.84 3.1572  0.035 * 
#   - Aspartic.Acid  1 497.02 3.3164  0.025 * 
#   - Glucose        1 497.22 3.4986  0.020 * 
#   - Glucosamine    1 498.42 4.6048  0.020 * 
#   - Fucose         1 497.50 3.7554  0.010 **
#   - Rhamnose       1 499.10 5.2462  0.005 **
#   - CSP.Area       1 499.68 5.7945  0.005 **
#   - TEP.Area       1 501.53 7.5781  0.005 **
# + Tyrosine           1 495.13 1.8209  0.065 .
# + Serine             1 495.24 1.7291  0.065 .
# + Glucuronic.Acid    1 495.35 1.6303  0.070 .
# + Arabinose          1 495.54 1.4717  0.090 .
# + Mannose_Xylose     1 495.55 1.4588  0.145  
# + Isoleucine         1 495.77 1.2701  0.175  
# + Threonine          1 495.99 1.0804  0.355  
# + Galactosamine      1 496.00 1.0724  0.355  
# + Phenylalanine      1 496.07 1.0130  0.355  
# + Arginine           1 496.00 1.0677  0.375  
# + Galacturonic.Acid  1 496.04 1.0367  0.375  
# + Valine             1 496.08 0.9997  0.430  
# + GABA               1 496.10 0.9853  0.465  
# + Leucine            1 496.34 0.7778  0.685  
# + Galactose          1 496.40 0.7332  0.805  
# + Glutamic.Acid      1 496.51 0.6336  0.910  
# + Alanine            1 496.59 0.5674  0.965    