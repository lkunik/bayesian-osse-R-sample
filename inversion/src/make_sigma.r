# make sigma (gridded prior uncertainty) vector for domain
# author: Lewis Kunik

## prerequisite scripts:
##  none (but prior_uncert.nc must be present in the include/ directory)
##
## output files:
##  sigma.rds - file containing a vector of length (#cells * #times)
##  with values of prior uncertainty for every timestep

# load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")



sprior_file <- paste(out_path, "sprior.rds", sep = "")
sprior_vec <- readRDS(sprior_file)


struth_file <- paste(out_path, "struth.rds", sep = "")
struth_vec <- readRDS(struth_file)

sigma <- abs(sprior_vec - struth_vec)


#assert a minimum uncertainty of 1 umol/(m^2/s)
ilt1 <- which(sigma < 1)
sigma[ilt1] = 1


# ~~~~~~~~~~~~~~~~~~~~ Save prior uncertainty file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Saving formatted prior emission uncertainty to sigma file")

filepath <- paste0(out_path, "sigma.rds")
saveRDS(sigma, filepath)
