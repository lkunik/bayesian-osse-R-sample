# convolve footprints with prior and true emissions, giving
# modeled enhancements at each receptor
# Author: Lewis Kunik

## prerequisite scripts:
##  make_receptors.r
##  make_sprior.r
##  make_struth.r
##  Hsplit.r
##
## output files:
##  Hsprior.rds - modeled prior enhancements (in order of receptors)
##  Hstruth.rds - modeled prior enhancements (in order of receptors)

# load package dependencies
library(ncdf4)

# run dependent scripts
source("config.r")


# ~~~~~~~~~~~~~~~ define a function to read sparse H file ~~~~~~~~~~~~~~~~ #
read_sparse_h <- function(timestep, nobs, ncells) {

    Hi_tmp <- readRDS(paste("H/H", formatC(timestep, width = filename_width, flag = "0"), ".rds", sep = ""))
    # Populate the H-slice matrix (nobs x ncells) with zeros
    Hi <- array(0, dim = c(nobs, ncells))

    # checking if Hi is empty
    if (length(Hi_tmp[, 1]) > 0) {
        # for every value in H file, locate where that exists and put a value there
        iobs <- Hi_tmp[, 1]
        icell <- Hi_tmp[, 2]
        Hi[cbind(iobs, icell)] <- Hi_tmp[, 3]
    }

    # return H timestep
    return(Hi)
}


# ~~~~~~~~~~~~~~~ Load receptor and grid info ~~~~~~~~~~~~~~~~#

lonlat_domain <- readRDS(lonlat_domain_file)
ncells <- nrow(lonlat_domain)
nsites <- length(sites)

recep_file <- paste0(out_path, "receptors.rds")
if(aggregate_obs)
  recep_file <- paste0(out_path, "receptors_aggr.rds")

receps <- readRDS(recep_file)
nobs <- nrow(receps)

# ~~~~~~~~~~~~~~~ Load hourly emissions files ~~~~~~~~~~~~~~~#

sprior_file <- paste(out_path, "sprior.rds", sep = "")
sprior_vec <- readRDS(sprior_file)
sprior_mat <- matrix(sprior_vec, nrow = ntimes, byrow = T)

struth_file <- paste(out_path, "struth.rds", sep = "")
struth_vec <- readRDS(struth_file)
struth_mat <- matrix(struth_vec, nrow = ntimes, byrow = T)

# set up the vectors
Hsprior <- rep(0, nobs)
Hstruth <- rep(0, nobs)

# convolve each H file and add to running total of Hsbio
for (ii in 1:ntimes) {

    Hi <- read_sparse_h(ii, nobs, ncells)

    # obtain biospheric contributions to observations
    Hstruth <- Hstruth + Hi %*% struth_mat[ii, ]
    Hsprior <- Hsprior + Hi %*% sprior_mat[ii, ]

}  #end ntimes for-loop

# save the receptor file as 2-D array of receptor times and files
# (times are saved as seconds since 1970-01-01Z)

filepath <- paste(out_path, "Hstruth.rds", sep = "")
saveRDS(Hstruth, filepath)

filepath <- paste(out_path, "Hsprior.rds", sep = "")
saveRDS(Hsprior, filepath)
