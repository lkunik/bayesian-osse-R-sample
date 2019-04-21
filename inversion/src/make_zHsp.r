# calculate observed minus prior enhancements
# author: Lewis Kunik

## prerequisite scripts:
##  make_sprior.r
##  Hsplit.r
##  make_z.r
##
## output files:
##  zHsp.rds - file containing a vector of values equal to (z - Hsp)

# run dependent scripts
source("config.r")

# load in z (obs - bg)
z_file <- paste0(out_path, "z.rds")
z <- readRDS(z_file)
nobs <- length(z)

Hsprior_file <- paste0(out_path, "Hsprior.rds")
Hsprior <- readRDS(Hsprior_file)

# subtract Hsprior from z
zHsp <- z - Hsprior

# save zHsp to file
print("saving z - Hsp to zHsp.rds")
filepath <- paste0(out_path, "zHsp.rds")
saveRDS(zHsp, filepath)
