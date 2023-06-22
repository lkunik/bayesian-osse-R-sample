# create sample domain lon/lat input file for Bayesian inversion OSSE

dest_dir <- "include/"

# Define the longitude and latitude values of the grid (center of cell)
lon <- round(seq(from = -112.295, to = -111.525, by = 0.01), 3)
lat <- round(seq(from = 40.395, to = 40.945, by = 0.01), 3)

# create an n x 2 array with all the lon/lat pairs
lonlat <- expand.grid(lon, lat)

# NOTE: It is important that all lat/lon pairs found in the lonlat_domain file are
# within the bounds of the prior, uncertainty (sigma) and truth NetCDF files

# NOTE: this sample grid is a rectangle covering the Salt Lake County area, including
# forested (non-urban) areas. For real inversion applications where you want to
# constrain only urban emissions, you may want to exclude all non-urban grid cells
# using a land cover product or personal knowledge of the domain. This is where
# you can define the specific grid cells you wish to use in the analysis. Because
# the lonlat_domain file is an n x 2 array, your domain area does not need to be
# a rectangular grid.

# ~~~~~~~~~~~~~~~~~~~~ Save lonlat_domain file ~~~~~~~~~~~~~~~~~~~~~~~#

print("Saving formatted domain grid cell coordinates to include/lonlat_domain.rds")

filepath <- paste0(dest_dir, "lonlat_domain.rds")
saveRDS(lonlat, filepath)
