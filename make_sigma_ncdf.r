# create sample gridded prior emissions uncertainty input file for Bayesian inversion OSSE

library(ncdf4)
dest_dir <- "include/"

# Define the dimension variables

lon_sigma <- round(seq(from = -112.295, to = -111.525, by = 0.01), 3)
lat_sigma <- round(seq(from = 40.395, to = 40.945, by = 0.01), 3)
time_sigma <- seq(from = ISOdatetime(2015, 1, 1, 7, 0, 0, tz = "UTC"),
                  to = ISOdatetime(2016, 1, 1, 6, 0, 0, tz = "UTC"), by = "hour")

nlon <- length(lon_sigma)
nlat <- length(lat_sigma)
ntime <- length(time_sigma)

#set values to the difference between truth and prior (in this case, 1 umol m-2 s-1)
prior_uncert <- array(1, dim = c(nlon, nlat, ntime))


time_dim <- ncdim_def("time", "seconds_since_1970_01_01", as.numeric(time_sigma),
    longname = "seconds since R epoch: 1970-01-01 00:00:00")
lat_dim <- ncdim_def("lat", "degrees_north", lat_sigma, longname = "latitude (center of cell)")
lon_dim <- ncdim_def("lon", "degrees_east", lon_sigma, longname = "longitude (center of cell)")
prior_uncert_var <- ncvar_def("uncertainty", "umol m-2 s-1", list(lon_dim, lat_dim, time_dim),
    longname = "prior emissions uncertainty")

prior_uncert_vars <- list(prior_uncert_var)

nc_filename <- paste0(dest_dir, "prior_uncert.nc")
nc_sigma <- nc_create(nc_filename, prior_uncert_vars)
ncvar_put(nc_sigma, prior_uncert_var, prior_uncert)

nc_close(nc_sigma)
