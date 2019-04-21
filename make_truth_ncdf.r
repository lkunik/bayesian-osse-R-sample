# create sample true emissions input file for Bayesian inversion OSSE

library(ncdf4)
dest_dir <- "include/"

# Define the dimension variables

lon_truth <- round(seq(from = -112.295, to = -111.525, by = 0.01), 3)
lat_truth <- round(seq(from = 40.395, to = 40.945, by = 0.01), 3)
time_truth <- seq(from = ISOdatetime(2015, 1, 1, 7, 0, 0, tz = "UTC"),
                  to = ISOdatetime(2016, 1, 1, 6, 0, 0, tz = "UTC"), by = "hour")

nlon <- length(lon_truth)
nlat <- length(lat_truth)
ntime <- length(time_truth)

# set values to 1 umol m-2 s-1
true_emiss <- array(1, dim = c(nlon, nlat, ntime))


time_dim <- ncdim_def("time", "seconds_since_1970_01_01", as.numeric(time_truth),
    longname = "seconds since R epoch: 1970-01-01 00:00:00")
lat_dim <- ncdim_def("lat", "degrees_north", lat_truth, longname = "latitude (center of cell)")
lon_dim <- ncdim_def("lon", "degrees_east", lon_truth, longname = "longitude (center of cell)")
true_emiss_var <- ncvar_def("emiss", "umol m-2 s-1", list(lon_dim, lat_dim, time_dim),
    longname = "true emissions")

true_vars <- list(true_emiss_var)

nc_filename <- paste0(dest_dir, "true_emiss.nc")
nc_truth <- nc_create(nc_filename, true_vars)
ncvar_put(nc_truth, true_emiss_var, true_emiss)

nc_close(nc_truth)
