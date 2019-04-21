# create sample prior emissions input file for Bayesian inversion OSSE

library(ncdf4)
dest_dir <- "include/"

# Define the dimension variables

lon_prior <- round(seq(from = -112.295, to = -111.525, by = 0.01), 3)
lat_prior <- round(seq(from = 40.395, to = 40.945, by = 0.01), 3)
time_prior <- seq(from = ISOdatetime(2015, 1, 1, 7, 0, 0, tz = "UTC"),
                  to = ISOdatetime(2016, 1, 1, 6, 0, 0, tz = "UTC"), by = "hour")

nlon <- length(lon_prior)
nlat <- length(lat_prior)
ntime <- length(time_prior)

# set values to 0 umol m-2 s-1
prior_emiss <- array(0, dim = c(nlon, nlat, ntime))


time_dim <- ncdim_def("time", "seconds_since_1970_01_01", as.numeric(time_prior),
    longname = "seconds since R epoch: 1970-01-01 00:00:00")
lat_dim <- ncdim_def("lat", "degrees_north", lat_prior, longname = "latitude (center of cell)")
lon_dim <- ncdim_def("lon", "degrees_east", lon_prior, longname = "longitude (center of cell)")
prior_emiss_var <- ncvar_def("emiss", "umol m-2 s-1", list(lon_dim, lat_dim, time_dim),
    longname = "prior emissions")

prior_vars <- list(prior_emiss_var)

nc_filename <- paste0(dest_dir, "prior_emiss.nc")
nc_prior <- nc_create(nc_filename, prior_vars)
ncvar_put(nc_prior, prior_emiss_var, prior_emiss)

nc_close(nc_prior)
