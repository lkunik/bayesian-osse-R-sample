# run a monte-carlo style set of ensembles of inversion runs to get overall
# expected value. 50 ensembles of 200 model iterations are run, totalling 10,000
# runs

# author: Lewis Kunik

## prerequisite scripts:
##  make_sprior.r
##  make_struth.r
##  make_sigma.r
##  convolve_emiss.r
##  make_sp_cov.r
##  make_tmp_cov.r
##  make_HQ.r
##  make_R.r
##
## output files:
##  mc_avg.rds - domain-averaged values for all iterations
##  mc_avg_subset.rds - time-subsetted (based on subset_hours_utc in config.r)
##                      domain-averaged values for all iterations
##  (mc_chi_sq.rds) - reduced chi squared values for each iteration
##                    (if run_chi_sq_TF = TRUE)

# load package dependencies
library(lubridate)

# run dependent scripts
source("config.r")

# start the clock to get a picture of how long this script runs
ptm1_all <- proc.time()

# ~~~~~~~~~~~~~~~~~~~~~~~ create time bins ~~~~~~~~~~~~~~~~~~~~~~~#

# list all desired obs times
flux_times <- seq(from = flux_start_POSIX + flux_t_res/2, to = flux_end_POSIX, by = flux_t_res)

times_hr <- hour(flux_times)  #this gives you the hour of timestamp at the CENTER of each timestep-bin
isub <- which(times_hr %in% subset_hours_utc)  #defines subsetted time indices (based on subset_hours_utc in config.r)

# load grid info
lonlat_domain <- readRDS(lonlat_domain_file)


# the difference between #ensembles and #iterations is arbitrary, except that a
# new random seed is set for each ensemble
n_ensembles <- 50
n_iter <- 200

# these will hold the total vals
shat_all_avgs <- array(NA, dim = c(n_ensembles, n_iter))
running_avgs <- array(NA, dim = c(n_ensembles, n_iter))

shat_sub_avgs <- array(NA, dim = c(n_ensembles, n_iter))
running_sub_avgs <- array(NA, dim = c(n_ensembles, n_iter))

chi_sq_arr <- array(NA, dim = c(n_ensembles, n_iter))

# if you desire to run an ensemble of chi-squared scripts, set this to TRUE but
# NOTE that you may wish to reduce the total number of runs because chi_sq.r is
# very slow for large inverse problems
run_chi_sq_TF <- F

# get a list of random seeds
random_seeds <- sample(1:10000, n_ensembles, replace = F)

# pre-define the aggregation operator so it can be used for each run
dlat <- 0.01  #0.01 is the spatial resolution of the state vector
dlon <- 0.01
finecellareas <- abs(2 * pi * (6371009^2) * (sin((lonlat_domain[, 2] - dlat/2) *
    pi/180) - sin((lonlat_domain[, 2] + dlat/2) * pi/180))/(360/dlon))
total_area <- sum(finecellareas)
W <- matrix(finecellareas/total_area, nrow = 1)  #W is the 'aggregation operator'

# loop through all model runs
for (iiii in 1:n_ensembles) {

    set.seed(random_seeds[iiii])

    for (jjjj in 1:n_iter) {
        print(paste("Ensemble", iiii, " Iteration", jjjj))
        # create new pseudo-obs
        source("src/make_z.r")
        source("src/make_zHsp.r")
        source("src/inversion.r")  #run the inversion

        if (run_chi_sq_TF) {
            source("src/chi_sq.r")
            chi_sq_arr[iiii, jjjj] <- Chi_sq_r
        }

        # aggregate posterior emissions on the time-step level.
        shat_mat <- matrix(s_hat, nrow = ntimes, byrow = T)
        shat_time_avgs <- apply(shat_mat, FUN = function(x) W %*% x, MARGIN = 1)

        # record the domain-averaged posterior value
        shat_all_avgs[iiii, jjjj] <- mean(shat_time_avgs)
        running_avgs[iiii, jjjj] <- mean(shat_all_avgs[iiii, 1:jjjj])  #this is the running mean of this ensemble

        # record the time-subsetted domain-averaged posterior value
        shat_sub_avgs[iiii, jjjj] <- mean(shat_time_avgs[isub])
        running_sub_avgs[iiii, jjjj] <- mean(shat_sub_avgs[iiii, 1:jjjj])  #this is the running mean of this ensemble

    }  #end n_iter

    # save the arrays of averages so you can analyze them later
    mc_avg_path <- paste(out_path, "mc_avg.rds", sep = "")
    saveRDS(shat_all_avgs, mc_avg_path)

    mc_avg_sub_path <- paste(out_path, "mc_avg_subset.rds", sep = "")
    saveRDS(shat_sub_avgs, mc_avg_sub_path)

    if (run_chi_sq_TF) {
        mc_chi_sq_path <- paste(out_path, "mc_chi_sq.rds", sep = "")
        saveRDS(chi_sq_arr, mc_chi_sq_path)
    }

}  #end n_ensembles


# print results
print("MONTE CARLO SIMULATION COMPLETE")
print(paste("Average posterior emissions =", round(mean(shat_all_avgs, na.rm = T), 5), flux_units))
print(paste(subset_hour_begin, "-", subset_hour_end, "UTC average posterior emissions =",
    round(mean(shat_sub_avgs, na.rm = T), 5), flux_units))


ptm2_all <- proc.time()
elapsed_seconds <- as.numeric(ptm2_all["elapsed"] - ptm1_all["elapsed"])
e_mins <- round(elapsed_seconds/60)
e_secs <- round(elapsed_seconds%%60, digits = 1)
print(paste("elapsed time: ", e_mins, " minutes, ", e_secs, " seconds", sep = ""))
