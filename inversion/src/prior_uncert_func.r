# Qsum_time_func.r author: Lewis Kunik

# This script defines a function that can be used by other scripts to calculate
# time-aggregated grid-scale prior uncertainty, given start/end timesteps

## input variables:
##  tstart - beginning timestep over which to aggregate Q
##  tstop - ending timestep over which to aggregate Q
##  sp_cov - spatial covariance, matrix of dims (#cells x #cells)
##  tmp_cov - temporal covariance, matrix of dims (#times x #times)
##  sigma - prior uncertainty, vector of length (#times * #cells)
##
## output variables:
##  Qsum - grid-scale aggregated prior error covariance (aggregated over
##  timesteps specified by tstart & tstop), matrix of dims (#cells x #cells)


get_prior_uncert_time <- function(tstart, tstop, sp_cov, tmp_cov, sigma) {

    # get spatial covariance matrix E
    E <- sp_cov
    # get temporal covariance matrix D
    D <- tmp_cov

    # often we calculate Qsum over just 1 timestep, so ensure that the D matrix is
    # conformable for matrix multiplication, i.e. 1x1
    if (ntimes == 1)
        D <- as.matrix(D)

    # sigma must be in format ntime x ncells
    sigma_mat <- matrix(sigma, nrow = ntimes, byrow = T)

    time_span <- tstop - tstart + 1

    S <- sigma_mat[tstart:tstop, ]  #in CT-Lagrange, this sub-matrix of sigma is denoted as S

    # ensure S is matrix-conformable, ntime x ncells even if ntime = 1
    if (time_span == 1)
        S <- matrix(S, nrow = time_span, byrow = T)


    # ~~~~~~~~~~~~~~~~~~~~~~~ calculate Qsum ~~~~~~~~~~~~~~~~~~~~~~~#

    # make Qsum - formula shown at this link:
    # https://www.esrl.noaa.gov/gmd/ccgg/carbontracker-lagrange/doc/inversion.html#posterior-uncertainty
    # (equation 5 shows Qsum = (t(S) * D * S) * E, substituting x = t(S) * D.
    x <- t(S) %*% D[tstart:tstop, tstart:tstop]
    Qsum <- (x %*% S) * E  #Multiplication with E is an element-by-element multiplication, not matrix multiplication.

    #return var
    Qsum

}
