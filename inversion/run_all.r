#!/usr/bin/env Rscript

# launch script to run all inversion parts in order
# author: Lewis Kunik

# start the clock to keep track of how long it takes for the script to run
ptm1_all <- proc.time()

# ~~~~~~~~~~~ Clear all existing data files before proceeding ~~~~~~~~~#

source("config.r")

# get all files in the outer directory
out_files <- paste0(out_path, list.files(out_path))

if(!dir.exists("H/"))
    dir.create("H/")

if(!dir.exists("HQ/"))
    dir.create("HQ/")

if(!dir.exists("out/"))
    dir.create("out/")


if (clear_H) {
    # H files
    if (length(list.files("H/")) > 0) {
        out_files <- c(out_files, paste0("H/", list.files("H/")))
    }
}

# HQ files
if (length(list.files("HQ/")) > 0) {
    out_files <- c(out_files, paste0("HQ/", list.files("HQ/")))
}

# determine which of these files exist (redundant but good)
iFiles <- which(file.exists(out_files))

# remove existing files to clean the directory
invisible(sapply(out_files[iFiles], FUN = function(x) system(paste("rm", x))))


# ~~~~~~~~~~~~~~~~~ run scripts in order ~~~~~~~~~~~~~~~~~#

# 1. make receptor list (list of all footprint files used for observations)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_receptors.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_receptors.r")

# 2. make sprior (prior emissions vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_sprior.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_sprior.r")

# 3. make struth (true emissions vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_struth.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_struth.r")

# 4. make sigma (prior uncertainty vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_sigma.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_sigma.r")


#clear_H=F
if (clear_H) {
    # 6. make H (footprint matrices)
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("running Hsplit.r")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    source("src/Hsplit.r")
}

#6. make modeled obs (Hsprior, Hstruth)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running convolve_emiss.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/convolve_emiss.r")


# 8. make spatial covariance matrix
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_sp_cov.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_sp_cov.r")

# 9. make temporal covariance matrix
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_tmp_cov.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_tmp_cov.r")

# 10. make HQ
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_HQ.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_HQ.r")

# 13. make R (model data mismatch)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_R.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_R.r")

# 12. make z (anthropogenic enhancement values)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_z.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_z.r")

# 14. make z - Hsp
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_zHsp.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_zHsp.r")

# 15. make s_hat (optimized emissions vector)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running inversion.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/inversion.r")

# 16. make_Vshat (posterior uncertainty - technically makes Vshat-bar, which is
# the grid-scale aggregated uncertainty)
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("running make_Vshat.r")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
source("src/make_Vshat.r")

if (compute_chi_sq) {
    # 18. calculate Chi-squared
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    print("running chi_sq.r")
    print("~~~~~~~~~~~~~~~~~~~~~~~~")
    source("src/chi_sq.r")
}


# ~~~~~~~~~~~~~~~~ get elapsed time data ~~~~~~~~~~~~~~~~#
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
print("~~~~~~~~~~~~~~~~~~~~~~~~")
ptm2_all <- proc.time()
elapsed_seconds <- as.numeric(ptm2_all["elapsed"] - ptm1_all["elapsed"])
e_mins <- round(elapsed_seconds/60)
e_secs <- round(elapsed_seconds%%60, digits = 1)
print(paste0("elapsed time: ", e_mins, " minutes, ", e_secs, " seconds"))
