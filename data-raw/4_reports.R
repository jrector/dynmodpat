library(dynmodpat)
library(stringr)
library(dplyr)
library(lme4)
library(lsmeans)
library(ggplot2)

# Set the projectroot variable
verify_dirs_exist()
projectroot <- paste0(getwd(), "/")

# current_prefix <- "NoBadTrials"

the_seed <- RAND_SEED()
if (is.null(the_seed) | is.na(the_seed)) {
  stop("In order to reproduce the dissertation analyses, the random seed must be set.\n")
}
set.seed(the_seed)
# Not strictly necessary to set the seed here unless any computations below use
# random number generator.

generate_boot_report_from_spec()

set.seed(the_seed)
compute_effects_from_spec()

compute_recency_interactions <- function() {
  # This could be improved by reading from a specification CSV as with
  # other similar analyses.
  uni_stubs <- c("mlefit_uni_a1v1_upright_mDBNMRiNxR", "mlefit_uni_a1v1_rev1_mDBNMRiBxNiNxR",
                 "mlefit_uni_a1v1_iden_mDBNMRiBxNiNxR", "mlefit_uni_a2v2_mDBNMRiNxR",
                 "mlefit_uni_a1v1_combined_mDBNMPRiBxNiBxPiNxR")
  uni_results <- sapply(uni_stubs, function(stub) {
    smry <- readRDS(paste0(projectroot, "cache/bootout_NoBadTrials_", stub, ".rds"))
    x <- smry$boot_out$t
    est_vec <- (x[,"maineff::acc::__nstim_fac_recency_fac__4__Last pos"] - x[,"maineff::acc::__nstim_fac_recency_fac__4__Other pos"]) - (x[,"maineff::acc::__nstim_fac_recency_fac__2__Last pos"] - x[,"maineff::acc::__nstim_fac_recency_fac__2__Other pos"])
    est_smry <- boot_res_compute_summary_stats(est_vec)
  })

  multi_stubs <- c("mlefit_multi2_a1v1_mDBSMRiBxSiBxM", "mlefit_multi2_a2v2_mDSMR")
  multi_results <- sapply(multi_stubs, function(stub) {
    smry <- readRDS(paste0(projectroot, "cache/bootout_NoBadTrials_", stub, ".rds"))
    x <- smry$boot_out$t
    est_vec <- x[,"maineff::acc::__recency_fac__Last pos"] - x[,"maineff::acc::__recency_fac__Other pos"]
    est_smry <- boot_res_compute_summary_stats(est_vec)
  })

  # write to file
  sink(file=str_c(projectroot, "reports/recency_interactions.txt", sep=""))
  print(uni_results)
  print(multi_results)
  sink()
}

set.seed(the_seed)
compute_recency_interactions()

# Analyze the lure trials in the multimodal Chapter 3 experiments to see if
# unattended stimuli influence responses.
set.seed(the_seed)
unattended_target_analysis()

#### Plots

# Generate plots for SDT parameters
set.seed(the_seed)
generate_plots_from_spec()

# Generate additional plots such as the summed similarity by experiment plot
set.seed(the_seed)
generate_custom_plots(dat_by_stub = wm_trials)
