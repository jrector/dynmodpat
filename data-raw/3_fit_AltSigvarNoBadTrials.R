library(dynmodpat)
library(stringr)
library(dplyr)
library(lme4)
library(lsmeans)

# Set the projectroot variable
verify_dirs_exist()
projectroot <- paste0(getwd(), "/")

current_prefix <- "AltSigvarNoBadTrials"

the_seed <- RAND_SEED()
if (is.null(the_seed) | is.na(the_seed)) {
  stop("In order to reproduce the dissertation analyses, the random seed must be set.\n")
}
set.seed(the_seed)

# fit the models to the sternberg task data with flagged bad trials excluded
fits <- process_experiments(experiment_stubs = names(wm_trials),
                             data_by_stub = wm_trials,
                             project_root = projectroot,
                             filter_f = filter_no_bad_trials,
                             analysis_prefix = current_prefix,
                             mle_batch_f = fit_mle_nchange_sequence,
                             fit_mle = TRUE,
                             fit_stan = FALSE
                             )

# explicitly set the seed again to ensure bootstrap results are stable
set.seed(the_seed)
bootstrap_cis_from_spec(ana_prefixes = c(current_prefix))
