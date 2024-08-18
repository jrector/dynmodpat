library(dynmodpat)
library(stringr)
library(dplyr)
library(lme4)
library(lsmeans)

# Set the projectroot variable
verify_dirs_exist()
projectroot <- paste0(getwd(), "/")

current_prefix <- "Acuity"

the_seed <- RAND_SEED()
if (is.null(the_seed) | is.na(the_seed)) {
  stop("In order to reproduce the dissertation analyses, the random seed must be set.\n")
}
set.seed(the_seed)

# fit the acuity data
fits <- process_experiments(experiment_stubs = names(acuity_trials),
                             data_by_stub = acuity_trials,
                             project_root = projectroot,
                             analysis_prefix = current_prefix,
                             fit_mle = TRUE,
                             fit_stan = FALSE,
                             mle_batch_f = fit_mle_acuity_sequence
                             )

# explicitly set the seed again to ensure bootstrap results are stable
set.seed(the_seed)
bootstrap_cis_from_spec(ana_prefixes = c(current_prefix))
