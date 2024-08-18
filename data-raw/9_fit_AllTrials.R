library(dynmodpat)
library(stringr)
library(dplyr)
library(lme4)
library(lsmeans)

# Set the projectroot variable
verify_dirs_exist()
projectroot <- paste0(getwd(), "/")

current_prefix <- NULL

the_seed <- RAND_SEED()
if (is.null(the_seed) | is.na(the_seed)) {
  stop("In order to reproduce the dissertation analyses, the random seed must be set.\n")
}
set.seed(the_seed)

# fit the models to the sternberg task data with all trial data
fits <- process_experiments(experiment_stubs = names(wm_trials),
                             data_by_stub = wm_trials,
                             project_root = projectroot,
                             fit_mle = TRUE,
                             fit_stan = FALSE
                             )

