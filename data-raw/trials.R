library(dplyr)

# Proprocessing script takes the full set of trials and does some pruning.
# Cut out the practice trials, and remove columns that are not used in any
# of the analysis. Save the results as a data set available tp package users.

# "trials" data frame is in a file in extdata
load(system.file("extdata", "alltrials.rda", package = "dynmodpat"))

# isolate the experiment trials
exp_trials <- trials[trials$is_experiment_block==1,]

#####
# prepare the sternberg task data
#####

# Define the mapping between the experiment name used in the original
# data and the abbreviation used in the analysis scripts.
wm_name_map <- list(
  ModFuncWM_a1v1="uni_a1v1_upright",
  ModFuncWM_a1v1_rev1_RotatedProbe="uni_a1v1_rev1",
  ModFuncWM_a1v1_iden="uni_a1v1_iden",
  ModFuncWM_a2v2="uni_a2v2",
  ModFuncWM_multi_a1_v1_short1="multi1_a1v1",
  ModFuncWM_multi_a1_v1_short2="multi2_a1v1",
  ModFuncWM_multi_a2_v2_short2="multi2_a2v2"
)

get_trials_by_experiment_factory <- function(filter_f) {
  function(experiment_name, data) {
    # Apply the filter function to the rows that are in the specified
    # experiment.
    filter_f(
      data[data$expName == experiment_name, ]
      )
  }
}

get_trials_for_wm_exp <- get_trials_by_experiment_factory(filter_keep_wm_task_cols)
get_trials_for_acuity_exp <- get_trials_by_experiment_factory(filter_keep_acuity_task_cols)

wm_trials <- lapply(names(wm_name_map), get_trials_for_wm_exp, exp_trials)
# assign the shorter experiment identifiers to each list item
names(wm_trials) <- wm_name_map

# Add a data set that combines the first 3 unimodal experiments, and another
# that includes the visual trials only.
wm_trials$uni_a1v1_combined <- dplyr::bind_rows(
  wm_trials$uni_a1v1_upright, wm_trials$uni_a1v1_rev1, wm_trials$uni_a1v1_iden
  )

wm_trials$a1v1_combo_v1only <- filter_keep_wm_task_cols(
  exp_trials %>%
    dplyr::filter(visual_probe_transform != "Irrelevant") %>%
    dplyr::filter(task_name %in% c("sternberg", "sternberg_wm")) %>%
    as.data.frame()
  )

#####
# prepare the acuity task data
#####
all_acuity_trials <- exp_trials[exp_trials$task_name=="2ifc_pair_acuity",]

acuity_name_map <- list(
  ModFuncWM_cross_a1_v1_2ifc="acuity_2ifcexp1",
  ModFuncWM_cross_a1_v1_2ifc2="acuity_2ifcexp2",
  ModFuncWM_cross_a1_v1_2ifc_iden="acuity_iden1",
  ModFuncWM_cross_a1_v1_2ifc3="acuity_acrossshape1"
)

acuity_trials <- lapply(
  names(acuity_name_map),
  get_trials_for_acuity_exp, all_acuity_trials
  )
# assign the shorter experiment identifiers to each list item
names(acuity_trials) <- acuity_name_map

# Add groupings of multiple experiments that will be analyzed together.
acuity_trials$acuity_transform1 <- dplyr::bind_rows(
  acuity_trials$acuity_2ifcexp1, acuity_trials$acuity_2ifcexp2
  )

acuity_trials$acuity_withinshape1 <- dplyr::bind_rows(
  acuity_trials$acuity_transform1, acuity_trials$acuity_iden1
  )

# Just the crossmodal trials
acuity_trials$acuity_crossmodal <- acuity_trials$acuity_withinshape1 %>%
  dplyr::filter(block_condition == "a1v1_crossmodal") %>% as.data.frame()

# Just auditory trials from first 3 acuity experiments
acuity_trials$acuity_auditory_transpose <- acuity_trials$acuity_withinshape1 %>%
  dplyr::filter(block_condition == "a1freq") %>% as.data.frame()

#####
# prepare the cross carrier direction task data
#####

directionmap_trials <- filter_keep_acuity_task_cols(
  exp_trials[exp_trials$task_name=="2ifc_pair",]
  )

# save for use in the package
devtools::use_data(wm_trials, acuity_trials, directionmap_trials, overwrite = TRUE)
