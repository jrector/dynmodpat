RAND_SEED <- function() {
  94023
}

NBOOTSTRAP_REPLICATES <- function() {
  1000
}

FILE_SUFFIX <- function() {
  paste0("seed", RAND_SEED())
}

# Columns used in the various models. Can be used to filter out unused columns
# to reduce the size of fitted models that include the source data.
MODEL_COLS <- function() {
  c(
    "participant_id", "is_bad_trial", "block_condition", "block_identifier",
    "stim_condition", "task_name", "probe_transform"
  )
} 

WM_MODEL_COLS <- function() {
  c(
  "is_old", "respondedsame", # essentials for SDT
  "trial_sumsim", "trial_sumsim_z", "trial_sumsim_trltype", "trial_sumsim_trltype_z",
  "stream_type", "nstim_fac", # experiment conditions
  "probe_is_monotonic", "probe_is_mono_fac",
  "target_last_pos_fac", "target_last_pos", "recency_effect", "recency_fac",
  "target_first_pos", "target_first_pos_fac", "target_serial_pos",
  "visual_probe_transform", "unimodal_probe_transform"
  )
}

ACUITY_MODEL_COLS <- function() {
  c(
  "resp_1st_interval", "same_in_1st_pair", "avg_time_difference_z",
  "freq_shift_type_fac", "ndirection_changes_fac", "block_first_stim"
  )
}

ALT_SIGVAR_MODEL_COLS <- function() {
  c(
  "is_old_nchange_and_direction",
  "is_old_unattended_stream"
  )
}

# Associate predictors to their abbreviation and logical indicating whether signal interaction is relevant.
# nested_within_sigvar indicates whether the predictor should be included in formulas
# nested under the signal indicator
PREDICTOR_INFO <- function() {
  data <- data.frame(
  abbrev=c("B", "N", "S", "M", "R", "T", "P", "C", "L", "V", "U", "F", "A", "D", "0WmDv", "0WmSig", "0AcuDv", "0AcuSig", "0Part"),
  predictor=c("block_condition", "nstim_fac", "stream_type", "probe_is_mono_fac",
    "recency_fac", "avg_time_difference_z", "probe_transform", "ndirection_changes_fac",
    "target_last_pos_fac", "visual_probe_transform", "unimodal_probe_transform",
    "block_first_stim", "freq_shift_type_fac", 
    "trial_sumsim_trltype_z",
    "is_old", "respondedsame", "resp_1st_interval", "same_in_1st_pair", "participant_id"),
  friendly_name=c("Modality", "ListLength", "PresentationType", "ProbeType",
    "Recency", "TimeDelta", "ProbeTransform", "DirectionChangeCount",
    "TargetLastPos", "VisualProbeTranspose", "UnimodalProbeTranspose",
    "FirstStimModality", "FreqShiftType",
    "ProbeItemSumSim",
    "IsTarget", "SaidTarget", "Is1stPair", "Said1stPair", "ParticipantId"),
  interaction_relevant=c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, NA, NA, NA, NA, NA),
  nested_within_sigvar=c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, NA, NA, NA, NA, NA),
  stringsAsFactors=FALSE
  )
  rownames(data) <- data$abbrev
  data
}

LSM_AT_OPTIONS <- function() {
  lsm_at_opts <- list()
  lsm_at_opts$binary <- c(0, 1)
  lsm_at_opts$contrast_neg_first <- c(-0.5, 0.5)
  lsm_at_opts$contrast_pos_first <- c(0.5, -0.5)
  lsm_at_opts$mean <- c(0.5)
  lsm_at_opts
}

