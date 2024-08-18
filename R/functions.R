filter_unmodeled_cols <- function(dat, keep_cols=MODEL_COLS()) {
  # filter out the columns not needed for any of the models specified in this file
  localdat <- select_(dat, .dots=unique(keep_cols)) %>% as.data.frame
  localdat
}

filter_keep_wm_task_cols <- purrr::partial(
  filter_unmodeled_cols,
  keep_cols = c(MODEL_COLS(), WM_MODEL_COLS(), ALT_SIGVAR_MODEL_COLS())
  )

filter_keep_acuity_task_cols <- purrr::partial(
  filter_unmodeled_cols,
  keep_cols = c(MODEL_COLS(), ACUITY_MODEL_COLS())
  )

filter_keep_task_cols <- purrr::partial(
  filter_unmodeled_cols,
  keep_cols = c(MODEL_COLS(), WM_MODEL_COLS(), ACUITY_MODEL_COLS(), ALT_SIGVAR_MODEL_COLS())
  )

# Optional filtering functions
filter_no_bad_trials <- function(dat) {
  out_dat <- dat %>% dplyr::filter(is_bad_trial==0) %>% as.data.frame()
  out_dat
}

filter_nonmono_only <- function(dat) {
  out_dat <- dat %>% dplyr::filter(probe_is_monotonic==0) %>% as.data.frame()
  out_dat
}

exclude_named_fits_from_list <- function(fit_list, suffixes_to_exclude = NULL) {
  # Make a mask with true values for fit names that end with any of the
  # suffixes provided in suffixes_to_exclude.
  msk <- sapply(names(fit_list), function(fitname, excl_suffixes) {
      any(sapply(suffixes_to_exclude, function(pat, fname) {
        str_detect(fname, paste(pat, "$", sep=""))
      }, fname=fitname))
    }, excl_suffixes=suffixes_to_exclude)
  # fit_list[! names(fit_list) %in% suffixes_to_exclude]
  fit_list[!msk]
}

verify_projectroot_exists <- function() {
  if (!exists("projectroot", inherits = FALSE)) {
    stop("Define the projectroot before running this file.")
  }
}

#' Takes a vector of directory names. Tests whether those directories
#' exist in the current working directory. All ust exist for this to pass.
verify_dirs_exist <- function(req_dirs = c("cache", "plots", "specs", "reports")) {
  all_dirs <- all(sapply(req_dirs, function(x) file.info(x)$isdir))
  if (is.na(all_dirs) || !all_dirs) {
    stop("Working directory must contain the following directories: ")
  }
}
