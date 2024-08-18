#' Placeholder replaces a complex function in the original analysis
#' that provides optional model summary. Not necessary for basic analysis.
#' 
process_stan_fits <- function(fit, base_vars = c(), data = NULL) {
	NULL
}

#' Placeholder replaces a complex function in the original analysis
#' that provides optional model summary. Not necessary for basic analysis.
#' Current version only computes in-sample error rate.
#' 
process_mle_fit <- function(fit, base_vars = c(), signal_var = NULL, data = NULL) {
	out_info = list()  # container for info this function will return

  # if the signal var was not supplied as an argument, extract it from the fit
  if (is.null(signal_var)) {
    signal_var <- get_signal_var_from_fit(fit)
  }
  # Compute in-sample error rate
  if (!is.null(data)) {
    out_info$insample_error_rate <- in_sample_error_rate(data, fit)
  }
  out_info
}

process_mle_fits <- function(fit_list, base_vars = c(), data = NULL) {
  # Consumes a list of fitted models.
  # Produces a revised list that has updated information in the info list element.

  # The elements of fit_list are lists where the element name is the model name and
  # the list contains an element name "fit" with the fitted model.
  for (mdl_name in names(fit_list)) {
    fit <- fit_list[[mdl_name]]$fit # convenience
    fit_list[[mdl_name]]$info <- process_mle_fit(fit, base_vars, data=data)
  }
  fit_list
} 

compare_mle_fits_from_list <- function(fit_list, filepath, outfile_append = FALSE) {
  # Consumes a list of lists where each entry in the top-level list
  # has named element "fit". Performs a comparison of each pair of
  # the model fits. 
  # Produces a report at the location filepath.

  # send output to file
  sink(file=filepath, append=outfile_append)

  # first report top 5 models by fit indicies
  report_top_n_by_aic(fit_list, n = 5)
  report_top_n_by_bic(fit_list, n = 5)
  report_models_by_aic(fit_list)
  report_models_by_bic(fit_list)
  report_models_by_error_rate(fit_list)
  
  sink(file=NULL)

}

sort_models_by_aic <- function(fit_list) {
  # uses AICcmodavg to produce a table of models sorted by AICc.
  fits <- sapply(fit_list, function(x) {x$fit})
  fittable <- AICcmodavg::aictab(fits, modnames=names(fits))
  fittable
}

report_models_by_aic <- function(fit_list) {
  # uses AICcmodavg to produce a table of models sorted by AICc.
  # fits <- sapply(fit_list, function(x) {x$fit})
  # fittable <- AICcmodavg::aictab(fits, modnames=names(fits))
  fittable <- sort_models_by_aic(fit_list)

  cat("Models sorted by AICc", sep="")
  cat("\n---------------\n")
  print(fittable)
  cat("---------------\n\n")
}

report_models_by_bic <- function(fit_list) {
  # uses AICcmodavg to produce a table of models sorted by BIC
  fits <- sapply(fit_list, function(x) {x$fit})
  fittable <- AICcmodavg::bictab(fits, modnames=names(fits))

  cat("Models sorted by BIC", sep="")
  cat("\n---------------\n")
  print(fittable)
  cat("---------------\n\n")
}

report_models_by_error_rate <- function(fit_list) {
  # Print the in-sample error rate
  mnames <- names(fit_list)
  err_list <- sapply(fit_list, function(x) {x$info$insample_error_rate})
  # transform the error list to a vector, sorted in ascending order
  errs <- round(sort(unlist(err_list)), digits=4)
  cat("In-sample error rate for each model", sep="")
  cat("\n---------------\n")
  for (entry_name in names(errs)) {
    if (!is.null(errs[entry_name])) {
      cat(errs[entry_name], "  ", entry_name, "\n")
    }
  }
  cat("---------------\n\n")
}

report_top_n_by_fit_criterion <- function(fit_list, n = length(fit_list), criterion_name = NULL, criterion_f = NULL) {
  # Consumes a list of model fits
  # Applies the criterion function to each fit, sorts the results by fit metric,
  # then reports the top n model names.
  scores <- sapply(fit_list, criterion_f)
  ranks <- rank(scores)
  df <- data.frame(model_name=names(ranks), rank=ranks, score=scores)
  ordered_df <- df[order(ranks),]

  cat("Top ", n, " by ", criterion_name, sep="")
  cat("\n---------------\n")
  for (i in seq_len(nrow(ordered_df))[1:n]) {
    item <- ordered_df[i,]
    cat(item$rank, ") ", as.character(item$model_name), " [", item$score, "]", "\n", sep="")
  }
  cat("---------------\n\n")
}

report_top_n_by_aic <- purrr::partial(
  report_top_n_by_fit_criterion,
  criterion_name = "AIC",
  criterion_f = function(x) {AIC(logLik(x$fit))}
  )

report_top_n_by_bic <- purrr::partial(
  report_top_n_by_fit_criterion,
  criterion_name = "BIC",
  criterion_f = function(x) {BIC(logLik(x$fit))}
  )


fit_model_if_relevant <- function(
  model_abbrev,  # string, each letter refers to a predictor to include
  model_flags,  # vector of logicals
  base_prefix,  # prefix assigned to the base model
  base_fit,  # optional fitted base model
  existing_fits,
  force_refit,
  formula_str = NULL, # optional string with the full model formula
  summary_f  # function used to summarize model fit
  ) {
  # build name for this model fit
  fitname <- str_c(base_prefix, model_abbrev)
  tm_fmt <- "%H:%M:%S"  # format for timestamps

  # Look for existing fit
  found_existing_fit <- FALSE
  if (exists(fitname, where=existing_fits)) {
    if (!is.null(existing_fits[[fitname]])) {
      found_existing_fit <- TRUE
    }
  }

  fit_info = existing_fits

  # Fit model only if needed
  if (model_flags & (!found_existing_fit | force_refit)) {
    cat(fitname, " started ", format(Sys.time(), tm_fmt), "...", "\n")
    # If formula string was not provided, build it using the model abbreviation
    if (missing(formula_str) | is.null(formula_str)) {
      signal_var <- get_signal_var_from_fit(base_fit)
      # build formula string with the predictors to append to base model
      formula_suffix <- build_formula_suffix_from_abbrev(abbrev_str=model_abbrev, signal_var=signal_var)
      # build complete formula string
      formula_str <- append_formula_from_abbrev(existing_fit=base_fit, formula_suffix=formula_suffix)
    }
    fit <- update(base_fit, formula. <- as.formula(formula_str))
    fit_info[[fitname]]$fit <- fit
    fit_info[[fitname]]$info <- summary_f(fit)
  }

  fit_info
}

fit_model_sequence <- function(
  dat,  # data frame
  stub,  # identifying name for this experiment
  estimation_lib = "mle",  # specify what estimation routine to use
  base_model_f = mD_mle,  # function that fits an appropriate base model
  base_model_fallback_f = mD_mle_noInteraction,  # alternate base model if other fails
  base_model_id = "_mD",  # identifier for variables corresponding to the provided base model
  modelfit_f = process_mle_fit,  # function used to compute various fit summaries
  reqenvir = .GlobalEnv,   # environment in which to save variables with output
  force_refit = FALSE,  # boolean, if true then all models refit even if previous analysis exists
  fit_info = list()  # fits from previous runs
  ) {
  # fits a series of models relevant to experiments with stream_type
  # stub is a name that will be used to prefix all of the variables here
  # so they can be saved in the global scope

  # Before fitting each model, we check and see if a variable with the
  # associated name already exists in the environment. If so, we do not
  # refit the model, unless the "force_refit" paramater is true.
  var_prefix = "mlefit_"  # default assumes we are using MLE
  if (estimation_lib == "stan") {
    var_prefix = "mfit_"
  }

  pfix <- str_c(var_prefix, stub, base_model_id)
  tm_fmt <- "%H:%M:%S"  # format for timestamps

  # local function to determine whether a particular prefix should have the model fit
  needs_fit <- function(stub, flags, existing_fits, force=force_refit) {
    found_existing_fit <- FALSE
    if (exists(stub, where=existing_fits)) {
      if (!is.null(existing_fits[[stub]])) {
        found_existing_fit <- TRUE
      }
    }
    # returns logical
    flags & (!found_existing_fit | force_refit)
  }

  # Set some flags based on the data to model. The idea here is to check for variability
  # along a few experiment condition variables to determine whether they are relevant for
  # the given dataset. By default these are all false.
  include_block = FALSE
  include_nstim = FALSE
  include_stream = FALSE
  include_mono = FALSE
  include_recency = FALSE
  include_timedelta = FALSE
  include_probetransform = FALSE
  include_nchanges = FALSE
  include_lastpos = FALSE
  include_visproberotate = FALSE
  include_uniprobetransform = FALSE
  include_firststim = FALSE
  include_freqshift = FALSE

  incl_B = FALSE
  incl_N = FALSE
  incl_S = FALSE
  incl_M = FALSE
  incl_R = FALSE
  incl_T = FALSE
  incl_P = FALSE
  incl_C = FALSE
  incl_L = FALSE
  incl_V = FALSE
  incl_U = FALSE
  incl_F = FALSE
  incl_A = FALSE

  # Flag to control fits to models that only apply to acuity data
  is_2ifc = any(dat$task_name %in% c("2ifc_pair_acuity"))

  # Flag to identify data sets with only one item per list, as these require
  # special handling. Not a very clean solution, but this is an unusual case. 
  is_only_nstim1 <- all(as.numeric(as.character(unique(dat$nstim_fac))) == 1)

  # Set flag to true if the relevant indicator variable has any variability in these data.
  if (length(unique(dat$block_condition)) > 1) {
    include_block = TRUE
    incl_B = TRUE
  }
  if (length(unique(dat$nstim_fac)) > 1) {
    include_nstim = TRUE
    incl_N = TRUE
  }
  if (length(unique(dat$stream_type)) > 1) {
    include_stream = TRUE
    incl_S = TRUE
  }
  if (length(unique(dat$probe_is_monotonic)) > 1) {
    include_mono = TRUE
    incl_M = TRUE
  }
  if (length(unique(dat$recency_fac)) > 1) {
    include_recency = TRUE
    incl_R = TRUE
  }
  if (length(unique(dat$avg_time_difference_z)) > 1) {
    include_timedelta = TRUE
    incl_T = TRUE
  }
  if (length(unique(dat$probe_transform)) > 1) {
    include_probetransform = TRUE
    incl_P = TRUE
  }
  if (length(unique(dat$ndirection_changes_fac)) > 1) {
    include_nchanges = TRUE
    incl_C = TRUE
  }
  if (length(unique(dat$target_last_pos_fac)) > 1) {
    include_lastpos = TRUE
    incl_L = TRUE
  }
  # Don't include the "Irrelevant" level of the visual probe transform in the count
  # of unique levels. Also restrict this to data with one modality.
  if ((length(setdiff(dat$visual_probe_transform, "Irrelevant")) > 1) & (!incl_B)) {
    include_visproberotate = TRUE
    incl_V = TRUE
  }
  # The unimodal probe transform indicator will be collinear with block_condition
  # unless there are more than 2 unique levels available. 
  if (length(setdiff(dat$unimodal_probe_transform, "Irrelevant")) > 2) {
    include_uniprobetransform = TRUE
    incl_U = TRUE
  }
  if (length(unique(dat$block_first_stim)) > 1) {
    include_firststim = TRUE
    incl_F = TRUE
  }
  if (length(unique(dat$freq_shift_type_fac)) > 1) {
    include_freqshift = TRUE
    incl_A = TRUE
  }

  ####
  # NOTE 20170828: disabling the 2-level recency factor. Estimates are typically
  # the same as with the 3-level factor, but the 3-level factor has the benefit
  # of being conceptually more clear.
  ####
  include_lastpos = FALSE
  incl_L = FALSE

  #### 
  # NOTE 20170912: disabling the U predictor which appears to generate estimates
  # that are misleading.
  include_uniprobetransform = FALSE
  incl_U = FALSE

  # We want to use a different base model for experiments with only nstim==1. 
  # The model to use is provided as a fallback function.
  if (is_only_nstim1) {
    base_model_f <- base_model_fallback_f
  }

  # store the results in a list of lists. Each entry in the top-level list contains
  # a name, fitted model, and optional additional data such as lsmeans. We use the
  # fit_info passed as a paramter in case any existing analyses are already done.

  # First fit the base model which includes any required covariates
  fitname <- pfix
  # We use try-catch here to handle cases where the base model issues a warning. This
  # happens in some cases where an interaction may be dropped because it is collinear.
  # In that scenario we want to use the supplied fallback base model (which we assume
  # is valid).

  if (needs_fit(stub=fitname, flags=TRUE, existing_fits=fit_info)) {
    cat(fitname, " started ", format(Sys.time(), tm_fmt), "...", "\n")
    # fit_info[[fitname]]$fit <- base_model_f(dat)
    fit_info[[fitname]]$fit <- tryCatch({
      fit <- base_model_f(dat)
    }, warning = function(w) {
      # There was a warning, so use the fallback model
      fit <- base_model_fallback_f(dat)
    })
  }

  # Define variable for base model for convenience
  base <- fit_info[[fitname]]$fit
  # Extract the predictor variable names from the base model. We will use this
  # information later to skip over these base variables when generating block-level
  # means, for example.
  base_model_preds <- all.vars(delete.response(terms(formula(base))))
  # There is one non-factor predictor we may use. If it is relevant, add it to the
  # vector of base predictors to ignore when computing lsmeans.
  # NOTE: If we end up with more non-factor predictors, this should be made more flexible.
  if (incl_T) {
    base_model_preds <- c(base_model_preds, pred_T)
  }

  # Extract the signal indicator for later use
  signal_var <- get_signal_var_from_fit(base)
  # Note that the default signal variable is used in some strings below, and is replaced by the local
  # signal_var in build_model_formula_str(). Another option would be to define a function here that is
  # applied to any explicitly defined predictor string that swaps in the correct signal variable name.

  # Perform any summaries on the base model
  fit_info[[fitname]]$info <- modelfit_f(fit_info[[fitname]]$fit, base_vars = base_model_preds, data = dat)

  # define function that computes model summaries given fixed base predictors
  func_model_summary <- purrr::partial(modelfit_f, base_vars = base_model_preds, data = dat)

  # define function to fit model using common arguments
  func_fit_abbrev <- purrr::partial(
    fit_model_if_relevant,
    base_prefix=pfix, base_fit=base, existing_fits=fit_info, 
    force_refit=force_refit, summary_f=func_model_summary
    )

  # add block condition
  fit_info <- func_fit_abbrev(
    model_abbrev="B", model_flags=(incl_B)
    )

  #####################################################################
  ## Models with stream type
  #####################################################################

  # add stream
  fit_info <- func_fit_abbrev(
    model_abbrev="S", model_flags=(incl_S)
    )

  # if modality and stream type both matter, this will be a good model
  fit_info <- func_fit_abbrev(
    model_abbrev="BS", model_flags=(incl_B & incl_S)
    )

  # fit a model with interaction between block and stream
  fit_info <- func_fit_abbrev(
    model_abbrev="BSiBxS", model_flags=(incl_B & incl_S)
    )

  #####################################################################
  ## Models with list length
  #####################################################################
  fit_info <- func_fit_abbrev(
    model_abbrev="N", model_flags=(incl_N)
    )

  # if modality and list length type both matter, this will be a good model
  fit_info <- func_fit_abbrev(
    model_abbrev="BN", model_flags=(incl_B & incl_N)
    )

  # fit a model with interaction between block and stream
  fit_info <- func_fit_abbrev(
    model_abbrev="BNiBxN", model_flags=(incl_B & incl_N)
    )

  #####################################################################
  ## Models with monotonic probe predictor
  #####################################################################

  # First try a basic model with just the monotonic probe predictor
  fit_info <- func_fit_abbrev(
    model_abbrev="M", model_flags=(incl_M)
    )

  # If modality is not important but stream type is, a model with stream type
  # and monotonic probe predictors this will be a good model
  fit_info <- func_fit_abbrev(
    model_abbrev="SM", model_flags=(incl_S & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NM", model_flags=(incl_N & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BM", model_flags=(incl_B & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BSM", model_flags=(incl_B & incl_S & incl_M)
    )

  # And one with the interaction and monotonic probes. This is a lot of terms
  # and may not be supported by the sample size.
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMiBxS", model_flags=(incl_B & incl_S & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNM", model_flags=(incl_B & incl_N & incl_M)
    )

  # And one with the interaction and monotonic probes. This is a lot of terms
  # and may not be supported by the sample size.
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMiBxN", model_flags=(incl_B & incl_N & incl_M)
    )

  #####################################################################
  ## Models with recency effect predictor
  ## Uses a 3-level factor that separates trials in which it is irrelevant
  ## from the target trials of interest. Note that one level of this factor
  ## is nearly collinear with the signal predictor because there is no
  ## variability on lure trials
  #####################################################################

  # First try a basic model with just the target in last position predictor
  fit_info <- func_fit_abbrev(
    model_abbrev="R", model_flags=(incl_R)
    )

  #### NOTE: we are missing many simpler models that include recency predictor, but
  #### will add in a couple of important ones here. If a complete set is needed, 
  #### several other combinations of predictors need to be added including interactions
  #### with other predictors like modality or stream type.

  # block + stream + monotonic + recency
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMR", model_flags=(incl_B & incl_S & incl_M & incl_R)
    )

  # stream + monotonic + recency but no effect of modality
  fit_info <- func_fit_abbrev(
    model_abbrev="SMR", model_flags=(incl_S & incl_M & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMR", model_flags=(incl_N & incl_M & incl_R)
    )

  # With recency and main factors but excluding monotonic. These should not be good
  # models, but useful as a comparison against others with similar degrees of freedom.
  fit_info <- func_fit_abbrev(
    model_abbrev="BNR", model_flags=(incl_B & incl_N & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BSR", model_flags=(incl_B & incl_S & incl_R)
    )

  # model suitable for unimodal experiments including monotonic and recency predictors
  
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMR", model_flags=(incl_B & incl_N & incl_M & incl_R)
    )

  # check for separate effect of recency and differnet levels of nstim
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMRiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R)
    )

  #######################################
  #### start of target_last_pos_fac models
  ## An alternative to the 3-level recency_fac predictor is a 2-level factor that
  ## lumps together lure trials, nstim=1 trials, and target trials where the target
  ## is not in last position. We will likely remove these models...

  # First try a basic model with just the target in last position predictor
  fit_info <- func_fit_abbrev(
    model_abbrev="L", model_flags=(incl_L)
    )

  # block + stream + monotonic + recency
  fit_info <- func_fit_abbrev(
    model_abbrev="BSML", model_flags=(incl_B & incl_S & incl_M & incl_L)
    )

  # stream + monotonic + recency but no effect of modality
  fit_info <- func_fit_abbrev(
    model_abbrev="SML", model_flags=(incl_S & incl_M & incl_L)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NML", model_flags=(incl_N & incl_M & incl_L)
    )

  # With recency and main factors but excluding monotonic. These should not be good
  # models, but useful as a comparison against others with similar degrees of freedom.
  fit_info <- func_fit_abbrev(
    model_abbrev="BNL", model_flags=(incl_B & incl_N & incl_L)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BSL", model_flags=(incl_B & incl_S & incl_L)
    )

  # model suitable for unimodal experiments including monotonic and recency predictors
  
  fit_info <- func_fit_abbrev(
    model_abbrev="BNML", model_flags=(incl_B & incl_N & incl_M & incl_L)
    )

  #### end of target_last_pos_fac models
  #######################################

  # Check for interaction between block_condition and nstim when both recency and 
  # monotonic probe predictors are in the model (for comparison with DBNMR)
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMRiBxN", model_flags=(incl_B & incl_N & incl_M & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMRiBxNiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMLiBxN", model_flags=(incl_B & incl_N & incl_M & incl_L)
    )

  # Basic interaction of block and nstim without recency. If this is a better fit than model without the interaction
  # then we need to consider the version with the recency predictor. 
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMiBxN", model_flags=(incl_B & incl_N & incl_M)
    )

  # Stream and block interaction, no recency
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMiBxS", model_flags=(incl_B & incl_S & incl_M)
    )

  # Stream and block interaction, with recency
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMRiBxS", model_flags=(incl_B & incl_S & incl_M & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BSMLiBxS", model_flags=(incl_B & incl_S & incl_M & incl_L)
    )

  ############
  # Does the monotonic probe predictor have interactions with the other factors?
  ############
  # includes block : nstim and block : monotonic interactions
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMiBxNiBxM", model_flags=(incl_B & incl_N & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMiBxNiNxM", model_flags=(incl_B & incl_N & incl_M)
    )

  # includes block : stream and block : monotonic interactions
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMiBxSiBxM", model_flags=(incl_B & incl_S & incl_M)
    )

  # includes block: stream and stream : monotonic interactions
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMiBxSiSxM", model_flags=(incl_B & incl_S & incl_M)
    )

  # Check for monotonic probe interactions in simpler model without block predictor
  fit_info <- func_fit_abbrev(
    model_abbrev="NMiNxM", model_flags=(incl_N & incl_M)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="SMiSxM", model_flags=(incl_S & incl_M)
    )

  # extend BSMiBxSiBxM by adding recency predictor
  fit_info <- func_fit_abbrev(
    model_abbrev="BSMRiBxSiBxM", model_flags=(incl_B & incl_S & incl_M & incl_R)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BSMRiBxSiSxM", model_flags=(incl_B & incl_S & incl_M & incl_R)
    )
  #### end the incomplete recency predictor section


  ############
  # Models for unimodal data sets that contain both iden and transformed
  # probe trials. We assume the N predictor is relevant in all cases.
  ############

  fit_info <- func_fit_abbrev(
    model_abbrev="NMP", model_flags=(incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMPR", model_flags=(incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNP", model_flags=(incl_B & incl_N & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNPR", model_flags=(incl_B & incl_N & incl_P & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMP", model_flags=(incl_B & incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPR", model_flags=(incl_B & incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  # interactions of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPiBxP", model_flags=(incl_B & incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxP", model_flags=(incl_B & incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  # interactions of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPiBxN", model_flags=(incl_B & incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxN", model_flags=(incl_B & incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  # interactions of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPiNxP", model_flags=(incl_B & incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiNxP", model_flags=(incl_B & incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  # interactions of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPiBxNiBxP", model_flags=(incl_B & incl_N & incl_M & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxNiBxP", model_flags=(incl_B & incl_N & incl_M & incl_P & incl_R & !is_2ifc)
    )

  # recency predictor is often important
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxNiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxNiNxPiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R & incl_P & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMPRiBxNiBxPiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R & incl_P & !is_2ifc)
    )

  ############
  # Models for unimodal data sets with a mixture of identical and transformed
  # probes. These models also distinguish between the different square probe
  # rotations (unlike the P predictor). Note that there is a lot of duplicated
  # code in this section, taken from the section above. All of P models were
  # copied and modified to use the U predictor.
  # Interactions BxU are not useful. Both modalities share one level of U and differ on the others.
  ############
  ## start U section

  fit_info <- func_fit_abbrev(
    model_abbrev="NMU", model_flags=(incl_N & incl_M & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMUR", model_flags=(incl_N & incl_M & incl_U & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNU", model_flags=(incl_B & incl_N & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNUR", model_flags=(incl_B & incl_N & incl_U & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMU", model_flags=(incl_B & incl_N & incl_M & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMUR", model_flags=(incl_B & incl_N & incl_M & incl_U & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMUiBxN", model_flags=(incl_B & incl_N & incl_M & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMURiBxN", model_flags=(incl_B & incl_N & incl_M & incl_U & incl_R & !is_2ifc)
    )

  # interactions of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMUiNxU", model_flags=(incl_B & incl_N & incl_M & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMURiNxU", model_flags=(incl_B & incl_N & incl_M & incl_U & incl_R & !is_2ifc)
    )

  # recency predictor is often important
  fit_info <- func_fit_abbrev(
    model_abbrev="BNMURiBxNiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R & incl_U & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="BNMURiBxNiNxUiNxR", model_flags=(incl_B & incl_N & incl_M & incl_R & incl_U & !is_2ifc)
    )
  ## end U section

  ############
  # Models that examine the contribution of probe rotation in the v1size condition.
  # Since this predictor is only relevant at one level of block condition, these are suitable
  # on data that includes only visual trials. The B and P predictors will not be relevant.
  # visual_probe_transform
  ############

  fit_info <- func_fit_abbrev(
    model_abbrev="NV", model_flags=(incl_N & incl_V & !is_2ifc)
    )

  # unlikely a model without N predictor will be a good model, so this could be excluded
  fit_info <- func_fit_abbrev(
    model_abbrev="MV", model_flags=(incl_M & incl_V & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMV", model_flags=(incl_N & incl_M & incl_V & !is_2ifc)
    )

  # interactions of list length and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="NMViNxV", model_flags=(incl_N & incl_M & incl_V & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMViMxV", model_flags=(incl_N & incl_M & incl_V & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMViMxViNxV", model_flags=(incl_N & incl_M & incl_V & !is_2ifc)
    )

  # Add recency predictor to all relevant models
  fit_info <- func_fit_abbrev(
    model_abbrev="NVR", model_flags=(incl_N & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="MVR", model_flags=(incl_M & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMVR", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  # models with interactions
  fit_info <- func_fit_abbrev(
    model_abbrev="NMVRiNxV", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMVRiMxV", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMVRiMxViNxV", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMVRiNxR", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="NMVRiNxViNxR", model_flags=(incl_N & incl_M & incl_V & incl_R & !is_2ifc)
    )

  # Note that we are not fitting a couple of variations, but could add them if needed. Specifically, we are not
  # adding the NxR interaction to models with MxV interaction (the MxV interaction seems unlikely).

  #####################################################################
  ## Models for acuity data
  #####################################################################

  # Average time difference between breakpoints in the "different" pair.
  # Models with the T predictor only apply to 2IFC acuity data
  # Models may have either B or F but not both since they are collinear for
  # auditory and visual blocks.
  fit_info <- func_fit_abbrev(
    model_abbrev="T", model_flags=(is_2ifc & incl_T)
    )

  # Block and time delta
  fit_info <- func_fit_abbrev(
    model_abbrev="TB", model_flags=(is_2ifc & incl_B & incl_T)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TF", model_flags=(is_2ifc & incl_F & incl_T)
    )

  # Probe type (identical or transformed)
  fit_info <- func_fit_abbrev(
    model_abbrev="P", model_flags=(is_2ifc & incl_P)
    )

  # The time difference almost certainly matters, so we include it in all basic models
  # moving forward.

  # Probe type (identical or transformed)
  fit_info <- func_fit_abbrev(
    model_abbrev="TP", model_flags=(is_2ifc & incl_T & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TBP", model_flags=(is_2ifc & incl_T & incl_B & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFP", model_flags=(is_2ifc & incl_T & incl_F & incl_P)
    )

  # Influence of time could be different by modality
  fit_info <- func_fit_abbrev(
    model_abbrev="TBiBxT", model_flags=(is_2ifc & incl_T & incl_B)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFiFxT", model_flags=(is_2ifc & incl_T & incl_F)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TBPiBxT", model_flags=(is_2ifc & incl_T & incl_B & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFPiFxT", model_flags=(is_2ifc & incl_T & incl_F & incl_P)
    )

  # Examine interaction of modality and probe transformation
  fit_info <- func_fit_abbrev(
    model_abbrev="TBPiBxP", model_flags=(is_2ifc & incl_T & incl_B & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFPiFxP", model_flags=(is_2ifc & incl_T & incl_F & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TBPiBxTiBxP", model_flags=(is_2ifc & incl_T & incl_B & incl_P)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFPiFxTiFxP", model_flags=(is_2ifc & incl_T & incl_F & incl_P)
    )

  # Number of direction changes in the patterns. Does not make sense for experiments
  # where the items within pairs can have different number of direction changes...
  fit_info <- func_fit_abbrev(
    model_abbrev="TBC", model_flags=(is_2ifc & incl_T & incl_B & incl_C)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TFC", model_flags=(is_2ifc & incl_T & incl_F & incl_C)
    )

  # Models appropriate when only auditory acuity trials are being considered. 
  # Assessing role of probe transposition distance.

  # Only include the A predictor if there is no variation in the block identifiers
  # to avoid collinearity.
  incl_A <- incl_A & !incl_B & !incl_F
  fit_info <- func_fit_abbrev(
    model_abbrev="A", model_flags=(is_2ifc & incl_A)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TA", model_flags=(is_2ifc & incl_T & incl_A)
    )

  fit_info <- func_fit_abbrev(
    model_abbrev="TAiAxT", model_flags=(is_2ifc & incl_T & incl_A)
    )

  # return the fit information
  fit_info
}

# define some convenience functions for fitting a sequence of stan models 
# and a sequence of MLE models
fit_stan_model_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="stan",
  base_model_f=mD_stan,
  base_model_fallback_f=mD_stan_zCovar_noInteraction,
  base_model_id="_mD",
  modelfit_f=process_stan_fits
  )
fit_mle_model_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="mle",
  base_model_f=mD_mle,
  base_model_fallback_f=mD_mle_noInteraction,
  base_model_id="_mD",
  modelfit_f=process_mle_fit
  )
fit_stan_acuity_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="stan",
  base_model_f=acuity_mBase_stan,
  base_model_fallback_f=acuity_mBase_stan,
  base_model_id="_m"
  )
fit_mle_acuity_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="mle",
  base_model_f=acuity_mBase_mle,
  base_model_fallback_f=acuity_mBase_mle,
  base_model_id="_m",
  modelfit_f=process_mle_fit
  )
# Fit models using an alternate signal variable for WM tasks.
# The nchange signal variable uses a match of initial direction and number of changes
# to define signal trials.
fit_mle_nchange_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="mle",
  base_model_f=nchange_mBase_mle,
  base_model_fallback_f=nchange_mBase_mle,
  base_model_id="_m",
  modelfit_f=process_mle_fit
  )

# Fit models to target trials to assess recency effect
# NOTE: Support for this function is not fully implemented and summary functions
# will fail. Do not use this...
fit_mle_recency_sequence <- purrr::partial(
  fit_model_sequence,
  estimation_lib="mle",
  base_model_f=nosig_mD_mle,
  base_model_fallback_f=nosig_mD_mle,
  base_model_id="_m",
  modelfit_f=process_mle_fit
  )

fit_nonsignal_model <- function(data, model_abbrev, formula_str = NULL, base_f) {
  # Consumes a data frame, an abbrevation indicating which predictors to include or
  # a formula string, and function specifying the base model to fit. 
  # Fits the base model, then fits an updated model with the requested predictors.
  # Returns the fit and a model summary.

  # First fit the base model
  base_fit <- base_f(data)

  # If a formula string was supplied, use it. Otherwise build it from abbreviation.
  if (missing(formula_str) | is.null(formula_str)) {
    
    formula_str <- build_model_formula_str(
      existing_fit = base_fit,
      predictor_str = predictor_vars_from_abbrev(abbrev=model_abbrev, signal_var=NULL)
      )
  }
  # produce a fit to this extended model
  fit <- update(base_fit, formula. <- as.formula(formula_str))

  fit
}

# filtering function that may be called before fitting recency model
filter_recency_only <- function(dat) {
  # The only levels we care about are Last pos and Other pos
  dat <- dat %>% dplyr::filter(recency_fac %in% c("Other pos", "Last pos")) %>%
    as.data.frame()
  dat
}

# Fit a model to the signal trials and summarize the effect of the recency predictor.
# We could build up the formulas from abbrevations as we do elsewhere, but for expediency
# we will just require the caller to provide model and lsmeans formula strings.
fit_and_report_recency_effect <- function(data, model_formula_suffix_str = NULL, lsm_formula_str = NULL, base_f) {
  # First fit the base model
  base_fit <- base_f(data)
  # extract formula string
  base_formula_str <- Reduce(paste, deparse(formula(base_fit)))
  # append the provided formula suffix
  model_formula_str <- str_c(base_formula_str, model_formula_suffix_str, sep=" + ")
  # produce a fit to this extended model
  fit <- update(base_fit, formula. <- as.formula(model_formula_str))

  # generate lsmeans
  lsm_obj <- lsmeans::lsmeans(fit, specs = as.formula(lsm_formula_str))

  # Transform the lsmeans values by computing the probability for each grouping
  # plus or minus one or two standard errors
  smry_df <- as.data.frame(summary(lsm_obj))
  smry_df <- smry_df %>% dplyr::rename(
    PrbtHit = lsmean
    ) %>% dplyr::mutate(
    Accuracy = pnorm(PrbtHit),
    one_se_lo = PrbtHit - SE,
    one_se_hi = PrbtHit + SE,
    two_se_lo = PrbtHit - 2 * SE,
    two_se_hi = PrbtHit + 2 * SE,
    prob_one_se_lo = pnorm(one_se_lo),
    prob_one_se_hi = pnorm(one_se_hi),
    prob_two_se_lo = pnorm(two_se_lo),
    prob_two_se_hi = pnorm(two_se_hi),
    asymp.LCL.resp = pnorm(asymp.LCL), # convert the confidence interval supplied by lsmeans
    asymp.UCL.resp = pnorm(asymp.UCL)
    )

  # 
  # return the fitted model and the summary
  list(
    fit=fit,
    info=smry_df
    )
}

recency_effect_from_lsm_str <- function(fit, lsm_formula_str) {
  # Consumes a fitted model and a string with lsmeans spec.
  # Produces data frame of results including the accuracy.
  lsm_obj <- lsmeans::lsmeans(fit, specs = as.formula(lsm_formula_str))
  smry_df <- as.data.frame(summary(lsm_obj))
  smry_df <- smry_df %>% dplyr::rename(
    PrbtHit = lsmean
    ) %>% dplyr::mutate(
    Accuracy = pnorm(PrbtHit)
    ) %>% as.data.frame
  smry_df
}

recency_accuracy_main_effects <- function(fit, preds) {
  y <- lapply(preds, function(pred, thefit) {
    recency_accuracy_main_effect(thefit, pred)
  }, thefit=fit)
  # y is a list of lists. Extract the estimates and assign name to each
  acc <- unlist(sapply(y, function (x) x$Accuracy))
  est_names <- unlist(sapply(y, function (x) x$effect_name))
  # apply effect name column as the row names and return accuracy
  setNames(acc, est_names)
}

recency_accuracy_main_effect <- function(fit, pred) {
  lsm_obj <- lsmeans(fit, pred, transform="response")
  smry_df <- as.data.frame(summary(lsm_obj))
  # rename the accuracy col
  colnames(smry_df)[colnames(smry_df) == "prob"] <- "Accuracy"
  # extract the accuracy column and name each value with the level of the predictor
  acc <- smry_df$Accuracy
  # build up effect names
  effect_names <- sapply(smry_df[, pred], as.character)
  # if there are multiple predictor columns we concatenate with double underscore
  if (!is.null(ncol(effect_names))) {
    effect_names <- apply(effect_names, 1, function(x) str_c(x, collapse ="__"))
  }

  names(acc) <- effect_names
  # add effect name as separate column
  smry_df$effect_val_str <- effect_names
  smry_df$effect_name <- str_c("maineff::acc::", str_c(pred, collapse="_"), smry_df$effect_val_str, sep="__")
  smry_df[,c("Accuracy", "effect_name")]
}

# Return median and 95%CI from a bootstrap analysis.
boot_summary <- function(merBoot) {
  return(
    data.frame(
      fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
      lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
      upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}

# boot_res_quant functions perform summaries of bootstrap results

quantile_factory <- function(prob) {
  function(x) {
    as.numeric(quantile(x, probs=prob, na.rm=TRUE))
  }
}

boot_res_quant_list <- function() {
  boot_res_quant_list <- list()
  boot_res_quant_list$median <- quantile_factory(.5)
  boot_res_quant_list$mean <- function(x) mean(x, na.rm=TRUE)
  boot_res_quant_list$lwr05 <- quantile_factory(.025) # lower CI for alpha .05
  boot_res_quant_list$upr05 <- quantile_factory(.975) # upper CI for alpha .05
  boot_res_quant_list$lwr01 <- quantile_factory(.005) # lower CI for alpha .01
  boot_res_quant_list$upr01 <- quantile_factory(.995) # upper CI for alpha .01
  boot_res_quant_list
}

# apply a list of summary functions to a vector of estimates
boot_res_compute_summary_stats <- function(vec, func_list = boot_res_quant_list()) {
  smry <- sapply(func_list, function(f) {f(vec)})
  names(smry) <- names(func_list)
  smry
}

# Given results from bootstrap and 2 effect names, summarize the difference
# between the estimates.
boot_res_effect_difference <- function(smry, minuend, subtrahend) {
  est_vec <- smry[,minuend] - smry[,subtrahend]
  est_smry <- boot_res_compute_summary_stats(est_vec)
  est_smry
}

