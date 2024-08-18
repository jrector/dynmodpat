# Functions to extract information from fitted models.

##########
# formula strings and variables

get_formula_str_from_fit <- function(fit) {
  gsub("  ", "", Reduce(paste, deparse(formula(fit))))
}

get_signal_var_from_fit <- function(fit) {
  # extract the signal variable by looking for a match between our vector
  # of possible signal vars and the variables observed in the supplied fit
  fit_vars <- all.vars(formula(fit))
  supported_signal_vars <- c("is_old", "same_in_1st_pair", "is_old_nchange_and_direction")
  signal_var <- intersect(supported_signal_vars, fit_vars)
  signal_var
}

get_response_var_from_fit <- function(fit) {
  # response variable is the second element in the terms list
  resp_var <- as.character(terms.formula(formula(fit))[[2]])
  # another option is the first element of all.vars
  # all.vars(formula(fit))[[1]]
  resp_var
}

get_predictors_from_fit <- function(fit) {
  all.vars(delete.response(terms(formula(fit))))
}

get_nonbase_predictors_from_fit <- function(fit) {
  # This contains some model-specific knowledge that should be moved out of the function.
  # Define the base predictors in relation to the response and signal variables.
  base_vars <- list()
  # sternberg task
  base_vars$respondedsame <- list()
  base_vars$respondedsame$is_old <- c("is_old", "trial_sumsim_trltype_z", "participant_id")
  base_vars$respondedsame$is_old_nchange_and_direction <- c("is_old_nchange_and_direction", "trial_sumsim_trltype_z", "participant_id")
  # acuity task
  base_vars$resp_1st_interval <- list()
  base_vars$resp_1st_interval$same_in_1st_pair <- c("same_in_1st_pair", "participant_id", "avg_time_difference_z")

  signal_var <- get_signal_var_from_fit(fit)
  resp_var <- get_response_var_from_fit(fit)
  # get all predictor variables
  pred_vars <- get_predictors_from_fit(fit)
  # look up the list of base predictors
  relevant_base_vars <- base_vars[[resp_var]][[signal_var]]
  if (is.null(relevant_base_vars)) {
    stop("Cannot handle model variables in get_nonbase_predictors_from_fit")
  } else {
    # filter the predictors to return only the non-base variables
    return(pred_vars[!pred_vars %in% relevant_base_vars])
  }
}

get_predictor_pairs <- function(preds) {
  # returns data frame with unique pairs of predictors
  if (length(preds) > 1) {
    pred_pairs = data.frame(combn(preds, 2), stringsAsFactors=FALSE)
  } else {
    pred_pairs <- NULL
  }
  pred_pairs
}

##########
# effect estimates

reformat_lsm_dprime_group <- function(lsm_obj) {
  # Consumes output of lsmeans contrast that computes dprime by group.
  # Produces revised data frame with unnecessary entries removed, and with
  # entries labeled with useful columns.

  # NOTE: Current implemention does not work on contrast of contrasts. 

  # Convert to data frame
  df <- as.data.frame(summary(lsm_obj))
  # contrast column is a string with hyphen separating the sides of the contrast.
  # This function currently works with either one or two contrasts
  
  # split the contrast string and reconstruct the 1st and 2nd parts of main contrast
  contrast_mtx <- str_split(df$contrast, " - ", simplify = TRUE)
  # If there are 2 elements, there is only a single contrast. If there are 4 elements,
  # this is a contrast of contrasts.
  if (ncol(contrast_mtx) == 2) {
    contrast_parts <- contrast_mtx
  } else if (ncol(contrast_mtx) == 4) {
    # concat the parts into a 2 column matrix
    thesep <- "__"
    contrast_parts <- matrix(c(
      str_c(contrast_mtx[,1], contrast_mtx[,2], sep=thesep),
      str_c(contrast_mtx[,3], contrast_mtx[,4], sep=thesep)
        ), ncol = 2)
  } else {
    contrast_parts <- NULL
    # issue warning or error
  }

  # contrast_parts <- str_split_fixed(df$contrast, " - ", 2)
  left_side <- contrast_parts[, 1]
  right_side <- contrast_parts[, 2]
  # Identify entries where the left side is evaluated at positive number and right
  # side at negative These are the contrasts of interest.
  left_positive <- str_detect(left_side, ",0.5$")
  right_negative <- str_detect(right_side, ",-0.5$")
  msk_desired_contrast_order <- left_positive & right_negative
  # Drop the contrast value from the contrast name
  left_name <- str_replace(left_side, ",-?0.5$", "")
  right_name <- str_replace(right_side, ",-?0.5$", "")
  # Identify contrasts where left and ride sides match
  msk_sides_match <- left_name == right_name
  # Construct a cleaned-up contrast name. For now will just use left side
  df$contrast_name <- left_name

  # Return the contrasts of interest
  df[msk_desired_contrast_order & msk_sides_match,]
}

generate_contrast_name <- function(lsm_summary_df, preds, separator=",") {
  contrast_name <- sapply(seq(nrow(lsm_summary_df)), function(i, dat) {
    # extract the predictor values as characters
    tmp = sapply(dat[i,preds], as.character)
    # concatenate into a string to match the output from lsmeans
    str_c(as.vector(tmp), collapse=separator)
  }, dat=lsm_summary_df)
}

contrast_codes_2nd_half_minus_1st <- function(df, preds, scaleby=1) {
  # Consumes a data frame that contains lsm summary output and a vector of
  # predictor variable names.
  # Produces a list of contrast codes suitable for subtracting the first
  # half of the data frame from the second half (by row).

  ncontrasts <- nrow(df) / 2
  # The contrast codes will subtract each row in the first half of the data
  # frame from its corresponding row in the second half.
  ccodes_list <- lapply(1:ncontrasts, function(i, ncontrasts, scale) {
    # make vector of zeros
    x <- rep(0, ncontrasts*2)
    # position i is row to be subtracted from row in position i+ncontrasts
    x[i] <- -1
    x[i+ncontrasts] <- 1
    # return x scaled by the supplied value
    x * scale
  }, ncontrasts=ncontrasts, scale=scaleby)

  # Get predictor values from the original data frame, and use these to
  # assign predictor columns and a contrast name to each row. The order of
  # the values should match between these data frames.
  pred_df <- df[1:ncontrasts, preds, drop = FALSE]

  names(ccodes_list) <- sapply(seq_len(ncontrasts), function (i, df) {
    str_c(sapply(pred_df[i,], as.character), collapse=",")
  }, df=pred_df)

  ccodes_list
}

lsm_accuracy_by_group <- function(lsm_obj, preds) {
  # Consumes lsm_obj with probit false alarm and hit rates for each grouping.
  # preds is a vector of predictor names.
  # Since accuracy is the average of the hit and correct rejection rates, we
  # will need to compute CR from FA and take an average. We use lsmeans to
  # generate the desired estimate so we get usable standard errors. This requires
  # a special set of contrast codes.
  df <- as.data.frame(summary(lsm_obj))

  # The supplied lsm_obj has the false alarm rates in the first half of the
  # data frame, hit rate in second half. We want to generate half as many
  # contrasts
  ncontrasts <- nrow(df) / 2
  # The contrast codes will subtract the center of the noise distribution
  # from the center of the signal distribution, giving us the center point.
  # Note that using custom contrast codes in this manner is probably overkill
  # as this could be done with suitable lsmeans syntax.
  ccodes_list <- contrast_codes_2nd_half_minus_1st(df, preds, scaleby=0.5)

  # Now we can compute the contrast. No offset is needed for codes in this form.
  acc_by_grp_lsm <- NULL
  try(
    acc_by_grp_lsm <- lsmeans::contrast(
      lsm_obj, ccodes_list, by = NULL, offset = 0
      )
    )

  # Get predictor values from the original data frame, and use these to
  # assign predictor columns and a contrast name to each row. The order of
  # the values should match between these data frames.
  pred_df <- df[1:ncontrasts, preds, drop = FALSE]
  # only proceed if the previous contrast was successful
  if (!is.null(acc_by_grp_lsm)) {
    # combo data frame with accuracy values and predictor values
    acc_by_grp <- cbind(pred_df, as.data.frame(summary(acc_by_grp_lsm)))

    # Apply tranformations to get estimates as probabilities, using standard
    # errors to get intervals for reporting
      acc_by_grp <- acc_by_grp %>%
      # The inferential stats columns are not useful so drop them. 
      dplyr::select(-c(df, z.ratio, p.value)) %>%
      dplyr::rename(
        PrbtAccuracy = estimate
      ) %>% dplyr::mutate(
      Accuracy = pnorm(PrbtAccuracy),
      one_se_lo = PrbtAccuracy - SE,
      one_se_hi = PrbtAccuracy + SE,
      two_se_lo = PrbtAccuracy - 2 * SE,
      two_se_hi = PrbtAccuracy + 2 * SE,
      prob_one_se_lo = pnorm(one_se_lo),
      prob_one_se_hi = pnorm(one_se_hi),
      prob_two_se_lo = pnorm(two_se_lo),
      prob_two_se_hi = pnorm(two_se_hi)
      ) %>% as.data.frame()
  } else {
    acc_by_grp <- data.frame(Accuracy=NA)
  }

  acc_by_grp
}

lsm_sdt_d_by_group <- function(lsm_obj, preds) {
  # NOTE: this function has considerable code overlap with lsm_accuracy_by_group()
  # Consumes lsm_obj with estimates for each grouping.
  # preds is a vector of predictor names.
  # Sensitivity is the difference between estimates at the signal variable
  # at 0.5 and -0.5, and these estimates have been provided. We generate
  # constrast codes to do the proper subtraction for each group, and this
  # gets proper standard errors since we use lsmeans. Pairwise contrasts
  # can then be generated from the results.
  df <- as.data.frame(summary(lsm_obj))

  # The supplied lsm_obj has the -0.5 estimate in the first half of the
  # data frame, 0.5 estimate in second half. We want to generate half as many
  # contrasts
  ncontrasts <- nrow(df) / 2
  # The contrast codes will subtract the -0.5 estimate from the 0.5 estimate
  ccodes_list <- contrast_codes_2nd_half_minus_1st(df, preds)

  # Now we can compute the contrast
  d_by_grp_lsm <- NULL
  try(
    d_by_grp_lsm <- lsmeans::contrast(
      lsm_obj, ccodes_list, by = NULL
      )
    )

  # only proceed if the previous contrast was successful
  if (!is.null(d_by_grp_lsm)) {
    # compute the pairwise contrast
    d_group_pairwise_lsm <- lsmeans::contrast(
      d_by_grp_lsm, by = NULL, method="pairwise"
    )
    d_group_pairwise <- as.data.frame(summary(d_group_pairwise_lsm))
  } else {
    d_group_pairwise <- data.frame()
    d_group_pairwise_lsm <- NULL
  }

  d_group_pairwise_lsm
}

lsm_main_effect_sdt_d <- function(
  fit, preds, signal_var, contrast_at, effect_names=lapply(preds, str_c, collapse="__")
  ) {
  if (length(preds) > 0) {
    y <- lapply(preds, function(x, thefit) {
      lsmeans::contrast(
        lsmeans(thefit, x, by=signal_var, at=contrast_at)
        , method = "pairwise", by = NULL)
    }, thefit=fit)
    names(y) <- effect_names
    # filter out the unnecessary rows in each list item
    y <- lapply(y, reformat_lsm_dprime_group)
    y
  } else {
    NULL
  }
}

lsm_single_effect <- function(
  fit, preds, signal_var, contrast_at, effect_names=lapply(preds, str_c, collapse="__")
  ) {
  # lsmeans::lsmeans(fit, nonbase, by=signal_var, at=lsm_binary_mean)
  if (length(preds) > 0) {
    y <- lapply(preds, function(x, thefit) {
      lsmeans(thefit, x, by=signal_var, at=contrast_at)
    }, thefit=fit)
    y <- mapply(function(lsm, predvec) {
      smry <- as.data.frame(summary(lsm))
      # add a column of effect names based on predictors
      smry$contrast_name <- generate_contrast_name(smry, predvec)
      smry
    }, y, preds, SIMPLIFY = FALSE)
    names(y) <- effect_names
    y
  } else {
    NULL
  }
}

lsm_main_effect_diff_sdt_d <- function(fit, preds, signal_var, contrast_at) {
  if (length(preds) > 0) {
    y <- lapply(preds, function(x, thefit) {
      lsmeans::contrast(
        lsmeans::contrast(
          lsmeans(thefit, x, by=signal_var, at=contrast_at)
          , method = "pairwise"),
        method = "pairwise", by = NULL)
      }, thefit=fit)
    names(y) <- lapply(preds, str_c, collapse="__")
    # filter out the unnecessary rows in each list item
    y <- lapply(y, reformat_lsm_dprime_group)
    y
  } else {
    NULL
  }
}

named_estimates_from_listoflist_apply <- function(lol) {
  x2 <- sapply(seq(length(lol)), function(idx) {
    mylist <- list(contrast_name=lol[[idx]]$contrast_name, est=lol[[idx]]$estimate)
    names(mylist$est) <- str_c(names(lol)[[idx]], mylist$contrast_name, sep="::")
    mylist$est
    })
  unlist(x2)
}

lsm_bias_estimates <- function(
  fit, preds, signal_var, contrast_at, effect_names=lapply(preds, str_c, collapse=":"), prefix, estimate_colname="lsmean"
  ) {
  # Get the estimate
  est <- lsm_single_effect(fit, preds=preds, signal_var=signal_var, contrast_at=contrast_at, effect_names=effect_names)
  named_est <- named_estimates_from_listoflists(est, prefix=prefix, estimate_colname=estimate_colname)
  # flip the sign to put this in SDT c units
  named_est <- sapply(named_est, function(x) {-1 * x})
  named_est
}

lsm_accuracy_from_fit <- function(
  fit, preds, signal_var, effect_names=lapply(preds, str_c, collapse=":"), prefix
  ) {
  # set up a binary contrast
  contrast_at = list()
  contrast_at[[signal_var]] <- c(0, 1)
  # get accuracy estimates for each grouping of predictors
  y <- lapply(preds, function(x, thefit) {
    groups <- lsm_accuracy_by_group(
      lsm_obj=lsmeans(thefit, x, by=signal_var, at=contrast_at),
      preds=x
      )
  }, thefit=fit)
  # extract only the Accuracy column
  named_est <- named_estimates_from_listoflists(y, prefix=prefix, estimate_colname="Accuracy", contrast_colname="contrast")

  named_est
}

named_estimates_from_listoflists <- function(lol, prefix="", estimate_colname="estimate", contrast_colname="contrast_name") {
  # extract estimates from the nested lists, and assign a meaningful name to each
  estimates <- unlist(sapply(lol, function(x) {x[[estimate_colname]]}))
  contrast_vals <- unlist(sapply(lol, function(x) {x[[contrast_colname]]}))
  contrast_names <- rep(names(lol), times=sapply(lol, nrow))
  if (!is.null(estimates)) {
    names(estimates) <- str_c(prefix, contrast_names, contrast_vals, sep="::")
  }
  estimates
}

boot_extract_stats <- function(fit, preds=get_nonbase_predictors_from_fit(fit)) {
  # get various useful values from the fitted object
  signal_var <- get_signal_var_from_fit(fit)
  # preds <- all.vars(delete.response(terms(formula(fit))))
  # ignore any variables that were part of the base model
  # nonbase <- preds[! preds %in% base_vars]
  # define the contrast for SDT sensitivity
  sensitivity_contrast_at <- list()
  sensitivity_contrast_at[[signal_var]] = c(0.5, -0.5)
  # get all pairs of predictors
  if (length(preds) > 1) {
    pred_pairs = data.frame(combn(preds, 2), stringsAsFactors=FALSE)
    # assign a meaningful name to each column
    colnames(pred_pairs) <- sapply(pred_pairs, str_c, collapse="__")
  } else {
    pred_pairs <- NULL
  }
  

  # extract the fixed effects
  x1 <- fixef(fit)
  # extract interactions
  x2 <- lsm_main_effect_sdt_d(fit, preds=preds, signal_var=signal_var, contrast_at=sensitivity_contrast_at)
  x2.est <- named_estimates_from_listoflists(x2, prefix="d_maineff")
  x3 <- lsm_main_effect_diff_sdt_d(fit, preds=preds, signal_var=signal_var, contrast_at=sensitivity_contrast_at)
  x3.est <- named_estimates_from_listoflists(x3, prefix="d_maindiff")
  x4 <- lsm_main_effect_sdt_d(fit, preds=pred_pairs, signal_var=signal_var, contrast_at=sensitivity_contrast_at)
  x4.est <- named_estimates_from_listoflists(x4, prefix="d_interact")
  # Get estimates for group means, where a group is a unique combination of all
  # of the predictors of interest.
  x5 <- lsm_main_effect_sdt_d(fit, preds=list(group=preds), signal_var=signal_var, contrast_at=sensitivity_contrast_at, effect_names="grp")
  x5.est <- named_estimates_from_listoflists(x5, prefix="d_groupeff")
  # Extract SDT bias estimates for each group
  est_mean_at = list()
  est_mean_at[[signal_var]] <- list()
  lsm_at_opts <- LSM_AT_OPTIONS()
  est_mean_at[[signal_var]] <- lsm_at_opts$mean
  # bias for each main effect
  maineff_bias.est <- lsm_bias_estimates(fit, preds=preds, signal_var=signal_var, contrast_at=est_mean_at, prefix="bias_maineff")
  # TODO: estimate bias for each element of pred_pairs for 2-way interactions
  group_bias.est <- lsm_bias_estimates(fit, preds=list(group=preds), signal_var=signal_var, contrast_at=est_mean_at, effect_names=c("c_grp"), prefix="bias_maineff")

  # get accuracy estimates for each main effect
  maineff_accuracy.est <- lsm_accuracy_from_fit(fit, preds=preds, signal_var=signal_var, prefix="accuracy_maineff")
  # get accuracy estimates for each group when all predictors are entered simultaneously
  group_accuracy.est <- lsm_accuracy_from_fit(fit, preds=list(group=preds), signal_var=signal_var, effect_names=c("acc_grp"), prefix="acc_maineff")
  # get accuracy estimates for two-way interaction groups if appropriate
  if (!is.null(pred_pairs)) {
    interact_accuracy.est <- lsm_accuracy_from_fit(fit, preds=pred_pairs, signal_var=signal_var, effect_names=lapply(pred_pairs, str_c, collapse=":"), prefix="acc_2wayeff")
  } else {
    interact_accuracy.est <- NULL
  }

  x <- c(x1, x2.est, x3.est, x4.est, x5.est, maineff_bias.est, group_bias.est, maineff_accuracy.est, group_accuracy.est, interact_accuracy.est)
  x
}

boot_extract_recency_factory <- function(lsm_formula_str = NULL) {
  function(fit) {
    smry_df <- recency_effect_from_lsm_str(fit, lsm_formula_str)
    # base predictors in recency models
    base_model_preds <- c("participant_id", "trial_sumsim_trltype_z")
    pred_vars <- get_predictors_from_fit(fit)
    nonbase <- pred_vars[!pred_vars %in% base_model_preds]
    # fixed effects
    fixed_eff <- fixef(fit)
    # main effect of each predictor
    main_eff <- recency_accuracy_main_effects(fit, nonbase)
    pred_pairs <- get_predictor_pairs(nonbase)
    if (!is.null(pred_pairs)) {
      twoway_eff <- recency_accuracy_main_effects(fit, pred_pairs)
    } else {
      twoway_eff <- NULL
    }
    # accuracy column with each estimate given a useful name
    groupacc <- setNames(smry_df$Accuracy, str_c("acc_maineff::group::", generate_contrast_name(smry_df, nonbase), sep=""))
    c(fixed_eff, main_eff, twoway_eff, groupacc)
  }
}

get_grouping_str_for_bootout <- function(bootout, preds) {
  # Each unique combination of non-base predictor values defines a grouping
  # in the data. Return a vector of strings with the groupings for the 
  # supplied bootstrap output.
  # extract the factor values as characters
  # Note that we use simplify=FALSE to support different lengths of preds
  factor_vals <- sapply(preds, function(x) unique(bootout$data[,x]), simplify=FALSE)
  factor_val_char <- sapply(factor_vals, as.character, simplify=FALSE)
  # data frame of unique combinations
  fac_val_cols <- expand.grid(factor_val_char, stringsAsFactors = FALSE)
  rownames(fac_val_cols) <- sapply(seq(nrow(fac_val_cols)), function(x) str_c(fac_val_cols[x,], collapse=","))
  fac_val_cols
}

in_sample_error_rate <- function(dat = NULL, fit = NULL) {
  # Quick check of the accuracy of the model. Tests predictions generated
  # by the model against the observed responses. Returns a proportion of
  # incorrect classifications.
  resp_var <- get_response_var_from_fit(fit)
  predictions <- predict(fit)
  # If a data frame with response values is not available, use the responses
  # that were provided when the model was fit
  if (is.null(dat)) {
      observed_vals <- model.frame(fit)[[resp_var]]
  } else {
      observed_vals <- dat[[resp_var]]
  }
  # It is an error if prediction is > 0.5 but observation is 0, 
  # or if prediction is < 0.5 but observervation is 1.
  err_rate <- mean((predictions > 0.5 & observed_vals==0) | (predictions <= 0.5 & observed_vals==1))
  err_rate
}

