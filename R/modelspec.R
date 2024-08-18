# Functions useful for parsing the simple model abbreviation language.

expand_abbrev_factory <- function(fieldname="predictor") {
  function(abbrev, defs) {
    defs[abbrev,fieldname]
  }
}

expand_abbrev <- expand_abbrev_factory(fieldname="predictor")
expand_abbrev_friendly <- expand_abbrev_factory(fieldname="friendly_name")

expand_interaction_factory <- function(separator="*", abbrev_func=expand_abbrev) {
  function(expr, defs) {
    # expr should be a 4-element string. We just want elements 2 and 4.
    parts <- unlist(strsplit(expr, ""))
    vars <- sapply(c(parts[2], parts[4]), abbrev_func, defs)
    paste(vars, collapse=separator)
  }
}

expand_interact_colon <- expand_interaction_factory(separator=":")
expand_interact_star <- expand_interaction_factory(separator="*")
expand_interact_by <- expand_interaction_factory(separator=" by ", abbrev_func=expand_abbrev_friendly)

separate_effects_from_abbrev_str <- function(expr, defs, x=list(fixed=NULL, interact=NULL), maineff_func=expand_abbrev, interact_func=expand_interact_colon) {
  # Separate the fixed effects from the interaction terms given a
  # model abbreviation string. The "ate" variable keeps track of
  # how many characters are consumed at each step.
  if (expr == "") {
    return(x)
  }
  token <- substr(expr, 1, 1)
  # is it uppercase?
  if (toupper(token) == token) {
    x$fixed <- c(x$fixed, maineff_func(token, defs))
    ate <- 1
  } else if (token == "i") {
    x$interact <- c(x$interact, interact_func(substr(expr, 1, 4), defs))
    ate <- 4
  } else if (token == "m") {
    # This is the start of each abbreviation string. Safe to ignore here.
    ate <- 1
  }
  separate_effects_from_abbrev_str(substr(expr, ate+1, nchar(expr)), defs, x=x, maineff_func=maineff_func, interact_func=interact_func)
}

extract_model_spec_from_fitname <- function(fitname) {
  # The fit name has several components separated by underscores.
  # We only care about the last component, and also want to strip
  # off the leading "m" from the abbreviation.

  # split by underscore
  parts <- unlist(str_split(fitname, "_"))
  # retain only the last part
  abbrev_str <- parts[length(parts)]
}

extract_formula_str_from_fit <- function(existing_fit) {
  # extract the formula from the existing fit
  existing_formula <- formula(existing_fit)
  # convert formula to a formula string
  src_formula_str <- Reduce(paste, deparse(existing_formula))
}

replace_default_signal_var <- function(x, sig_var, default_sig_var="is_old") {
  # Only want to perform a replacement if the default_sig_var is present with a terminating
  # character such as space, colon, etc. This is necessary to avoid replacing "is_old" in
  # signal variables that start with that string. This does not handle a possible case where
  # the new signal var ends with the default_sig_var string.
  if (str_detect(x, str_c(default_sig_var, "[:\\s+*]"))) {
    stringr::str_replace_all(x, default_sig_var, sig_var)
  } else {
    x
  }
}

build_predictor_str <- function(predictors, interacts_with_signal_var, nested_under_signal_var, signal_var) {
  if (!is.null(signal_var)) {
    pred_stubs <- mapply(function(pred, include_interact, is_nested, sig_var) {
      if (include_interact) {
        str_c(pred, str_c(pred, sig_var, sep=":"), sep=" + ")
      } else if (is_nested) {
        # Option to include only interaction, no main effect
        str_c(sig_var, pred, sep=":")
      } else {
        pred
      }
    }, predictors, interacts_with_signal_var, nested_under_signal_var, signal_var)
  } else {
    # No signal variable, so nothing for predictors to interact with. Pass along
    # the supplied predictors.
    pred_stubs <- predictors
  }
  predictor_str <- str_c(pred_stubs, collapse=" + ")
  predictor_str
}

predictor_vec_from_abbrev <- function(abbrev_str, pred_info=PREDICTOR_INFO()) {
  # Consumes abbreviation string and predictor information data frame.
  # Returns a list of all predictor variables in the main effects and
  # interactions.

  effects_list <- separate_effects_from_abbrev_str(abbrev_str, defs=pred_info)
  pred_vec <- effects_list %>%
    unlist %>% # collapse main effects and interactions to one vector
    stringr::str_split(pattern="\\:") %>% # split interaction strings
    unlist %>%
    unique
  pred_vec
}

append_formula_from_abbrev <- function(existing_fit, formula_suffix) {
  # Consumes an existing fit and a fragment of a formula as a string.
  # Extracts the formula from the existing fit and appends the supplied
  # suffix. Returns a string.

  src_formula_str <- extract_formula_str_from_fit(existing_fit)
  # build a string that includes the original formula with the suffix appended
  str_c(src_formula_str, formula_suffix, sep=" + ")
}

build_formula_suffix_from_abbrev <- function(abbrev_str, pred_info=PREDICTOR_INFO(), signal_var=NULL) {
  # Consumes an abbreviation string and data frame with predictor information.
  # Parses the abbreviation string and separates out main effects and
  # interactions. Uses the supplied predictor information to determine whether
  # or not the signal variable should be included in each term based on the
  # predictor.

  # identify the main effect predictors and interactions
  effects_list <- separate_effects_from_abbrev_str(abbrev_str, defs=pred_info)

  # Create string for the main effects
  main_preds_msk <- pred_info$predictor %in% effects_list$fixed
  main_info <- pred_info[main_preds_msk,]
  main_effect_str <- build_predictor_str(main_info$predictor,
                                       interacts_with_signal_var=main_info$interaction_relevant,
                                       nested_under_signal_var=main_info$nested_within_sigvar,
                                       signal_var=signal_var)

  # The interaction terms are a bit more complex. For each term, break
  # the interaction into the component predictors. Note that interactions
  # greater than 2-way are not supported.
  inter_str_list <- lapply(effects_list$interact, function(x) {
                             # split into 2 predictors
                             preds <- str_split(x, pattern="\\*")
                             msk <- pred_info$predictor %in% preds
                             both_interact <- all(pred_info[msk,"interaction_relevant"])
                             both_nested <- all(pred_info[msk,"nested_within_sigvar"])
                             # Now create the formula string for this term
                             y <- build_predictor_str(x, interacts_with_signal_var=both_interact, nested_under_signal_var=both_nested, signal_var=signal_var)
                             y
                                       }
  )
  # Collapse the list into a single string
  inter_effect_str <- str_c(inter_str_list, collapse=" + ")
  # Combine main effects and interaction strings
  str_c(main_effect_str, inter_effect_str, sep=" + ")
}

build_model_formula_str <- function(existing_fit = NULL, predictors, 
  interacts_with_signal_var = rep_len(TRUE, length(predictors)), 
  nested_under_signal_var = rep_len(FALSE, length(predictors)), 
  signal_var=NULL, predictor_str=NULL) {
  # Consumes an existing model fit. Extracts the formula, adds the requested predictors.
  # Produces a formula string.
  # The interacts_with_signal_var is a vector of logicals with the same length as the predictors. Indicates
  # whether the predictor should include interaction with the signal variable.

  # TODO: call extract_formula_str_from_fit() here
  # extract the formula from the existing fit
  existing_formula <- formula(existing_fit)
  
  # convert formula to a formula string
  src_formula_str <- Reduce(paste, deparse(existing_formula))

  # determine the signal variable if needed
  if (is.null(signal_var)) {
    signal_var <- get_signal_var_from_fit(existing_fit)
  }

  # If no predictor string provided, build one from the other arguments
  if (is.null(predictor_str)) {
    # prepare formula stubs for each predictor
    predictor_str <- build_predictor_str(predictors, interacts_with_signal_var, nested_under_signal_var, signal_var)
  } else {
    # A predictor string was provided, but it might contain a default signal variable. Replace any
    # instances of the default with the current signal variable.
    predictor_str <- replace_default_signal_var(predictor_str, sig_var=signal_var)
  }
  
  # build a string that includes the original formula with the new predictors appended
  str_c(src_formula_str, predictor_str, sep=" + ")
}

predictor_vars_from_abbrev <- function(abbrev, pred_info=PREDICTOR_INFO(), signal_var=NULL) {
  # Consumes a string, each letter refers to a predictor to include.
  # Produces a string with all predictors suitable for a formula

  # Mask indicates which predictors are requested
  pred_msk <- sapply(pred_info$abbrev, function(x, abbrev) {str_detect(abbrev, x)}, abbrev=abbrev)

  # Build predictor string to include the reqested predictors and possible interactions
  pred_str <- build_predictor_str(
    predictors=pred_info[pred_msk, "predictor"], 
    interacts_with_signal_var=pred_info[pred_msk, "interaction_relevant"],
    nested_under_signal_var=pred_info[pred_msk, "nested_within_sigvar"],
    signal_var=signal_var
    )
  pred_str
}

