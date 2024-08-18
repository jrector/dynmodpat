save_mincachefile <- function(outfile=minimaldatafile, dat=list(dat_by_stub=dat_by_stub, acuity_by_stub=acuity_by_stub)) {
  saveRDS(object = dat, file = outfile)
}

get_cache_filename <- function(project_root=projectroot, analysis_prefix=NULL, outfile_stub=FILE_SUFFIX()) {
  outputfile <- str_c(project_root, "cache/SdtModelFits_", str_c(analysis_prefix, outfile_stub, sep="_"), ".rds")
  outputfile
}

read_cached_fits <- function(outputfile=get_cache_filename()) {
  if (file.exists(outputfile)) {
    cat("loading cached model fits...\n")
    fits <- readRDS(outputfile)
  } else {
    fits <- list()
  }
  fits
}

get_boot_cache_filename <- function(project_root=projectroot, analysis_prefix=NULL, exp_stub=NULL, fitname=NULL) {
  bootcachefile <- str_c(project_root, "cache/bootout_", str_c(analysis_prefix, exp_stub, fitname, sep="_"), ".rds")
}

write_results_file <- function(results, outputfile) {
  cat("saving to ", outputfile, "\n")
  saveRDS(results, outputfile)
}

get_mle_fit_by_name <- function(fits, exp_stub, fitname) {
  # return the model fit given a list of fits and the name of the desired model
  fit <- fits[[exp_stub]]$mle[[fitname]]$fit
  fit
}

process_experiments <- function(
  experiment_stubs = exp_stubs,
  data_by_stub = dat_by_stub,
  outfile_stub = FILE_SUFFIX(),
  project_root = projectroot,
  filter_f = NULL,
  analysis_prefix = NULL,
  excluded_names = NULL,
  fit_mle = TRUE,
  fit_stan = FALSE,
  mle_batch_f = fit_mle_model_sequence,
  stan_batch_f = fit_stan_model_sequence,
  save_output = FALSE,
  force_model_refit = FALSE
  ) {
  # experiment_stubs is a vector of experiment identifiers.
  # data_by_stub is a list associating the identifiers to data frames with
  # experiment data.
  # filter_f is an optional function that applies a filter to the experiment
  # data, useful for fitting models to a subset of the data
  # analysis_prefix is an optional prefix to include in the output, typically
  # used in conjunction with filter_f to work on a subset of data
  # excluded_names can be used to remove some analyses from the results. Note
  # that this is not efficient in the current implementation, because the fit
  # is performed, and then excluded (instead of skipped entirely).

  # Define the file that will contain the analysis results.
  outputfile <- get_cache_filename(project_root, analysis_prefix, outfile_stub)
  # outputfile <- str_c(project_root, "cache/SdtModelFits_", str_c(analysis_prefix, outfile_stub, sep="_"), ".rds")
  if (save_output | !file.exists(outputfile)) {
    cat("output file will be ", outputfile, "\n")
  } else {
    cat("will not overwrite existing cache file ", outputfile, "\n")
  }

  # Container to store all of the fitted models and summaries
  fits <- read_cached_fits(outputfile=outputfile)
  # If any of the provided experiment stubs do not have corresponding
  # entries in "fits", create placeholders for them.
  new_stubs <- setdiff(experiment_stubs, names(fits))
  fits[new_stubs] <- NULL

  # Iterate over the experiment stubs to process each requested data set.
  # No harm in using for-loop since we are not growing a data structure.
  for (stub in experiment_stubs) {
    # filter the data if needed
    if (is.null(filter_f)) {
      # no filter, so retain all raw data
      filtered_data <- data_by_stub[[stub]]
    } else {
      # apply the filter function
      filtered_data <- filter_f(data_by_stub[[stub]])
    }

    # define stub to use for report names for this analysis
    report_stub <- str_c(analysis_prefix, stub, sep="_")

    if (fit_mle) {
      if (is.null(fits[[stub]]) | is.null(fits[[stub]]$mle)) {
        fits[[stub]]$mle <- list()
      }
      fits[[stub]]$mle <- mle_batch_f(
        dat=filtered_data,
        stub=stub,
        fit_info=fits[[stub]]$mle,
        force_refit=force_model_refit
        )
      fits[[stub]]$mle <- exclude_named_fits_from_list(fits[[stub]]$mle, excluded_names)
      fname <- str_c(projectroot, "reports/mlecompare_", report_stub, ".txt")
      compare_mle_fits_from_list(fits[[stub]]$mle, fname)
    }
  }

  #######
  # save results
  #######
  # Only save if flag is set or if file does not exist.
  if (save_output | !file.exists(outputfile)) {
    cat("saving to ", outputfile, "\n")
    saveRDS(fits, outputfile)
  }

  # return the fit information
  fits
}

update_sternberg_mle_reports <- function(fits, data, analysis_prefix = NULL) {
  # Consumes a list of fits organized by experiment.
  # Runs the function that processes the fitted model to compute summaries.
  # Products a revised list of fits with updated reports attached
  # Note that this does not touch the cache file of model fits, so if you want
  # to retain the results, save the cache file manually after running this.
  stubs <- names(fits)
  base_vars <- c("respondedsame", "is_old", "participant_id", "trial.sumsim.trltype.z")
  for (stub in stubs) {
    # define stub to use for report names for this analysis
    report_stub <- str_c(analysis_prefix, stub, sep="_")
    # Process the existing fits. Assumes a fit exists for each experiment stub.
    fits[[stub]]$mle <- process_mle_fits(
      fit_list = fits[[stub]]$mle,
      base_vars = base_vars,
      data = data[[stub]]
      )

    fname <- str_c(projectroot, "reports/mlecompare_", report_stub, ".txt")
    compare_mle_fits_from_list(fits[[stub]]$mle, fname)
  }
  fits
}

process_bootstrap_for_named_fit <- function(analysis_prefix, fit, exp_stub, fitname, row_spec=NULL) {
  # exp_stub is redundant because exp name is part of the model fitname
  outputfile <- get_boot_cache_filename(analysis_prefix=analysis_prefix, exp_stub=NULL, fitname=fitname)
  if (file.exists(outputfile)) {
    cat("Skipping ", fitname, "...\n")
  } else {
    cat("Starting bootstrap ", fitname, "...\n")
    # run and save the bootstrap results
    results <- run_bootstrap_for_fit(fit)
    write_results_file(results, outputfile)
  }
}

run_bootstrap_for_fit <- function(fit, nsim=NBOOTSTRAP_REPLICATES(), func=boot_extract_stats) {
  boot_out <- bootMer(fit, FUN=func, nsim=nsim)
  boot_smry <- boot_summary(boot_out)
  results <- list(boot_out=boot_out, boot_summary=boot_smry)
}

# create functions that return mappings between func values in the plot
# spec CSV and a target variable
plot_value_from_name <- function(values) {
  function(func_name) {
    switch(func_name,
      plot_dprime_by_block = values[1],
      plot_bias_by_block = values[2],
      plot_accuracy_by_block = values[3],
      plot_accuracy_density = values[4]
      )
  }
}

plot_func_from_name <- plot_value_from_name(c(group_plot_sdt_d, group_plot_sdt_c, group_plot_sdt_acc, group_plot_density))
# Note that for the density plot (position 4) any valid effect name is suitable here
plot_effectprefix_from_name <- plot_value_from_name(c("d_groupeff::grp::", "bias_maineff::c_grp::", "acc_maineff::group::", "acc_maineff::group::"))

plot_bootstrap_for_named_fit <- function(analysis_prefix, fit, exp_stub, fitname, row_spec) {
  bootcachefile <- get_boot_cache_filename(analysis_prefix=analysis_prefix, exp_stub=NULL, fitname=fitname)
  if (file.exists(bootcachefile)) {
    boot_info <- readRDS(bootcachefile)
    nonbase <- get_nonbase_predictors_from_fit(fit)
    resp_var <- get_response_var_from_fit(fit)
    sig_var <- get_signal_var_from_fit(fit)
    grouping_id_info <- get_grouping_str_for_bootout(boot_info$boot_out, preds=nonbase)
    # build names that match the rownames in the bootstrap output
    eff_prefix <- plot_effectprefix_from_name(row_spec$func)
    grouping_id_info$group_symbol <- str_c(eff_prefix, rownames(grouping_id_info), sep="")
    # extract data frame of estimated marginal mean summary for each group effect from 
    # the boostrap summary (median, 95%CI)
    msk <- rownames(boot_info$boot_summary) %in% grouping_id_info$group_symbol
    group_emms <- boot_info$boot_summary[msk,, drop = FALSE]
    # pull out the rowname as a column to make it easier to join
    group_emms$group_symbol <- rownames(group_emms)
    # attach the separate predictor values to each EMM entry
    group_emms <- left_join(group_emms, grouping_id_info, by="group_symbol")
    bootstrap_data <- boot_info$boot_out$data
    # If the signal variable is not present in the data that was used for bootstrapping,
    # it should be excluded from the summary function by setting it to null here.
    if (!sig_var %in% colnames(boot_info$boot_out$data)) {
      sig_var <- NULL
    }
    # summarize the trial-level data for possible inclusion in the plot
    trial_smry <- grp_data_summary(raw_data=bootstrap_data, preds=nonbase, resp_var=resp_var, sig_var=sig_var)
    plt_func <- plot_func_from_name(row_spec$func)[[1]]
    plt <- plt_func(model_smry=group_emms, trial_smry=trial_smry, preds=nonbase, plot_spec=row_spec)
    # save to file
    pdf_name <- str_c(projectroot, "plots/", analysis_prefix, "_", row_spec$plot_id, ".pdf")
    pdf(pdf_name, width=row_spec$width, height=row_spec$height)
    print(plt)
    dev.off()
  } else {
    cat ("skipping ", fitname, "\n")
    plt <- NULL
  }
  plt
}

report_boot_summary <- function(analysis_prefix, fit, exp_stub, fitname, row_spec) {
  bootcachefile <- get_boot_cache_filename(analysis_prefix=analysis_prefix, exp_stub=NULL, fitname=fitname)
  if (file.exists(bootcachefile)) {
    boot_info <- readRDS(bootcachefile)

    # increase width of output
    wid <- options()$width
    options(width = 120)

    txt <- sprintf("Report for %s", fitname)
    txt <- c(txt, sprintf("Analysis %s", analysis_prefix))
    txt <- c(txt, "----------")
    txt <- c(txt, capture.output(print(boot_info$boot_summary)))
    txt <- c(txt, "----------")

    options(width = wid)

    # return a vector of character strings
    txt
  } else {
    cat ("skipping ", fitname, "\n")
  }
}

report_model_prettyprint <- function(analysis_prefix, fit, exp_stub, fitname, row_spec) {
  pred_info <- PREDICTOR_INFO()
  # Given a model fit, extract the formula as a string. Then go through each possible
  # predictor and do a find / replace to give each predictor a friendly name for print.
  formula_str <- get_formula_str_from_fit(fit=fit)
  txt <- formula_str
  # Iterate over the various predictor variables
  for(i in seq_len(nrow(pred_info)))
  {
      txt <- gsub(pred_info[i,"predictor"], pred_info[i, "friendly_name"], txt)
  }
  # Spaces are inconsistent. Clean them up by first removing all existing spaces,
  # then inserting spaces around the relevant operators.
  txt <- gsub("\\s+", "", txt)
  # Move random effects statement to the end of the expression.
  # NOTE: this only addresses expressions with a single random effect paren grouping
  # that appears somewhere in the middle of RHS of the expression.
  # \\+\\((.+)\\) matches a plus sign followed by statement in parens.
  txt <- gsub("^(.+)\\+\\((.+)\\)(.+)$", "\\1\\3+(\\2)", txt)
  # spaces around any of the following: plus tilde pipe
  txt <- gsub("(\\+|~|\\|)", " \\1 ", txt)
  # prepend the analysis and experiment identifiers
  paste(analysis_prefix, exp_stub, fitname, txt, "", "", sep="\n")
}

report_predictor_info <- function(pred_info=PREDICTOR_INFO()) {
  # Consumes data frame with predictor information.
  # Produces a tab-delimited string with a report-friendly table of
  # key information about each predictor.

  # load the full data file so we have access to levels and type information
  # for each predictor
  load(system.file("extdata", "alltrials.rda", package = "dynmodpat"))

  # build data with predictor info of interest
  report_info <- pred_info[,c("predictor", "friendly_name")]
  # Special variables that should be treated separately have abbreviations that
  # start with a zero.
  msk_special <- str_detect(pred_info$abbrev, "^0")
  # add columns for type and levels for each predictor
  report_info$rtype <- sapply(report_info$predictor, function(x) typeof(trials[[x]]))
  report_info$is_factor <- sapply(report_info$predictor, function(x) is.factor(trials[[x]]))
  report_info$type <- ""
  report_info$levels <- ""
  report_info$notes <- ""  # placeholder for any comments about the variables
  report_info[!msk_special,]$type <- sapply(report_info[!msk_special,]$predictor, function(x) ifelse(is.factor(trials[[x]]), "Categorical", "Quantitative"))
  report_info[msk_special,]$type <- sapply(report_info[msk_special,]$predictor, function(x) ifelse(is.factor(trials[[x]]), "Categorical", "Binary"))
  report_info[!msk_special,]$levels <- sapply(report_info[!msk_special,]$predictor, function (x) ifelse(is.factor(trials[[x]]), paste(levels(trials[[x]]), collapse="; "), ""))
  # return a tab-delimited string of key columns
  write.table(x=report_info[,c("friendly_name", "type", "levels", "notes")], file="", sep="\t", quote=FALSE, row.names=FALSE)
}

fit_recency_model_for_rowspec <- function(dat, row_spec) {
  model_func <- get(row_spec$model_func)
  # Get filtering function from spreadsheet
  filter_func <- get(row_spec$filter_func)
  # apply it
  data_to_fit <- filter_func(dat)

  rec_fit <- fit_and_report_recency_effect(
    data = data_to_fit, model_formula_suffix_str = row_spec$predictor_str,
    lsm_formula_str = row_spec$lsm_str, base_f = model_func)
}

analyze_recency_for_named_fit <- function(analysis_prefix, fit, exp_stub, fitname, row_spec=NULL) {
  # exp_stub is redundant because exp name is part of the model fitname
  outputfile <- get_boot_cache_filename(analysis_prefix=analysis_prefix, exp_stub=NULL, fitname=fitname)
  if (file.exists(outputfile)) {
    cat("Skipping ", fitname, "...\n")
  } else {
    # extract full experiment data from existing fit
    dat <- fit@frame
    # fit the model and capture output
    rec_fit <- fit_recency_model_for_rowspec(dat, row_spec)
    # define a summary model for this bootstrap
    smry_func <- boot_extract_recency_factory(lsm_formula_str=row_spec$lsm_str)

    cat("Starting bootstrap ", fitname, "...\n")
    # run and save the bootstrap results
    results <- run_bootstrap_for_fit(rec_fit$fit, func=smry_func)
    write_results_file(results, outputfile)
  }
}

analyze_effect_diffs_from_boot <- function(analysis_prefix, exp_stub, fitname, row_spec=NULL) {
  # For a given row of the CSV, compute the requested effect
  bootcachefile <- get_boot_cache_filename(analysis_prefix=analysis_prefix, exp_stub=NULL, fitname=fitname)
  if (file.exists(bootcachefile)) {
    boot_info <- readRDS(bootcachefile)
    # In the future, the function to apply here could be specified in the CSV if
    # additional flexibility is needed.
    smry <- boot_res_effect_difference(boot_info$boot_out$t, minuend=row_spec$minuend, subtrahend=row_spec$subtrahend)
    list(summary=smry, effect=str_c(row_spec$ana_prefix, row_spec$model, row_spec$effect_id, sep=":"), title=row_spec$title)
  } else {
    NULL
  }
}

process_fits_factory <- function(func) {
  function(ana_prefix, to_process) {
    # Consumes the analysis prefix as a string, and a data frame that
    # contains the experiment stub and fitname to bootstrap.
    # Runs the bootstrap and saves the results to a file.

    # annoying workaround
    if (ana_prefix=="" | is.na(ana_prefix)) {
      ana_prefix <- NULL
    }

    # read in the cached model fits for the selected analysis
    fits <- read_cached_fits(outputfile=get_cache_filename(analysis_prefix=ana_prefix))

    # apply the function to each row
    lapply(seq(nrow(to_process)), function (idx) {
      row <- to_process[idx,]
      fit <- get_mle_fit_by_name(fits, exp_stub=row$exp_stub, fitname=row$model)
      # the row_spec contains other information from the associated row in the plot specs
      func(analysis_prefix=ana_prefix, fit=fit, exp_stub=row$exp_stub, fitname=row$model, row_spec=row)
    })
  }
}

# Note that this could be combined with process_fits_factory if that included
# an option to only load the fitted model when required. There is no need to
# iterate over the analysis prefixes here, other than to retain compatibility
# with process_fits_factory.
process_boot_factory <- function(func) {
  function(ana_prefix, to_process) {
    # annoying workaround
    if (ana_prefix=="" | is.na(ana_prefix)) {
      ana_prefix <- NULL
    }
    # apply the function to each row
    lapply(seq(nrow(to_process)), function (idx) {
      row <- to_process[idx,]
      # the row_spec contains other information from the associated row in the plot specs
      func(analysis_prefix=ana_prefix, exp_stub=row$exp_stub, fitname=row$model, row_spec=row)
    })
  }
}

process_fits_for_boot <- process_fits_factory(process_bootstrap_for_named_fit)
process_fits_for_plots <- process_fits_factory(plot_bootstrap_for_named_fit)
process_fits_for_boot_report <- process_fits_factory(report_boot_summary)
process_fits_for_model_prettyprint <- process_fits_factory(report_model_prettyprint)
process_fits_for_recency_fit <- process_fits_factory(analyze_recency_for_named_fit)
process_boot_compute_diffs <- process_boot_factory(analyze_effect_diffs_from_boot)

read_csv_factory <- function(file_path) {
  function() {
    fcontents <- read.csv(file_path, stringsAsFactors = FALSE)
    fcontents
  }
}

read_effect_csv <- read_csv_factory(file_path=str_c(projectroot, "specs/effects_to_compute.csv"))
read_plot_csv <- read_csv_factory(file_path=str_c(projectroot, "specs/dissertation.csv"))

get_plot_specs <- function(keepStandardRowsOnly=TRUE, keepMinCols=FALSE) {
  pinfo <- read_plot_csv()
  if (keepStandardRowsOnly) {
    # filter out the non-standard rows
    msk <- pinfo$model_func == ""
    pinfo <- pinfo[msk,]
  }
  if (keepMinCols) {
  # only keep a few key columns
  pinfo <- pinfo[,c("ana_prefix", "exp_stub", "model")]
  }
  pinfo
}

remove_dupes_multiple_cols <- function(df, cols) {
  # Returns a data frame that has removed any rows that are duplicates based
  # on the values in the columns cols
  dupes <- duplicated(df[, cols])
  df[!dupes,]
}

####################
# Functions to run after the models have been fit.
####################

bootstrap_cis_from_spec <- function(plot_info = read_plot_csv(), ana_prefixes = NULL) {
  # bootstrap models found in the plots file
  boot_info <- plot_info[,c("ana_prefix", "exp_stub", "model")]
  # skip any models that have an model function associated with them
  # (these will be handled separately)
  msk_ready_to_boot <- plot_info$model_func == ""
  boot_info <- boot_info[msk_ready_to_boot,]

  if (is.null(ana_prefixes)) {
    ana_prefixes <- unique(boot_info$ana_prefix)
  }

  # apply a bootstrap function to each of the fits
  lapply(ana_prefixes, function (prefix) {
    process_fits_for_boot(ana_prefix=prefix, to_process=boot_info[boot_info$ana_prefix==prefix,])
  })

}

generate_plots_from_spec <- function(ana_prefixes = NULL) {
  plot_info <- get_plot_specs(keepStandardRowsOnly = FALSE)

  if (is.null(ana_prefixes)) {
    ana_prefixes <- unique(plot_info$ana_prefix)
  }

  plts <- lapply(ana_prefixes, function (prefix) {
    process_fits_for_plots(ana_prefix=prefix, to_process=plot_info[plot_info$ana_prefix==prefix,])
  })

}

generate_custom_plots <- function(dat_by_stub = NULL) {
  # Generate grayscale unimodal by carrier plot
  plt <- plot_similarity_glm_smoothed(
    dat=rbind(dat_by_stub$uni_a1v1_combined, dat_by_stub$uni_a2v2),
    the.fill = "nstim_fac", the.facet.var = c("is_old", "block_condition"), xlabel="Summed similarity",
    in_color = FALSE, ncol=4
  )
  pdf_name <- str_c(projectroot, "plots/", "UnimodalGraySummedSimilarityGLM-BlockCond.pdf")
  pdf(pdf_name, width=8, height=12)
  print(plt)
  dev.off()
}

generate_report_from_fits_factory <- function(fit_process_func, report_filename, keepStandardRowsOnly=TRUE) {
  function() {
    plot_info <- get_plot_specs(keepStandardRowsOnly = keepStandardRowsOnly)

    # We don't want to report any analysis+model combos more than once
    plot_info <- remove_dupes_multiple_cols(plot_info, c("ana_prefix", "model"))

    ana_prefixes <- unique(plot_info$ana_prefix)
    # apply a report function to each of the fits
    reports <- lapply(ana_prefixes, function (prefix) {
      fit_process_func(ana_prefix=prefix, to_process=plot_info[plot_info$ana_prefix==prefix,])
    })

    # dump the output to file
    sink(file=str_c(projectroot, "reports/", report_filename, sep=""))
    sapply(reports, function(x) {
      y <- unlist(x)
      cat(y, sep="\n")
    })
    sink()
  }
}

generate_boot_report_from_spec <- generate_report_from_fits_factory(
  fit_process_func=process_fits_for_boot_report,
  report_filename="boot_reports.txt"
  )

generate_model_prettyprint_from_spec <- generate_report_from_fits_factory(
  fit_process_func=process_fits_for_model_prettyprint,
  report_filename="model_formulas.txt",
  keepStandardRowsOnly=FALSE
  )

fit_recency_models_from_spec <- function(ana_prefixes = NULL) {
  # These models are specified in the plot_specs file.
  plot_info <- get_plot_specs(keepStandardRowsOnly=FALSE)
  # Filter to retain only the rows with an explicit model spec
  msk <- plot_info$model_func != ""
  plot_info <- plot_info[msk,]

  if (is.null(ana_prefixes)) {
    ana_prefixes <- unique(boot_info$ana_prefix)
  }

  rec_fits <- lapply(ana_prefixes, function(prefix) {
    process_fits_for_recency_fit(ana_prefix=prefix, to_process=plot_info[plot_info$ana_prefix==prefix,])
  })

}

compute_effects_from_spec <- function(eff_info = read_effect_csv(), ana_prefixes = NULL) {
  # drop any rows without a minuend for now
  eff_info <- eff_info[eff_info$minuend != "",]

  # for consistency with other methods, iterate over the analysis prefixes
  if (is.null(ana_prefixes)) {
    ana_prefixes <- unique(eff_info$ana_prefix)
  }

  eff_ests <- lapply(ana_prefixes, function(prefix) {
    process_boot_compute_diffs(ana_prefix=prefix, to_process=eff_info[eff_info$ana_prefix==prefix,])
  })

  # increase width of output
  wid <- options()$width
  options(width = 120)

  # write report to file
  sink(file=str_c(projectroot, "reports/boot_computed_effects.txt", sep=""))
  lapply(eff_ests, function (eff_est_list) {
    lapply(eff_est_list, function (x) {
      cat(x$effect, x$title, sep="\n")
      y <- capture.output(print(x$summary))
      cat(y, sep="\n")
      # determine if zero is in the interval
      lwr <- min(x$summary['upr05'], x$summary['lwr05'])
      upr <- max(x$summary['upr05'], x$summary['lwr05'])
      if (lwr < 0 & upr > 0) {
        cat("n.s.\n")
      } else {
        cat("* Significant at alpha .05\n")
      }
      cat("-----\n\n")
    })
  })
  sink()

  options(width = wid)
}

unattended_target_analysis <- function(filter_f = filter_no_bad_trials) {
  # Load the full data set including all trials and columns. This will be available
  # in the workspace as the "trials" data frame.
  load(system.file("extdata", "alltrials.rda", package = "dynmodpat"))

  # define the experiments to use by dissertation name
  stubs <- c("Exp 5", "Exp 6", "Exp 7")

  unatt_filt_f <- function(dat) {
    filtereddat <- dat %>%
      dplyr::filter(!is.na(is_old_unattended_stream)) %>%
      dplyr::filter(stream_type=="async") %>%
      as.data.frame
    filtereddat
  }

  unatt_fit_f <- function(dat) {
    glmer(respondedsame ~ is_old_unattended_stream + (1 | participant_id), 
          data=dat, 
          family=binomial(link="logit"), 
          control = glmerControl(
            optimizer = "nloptwrap",
            calc.derivs = FALSE,
            optCtrl = list(maxfun = 100000)
            )
          )
  }

  unatt_summarize_f <- function(fit) {
    # get contrast between 1 and 0 and return the estimate
    smry <- as.data.frame(summary(
      lsmeans::contrast(
        lsmeans(fit, ~ is_old_unattended_stream, at=list(is_old_unattended_stream=c(1, 0)
          ), transform="response"), method="pairwise")
      ))
    smry$estimate
  }

  unatt_run_f <- function(stub, filter_f = NULL) {
    # load the data based on the stub
    dat <- trials %>%
      dplyr::filter(exp_diss_name_fac == stub) %>%
      as.data.frame()
    if (!is.null(filter_f)) {
      # apply the optional filtering function
      dat <- filter_f(dat)
    }
    # filter, fit the model, extract the estimates
    filtered_dat <- unatt_filt_f(dat)
    fit <- unatt_fit_f(filtered_dat)
    run_bootstrap_for_fit(fit=fit, func=unatt_summarize_f)
  }


  unattcachefile <- str_c(projectroot, "cache/unattended_target_analysis", ".rds")

  if (file.exists(unattcachefile)) {
    unatt_targ_analysis <- readRDS(unattcachefile)
  } else {
    # process each data set
    unatt_targ_analysis <- lapply(stubs, unatt_run_f, filter_f = filter_f)
    names(unatt_targ_analysis) <- stubs

    # save the results
    write_results_file(unatt_targ_analysis, unattcachefile)
  }

  # write to file
  sink(file=str_c(projectroot, "reports/unattended_target_analysis.txt", sep=""))
  y <- lapply(unatt_targ_analysis, function(x) x$boot_summary)
  print(y)
  sink()

  unatt_targ_analysis
}

process_exps_for_analysis_prefix <- function(ana_prefix, func) {
  fits <- read_cached_fits(outputfile=get_cache_filename(analysis_prefix=ana_prefix))
  exp_stubs <- names(fits)
  sapply(exp_stubs, model_spec_report_for_exp, fits=fits, prefix=ana_prefix)
}

predictor_separation_from_fitname <- function(fitname) {
  pred_info <- separate_effects_from_abbrev_str(expr=extract_model_spec_from_fitname(fitname), defs=PREDICTOR_INFO, maineff_func = expand_abbrev_friendly, interact_func = expand_interact_by)
}

model_spec_report_for_exp <- function(exp_stub, fits, prefix=NULL) {
  aic_info <- sort_models_by_aic(fits[[exp_stub]]$mle)
  fit_stubs <- sapply(aic_info[[1]], as.character)
  pred_info <- sapply(fit_stubs, predictor_separation_from_fitname)
  # Build data frame with the AIC info and with print-friendly names
  # for the fixed effects and interactions in the model based on the
  # model abbreviation string.
  model_info <- data.frame(aic_info[, c("Modnames", "AICc", "Delta_AICc", "LL", "K")], fixed=NA_character_, interact=NA_character_, fitname=fit_stubs)
  model_info$fixed <- sapply(fit_stubs, function(x) {paste(pred_info[, x]$fixed, collapse=", ")})
  model_info$interact <- sapply(fit_stubs, function(x) {paste(pred_info[, x]$interact, collapse=", ")})
  # reorder columns for report
  model_report_data <- model_info[, c("fixed", "interact", "K", "AICc", "Delta_AICc", "LL")]
  colnames(model_report_data)[colnames(model_report_data)=="K"] <- "nParams"
  # replace numbers with 2-digit versions
  model_report_data$AICc <- sprintf("%.2f", model_report_data$AICc)
  model_report_data$Delta_AICc <- sprintf("%.2f", model_report_data$Delta_AICc)
  model_report_data$LL <- sprintf("%.2f", model_report_data$LL)
  cat("Model report for", exp_stub, " [", prefix, "]", "\n")
  write.table(x=model_report_data, file="", sep="\t", quote=FALSE)
  cat("\n\n")
}

run_model_spec_report <- function(ana_prefixes) {
  sapply(ana_prefixes, process_exps_for_analysis_prefix)
  cat("\nDone\n")
}
