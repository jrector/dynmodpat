get_text_margin_pad <- function() {8}
get_errorbar_size <- function() {0.5}
get_lineplot_size <- function() {0.75}
get_theme_basesize <- function() {10}

theme_jr <- ggplot2::theme(
  aspect.ratio = 1.0,
  panel.border = ggplot2::element_rect(fill=NA, color="black"),
  panel.grid.major = ggplot2::element_line(color="#CDCDCD", linetype="dotted"),
  panel.background = ggplot2::element_blank(),
  legend.background = ggplot2::element_blank(),
  legend.key = ggplot2::element_blank(),
  legend.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
  plot.title = ggplot2::element_text(margin=ggplot2::margin(0,0,get_text_margin_pad(),0)),
  axis.title.x = ggplot2::element_text(margin=ggplot2::margin(get_text_margin_pad(),0,0,0)),
  axis.title.y = ggplot2::element_text(margin=ggplot2::margin(0,get_text_margin_pad(),0,0)),
  text = ggplot2::element_text(size=get_theme_basesize(), family="Helvetica")
  )

####################
# The code below needs significant cleanup
####################

get_custom_palette <- function() {
  custom_palette <- list(wm = list(), acuity = list())
  custom_palette$wm$block_condition <- grDevices::gray.colors(4)
  custom_palette$wm$stream_type <- grDevices::gray.colors(3)
  custom_palette$wm$nstim_fac <- grDevices::gray.colors(3)
  custom_palette$wm$recency_fac <- grDevices::gray.colors(4)
  custom_palette$acuity$block_condition <- grDevices::gray.colors(5)
  custom_palette$acuity$ndirection_changes_fac <- grDevices::gray.colors(4)
  custom_palette$acuity$probe_transform <- grDevices::gray.colors(3)
  custom_palette
}

get_eps <- function() {
  1e-5
}

# Define labels and levels for variables
get_label_info_all_variables <- function() {
  label_info_by_var = list()
  label_info_by_var$block_condition = list()
  label_info_by_var$block_condition$original_names <- c("a1freq", "v1size", "a2gain", "v2freq", "a3spec", "a4harmo")
  label_info_by_var$block_condition$abbrevs <- c(
    "A: pitch", "V: size", "A: gain", "V: spatial freq", "A: spectral env", "A: Harmonicity"
    )
  label_info_by_var$block_condition$fullnames <- c(
    "Auditory: Pitch", "Visual: Size", "Auditory: Gain", "Visual: Spatial Frequency",
    "Auditory: Spectral Envelope", "Auditory: Harmonicity"
    )

  label_info_by_var$nstim_fac = list()
  label_info_by_var$nstim_fac$original_names <- c(1, 2, 4)
  label_info_by_var$nstim_fac$abbrevs <- c(
    "List 1", "List 2", "List 4"
    )
  label_info_by_var$nstim_fac$fullnames <- c(
    "List length 1", "List length 2", "List length 4"
    )

  label_info_by_var$is_old_fac = list()
  label_info_by_var$is_old_fac$original_names <- c(0, 1)
  label_info_by_var$is_old_fac$abbrevs <- c(
    "Lure", "Target"
    )
  label_info_by_var$is_old_fac$fullnames <- c(
    "Lure trial", "Target trial"
    )

  label_info_by_var
}

get_label_info_by_var <- function(var) {
  label_info <- get_label_info_all_variables()
  label_info[[var]]
}

get_manual_legend_for_var <- function(var, scale_type="fill", display_name) {
  custom_palette <- get_custom_palette()
  lbl_info <- get_label_info_by_var(var)
  if (scale_type=="fill") {
    scl <- ggplot2::scale_fill_manual(
      values = custom_palette$wm[[var]],
      name = display_name,
      breaks = lbl_info$abbrevs,
      labels = lbl_info$fullnames
      )
  }
  scl
}

recode_factor <- function(var, old_names, new_names, reorder_levels = TRUE) {
  # Consumes a factor variable and vectors mapping the old names of factor
  # levels to new names.
  # Produces a factor variable using the new level names.

  # Initialize the new variable
  revised_var <- var
  # Assuming old_names and new_names are the same length, replace level names
  for (i in seq_along(old_names)) {
    levels(revised_var)[levels(revised_var)==old_names[i]] <- new_names[i]
  }

  # Use the order of names in new_names to order the labels in the new variable
  if (reorder_levels) {
    revised_var <- factor(revised_var, levels = new_names)
  }

  revised_var
}

recode_stim_condition <- function(var) {
  recode_factor(var,
    old_names=c(
      "a1freq", "a1freq_iden", "a1freq_multi_sync", "a1freq_multi_async",
      "v1size", "v1size_sqr_stndrd", "v1size_sqr_rot", "v1size_iden", "v1size_multi_sync", "v1size_multi_async",
      "a2gain", "a2gain_multi_sync", "a2gain_multi_async",
      "v2freq", "v2freq_multi_sync", "v2freq_multi_async"
      ),
    new_names=c(
      "Pitch, with transpose",
      "Pitch, no transpose",
      "Pitch, multimodal sync",
      "Pitch, multimodal async",
      "Size, with transpose",
      "Size, Exp 1",
      "Size, Exp 2",
      "Size, no transpose",
      "Size, multimodal sync",
      "Size, multimodal async",
      "Auditory gain",
      "Auditory gain, multimodal sync",
      "Auditory gain, multimodal async",
      "Spatial frequency",
      "Spatial frequency, multimodal sync",
      "Spatial frequency, multimodal async"
      )
    )
}

cleanup_factor_labels <- function(the_dat, the.x, preds=NULL, first_col_to_ignore="PrbtFalseAlarm") {
  # the_dat contains an enhanced summary data frame based on an lsm contrast.
  # (Previous statement may not be true...)

  # If predictors were not supplied, infer them from the data
  if (is.null(preds)) {
    # The predictors are the first columns in the data frame; all columns from
    # PrbtFalseAlarm to the end are summary data.

    # Determine the position of the first_col_to_ignore column
    prbt_col_id <- which(str_detect(colnames(the_dat), first_col_to_ignore))
    # Take a subset that contains only the predictor factors
    pred_df <- the_dat[, -c(prbt_col_id:ncol(the_dat)), drop = FALSE]

    preds <- colnames(pred_df)
  }

  # make sure all predictors are factors
  the_dat[, preds] <- lapply(the_dat[, preds, drop = FALSE], as.factor)

  has_B <- any(str_detect(preds, "block_condition"))
  has_S <- any(str_detect(preds, "stream_type"))
  has_N <- any(str_detect(preds, "nstim_fac"))
  has_M <- any(str_detect(preds, "probe_is_mono_fac"))
  has_L <- any(str_detect(preds, "target_last_pos_fac"))
  has_R <- any(str_detect(preds, "recency_fac"))
  has_C <- any(str_detect(preds, "ndirection_changes_fac"))
  has_F <- any(str_detect(preds, "block_first_stim"))

  custom_palette <- get_custom_palette()

  # container list for ggplot legends
  legend_scales <- list()
  legend_scales$block_condition = list()
  legend_scales$stream_type = list()
  legend_scales$nstim_fac = list()
  legend_scales$probe_is_mono_fac = list()
  legend_scales$target_last_pos_fac = list()
  legend_scales$recency_fac = list()
  legend_scales$ndirection_changes_fac = list()
  legend_scales$block_first_stim = list()

  xlabel <- the.x

  # Labels may differ for acuity task (e.g. block_condition will use different values) so
  # we need to determine if this is the acuity task by looking for a crossmodal condition
  is_acuity_task <- FALSE
  if (has_B | has_F) {
    # Looking for a condition that contains "crossmodal", or just using the presence of the 
    # "F" predictor to indicate that this is an acuity task
    is_acuity_task <- any(str_detect(the_dat$block_condition, "crossmodal")) | has_F
  }

  if (has_B & has_N & !has_S & !is_acuity_task) {
    # Unimodal experiments have block_condition and nstim_fac.
    block_oldnames <- c("a1freq", "v1size", "a2gain", "v2freq", "a3spec", "a4harmo")
    block_abbrevs <- c("A: pitch", "V: size", "A: gain", "V: spatial freq", "A: spectral env", "A: Harmonicity")
    block_fullnames <- c("Auditory: Pitch", "Visual: Size", "Auditory: Gain", "Visual: Spatial Frequency", "Auditory: Spectral Envelope", "Auditory: Harmonicity")
    the_dat$block_condition <- recode_factor(
      the_dat$block_condition,
      old_names = block_oldnames,
      new_names = block_abbrevs
      )

    legend_scales$block_condition$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$wm$block_condition,
      name = "Modality",
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$block_condition$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$block_condition,
      name = "Modality",
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$nstim_fac$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$wm$nstim_fac,
      name = "List Length"
      )
    legend_scales$nstim_fac$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$nstim_fac,
      name = "List Length"
      )

    # determine label for x-axis
    if (the.x == "nstim_fac") {
      xlabel <- "List Length"
    } else {
      xlabel <- "Modality"
    }
  }

  if (has_B & has_S) {
    # Multimodal experiments have block_condition and nstim_fac.
    block_oldnames <- c("avShortAa1v1", "avShortVa1v1", "avShortAa2v2", "avShortVa2v2")
    block_abbrevs <- c("A: pitch", "V: size", "A: gain", "V: spatial freq")
    block_fullnames <- c("Auditory: Pitch", "Visual: Size", "Auditory: Gain", "Visual: Spatial Frequency")
    the_dat$block_condition <- recode_factor(
      the_dat$block_condition,
      old_names = block_oldnames,
      new_names = block_abbrevs
      )

    legend_scales$block_condition$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$wm$block_condition,
      name = "Attended Modality",
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$block_condition$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$block_condition,
      name = "Attended Modality",
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$stream_type$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$wm$stream_type,
      name = "List Length"
      )
    legend_scales$stream_type$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$stream_type,
      name = "List Length"
      )

    stream_oldnames <- c("unimodal", "sync", "async")
    stream_abbrevs <- c("Unimodal", "Sync", "Async")
    stream_fullnames <- c("Unimodal", "Synchronous (congruent)", "Asynchronous (incongruent)")
    the_dat$stream_type <- recode_factor(
      the_dat$stream_type,
      old_names = stream_oldnames,
      new_names = stream_abbrevs
      )

    # determine label for x-axis
    if (the.x == "stream_type") {
      xlabel <- "Presentation Type"
    } else {
      xlabel <- "Attended Modality"
    }
  }

  if (has_S & !has_B) {
    stream_oldnames <- c("unimodal", "sync", "async")
    stream_abbrevs <- c("Unimodal", "Sync", "Async")
    stream_fullnames <- c("Unimodal", "Synchronous (congruent)", "Asynchronous (incongruent)")
    the_dat$stream_type <- recode_factor(
      the_dat$stream_type,
      old_names = stream_oldnames,
      new_names = stream_abbrevs
      )

    legend_scales$stream_type$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$wm$stream_type,
      name = "Presentation Type"
      )
    legend_scales$stream_type$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$stream_type,
      name = "Presentation Type"
      )

    # determine label for x-axis
    if (the.x == "stream_type") {
      xlabel <- "Presentation Type"
    } 
  }


  if (has_M) {
    mono_oldnames <- c("Nonmono", "Mono")
    mono_abbrevs <- c("Multidirectional probe", "Unidirectional probe")
    the_dat$probe_is_mono_fac <- recode_factor(
      the_dat$probe_is_mono_fac,
      old_names = mono_oldnames,
      new_names = mono_abbrevs
      )

  }

  if (has_L) {
    legend_scales$target_last_pos_fac$fill <- ggplot2::scale_fill_brewer(
      palette = "Dark2",
      name = "Target Serial Position"
      )
  }

  if (has_C) {
    legend_scales$ndirection_changes_fac$fill <- ggplot2::scale_fill_brewer(
      palette = "Dark2",
      name = "Direction Change Count"
      )
  }

  if ((has_B | has_F) & is_acuity_task) {
    # Acuity experiments with block condition and number of direction changes
    block_oldnames <- c("a1freq", "v1size", "a1v1_crossmodal")
    block_abbrevs <- c("Auditory", "Visual", "Crossmodal")
    block_fullnames <- c("Auditory: Pitch", "Visual: Size", "Crossmodal: Pitch and Size")
    block_first_stim_names <- c("Auditory", "Visual", "Crossmodal A V", "Crossmodal V A")
    block_first_stim_newnames <- c("Aud-Aud", "Vis-Vis", "Aud-Vis", "Vis-Aud")
    xlabel_block_condition <- "Comparison Modality"
    xlabel_block_first_stim <- "Comparison Type"
    if (has_B) {
      the_dat$block_condition <- recode_factor(
        the_dat$block_condition,
        old_names = block_oldnames,
        new_names = block_abbrevs
        )
    }

    if (has_F) {
      # The names are fine, but not in a good order. Change so we get unimodal first.
      the_dat$block_first_stim <- recode_factor(
          the_dat$block_first_stim,
          old_names = block_first_stim_names,
          new_names = block_first_stim_newnames
        )
    }

    legend_scales$block_condition$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$acuity$block_condition,
      name = xlabel_block_condition,
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$block_condition$color <- ggplot2::scale_color_manual(
      values = custom_palette$acuity$block_condition,
      name = xlabel_block_condition,
      breaks = block_abbrevs,
      labels = block_fullnames
      )
    legend_scales$block_condition$linetype <- ggplot2::scale_linetype_discrete(
      name = xlabel_block_condition,
      breaks = block_abbrevs,
      labels = block_fullnames
      )

    legend_scales$ndirection_changes_fac$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$acuity$ndirection_changes_fac,
      name = "Direction Change Count"
      )
    legend_scales$ndirection_changes_fac$color <- ggplot2::scale_color_manual(
      values = custom_palette$acuity$ndirection_changes_fac,
      name = "Direction Change Count"
      )

    legend_scales$probe_transform$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$acuity$probe_transform,
      name = "Probe Transformation"
      )
    legend_scales$probe_transform$color <- ggplot2::scale_color_manual(
      values = custom_palette$acuity$probe_transform,
      name = "Probe Transformation"
      )

    legend_scales$block_first_stim$fill <- ggplot2::scale_fill_manual(
      values = custom_palette$acuity$block_condition,
      name = xlabel_block_first_stim
      )
    legend_scales$block_first_stim$color <- ggplot2::scale_color_manual(
      values = custom_palette$acuity$block_condition,
      name = xlabel_block_first_stim
      )
    legend_scales$block_first_stim$linetype <- ggplot2::scale_linetype_discrete(
      name = xlabel_block_first_stim
      )

    # determine label for x-axis
    if (the.x == "ndirection_changes_fac") {
      xlabel <- "Direction Change Count"
    } else if (the.x == "block_first_stim") {
      xlabel <- xlabel_block_first_stim
    } else if (the.x == "block_condition") {
      xlabel <- xlabel_block_condition
    } else {
      xlabel <- "Comparison Modality"
    }
  }

  # recency effect plots
  if (has_R) {
    # Force the order to "other" followed by "last".
    # This is a HACK to ensure that jittered scatterplots are aligned with their
    # corresponding error bars.
    the_dat$recency_fac <- recode_factor(
      the_dat$recency_fac,
      old_names = c("Other pos", "Last pos"),
      new_names = c("Other pos", "Last pos")
    )
    legend_scales$recency_fac$color <- ggplot2::scale_color_manual(
      values = custom_palette$wm$recency_fac,
      name = "Target position"
      )
  }

  list(data=the_dat, legends=legend_scales, xlabel=xlabel)
}

adjust_rate_for_sdt <- function(rate, eps=get_eps()) {
  # Avoid values of 0 or 1 for rates associated with SDT values. 
  # This prevents getting +-Inf for an SDT parameter.
  max(min(rate, 1-eps), eps)
}

grp_summarize_incl_signal_var <- function(grouped_dat, signal, eps = get_eps()) {
  # prepare summary data frame for data that includes signal variable
  dat_smry <- grouped_dat %>%
    dplyr::mutate(signal_var=get(signal)) %>%
    dplyr::summarize(
      n_hit = sum(response_var==1 & signal_var==1),
      n_cr = sum(response_var==0 & signal_var==0),
      # compute hit and fa rates, but adjust to avoid 0 or 1
      hit_rate = adjust_rate_for_sdt(sum(response_var==1 & signal_var==1) / (sum(signal_var==1) + eps)),
      fa_rate = adjust_rate_for_sdt(sum(response_var==1 & signal_var==0) / (sum(signal_var==0) + eps)),
      n_trials = n(),
      accuracy = (n_hit+n_cr) / n_trials,
      sdt_d = qnorm(hit_rate) - qnorm(fa_rate),
      sdt_c = -0.5 * (qnorm(hit_rate) + qnorm(fa_rate))
    )
  dat_smry
}

grp_summarize_no_signal_var <- function(grouped_dat) {
  # prepare summary data frame for data without a signal variable
  dat_smry <- grouped_dat %>%
    dplyr::summarize(
      n_hit = sum(response_var==1),
      n_cr = 0,
      hit_rate = n_hit / n(),
      fa_rate = 0,
      n_trials = n(),
      accuracy = hit_rate,
      sdt_d = NA,
      sdt_c = NA
    )
  dat_smry
}

grp_data_summary <- function(raw_data, preds, resp_var, sig_var=NULL) {
  # The provided predictors define the groups to use for summary. 
  # Computes accuracy and SDT measures. Returns a data frame suitable
  # for some of the barplot and errorbar plot functions.
  grouping_cols <- c("participant_id", preds)
  report_cols <- c("n_hit", "n_cr", "hit_rate", "fa_rate", "n_trials", "accuracy", "sdt_d", "sdt_c")
  output_cols <- c(grouping_cols, report_cols)
  grouped_dat <- raw_data %>%
    dplyr::group_by_(.dots=grouping_cols) %>%
    dplyr::mutate(response_var=get(resp_var)
      )
  # Produce different summaries depending on whether there is a signal variable.
  # If there is a sig var, compute SDT parameter. If not, compute accuracy only.
  if (is.null(sig_var)) {
    dat_smry <- grp_summarize_no_signal_var(grouped_dat)
  } else {
    dat_smry <- grp_summarize_incl_signal_var(grouped_dat, signal=sig_var)
  }

  dat_smry %>% dplyr::select_(.dots=output_cols) %>% as.data.frame
}

merge_analyses <- function(ana_list, ids) {
  # Consumes a list of analysis data frames and a list of identifiers to match
  # those analysis dfs. 
  # Produces a single data frame that contains the analyses bound to a single 
  # data frame, with an "identifer" column with the supplied ID values.

  # This is primarily useful for generating a plot that combines results from
  # multiple experiments so they can be shown together on the same plot, or
  # using facets to put them in the same figure.

  # Currently does not do any checking that it is sensible to combine the data frames,
  # and assumes ana_list and ids are both lists and have the number of entries.

  # Attach the identifier to each analysis
  rev_analysis_list <- lapply(seq_along(ana_list), function (i, ana, id) {
    ana[[i]]$identifier <- id[[i]]
    ana[[i]]
  }, ana=ana_list, id=ids)
  # rev_analysis_list <- ana_list
  # for (i in seq_along(ana_list)) {
  #   rev_analysis_list[[i]]$identifier <- ids[[i]]
  # }
  # Combine analyses into a single df
  combo_df <- dplyr::bind_rows(rev_analysis_list)
  combo_df$identifier <- factor(combo_df$identifier)
  combo_df
}

modfunc_result_barplot <- function(data, the.x, the.y, the.fill, the.facet.var=NULL,
  the.min.y=0, the.max.y, legends, xlabel, fill_adjust = NULL,
  errbar_lower, errbar_upper, the.title=NULL) {
  plt <- ggplot(
    data=data,
    aes_string(x=the.x, y=the.y, fill=the.fill)
    ) +
    geom_bar(stat="identity", position=position_dodge(), colour="black") +
    ylim(the.min.y, the.max.y) +
    theme_jr +
    legends[[the.fill]]$fill +
    xlab(xlabel) +
    geom_errorbar(aes(ymin=errbar_lower, ymax=errbar_upper), position=position_dodge(), size=get_errorbar_size())
  # add facet wrapping if appropriate
  plt <- add_facet_wrap_to_plot(plt, the.facet.var)
  if (!is.null(fill_adjust)) {
    # This is an optional adjustment that can be supplied to change the palette
    # by luminance, for example. Can be useful when plotting different parameters
    # for the same group of variables (to distinguish the related parameters from
    # each other). Not useful when using a custom or brewer palette.
    plt <- plt + fill_adjust
  }
  if (!is.null(the.title)) {
    plt <- plt + ggtitle(the.title)
  }

  plt
}

modfunc_result_scatter_with_model <- function(model_smry, trial_smry, the.x, the.y, the.fill,
  the.facet.var=NULL, the.min.y=0, the.max.y, legends, xlabel,
  errbar_lower, errbar_upper, the.title=NULL) {
  # bar width is different on plots where 2 factors are plotted compared to those
  # with a single factor
  ebar_width <- 0.6
  if (the.x == the.fill) {
    ebar_width <- 0.4
  }
  plt <- ggplot(
    data=trial_smry,
    aes_string(x=the.x, y=the.y)
    ) +
  geom_point(aes_string(x=the.x, y=the.y, color=the.fill), size=0.7,
               position = position_jitterdodge(dodge.width=0.8, jitter.width=0.6)) +
    ylim(the.min.y, the.max.y) +
    theme_jr +
    legends[[the.fill]]$color +
    xlab(xlabel) +
    # errorbar based on the model estimates
    geom_errorbar(
      data=model_smry,
      aes(ymin=errbar_lower, ymax=errbar_upper, color=get(the.fill)),
      # aes_(x=as.name(the.x), ymin=~errbar_lower, ymax=~errbar_upper,
      #      color=as.name(the.fill), group=as.name(the.fill)),
      size=get_errorbar_size(),
      width=ebar_width,
      position = position_dodge(width=0.8)
      )
  # add facet wrapping if appropriate
  plt <- add_facet_wrap_to_plot(plt, the.facet.var)
  if (!is.null(the.title)) {
    plt <- plt + ggtitle(the.title)
  }

  plt
}

add_facet_wrap_to_plot <- function(plt, facet_var, ncol=NULL) {
  if (!is.null(facet_var)) {
    if (length(facet_var) == 1) {
      plt <- plt + facet_wrap(as.formula(paste("~", facet_var)))
    } else {
      # assume first element in facet_var goes on left side of the tilde
      first_var = facet_var[1]
      remaining_var = facet_var[2:length(facet_var)]
      if (is.null(ncol)) {
        plt <- plt + facet_wrap(as.formula(paste(first_var, "~", str_c(remaining_var, sep="+"))))
      } else {
        plt <- plt + facet_wrap(as.formula(paste(first_var, "~", str_c(remaining_var, sep="+"))), ncol=ncol)
      }
    }
  }
  plt
}

modfunc_plot_density <- function(trial_smry, the.x, the.y, the.fill, the.facet.var=NULL,
  the.min.y, the.max.y, legends, xlabel, the.title=NULL) {
  # TODO: remove the.x as a param if it is unused
  plt <- ggplot(
    data=trial_smry,
    aes_string(x=the.y)
    ) +
    geom_density(aes_string(color=the.fill, linetype=the.fill), size=get_lineplot_size()) +
    xlim(0.25, 1) +
    ylim(the.min.y, the.max.y) +
    legends[[the.fill]]$color +
    legends[[the.fill]]$linetype +
    xlab(xlabel) +
    theme_jr +
    ylab("Accuracy density")
  plt <- add_facet_wrap_to_plot(plt, the.facet.var)
  if (!is.null(the.title)) {
    plt <- plt + ggtitle(the.title)
  }
  plt
}

modfunc_plot_density_with_chance <- function(trial_smry, the.x, the.y, the.fill, the.facet.var=NULL,
  the.min.y, the.max.y, legends, xlabel, the.title=NULL) {
  # TODO: add option to show chance level distribution
  ntrials <- 64
  probsuccess <- 0.5
  xbinom <- seq(0, ntrials, length=(ntrials+1))
  ybinom <- dbinom(xbinom, ntrials, probsuccess)
  plt <- ggplot(
    data=trial_smry,
    aes_string(x=the.y)
    ) +
    geom_density(aes_string(color=the.fill, linetype=the.fill), size=get_lineplot_size()) +
    xlim(0.25, 1) +
    ylim(the.min.y, the.max.y) +
    legends[[the.fill]]$color +
    legends[[the.fill]]$linetype +
    xlab(xlabel) +
    geom_line(aes(xbinom, y=ybinom), linetype="F1") +
    theme_jr +
    ylab("Accuracy density")
  plt <- add_facet_wrap_to_plot(plt, the.facet.var)
  if (!is.null(the.title)) {
    plt <- plt + ggtitle(the.title)
  }
  plt
}

extract_facets_from_spec <- function(entry) {
  # There could be multiple facet variables, so split the string to get a
  # vector of variable names
  facet_vars <- unlist(sapply(entry$facet, function(x) {
    strsplit(x, split = "|", fixed = TRUE)
  }))
  if (length(facet_vars) < 1) {
    # No facet vars were found, so set to NULL so they are skipped in the plot
    facet_vars <- NULL
  }
  facet_vars
}

extract_maxy_from_spec <- function(entry) {
  the.max.y <- NULL
  the.min.y <- NULL
  if (!is.null(entry$max_y)) {
    # if there is a pipe character in max y split this into min y and max y
    if (str_detect(entry$max_y, "\\|")) {
      y_parts <- str_split_fixed(entry$max_y, "\\|", 2)
      # first part is min y, second part is max y
      the.min.y <- as.numeric(y_parts[[1]])
      the.max.y <- as.numeric(y_parts[[2]])
    } else {
      # if only one value, it is max y
      the.max.y <- as.numeric(entry$max_y)
    }
  }
  list(the.max.y=the.max.y, the.min.y=the.min.y)
}

group_plot_sdt_model <- function(the.y="fit", y.axis.label=NULL, trial_y=NULL, plot_type="bar") {
  function(model_smry, trial_smry = NULL, preds = NULL, plot_spec = NULL, the.x, the.fill = the.x, the.facet.var = NULL, 
    the.max.y = NULL, fill_adjust = NULL, the.title=NULL) {

    # If the plot spec is provided, use this to set several of the plot parameters. Note that
    # the plot spec could override other parameters such as "the.x"
    the.min.y <- NULL  # TODO: should move this to function argument or change existing 
                       # max y argument to be more general
    if (!is.null(plot_spec)) {
      the.facet.var <- extract_facets_from_spec(plot_spec)
      # extract ymin and ymax if available
      y_parts <- extract_maxy_from_spec(plot_spec)
      the.max.y <- y_parts$the.max.y
      the.min.y <- y_parts$the.min.y
      the.x <- plot_spec$x
      the.fill <- plot_spec$fill
      the.title <- plot_spec$title
    }
  
    if (is.null(preds)) {
      # this should be an error
    }

    # Change the names of the dependent measure
    names(model_smry)[names(model_smry)==the.y] <- y.axis.label
    if (!is.null(trial_y) & !is.null(trial_smry)) {
      names(trial_smry)[names(trial_smry)==trial_y] <- y.axis.label
    }
    tmp <- cleanup_factor_labels(model_smry, the.x, preds=preds)

    # extract the revised data and the legends from the list returned by the cleanup function
    model_smry <- tmp$data
    legends <- tmp$legends
    xlabel <- tmp$xlabel

    if (!is.null(trial_smry)) {
      tmp <- cleanup_factor_labels(trial_smry, the.x, preds=preds)
      trial_smry <- tmp$data
    }

    # choose sensible defaults for min and max y if they were not specified
    if (is.null(the.max.y)) {
      the.max.y <- ceiling(max(model_smry$upr)) # FIX THIS
    }

    if (is.null(the.min.y)) {
      the.min.y <- floor(min(model_smry$lwr))
      the.min.y <- min(the.min.y, 0)
    }

    model_smry$errbar_lower <- model_smry$lwr
    model_smry$errbar_upper <- model_smry$upr

    # Plot it
    if (plot_type=="bar") {
      plt <- modfunc_result_barplot(data=model_smry, the.x=the.x, the.y=y.axis.label, the.fill=the.fill,
        the.facet.var=the.facet.var, the.min.y=the.min.y, the.max.y=the.max.y, legends=legends,
        xlabel=xlabel, fill_adjust=fill_adjust, errbar_lower=errbar_lower, errbar_upper=errbar_upper, the.title=the.title)
    } else if (plot_type=="scatter") {
      plt <- modfunc_result_scatter_with_model(
        model_smry=model_smry, trial_smry=trial_smry, the.x=the.x, the.y=y.axis.label, the.fill=the.fill,
        the.facet.var=the.facet.var, the.min.y=the.min.y, the.max.y=the.max.y, legends=legends,
        xlabel=xlabel, errbar_lower=errbar_lower, errbar_upper=errbar_upper, the.title=the.title
        )
    } else if (plot_type=="density") {
      plt <- modfunc_plot_density(
        trial_smry=trial_smry, the.x=the.x, the.y=y.axis.label, the.fill=the.fill,
        the.facet.var=the.facet.var, the.min.y=the.min.y, the.max.y=the.max.y, legends=legends,
        xlabel=xlabel, the.title=the.title
        )
    } else {
      plt <- NULL
    }

    plt
  }
}

group_plot_sdt_d <- group_plot_sdt_model(the.y="fit", y.axis.label="Sensitivity")
group_plot_sdt_c <- group_plot_sdt_model(the.y="fit", y.axis.label="Bias")
group_plot_sdt_acc <- group_plot_sdt_model(the.y="fit", y.axis.label="Accuracy", plot_type="scatter", trial_y="accuracy")
group_plot_density <- group_plot_sdt_model(the.y="Accuracy", y.axis.label="Accuracy", plot_type="density", trial_y="accuracy")

extract_predictor_names <- function(data, first_col_to_ignore=NULL) {
  non_pred_idx <- which(stringr::str_detect(names(data), first_col_to_ignore))
  if (non_pred_idx > 1) {
    preds <- names(data)[1:(non_pred_idx-1)]
  }
  preds
}

# NOTE:
# plot_sdt_parameters_by_block is used in the "plots_for_2017committee.R" script but
# should not be used for future plots. Use group_plot_sdt_model() instead.
plot_sdt_parameters_by_block <- function(sdt_dat, preds = NULL, the.x, the.y="sdt_d_computed_with_se",
  the.fill = the.x, the.facet.var = NULL, the.max.y = NULL, fill_adjust = NULL, the.title=NULL) {
  # sdt_dat is a data frame with estimates of SDT parameters by experiment
  # condtion. 

  # If the predictor variables are not provided, assume they are the first
  # few columns in the data frame. The first non-predictor column is 
  # "PrbtFalseAlarm".
  if (is.null(preds)) {
    non_pred_idx <- which(stringr::str_detect(names(sdt_dat), "PrbtFalseAlarm"))
    if (non_pred_idx > 1) {
      preds <- names(sdt_dat)[1:(non_pred_idx-1)]
    }
  }

  # Change the names of the SDT columns
  names(sdt_dat)[names(sdt_dat)=="sdt_d_computed_with_se"] <- "Sensitivity"
  names(sdt_dat)[names(sdt_dat)=="sdt_c_computed_with_se"] <- "Bias"

  if (the.y == "sdt_d_computed_with_se") {
    the.y = "Sensitivity"
    the.y.se = "sdt_d_se"
  } else if (str_detect(the.y, "^sdt_c")) {
    # Start with sdt_c, so must be plotting bias
    the.y = "Bias"
    the.y.se = "sdt_c_se"
  }

  tmp <- cleanup_factor_labels(sdt_dat, the.x)
  # extract the revised data and the legends from the list returned by the cleanup function
  sdt_dat <- tmp$data
  legends <- tmp$legends
  xlabel <- tmp$xlabel

  if (is.null(the.max.y)) {
    the.max.y <- ceiling(max(sdt_dat[, the.y] + sdt_dat[, the.y.se])) # FIX THIS
  }
  the.min.y <- floor(min(sdt_dat[, the.y] + sdt_dat[, the.y.se]))
  the.min.y <- min(the.min.y, 0)
  
  # compute error bar bounds
  sdt_dat$errbar_lower <- sdt_dat[, the.y] - sdt_dat[, the.y.se]
  sdt_dat$errbar_upper <- sdt_dat[, the.y] + sdt_dat[, the.y.se]

  # Plot it
  plt <- modfunc_result_barplot(data=sdt_dat, the.x=the.x, the.y=the.y, the.fill=the.fill,
    the.facet.var=the.facet.var, the.min.y=the.min.y, the.max.y=the.max.y, legends=legends,
    xlabel=xlabel, fill_adjust=fill_adjust, errbar_lower=errbar_lower, errbar_upper=errbar_upper, the.title=the.title)

  plt
}

# Define aliases for convenience
plot_dprime_by_block <- purrr::partial(
  plot_sdt_parameters_by_block,
  the.y="sdt_d_computed_with_se"
  )
plot_bias_by_block <- purrr::partial(
  plot_sdt_parameters_by_block,
  the.y="sdt_c_computed_with_se"
  )
plot_criterion_by_block <- purrr::partial(
  plot_sdt_parameters_by_block,
  the.y="sdt_c_computed_with_se"
  )

plot_accuracy_by_block <- function(acc_dat, preds = NULL, the.x, the.y="Accuracy",
  the.fill = the.x, the.facet.var = NULL, the.max.y = NULL, errbar_col_stub="prob_one_se",
  the.title = NULL) {
  # acc_dat is a data frame with accuracy estimates including upper and lower bounds

  # If the predictor variables are not provided, assume they are the first
  # few columns in the data frame. The first non-predictor column is 
  # "PrbtHit".
  first_nonpred_col_name <- colnames(acc_dat)[colnames(acc_dat) %in% c("PrbtHit", "PrbtAccuracy")]
  if (is.null(preds)) {
    preds <- extract_predictor_names(acc_dat, first_nonpred_col_name)
    # non_pred_idx <- which(stringr::str_detect(names(acc_dat), "PrbtHit"))
    # if (non_pred_idx > 1) {
    #   preds <- names(acc_dat)[1:(non_pred_idx-1)]
    # }
  }

  tmp <- cleanup_factor_labels(acc_dat, the.x, first_col_to_ignore=first_nonpred_col_name)
  # extract the revised data and the legends from the list returned by the cleanup function
  acc_dat <- tmp$data
  legends <- tmp$legends
  xlabel <- tmp$xlabel

  if (is.null(the.max.y)) {
    the.max.y <- ceiling(max(acc_dat[, str_c(errbar_col_stub, "_hi")]))
  }
  
  # The error bar boundaries are determined by taking the provided stub and appending
  # "lo" and "hi" to get lower and upper bounds
  acc_dat$errbar_lower <- acc_dat[, str_c(errbar_col_stub, "_lo")]
  acc_dat$errbar_upper <- acc_dat[, str_c(errbar_col_stub, "_hi")]

  # Plot it
  plt <- modfunc_result_barplot(data=acc_dat, the.x=the.x, the.y=the.y, the.fill=the.fill,
    the.facet.var=the.facet.var, the.max.y=the.max.y, legends=legends,
    xlabel=xlabel, errbar_lower=errbar_lower, errbar_upper=errbar_upper, the.title=the.title)

  plt
}

###############################
# diagnostic plots
###############################

plot_similarity_glm_smoothed <- function(dat, 
  the.x = "trial_sumsim_trltype_z", the.y = "respondedsame", the.fill = "block_condition", the.facet.var = NULL,
  legends = NULL, xlabel = NULL, ylabel = "P(old)", in_color=TRUE, ncol=NULL, the.title = NULL
  )
{
  if (the.fill == "stim_condition") {
    dat$stim_condition <- recode_stim_condition(dat$stim_condition)
  }
  plt <- ggplot(dat)
  if (in_color) {
    plt <- plt + aes_string(x=the.x, y=the.y, color=the.fill)
    plt <- plt + stat_smooth(method="glm", method.args = list(family = "binomial"), formula=y~x, alpha=0.2, size=get_lineplot_size(), aes_string(fill=the.fill))
  } else {
    plt <- plt + aes_string(x=the.x, y=the.y, linetype=the.fill)
    plt <- plt + stat_smooth(method="glm", method.args = list(family = "binomial"), formula=y~x, alpha=0.2, size=get_lineplot_size(), color="black")
  }
  plt <- plt + xlab(xlabel) + ylab(ylabel) +
    ylim(0., 1.0)
  plt <- plt + theme_jr
  # add facet wrapping if appropriate
  plt <- add_facet_wrap_to_plot(plt, the.facet.var, ncol=ncol)

  if (!is.null(the.title)) {
    plt <- plt + ggtitle(the.title)
  }

  plt

}
