## stan versions
# Removed for simplicity
mBase_stan <- NULL

mD_stan_zCovar_noInteraction <- NULL

mD_stan_zCovar <- NULL

# Use the basic SDT model with z-scored difficulty covariate as the base model.
mD_stan <- mD_stan_zCovar

# Basic SDT model with no predictors
mBase_mle <- function(dat) {
  fit <- glmer(
    respondedsame ~
      is_old +
      (is_old | participant_id),
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
}

# Simple SDT model with z-scored difficulty covariate, useful as a base model.
# The interaction of the signal indicator and the difficulty covariate is included
# in this model, although it is dropped automatically in experiments that only have
# 1 stim per trial. Those experiments have no variation when the signal var is 1
# (i.e. on target trials) so the interaction becomes collinear with the main effect.
mD_mle <- function(dat) {
  fit <- glmer(
    respondedsame ~
      is_old +
      trial_sumsim_trltype_z + is_old:trial_sumsim_trltype_z + # covariate
      (is_old | participant_id),  # participant_id gets random criterion, is_old gets random discriminability
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
}

# SDT model without the interaction of the signal variable and covariate
mD_mle_noInteraction <- function(dat) {
  fit <- glmer(
    respondedsame ~
      is_old +
      trial_sumsim_trltype_z + # covariate
      (is_old | participant_id),  # participant_id gets random criterion, is_old gets random discriminability
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
}

# Models using an alternate signal variable that does not include a difficulty covariate. If any study item
# has the same initial direction and number of changes as the probe, that is considered a signal trial.
nchange_mBase_mle <- function(dat) {
  fit <- glmer(
    respondedsame ~
      is_old_nchange_and_direction +
      (is_old_nchange_and_direction | participant_id),  # participant_id gets random criterion, is_old gets random discriminability
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
}


########
# SDT models for 2IFC acuity task
########

acuity_mBase_stan <- NULL

acuity_mBase_mle <- function(dat) {
  fit <- glmer(
    resp_1st_interval ~ 
      same_in_1st_pair +
      (same_in_1st_pair | participant_id),
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
}

########
# Models with no signal indicator, useful for predictors such as recency effect
########
nosig_mD_stan <- NULL

nosig_mD_mle <- function(dat) {
  fit <- glmer(
    respondedsame ~
      trial_sumsim_trltype_z +
      (1 | participant_id),
    data=dat,
    family=binomial(link="probit"),
    control = glmerControl(
      optimizer = "nloptwrap",
      calc.derivs = FALSE,
      optCtrl = list(maxfun = 100000)
      )
  )
  fit
}
