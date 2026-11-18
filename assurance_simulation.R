#-------------------------------------------------------------------------------
# SIMULATION BASED SAMPLE SIZE AND POWER ANALYSIS
#-------------------------------------------------------------------------------
set.seed(13181913)
options(brms.backend = "cmdstanr")

#Dependencies
library(tidyverse)
library(tidybayes)
library(readr)
library(brms)
library(bayestestR)
library(posterior)
library(bayesplot)
library(emmeans)
library(cmdstanr)
source("prereg_code/util.R")

#-------------------------------------------------------------------------------
# CREATING BASELINE WITH PREVIOUS RESEARCH RESULTS
#-------------------------------------------------------------------------------
dat_old = readRDS("dat_long.RDS")
mean(log(dat_old$abs_error))
sd(log(dat_old$abs_error))
priors_old = c(
  prior(normal(0,1), class = "b"),
  prior(normal(1.43, 2), class = "Intercept"), #mean = mean of log abs error from data sd = 2* sd of log abs error
  prior(exponential(1), class = "sd"),
  prior(exponential(1), class = "shape")
)


fit_old = readRDS("material/fit_old.RDS")
if(is.null(fit_old)){
  fit_old = brm(
    abs_error ~
      act_pass * abs_dist * direction * speed +
      act_pass * start * abs_dist * speed +
      time +
      (1 | ID),
    data = dat_old,
    family = Gamma(link = "log"),
    prior = priors_old,
    chains = 4,
    iter = 7500,
    warmup = 3000,
    cores = 4
  )
  saveRDS(fit_old, "material/fit_old.RDS")
}

#posterior draws
posterior_old = posterior_summary(fit_old)
draws = as_draws_df(fit_old)

#ellicit shape
mcmc_dens(fit_old, pars = "shape")
mu_shape = mean(draws$shape)
var_shape = var(draws$shape)
gamma_shape_pars = gamma_moment_match(mu_shape, var_shape, 1) #no increase in uncertainty
q025 = qgamma(.025,  gamma_shape_pars[1], gamma_shape_pars[2])
q975 = qgamma(.975, gamma_shape_pars[1], gamma_shape_pars[2])
shape_mean = gamma_shape_pars[1] /  gamma_shape_pars[2]
hist(rgamma(1000, gamma_shape_pars[1], gamma_shape_pars[2]))

#ellicit sd_ID__Intercept
mcmc_dens(fit_old, pars = "sd_ID__Intercept")
mu_rIntercept = mean(draws$sd_ID__Intercept)
sd_rIntercept = sqrt(var(draws$sd_ID__Intercept)) #no increase in uncertainty
hist(rnorm(1000, mu_rIntercept, sd_rIntercept))

#ellicit b_Intercept
mcmc_dens(fit_old, pars = "b_Intercept")
mu_Intercept = mean(draws$b_Intercept)
sd_Intercept = sd(draws$b_Intercept)
hist(rnorm(1000, mu_Intercept, sd_Intercept))

#ellicit group effect
draws_actpass = emmeans(fit_old, ~ act_pass, type = "response", epred = TRUE) %>%
  gather_emmeans_draws()
rr_actpass = log(draws_actpass[draws_actpass[,1] == "PASS", 5] / draws_actpass[draws_actpass[,1] == "ACT", 5])
mu_b_act_pass = mean(rr_actpass$.value)
sd_b_act_pass = sqrt(var(rr_actpass$.value) * 2)


derived_params = list(
  b_Intercept = c(mu_Intercept, sd_Intercept),
  sd_ID__Intercept = c(mu_rIntercept, sd_rIntercept),
  shape = c(gamma_shape_pars[1], gamma_shape_pars[2]),
  b_actpassPASS = c(mu_b_act_pass, sd_b_act_pass)
)

#-------------------------------------------------------------------------------
# GENERATING PRIOR
#-------------------------------------------------------------------------------
#...............................................................................
#confirmatory directional hypotheses: beta_main + 1/2 * beta_mainXtiming " 1/2 beta_mainXmovement_type > 0
# small to moderate effect
# mean according to d = 0.35, sd according to 95CI[d = 0.2, d = 0.5]
mu_confirmatory = approx_beta_cohend(0.35, mu_shape)
ci_confirmatory = approx_beta_cohend(c(0.2, 0.5), mu_shape)
sigma_confirmatory = approx_sigma(mu_confirmatory, ci_confirmatory)

#validity check
theta = rnorm(
  1e6,
  mu_confirmatory,
  sigma_confirmatory
)
quantile(theta, c(0.025, 0.5, 0.975))
hist(theta)

d_theta = approx_cohend_beta(theta, mu_shape)
quantile(d_theta, c(0.025, 0.5, 0.975))
#...............................................................................
#confirmatory null-effect hypotheses: beta_main + 1/2 * beta_mainXtiming " 1/2 beta_mainXmovement_type = 0
# ROPE according to +/- d = 0.1 
theta_null = 0

rope_d = c(-0.1, 0.1)
rope_beta = approx_beta_cohend(rope_d, mu_shape)
#...............................................................................
#exploratory hypotheses (mainly interactions with timing and movement_type)
# mean = 0, no directionality assumption

interaction_ci = approx_beta_cohend(c(-0.15,0.15), mu_shape)
sigma_interaction = approx_sigma(0,interaction_ci)
#...............................................................................
#main effect priors derived from confirmatory effect and interactions
var_main_confirmatory = sigma_confirmatory^2 - 0.25 * sigma_interaction^2 - 0.25 * sigma_interaction^2
sigma_main_confirmatory = sqrt(var_main_confirmatory)

beta_params_confirmatory = c(mu_confirmatory, sigma_main_confirmatory)
#verification that the prior choice results in approximate cohen d-s defined
b_main = rnorm(1e6, mu_confirmatory, sigma_main_confirmatory)
b_timing = b_movement = rnorm(1e6, 0, sigma_interaction)
theta = b_main + 0.5 * b_timing + 0.5 * b_movement
d_theta = approx_cohend_beta(theta, mu_shape)
quantile(d_theta, c(.025, .5, .975))

#...............................................................................
#main effect of visual and verbal derived from the assumption that the overall marginal effect is 0
# FIXED TO 0

#...............................................................................
#ellicit random slopes 
#population mean based on beta_params_confirmatory and from sigma_main_null, sigma_interaction null etc
# at 95th percentile we set RR=1.10 (10 percent variation in error around the fixed effect)
upper_rr = 1.10
upper_beta = log(upper_rr)
sd_slope_scale = upper_beta / qnorm(0.975)

#verify the implied distribution
median_slope_sd = sd_slope_scale * qnorm(0.75)
exp(median_slope_sd)
q95_slope_sd = sd_slope_scale * qnorm(0.975)
exp(q95_slope_sd)

#...............................................................................
#ellicit random intercept for trial error 
upper_rr_trial = 1.10
upper_beta_trial = log(upper_rr_trial)
sd_trial_scale = upper_beta_trial / qnorm(0.975)
#...............................................................................
# duration effect:
# assumed to be positive (with only unlikely <0 values), CI95 based on d = 0 - 0.2
beta_duration_ci = approx_beta_cohend(c(0,0.2), mu_shape)
beta_duration_mean = approx_beta_cohend(0.1, mu_shape)
beta_duration_sigma = approx_sigma(beta_duration_mean, beta_duration_ci)

#...............................................................................
# timing main effect:
beta_timing_ci = approx_beta_cohend(c(-0.1, 0.1), mu_shape)
beta_timing_sigma = approx_sigma(0, beta_timing_ci)

gen_prior = c(
  prior(normal(0.299, 0.051), class = "b", coef = "interferencemotor"), #beta_params_confirmatory
  prior(normal(0.299, 0.051), class = "b", coef = "interferencespatial"), #beta_params_confirmatory
  prior(constant(0), class = "b", coef = "interferenceverbal"), #beta_params_null 
  prior(constant(0), class = "b", coef = "interferencevisual"),
  prior(normal(0.059,0.03), class = "b", coef = "movement_typepassive"), #moment mathed empirical
  prior(normal(0, 0.064), class = "b", coef = "interferencemotor:timingretention"), #sigma interaction
  prior(normal(0, 0.064), class = "b", coef = "interferencespatial:timingretention"), #sigma interaction
  prior(constant(0), class = "b", coef = "interferenceverbal:timingretention"),
  prior(constant(0), class = "b", coef = "interferencevisual:timingretention"),
  prior(normal(0, 0.064), class = "b", coef = "interferencemotor:movement_typepassive"), #sigma interaction
  prior(normal(0, 0.064), class = "b", coef = "interferencespatial:movement_typepassive"), #sigma interaction
  prior(constant(0), class = "b", coef = "interferenceverbal:movement_typepassive"), #beta null sigma
  prior(constant(0), class = "b", coef = "interferencevisual:movement_typepassive"), #beta null sigma
  prior(normal(0, 0.064), class = "b", coef = "timingretention:movement_typepassive"), #sigma interaction
  prior(normal(0.084,0.043), class = "b", coef = "duration"), #beta duration mean and sigma
  prior(normal(0, 0.042), class = "b", coef = "timingretention"), #beta timing 
  prior(normal(1.623, 0.055), class = "Intercept"), #moment matched empirical
  prior(normal(0.214, 0.017), class = "sd", coef = "Intercept", group = "participant"), #moment matched empirical
  prior(normal(0,0.049), class = "sd", group = "participant"), #sd_slope_scale
  prior(lkj(2), class = "cor", group = "participant"),
  prior(normal(0, 0.049), class = "sd", coef = "Intercept", group = "movement_trial_id"),
  prior(gamma(3870.081,2719.264), class = "shape") #moment matched empirical
)

#...............................................................................
# validate generator
trials = read_csv("material/trial_params.csv")

design_data = create_design(50,trials)
gen_fit = brm(data = design_data,
              formula = abs_error ~ interference + 
                timing + 
                movement_type +
                interference*timing + 
                interference*movement_type +
                movement_type*timing +
                duration + 
                (1 + interference | participant) +
                (1 | movement_trial_id),
              family = Gamma(link = "log"),
              prior = gen_prior,
              sample_prior = "only",
              iter = 2000,
              warmup = 500,
              chains = 4,
              cores = 4)
prior_summary(gen_fit)

mu_rep = posterior_epred(
  gen_fit,
  ndraws = 1000
)

mu_check = tibble(
  draw = 1:nrow(mu_rep),
  mean_mu = apply(mu_rep, 1, mean),
  median_mu = apply(mu_rep, 1, median),
  p95_mu = apply(mu_rep, 1, quantile, .95),
  max_mu = apply(mu_rep, 1, max)
)

summary(mu_check)
pp_generator = check_generator(gen_fit)
summary(pp_generator)

yrep = posterior_predict(gen_fit, ndraws = 1, re_formula = NULL)
genprior_data = design_data
genprior_data$abs_error = as.numeric(yrep[1,])

genprior_data %>%
  group_by(participant) %>%
  summarise(mean_abserror = mean(abs_error), .groups = "drop") %>%
  ggplot(aes(x = mean_abserror)) +
  geom_histogram(bins = 20)

dat_old %>%
  group_by(ID) %>%
  summarise(mean_abserror = mean(abs_error), .groups = "drop") %>%
  ggplot(aes(x = mean_abserror)) +
  geom_histogram(bins = 20)


#-------------------------------------------------------------------------------
# ANALYSIS PRIOR
#-------------------------------------------------------------------------------
#...............................................................................
#analysis priors for beta parameters such that the 95CI is generated from d = +/- 1
analysis_beta_ci = approx_beta_cohend(c(-1,1), mu_shape)
analysis_beta_sigma = approx_sigma(0, analysis_beta_ci)

#...............................................................................
#shape parameter from empirical but with seriously inflated variance

analysis_gammas = gamma_moment_match(mu_shape, 0.5^2, 1)
hist(rnorm(1000, mu_shape, 0.5^2))
hist(rgamma(1000, analysis_gammas[1], analysis_gammas[2]))

analysis_prior = c(
  prior(normal(0,0.510), class = "b"), #+/- 1 cohen d
  prior(normal(1.623, 0.11), class = "Intercept"), #double sd, empirical
  prior(exponential(3), class = "sd"), 
  prior(lkj(2), class = "cor"),
  prior(gamma(8.102, 5.693), class = "shape")
)

#-------------------------------------------------------------------------------
# ASSURANCE FOR NONOPTIMAL DESIGN (N = 80, 100, 110,120, 140)
#-------------------------------------------------------------------------------
criteria = matrix(c(
  c(0.95, 0.95, 0.95, 0.95, 0.95, 0.95, 0.05, 0.05, 0.95, 0.95, 0.95, 0.95, 0.05, 0.05,
    1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1)
), ncol = 14, nrow = 2, byrow = TRUE)


assurance_80 = one_assurance(40, trials, 100, 
                              gen_priors = gen_prior, fit_priors = analysis_prior, criteria = criteria,
                              warmup = 750, iter = 2750, chains = 1, cores = 6)
assurance_80 %>%
  group_by(variable) %>%
  summarise(across(c("mean", 
                     "sd", 
                     "ess_bulk", 
                     "mcse_mean", 
                     "decision"), ~ mean(.x, na.rm = TRUE)))

saveRDS(assurance_80, file = "material/assurance_80.RDS")

assurance_100 = one_assurance(50, trials, 100, 
                             gen_priors = gen_prior, fit_priors = analysis_prior, criteria = criteria,
                             warmup = 750, iter = 2750, chains = 1, cores = 6)
assurance_100 %>%
  group_by(variable) %>%
  summarise(across(c("mean", 
                     "sd", 
                     "ess_bulk", 
                     "mcse_mean", 
                     "decision"), ~ mean(.x, na.rm = TRUE)))

saveRDS(assurance_100, file = "material/assurance_100.RDS")

assurance_110= one_assurance(55, trials, 100, 
                              gen_priors = gen_prior, fit_priors = analysis_prior, criteria = criteria,
                              warmup = 750, iter = 2750, chains = 1, cores = 6)
assurance_110 %>%
  group_by(variable) %>%
  summarise(across(c("mean", 
                     "sd", 
                     "ess_bulk", 
                     "mcse_mean", 
                     "decision"), ~ mean(.x, na.rm = TRUE)))

saveRDS(assurance_110, file = "material/assurance_110.RDS")


assurance_120 = one_assurance(60, trials, 100, 
                              gen_priors = gen_prior, fit_priors = analysis_prior, criteria = criteria,
                              warmup = 750, iter = 2750, chains = 1, cores = 6)
assurance_120 %>%
  group_by(variable) %>%
  summarise(across(c("mean", 
                     "sd", 
                     "ess_bulk", 
                     "mcse_mean", 
                     "decision"), ~ mean(.x, na.rm = TRUE)))

saveRDS(assurance_120, file = "material/assurance_120.RDS")

assurance_140 = one_assurance(70, trials, 100, 
                              gen_priors = gen_prior, fit_priors = analysis_prior, criteria = criteria,
                              warmup = 750, iter = 2750, chains = 1, cores = 6)
assurance_140 %>%
  group_by(variable) %>%
  summarise(across(c("mean", 
                     "sd", 
                     "ess_bulk", 
                     "mcse_mean", 
                     "decision"), ~ mean(.x, na.rm = TRUE)))

saveRDS(assurance_140, file = "material/assurance_140.RDS")


