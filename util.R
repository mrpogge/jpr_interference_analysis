################################################################################
# FUNCTIONS OF THE PREREG CODE
################################################################################

#-------------------------------------------------------------------------------
# ASSURANCE ANALYSIS FUNCTIONS
#-------------------------------------------------------------------------------
#based on a posterior moments find gamma pars (k governs uncertainty)
gamma_moment_match = function(mu, sigma2, k){
  sigma2_k = k * sigma2
 return(c(mu^2/sigma2_k, mu/sigma2_k)) 
}

rescale_01 = function(X, minmax = NULL){
  if(!is.null(minmax)){
    return((X - minmax[1]) / (minmax[2]-minmax[1]))
  } else {
    return((X - min(X)) / (max(X) - min(X)))
  }
  
}

create_design = function(n_sample_half, trials){
  
  #---------------------------------------------------------------------------
  # CREATE PARTICIPANT STRUCTURE
  #---------------------------------------------------------------------------
  
  dat = expand.grid(
    participant = factor(1:(n_sample_half * 2)),
    timing = c("encoding", "retention"),
    interference = c("control", "verbal", "visual", "spatial", "motor"),
    repetition = 1:5
  ) %>%
    mutate(
      movement_type = ifelse(
        as.numeric(as.character(participant)) <= n_sample_half,
        "passive",
        "active"
      )
    )
  
  #---------------------------------------------------------------------------
  # CREATE ONE TRIAL-INTERFERENCE MAPPING PER PARTICIPANT
  #---------------------------------------------------------------------------
  
  trial_map = expand.grid(
    participant = factor(1:(n_sample_half * 2)),
    interference = c("control", "verbal", "visual", "spatial", "motor"),
    repetition = 1:5
  ) %>%
    group_by(participant) %>%
    mutate(
      movement_trial_id = sample(
        1:25,
        25,
        replace = FALSE
      )
    ) %>%
    ungroup()
  
  #---------------------------------------------------------------------------
  # APPLY THE SAME MAPPING TO BOTH TIMING CONDITIONS
  #---------------------------------------------------------------------------
  
  dat = dat %>%
    left_join(
      trial_map,
      by = c(
        "participant",
        "interference",
        "repetition"
      )
    )
  
  #---------------------------------------------------------------------------
  # ADD TRIAL PARAMETERS
  #---------------------------------------------------------------------------
  
  dat$abs_error = runif(nrow(dat), 0.5, 2)
  
  dat = dat %>%
    left_join(
      trials,
      by = "movement_trial_id"
    ) %>%
    mutate(
      duration = 
        calculate_interference_time(
          amplitude,
          speed
        )
    ) %>%
    mutate(
      duration =
        (duration - mean(duration)) /
        sd(duration)
    )
  
  return(dat)
}

duration_calibration = readRDS(
  "duration_calibration.rds"
)

calculate_interference_time = function(amplitude, speed){
  
  new_forward = data.frame(
    amplitude = amplitude,
    speed = speed
  )
  
  new_return = data.frame(
    amplitude = amplitude
  )
  
  forward_time = predict(
    duration_calibration$fit_forward,
    newdata = new_forward
  )
  
  return_time = predict(
    duration_calibration$fit_return,
    newdata = new_return
  )
  
  as.numeric(
    forward_time +
      4 +
      return_time + duration_calibration$duration_overhead
  )
}

approx_beta_cohend_tay = function(d, shape){
  return(log(1 + d/sqrt(shape)))
}

approx_beta_cohend = function(d, shape){
  return(2*atanh(d/(sqrt(4*shape - d^2))))
}

approx_cohend_beta = function(beta, shape){
  return(sqrt(2*shape) * (exp(beta) - 1) / sqrt(1 + exp(2 * beta)))
}

approx_sigma = function(mu, ci){
  sigma = mean(c(mu - ci[1], ci[2]-mu)) / qnorm(0.975)
}


check_generator = function(gen_fit, n_sim = 100){
  yrep = posterior_predict(gen_fit, ndraws = n_sim)
  return(
    tibble(
      draw = 1:nrow(yrep),
      mean_y = apply(yrep, 1, mean),
      median_y = apply(yrep, 1, median),
      p05_y    = apply(yrep, 1, quantile, .05),
      p95_y    = apply(yrep, 1, quantile, .95),
      max_y    = apply(yrep, 1, max)
    )
  )
}


#------------------
init_generator_template = function(n_sample_half, data = NULL, trials, priors, type, warmup = 500, iter = 2000, chains = 4, cores = 4){
  
  if(is.null(data)){
    design_data = create_design(n_sample_half, trials)
  } else {
    design_data = data
  }
  
  if(type == "gen"){
    fit = brm(data = design_data,
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
              prior = priors,
              sample_prior = "only",
              iter = iter,
              warmup = warmup,
              chains = chains,
              cores = cores)
  } else if(type == "template"){
    fit =  brm(data = design_data,
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
               prior = priors,
               sample_prior = "no",
               iter = iter,
               warmup = warmup,
               chains = chains,
               cores = cores)
  }
  
  return(list(fit = fit, design_data = design_data))
}

decide_crit = function(x, criteria, direction){
  return(ifelse(direction == 1, as.numeric(x > criteria), as.numeric(x < criteria)))
}


sample_assurance = function(design_data, gen_fit, template_fit, criteria,
                            warmup = 500, iter = 2000, chains = 4, cores = 4){
  results = numeric(6)
  #sample one set of parameters
  yrep = posterior_predict(
    gen_fit,
    ndraws = 1
  )
  design_data$abs_error = as.numeric(yrep)
  
  test_fit = update(
    template_fit,
    newdata = design_data,
    recompile = FALSE,
    backend = "cmdstanr",
    refresh = 0
  )
  posterior_draw = as_draws_df(test_fit)
  mu_shape = mean(posterior_draw$shape)
  posterior_test = mutate_variables(
    posterior_draw,
    theta_1 =
      b_interferencemotor +
      0.5 * `b_interferencemotor:timingretention` +
      0.5 * `b_interferencemotor:movement_typepassive`,
    theta_2 =
      b_interferencespatial +
      0.5 * `b_interferencespatial:timingretention` +
      0.5 * `b_interferencespatial:movement_typepassive`,
    theta_3 =
      b_interferenceverbal +
      0.5 * `b_interferenceverbal:timingretention` +
      0.5 * `b_interferenceverbal:movement_typepassive`,
    theta_4 =
      b_interferencevisual +
      0.5 * `b_interferencevisual:timingretention` +
      0.5 * `b_interferencevisual:movement_typepassive`,
    rope_upper = approx_beta_cohend(0.1, mu_shape),
    rope_lower = approx_beta_cohend(-0.1, mu_shape),
    rope_upper2 = approx_beta_cohend(0.2, mu_shape),
    rope_lower2 = approx_beta_cohend(-0.2, mu_shape),
    H1_0 = as.numeric(theta_1 > 0),
    H2_0 = as.numeric(theta_2 > 0),
    H1_d = as.numeric(theta_1 > rope_upper),
    H2_d = as.numeric(theta_2 > rope_upper),
    H3_rope = as.numeric(theta_3 > rope_lower & theta_3 < rope_upper),
    H4_rope = as.numeric(theta_4 > rope_lower & theta_4 < rope_upper),
    H3_d = as.numeric(theta_3 > rope_upper),
    H4_d = as.numeric(theta_4 > rope_upper),
    H1_d2 = as.numeric(theta_1 > rope_upper2),
    H2_d2 = as.numeric(theta_2 > rope_upper2),
    H3_rope2 = as.numeric(theta_3 > rope_lower2 & theta_3 < rope_upper2),
    H4_rope2 = as.numeric(theta_4 > rope_lower2 & theta_4 < rope_upper2),
    H3_d2 = as.numeric(theta_3 > rope_upper2),
    H4_d2 = as.numeric(theta_4 > rope_upper2),
  )
  
  vars_save = vars_save = c(
    "H1_0",
    "H2_0",
    "H1_d",
    "H2_d",
    "H3_rope",
    "H4_rope",
    "H3_d",
    "H4_d",
    "H1_d2",
    "H2_d2",
    "H3_rope2",
    "H4_rope2",
    "H3_d2",
    "H4_d2"
  )
  
  diag = posterior_test %>%
    subset_draws(vars_save) %>%
    summarise_draws(
      "mean",
      "sd",
      "ess_bulk",
      "mcse_mean"
    ) %>%
    mutate(val = as.numeric(criteria[1,]),
           dir = as.numeric(criteria[2,]),
           decision = decide_crit(mean, val, dir)) %>%
    select(-val, -dir)
  
  
  return(diag)
}

eval_assurance = function(probs, threshold = 0.95){
  return(as.integer(probs > threshold))
}

one_assurance = function(n_sample_half, trials, n_sim = 100, gen_priors, fit_priors, 
                         criteria,
              seed = 13182813, warmup = 500, iter = 2000, chains = 4, cores = 4,
              data = NULL){
  set.seed(seed)
  container = as.data.frame(matrix(0, ncol = 8, nrow = 14*n_sim))
  colnames(container) = c("hip","n_sim", "variable", "mean", "sd", "ess_bulk", 
                          "mcse_mean", "decision")
  
  gen_obj = init_generator_template(n_sample_half, data, trials, gen_priors, type = "gen", 
                                    warmup = warmup, iter = iter, chains = chains, cores = cores)
  template_fit = init_generator_template(n_sample_half, data, trials, fit_priors, type = "template", 
                                         warmup = warmup, iter = iter, chains = chains, cores = cores)
  gen_fit = gen_obj[["fit"]]
  template_fit = template_fit[["fit"]]
  design_data = gen_obj$design_data
  
  for(i in 1:n_sim){
    assurance = sample_assurance(design_data = design_data, 
                                 gen_fit = gen_fit, template_fit = template_fit, criteria = criteria,
                                 warmup = warmup, iter = iter, chains = chains, cores = cores)
    assurance_dat = as.data.frame(
      cbind(1:14, i, assurance)
    )
    container[((i-1)*14+1):(i*14),] = assurance_dat
  }
  
  return(container)
}

#-------------------------------------------------------------------------------
# COUNTERBALANCING FUNCTIONS
#-------------------------------------------------------------------------------

#interference penalities
interference_pos = function(dat_tidy, interference_levels){
  return(
    dat_tidy %>%
      count(timing,interference,run_inside_timing,name = "n" ) %>%
      complete(timing,interference = interference_levels,run_inside_timing = 1:25, fill = list(n = 0)) %>%
      group_by(timing) %>%
      mutate(
        count = n / sum(n)
      ) %>%
      ungroup()
  )
}  

p_interference_pos = function(dat_tidy, interference_levels){
  return(
    interference_pos(dat_tidy, interference_levels) %>%
      group_by(timing) %>%
      mutate(expected =1 / (5 * 25),
             penalty = (count - expected)^2) %>%
      summarise(m_penalty = mean(penalty), .groups = "drop") %>%
      summarise(sum_penalty = sum(m_penalty)) %>%
      as.numeric()
  )
}

interference_trans = function(dat_tidy, interference_levels){
  
  all_trans = expand_grid(
    interference = interference_levels,
    next_interference = interference_levels
  ) %>%
    mutate(
      trans = paste0(interference, "_", next_interference)
    ) %>%
    pull(trans)
  
  dat_tidy %>%
    group_by(design_id, timing) %>%
    arrange(run_inside_timing, .by_group = TRUE) %>%
    mutate(
      next_interference = lead(interference)
    ) %>%
    filter(!is.na(next_interference)) %>%
    mutate(
      trans = paste0(interference, "_", next_interference)
    ) %>%
    ungroup() %>%
    count(timing, trans, name = "n") %>%
    complete(
      timing,
      trans = all_trans,
      fill = list(n = 0)
    ) %>%
    group_by(timing) %>%
    mutate(
      count = n / sum(n)
    ) %>%
    ungroup()
}

p_interference_trans = function(dat_tidy, interference_levels){
  return(
    interference_trans(dat_tidy, interference_levels) %>%
      group_by(timing) %>%
      mutate(expected =1 / (5 * 5),
             penalty = (count - expected)^2) %>%
      summarise(m_penalty = mean(penalty), .groups = "drop") %>%
      summarise(sum_penalty = sum(m_penalty)) %>%
      as.numeric()
  )
}

p_interference_total = function(dat_tidy, weights, interference_levels){
  return(weights[1] * p_interference_pos(dat_tidy, interference_levels) + 
           weights[2] * p_interference_trans(dat_tidy, interference_levels))
}

# montecarlo optimisation
monte_carlo_opt = function(dat_balanced_t, interference_levels, monte_carlo_n, seeds, weights){
  interference_penalties = matrix(0, nrow = monte_carlo_n, ncol = 3)
  
  for(i in 1:monte_carlo_n){
    if(i %% 100 == 0){
      print(paste0("Iter: ",i))
    }
    set.seed(seeds[i])
    temp = dat_balanced_t %>%
      group_by(design_id, timing) %>%
      mutate(interference = sample(rep(interference_levels, each = 5), 25, replace = FALSE)) %>%
      ungroup()
    
    temp_pos = p_interference_pos(temp, interference_levels)
    temp_trans = p_interference_trans(temp, interference_levels)
    interference_penalties[i, ] = c(temp_pos, 
                                    temp_trans,
                                    weights[1]*temp_pos + weights[2]*temp_trans)
    
  }
  return(interference_penalties)
}

select_best_designs = function(top_n, interference_penalties, seeds, weights, dat_balanced_t){
  top_idx = order(interference_penalties[,3])[1:top_n]
  top_seeds = seeds[top_idx]
  top_score = interference_penalties[top_idx, 3]
  
  top_designs = vector("list", top_n)
  best = 999
  top_design_id = -999
  
  for(i in 1:top_n){
    set.seed(top_seeds[i])
    curr_balanced_t = dat_balanced_t %>%
      group_by(design_id, timing) %>%
      mutate(interference = sample(rep(interference_levels, each = 5), 25, replace = FALSE)) %>%
      ungroup()
    
    top_designs[[i]] = curr_balanced_t
    
    curr_p = p_interference_total(curr_balanced_t, weights, interference_levels)
    if(curr_p < best){
      top_design_id = i
      best = curr_p
    } 
  }
  return(list(top_designs = top_designs,
              top_idx = top_idx, 
              top_seeds = top_seeds, 
              top_score = top_score, 
              top_design_id = top_design_id))
}

hamming_distance = function(dat_tidy1, dat_tidy2, is.interference = TRUE){
  if(is.interference){
    sum(dat_tidy1$interference != dat_tidy2$interference)
  } else {
    sum(dat_tidy1$movement_trial_id != dat_tidy2$movement_trial_id)
  }
}

select_maximin_greedy = function(dist_mat, target_n = 10, start_id = NULL) {
  
  n_designs = nrow(dist_mat)
  if(is.null(start_id)){
    dist_temp = dist_mat
    diag(dist_temp) = -Inf
    max_pair = which(dist_temp == max(dist_temp),arr.ind = TRUE)[1, ]
    selected = max_pair
    
  } else {
    selected = start_id
  }
  
  while(length(selected) < target_n){
    candidates = setdiff(1:n_designs, selected)
    candidate_scores = tibble(
      id = candidates,
      min_dist = sapply(
        candidates,
        function(i) min(dist_mat[i, selected])
      ),
      mean_dist = sapply(
        candidates,
        function(i) mean(dist_mat[i, selected])
      )
    )
    
    # Maximize minimum distance first;
    # maximize mean distance only as tie-breaker
    next_id = candidate_scores %>%
      arrange(desc(min_dist), desc(mean_dist)) %>%
      slice(1) %>%
      pull(id)
    
    selected = c(selected, next_id)
  }
  
  return(selected)
}

#--------------------------
# simulated annealing
neighbour_i = function(dat_tidy){
  selected_design = sample(unique(dat_tidy$design_id),1)
  selected_timing = sample(unique(dat_tidy$timing), 1)
  
  idx = which(
    dat_tidy$design_id == selected_design &
      dat_tidy$timing == selected_timing
  )
  swap_idx = sample(idx,2)
  dat_tidy$interference[swap_idx] = dat_tidy$interference[rev(swap_idx)]
  
  return(dat_tidy)
}


sim_annealing = function(dat_start, 
                         obj_fun, 
                         weights,
                         neighbour_fun, 
                         levels,
                         n_iter = 10000, 
                         alpha = 0.9995,
                         target_accept = 0.8, 
                         calib_iter = 500
){
  #calibrate Temperature
  delta = replicate(calib_iter, {
    proposal = neighbour_fun(dat_start)
    obj_fun(proposal, weights, levels) - obj_fun(dat_start, weights, levels)
  })
  delta_pos = delta[delta>0]
  temp_start = median(delta_pos) / (-log(target_accept))
  
  #init
  current = dat_start
  current_score = obj_fun(current, weights, levels)
  best = current
  best_score = current_score
  temp = temp_start
  
  container = data.frame(
    iter = 0:n_iter,
    current = NA_real_,
    best = NA_real_,
    temperature = NA_real_,
    alpha = NA_real_
  )
  
  container$current = current_score
  container$best = best_score
  container$temperature = temp
  container$alpha = alpha
  
  #SA loop
  for(i in 1:n_iter){
    if(i %% 100 == 0){
      print(paste0("Iter: ", i))
    }
    proposal = neighbour_fun(current)
    proposal_score = obj_fun(proposal, weights, levels)
    delta = proposal_score - current_score
    
    accept = FALSE
    if(delta <=0){
      accept = TRUE
    } else {
      p_accept = exp(-delta/temp)
      accept = runif(1) < p_accept
    }
    
    if(accept){
      current = proposal
      current_score = proposal_score
      if(current_score < best_score){
        best = current
        best_score = current_score
      }
    }
    
    temp = alpha * temp
    
    container$current[i+1] = current_score
    container$best[i+1] = best_score
    container$temperature[i+1] = temp
    
  }
  
  list(
    best_design = best,
    best_score = best_score,
    container = container,
    temp_start = temp_start
  )
}

#trial position functions 
trial_pos = function(dat_tidy){
  
  dat_tidy %>%
    count(
      timing,
      movement_trial_id,
      run_inside_timing,
      name = "n"
    ) %>%
    complete(
      timing,
      movement_trial_id = 1:25,
      run_inside_timing = 1:25,
      fill = list(n = 0)
    ) %>%
    group_by(timing) %>%
    mutate(
      count = n / sum(n)
    ) %>%
    ungroup()
}


p_trial_pos = function(dat_tidy){
  
  trial_pos(dat_tidy) %>%
    mutate(
      expected = 1 / (25 * 25),
      penalty = (count - expected)^2
    ) %>%
    group_by(timing) %>%
    summarise(
      m_penalty = mean(penalty),
      .groups = "drop"
    ) %>%
    summarise(
      sum_penalty = sum(m_penalty)
    ) %>%
    pull(sum_penalty)
}

trial_pos_interference = function(dat_tidy, interference_levels){
  
  dat_tidy %>%
    count(
      timing,
      interference,
      movement_trial_id,
      run_inside_timing,
      name = "n"
    ) %>%
    complete(
      timing,
      interference = interference_levels,
      movement_trial_id = 1:25,
      run_inside_timing = 1:25,
      fill = list(n = 0)
    ) %>%
    group_by(timing, interference) %>%
    mutate(
      count = n / sum(n)
    ) %>%
    ungroup()
}


p_trial_pos_interference = function(dat_tidy, interference_levels){
  
  trial_pos_interference(dat_tidy, interference_levels) %>%
    mutate(
      expected = 1 / (25 * 25),
      penalty = (count - expected)^2
    ) %>%
    group_by(timing, interference) %>%
    summarise(
      m_penalty = mean(penalty),
      .groups = "drop"
    ) %>%
    summarise(
      sum_penalty = sum(m_penalty)
    ) %>%
    pull(sum_penalty)
}

p_trial_total = function(dat_tidy, weights, interference_levels){
  
  return(
    weights[1] * p_trial_pos(dat_tidy) +
      weights[2] * p_trial_pos_interference(dat_tidy, interference_levels)
  )
}

assign_paired_trials = function(dat, interference_levels){
  
  # Randomly assign the 25 trials to the five interference conditions
  trial_map = data.frame(
    design_id = rep(unique(dat$design_id), each = 25),
    interference = rep(
      rep(interference_levels, each = 5),
      times = length(unique(dat$design_id))
    )
  ) %>%
    group_by(design_id) %>%
    mutate(
      movement_trial_id = sample(1:25, 25, replace = FALSE)
    ) %>%
    group_by(design_id, interference) %>%
    mutate(slot = 1:5) %>%
    ungroup()
  
  # Independently order those same five trials within encoding/retention
  dat %>%
    group_by(design_id, timing, interference) %>%
    mutate(slot = sample(1:5, 5, replace = FALSE)) %>%
    ungroup() %>%
    left_join(
      trial_map,
      by = c("design_id", "interference", "slot"),
      suffix = c("", ".new")
    ) %>%
    mutate(
      movement_trial_id = movement_trial_id.new
    ) %>%
    select(-movement_trial_id.new, -slot)
}

monte_carlo_opt_trial = function(
    dat_balanced_i,
    interference_levels,
    monte_carlo_n,
    seeds,
    weights
){
  
  trial_penalties = matrix(
    0,
    nrow = monte_carlo_n,
    ncol = 3
  )
  
  for(i in 1:monte_carlo_n){
    
    if(i %% 100 == 0){
      print(paste0("Iter: ", i))
    }
    
    set.seed(seeds[i])
    
    temp =  assign_paired_trials(
      dat_balanced_i,
      interference_levels
    )
    
    temp_pos = p_trial_pos(temp)
    
    temp_pos_interference =
      p_trial_pos_interference(
        temp,
        interference_levels
      )
    
    trial_penalties[i, ] = c(
      temp_pos,
      temp_pos_interference,
      weights[1] * temp_pos +
        weights[2] * temp_pos_interference
    )
  }
  
  return(trial_penalties)
}

select_best_trial_designs = function(
    top_n,
    trial_penalties,
    seeds,
    weights,
    dat_balanced_i,
    interference_levels
){
  
  top_idx = order(trial_penalties[,3])[1:top_n]
  top_seeds = seeds[top_idx]
  top_score = trial_penalties[top_idx,3]
  
  top_designs = vector("list", top_n)
  
  best = Inf
  top_design_id = NA
  
  for(i in 1:top_n){
    
    set.seed(top_seeds[i])
    
    curr_design = assign_paired_trials(
      dat_balanced_i,
      interference_levels
    )
    
    top_designs[[i]] = curr_design
    
    curr_p = p_trial_total(
      curr_design,
      weights,
      interference_levels
    )
    
    if(curr_p < best){
      top_design_id = i
      best = curr_p
    }
  }
  
  return(
    list(
      top_designs = top_designs,
      top_idx = top_idx,
      top_seeds = top_seeds,
      top_score = top_score,
      top_design_id = top_design_id
    )
  )
}

neighbour_trial = function(dat_tidy){
  
  selected_design = sample(
    unique(dat_tidy$design_id),
    1
  )
  
  trial_ids = sample(
    unique(
      dat_tidy$movement_trial_id[
        dat_tidy$design_id == selected_design
      ]
    ),
    2
  )
  
  idx1 = which(
    dat_tidy$design_id == selected_design &
      dat_tidy$movement_trial_id == trial_ids[1]
  )
  
  idx2 = which(
    dat_tidy$design_id == selected_design &
      dat_tidy$movement_trial_id == trial_ids[2]
  )
  
  dat_tidy$movement_trial_id[idx1] = trial_ids[2]
  dat_tidy$movement_trial_id[idx2] = trial_ids[1]
  
  return(dat_tidy)
}