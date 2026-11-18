#-------------------------------------------------------------------------------
# SELECT BALANCED POSITIONS
#-------------------------------------------------------------------------------

set.seed(13181913)
#-------------------------------------------------------------------------------
# DEPENDENCIES
#-------------------------------------------------------------------------------
library(gtools)
library(tidyverse)
source("prereg_code/util.R")

monte_carlo_n = 10000
n_annealig = 10000

#-------------------------------------------------------------------------------
# DATAFORMAT
#-------------------------------------------------------------------------------
n_samp_half = 50

dat_balanced = data.frame(
  design_id = rep(1:n_samp_half, each = 50),
  run_order = rep(1:50, times = n_samp_half),
  run_inside_timing = rep(1:25, times = n_samp_half*2),
  timing = 0,
  interference = 0,
  movement_trial_id = 0
) 

#-------------------------------------------------------------------------------
#TIMING ORDER BALANCE
#-------------------------------------------------------------------------------
# ER - 0, RE -1
timing_assignment = data.frame(
  design_id = 1:50,
  timing_assignment = sample(rep(c(0,1), each = n_samp_half/2), n_samp_half, replace = FALSE)
)
  
dat_balanced_t = dat_balanced %>%
  left_join(timing_assignment, by = "design_id") %>%
  mutate(timing =as.factor(case_when(
    timing_assignment == 1 & run_order <= 25 ~ "retention",
    timing_assignment == 1 & run_order > 25  ~ "encoding",
    timing_assignment == 0 & run_order <= 25 ~ "encoding",
    timing_assignment == 0 & run_order > 25  ~ "retention"
  ))) 

#timing balance check:
dat_balanced_t %>%
  distinct(design_id, timing_assignment) %>%
  count(timing_assignment) #TIMING POS BALANCE SATISFIED!!

#-------------------------------------------------------------------------------
#INTERFERENCE WITHIN TIMING ORDER ANFD TRANSFER BALANCE
#-------------------------------------------------------------------------------

interference_levels = c(
  "control",
  "motor",
  "spatial",
  "visual",
  "verbal"
)

#-------------------------------------------------------------------------------
#FIND WEIGHTS
#-------------------------------------------------------------------------------
emp_weights = matrix(0, nrow = 1000, ncol = 2)

for(i in 1:1000){
  if(i %% 100 == 0){
    print(paste0("Iter: ",i))
    }
  temp = dat_balanced_t %>%
    group_by(design_id, timing) %>%
    mutate(interference = sample(rep(interference_levels, each = 5), 25, replace = FALSE)) %>%
    ungroup()
  
  emp_weights[i, ] = c(p_interference_pos(temp, interference_levels), 
                       p_interference_trans(temp, interference_levels))
}

weights_interference = c(1/sqrt(sd(emp_weights[,1])), 1/sqrt(sd(emp_weights[,2])))
w = weights_interference / sum(weights_interference)

#-------------------------------------------------------------------------------
#MONTE CARLO SEARCH
#-------------------------------------------------------------------------------

seeds = sample.int(.Machine$integer.max, monte_carlo_n)
interference_penalties = monte_carlo_opt(dat_balanced_t, interference_levels, monte_carlo_n, seeds, w)

#best selection
set.seed(seeds[which.min(interference_penalties[,3])])
dat_balanced_i = dat_balanced_t %>%
  group_by(design_id, timing) %>%
  mutate(interference = sample(rep(interference_levels, each = 5), 25, replace = FALSE)) %>%
  ungroup()

#validate
interference_pos(dat_balanced_i, interference_levels)
interference_trans(dat_balanced_i,interference_levels)
p_interference_total(dat_balanced_i, w, interference_levels)


#best 100 selection
top_ns = 100
best_designs = select_best_designs(top_ns, interference_penalties, seeds, w, dat_balanced_t)
top_designs = best_designs$top_designs
best_id = best_designs$top_design_id
dist_mat_i = matrix(0, nrow = top_ns, ncol = top_ns)

for(i in 1:(top_ns-1)){
  for(j in (i+1):top_ns){
    d = hamming_distance(top_designs[[i]], top_designs[[j]])
    dist_mat_i[i,j] = d
    dist_mat_i[j,i] = d
  }
}

selected_ids = select_maximin_greedy(dist_mat_i, 10, best_id)

#-------------------------------------------------------------------------------
#SIMULATED ANNEALING
#-------------------------------------------------------------------------------
sa_results = lapply(
  selected_ids,
  function(i){
    print(paste0("Starting point no.:", i))
    sim_annealing(
      dat_start = top_designs[[i]],
      obj_fun = p_interference_total,
      levels = interference_levels,
      weights = w,
      neighbour_fun = neighbour_i,
      n_iter = n_annealig
    )
  }
)

sa_scores = sapply(sa_results, function(x) x$best_score)
best_sa = sa_results[[which.min(sa_scores)]]
dat_balanced_i = best_sa$best_design

interference_pos(dat_balanced_i, interference_levels)
interference_trans(dat_balanced_i,interference_levels)
p_interference_total(dat_balanced_i, w, interference_levels)

#-------------------------------------------------------------------------------
# TRIAL POSITION OPTIMISATION
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# FIND WEIGHTS
#-------------------------------------------------------------------------------

emp_weights_trial = matrix(
  0,
  nrow = 1000,
  ncol = 2
)

for(i in 1:1000){
  
  if(i %% 100 == 0){
    print(paste0("Iter: ", i))
  }
  
  temp = assign_paired_trials(dat_balanced_i, interference_levels)
  
  emp_weights_trial[i, ] = c(
    p_trial_pos(temp),
    p_trial_pos_interference(
      temp,
      interference_levels
    )
  )
}

weights_trial = c(
  1 / sqrt(sd(emp_weights_trial[,1])),
  1 / sqrt(sd(emp_weights_trial[,2]))
)

w_trial = weights_trial / sum(weights_trial)

#-------------------------------------------------------------------------------
# MONTE CARLO SEARCH
#-------------------------------------------------------------------------------
seeds_trial = sample.int(
  .Machine$integer.max,
  monte_carlo_n
)

trial_penalties = monte_carlo_opt_trial(
  dat_balanced_i = dat_balanced_i,
  interference_levels = interference_levels,
  monte_carlo_n = monte_carlo_n,
  seeds = seeds_trial,
  weights = w_trial)
  

set.seed(
  seeds_trial[
    which.min(trial_penalties[,3])
  ]
)
  
dat_balanced_trial_mc = assign_paired_trials(
  dat_balanced_i,
  interference_levels
)

trial_pos(dat_balanced_trial_mc)

trial_pos_interference(
  dat_balanced_trial_mc,
  interference_levels
)

p_trial_total(
  dat_balanced_trial_mc,
  w_trial,
  interference_levels
)

#-------------------------------------------------------------------------------
# BEST 100
#-------------------------------------------------------------------------------

best_trial_designs = select_best_trial_designs(
  top_n = 100,
  trial_penalties = trial_penalties,
  seeds = seeds_trial,
  weights = w_trial,
  dat_balanced_i = dat_balanced_i,
  interference_levels = interference_levels
)

top_trial_designs = best_trial_designs$top_designs

best_trial_id = best_trial_designs$top_design_id

top_n = length(top_trial_designs)

dist_mat_trial = matrix(
  0,
  nrow = top_n,
  ncol = top_n
)

for(i in 1:(top_n - 1)){
  
  for(j in (i + 1):top_n){
    
    d = hamming_distance(
      top_trial_designs[[i]],
      top_trial_designs[[j]],
      is.interference = FALSE
    )
    
    dist_mat_trial[i,j] = d
    dist_mat_trial[j,i] = d
  }
}

selected_trial_ids = select_maximin_greedy(
  dist_mat = dist_mat_trial,
  target_n = 10,
  start_id = best_trial_id
)

#-------------------------------------------------------------------------------
# SIMULATED ANNEALING
#-------------------------------------------------------------------------------

sa_trial_results = lapply(
  selected_trial_ids,
  function(i){
    
    print(
      paste0(
        "Starting point no.: ",
        i
      )
    )
    
    sim_annealing(
      dat_start = top_trial_designs[[i]],
      obj_fun = p_trial_total,
      weights = w_trial,
      neighbour_fun = neighbour_trial,
      levels = interference_levels,
      n_iter = n_annealig
    )
  }
)

sa_trial_scores = sapply(
  sa_trial_results,
  function(x) x$best_score
)

best_sa_trial =
  sa_trial_results[[which.min(sa_trial_scores)]]

dat_balanced_final =
  best_sa_trial$best_design

trial_pos(
  dat_balanced_final
)

trial_pos_interference(
  dat_balanced_final,
  interference_levels
)

p_trial_total(
  dat_balanced_final,
  w_trial,
  interference_levels
)

saveRDS(dat_balanced_final, "material/dat_balanced_final.RDS")


#-------------------------------------------------------------------------------
# EXPORT DESIGN FOR PYTHON ALLOCATOR APP
#-------------------------------------------------------------------------------

design_export = dat_balanced_final %>%
  ungroup() %>%
  mutate(
    # 1 = first 25-trial timing period
    # 2 = second 25-trial timing period
    timing_id = ifelse(run_order <= 25, 1L, 2L),
    
    # Positional five-trial chunks within each timing period
    block = ceiling(run_inside_timing / 5),
    
    # Position 1-5 within positional block
    trial_in_block = ((run_inside_timing - 1) %% 5) + 1
  ) %>%
  select(
    design_id,
    run_order,
    timing_id,
    timing,
    block,
    trial_in_block,
    interference,
    movement_trial_id
  ) %>%
  arrange(
    design_id,
    run_order
  )

write_csv(
  design_export,
  "material/design_trials.csv"
)



#-------------------------------------------------------------------------------
# EXPORT MOVEMENT TRIAL PARAMETERS
#-------------------------------------------------------------------------------

trial_params_export = read_csv(
  "material/selected_trials.csv"
) %>%
  mutate(
    movement_trial_id = row_number()
  ) %>%
  select(
    movement_trial_id,
    start,
    target,
    amplitude,
    speed,
    direction
  )

write_csv(
  trial_params_export,
  "material/trial_params.csv"
)


#-------------------------------------------------------------------------------
# FINAL DESIGN VALIDATION
#-------------------------------------------------------------------------------

# Correct total size
stopifnot(nrow(dat_balanced_final) == n_samp_half * 50)

# Each design has 25 encoding + 25 retention trials
stopifnot(
  all(
    dat_balanced_final %>%
      count(design_id, timing) %>%
      pull(n) == 25
  )
)

# Each interference occurs exactly 5 times per timing/design
stopifnot(
  all(
    dat_balanced_final %>%
      count(design_id, timing, interference) %>%
      pull(n) == 5
  )
)

# Each movement trial occurs exactly once per timing/design
stopifnot(
  all(
    dat_balanced_final %>%
      count(design_id, timing, movement_trial_id) %>%
      pull(n) == 1
  )
)

# Crucial pairing check:
# each movement_trial_id is assigned to the same interference
# in encoding and retention within a design
stopifnot(
  all(
    dat_balanced_final %>%
      distinct(design_id, timing, interference, movement_trial_id) %>%
      group_by(design_id, movement_trial_id) %>%
      summarise(
        n_timing = n_distinct(timing),
        n_interference = n_distinct(interference),
        .groups = "drop"
      ) %>%
      transmute(ok = n_timing == 2 & n_interference == 1) %>%
      pull(ok)
  )
)

# Timing-order balance: 25 ER and 25 RE designs
stopifnot(
  all(
    dat_balanced_final %>%
      distinct(design_id, timing_assignment) %>%
      count(timing_assignment) %>%
      pull(n) == n_samp_half / 2
  )
)

message("All final design validation checks passed.")

