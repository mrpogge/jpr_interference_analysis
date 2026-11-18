#-------------------------------------------------------------------------------
# FINDING TRIAL POSITION SETTINGS
#-------------------------------------------------------------------------------

set.seed(13181913)

#-------------------------------------------------------------------------------
# DEPENDENCIES
#-------------------------------------------------------------------------------
library(tidyverse)
library(plotly)
library(lme4)
library(cluster)
library(maximin)

#-------------------------------------------------------------------------------
# SPECIFYING PARAMETER SPACE
#-------------------------------------------------------------------------------
all_trials = expand.grid(
  amplitude = 30:120,
  target = 30:150,
  speed= 10:30,
  direction=c(-1,1)
) %>%
  mutate(start = target + direction * amplitude,
         is_valid = ifelse(start <= 150 & start >= 30, 1,0))

valid_trials = all_trials %>% 
  filter(is_valid == 1) %>%
  select(start, target, amplitude, speed, direction) %>%
  mutate(start_c = (start-30)/120,
         target_c = (target-30)/120,
         amplitude_c = (amplitude-30)/90,
         speed_c = (speed-10)/20,
         direction_c = ifelse(direction == -1, 0, 1))

X_all = valid_trials %>%
  select(
    target_c,
    amplitude_c,
    speed_c,
    direction_c
  ) %>%
  as.matrix()

set.seed(13181913)

all_mm = maximin.cand(
  n = 25,
  Xcand = X_all
)

final_trials = valid_trials[all_mm$inds, ] %>%
  select(
    start,
    target,
    amplitude,
    speed,
    direction
  )

plot_ly(
  data = final_trials,
  x = ~start,
  y = ~target,
  z = ~speed,
  type = "scatter3d"
)
final_trials %>%
  ggplot(aes(x = start, y = amplitude, colour = speed)) +
  geom_point() 

final_trials %>%
  group_by(direction) %>%
  summarise(
    n = n(),
    amp_min = min(amplitude),
    amp_mean = mean(amplitude),
    amp_sd = sd(amplitude),
    amp_max = max(amplitude),
    
    speed_min = min(speed),
    speed_mean = mean(speed),
    speed_sd = sd(speed),
    speed_max = max(speed)
  )

  
write_csv(final_trials, "material/selected_trials.csv")
final_trials = read_csv("material/selected_trials.csv")
