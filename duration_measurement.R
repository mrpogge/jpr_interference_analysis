#-------------------------------------------------------------------------------
# DEPENDENCIES
#-------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lubridate)
library(av)
library(tuneR)


dat_raw_000 = read_excel("P000_raw.xlsx")
trial_params = read.csv("material/trial_params.csv")
design_trials = read.csv("material/design_trials.csv")

curr_design = design_trials %>%
  filter(design_id == 50) %>%
  left_join(trial_params, by = "movement_trial_id")

dat = dat_raw_000 %>%
  filter(
    !is.na(Feladat),
    !is.na(Lepes),
    !is.na(`BAL KAR`)
  ) %>%
  mutate(
    run_order = as.integer(Feladat),
    Lepes = as.integer(Lepes),
    time = period_to_seconds(hms(Ido)),
    timing = if_else(
      Feladat <= 25,
      "encoding",
      "retention"
    )
  )

#-------------------------------------------------------------------------------
# EXTRACT ROTATIONS
#-------------------------------------------------------------------------------

dat_rotations = dat %>%
  mutate(rotation_time = case_when(timing == "retention" & Lepes >= 3 & Lepes <= 6 ~ "rotate",
                              timing == "encoding" & Lepes >= 4 & Lepes <= 7 ~ "rotate")) %>%
  group_by(run_order,timing) %>%
  mutate(time_rel = time - first(time), angle = `BAL KAR`) %>%
  ungroup()


rotation_summary = dat_rotations[73:nrow(dat_rotations),] %>% #the first 73 was a test which was logged by the proprioceptor so we dont need it
  group_by(run_order) %>%
  filter(rotation_time == "rotate") %>%
  summarise(rotation_time = last(time_rel) - first(time_rel))


curr_design = curr_design %>%
  left_join(rotation_summary, "run_order")

duration_to_trial = curr_design %>%
  select(movement_trial_id, rotation_time) %>%
  group_by(movement_trial_id) %>%
  summarise(duration = mean(rotation_time))

write.csv(duration_to_trial, "material/duration_to_trial.csv")

#-------------------------------------------------------------------------------
# Generate silence for the exact retention period
#-------------------------------------------------------------------------------
make_silence_wav <- function(duration_sec,
                             filename,
                             sample_rate = 44100) {
  
  n_samples <- round(duration_sec * sample_rate)
  
  silence <- Wave(
    left = rep(0, n_samples),
    samp.rate = sample_rate,
    bit = 16
  )
  
  writeWave(silence, filename)
}


for(i in seq_len(nrow(duration_to_trial))) {
  
  make_silence_wav(
    duration_sec = duration_to_trial$duration[i],
    filename = file.path(
      "material",
      sprintf(
        "trial_%02d.wav",
        duration_to_trial$movement_trial_id[i]
      )
    )
  )
}


