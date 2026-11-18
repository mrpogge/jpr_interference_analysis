#-------------------------------------------------------------------------------
# DURATION CALIBRATION
#-------------------------------------------------------------------------------
# Purpose:
#   Estimate a deterministic duration-matching rule from the P004 pilot logs.
#   The rule predicts forward and return rotation durations from programmed
#   movement parameters and can be used prospectively when generating sequences.
#
#   Forward model: duration ~ amplitude / speed + speed
#   Return model:  duration ~ amplitude

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# DEPENDENCIES
#-------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(lubridate)
library(av)

audio_length1 = av_media_info("material/target_position.wav")$duration
audio_length2 = av_media_info("material/interference_end.wav")$duration
duration_overhead = audio_length1
#-------------------------------------------------------------------------------
# INPUT FILES
#-------------------------------------------------------------------------------

raw_file = "P004_raw.xlsx"
processed_file = "feldolgozott_regisztratum.xlsx"

#-------------------------------------------------------------------------------
# READ DATA
#-------------------------------------------------------------------------------

dat_raw = read_excel(raw_file)
dat_processed = read_excel(processed_file)

#-------------------------------------------------------------------------------
# PREPARE RAW LOG
#-------------------------------------------------------------------------------
# For P004, retention was administered first (Feladat 1-25) and encoding second
# (Feladat 26-50).
#
# Relevant proprioceptor rotations:
#   retention: forward = step 3, return = step 6
#   encoding:  forward = step 4, return = step 7
#-------------------------------------------------------------------------------

dat = dat_raw %>%
  filter(
    !is.na(Feladat),
    !is.na(Lepes),
    !is.na(`BAL KAR`)
  ) %>%
  mutate(
    Feladat = as.integer(Feladat),
    Lepes = as.integer(Lepes),
    time = period_to_seconds(hms(Ido)),
    timing = if_else(
      Feladat <= 25,
      "retention",
      "encoding"
    )
  )

#-------------------------------------------------------------------------------
# EXTRACT ROTATIONS
#-------------------------------------------------------------------------------

dat_rotations = dat %>%
  mutate(
    rotation = case_when(
      timing == "retention" & Lepes == 3 ~ "forward",
      timing == "retention" & Lepes == 6 ~ "return",
      timing == "encoding"  & Lepes == 4 ~ "forward",
      timing == "encoding"  & Lepes == 7 ~ "return",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(rotation)) %>%
  group_by(
    Feladat,
    timing,
    rotation
  ) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    rotation_id = paste(
      Feladat,
      timing,
      rotation,
      sep = "_"
    ),
    time_rel = time - first(time),
    angle = `BAL KAR`
  ) %>%
  ungroup()

#-------------------------------------------------------------------------------
# SUMMARISE OBSERVED ROTATION DURATIONS
#-------------------------------------------------------------------------------

rotation_kinematics = dat_rotations %>%
  group_by(
    Feladat,
    timing,
    rotation
  ) %>%
  summarise(
    observed_start = first(angle),
    observed_end = last(angle),
    observed_amplitude = abs(observed_end - observed_start),
    duration = max(time_rel),
    n_samples = n(),
    .groups = "drop"
  ) %>%
  left_join(
    dat_processed %>%
      select(
        Feladat,
        movement_trial_id,
        start,
        target,
        amplitude,
        speed
      ),
    by = "Feladat"
  )

# Basic extraction checks
stopifnot(nrow(rotation_kinematics) == 100)
stopifnot(all(rotation_kinematics$duration > 0))
stopifnot(all(!is.na(rotation_kinematics$amplitude)))
stopifnot(all(!is.na(rotation_kinematics$speed)))

print(
  rotation_kinematics %>%
    count(timing, rotation)
)

#-------------------------------------------------------------------------------
# CALIBRATION MODELS
#-------------------------------------------------------------------------------

forward_dat = rotation_kinematics %>%
  filter(rotation == "forward")

return_dat = rotation_kinematics %>%
  filter(rotation == "return")

fit_forward = lm(
  duration ~ I(amplitude / speed) + speed,
  data = forward_dat
)

fit_return = lm(
  duration ~ amplitude,
  data = return_dat
)

print(summary(fit_forward))
print(summary(fit_return))


#-------------------------------------------------------------------------------
# PROSPECTIVE DURATION FUNCTION
#-------------------------------------------------------------------------------

calculate_interference_time = function(amplitude, speed){

  new_forward = data.frame(
    amplitude = amplitude,
    speed = speed
  )

  new_return = data.frame(
    amplitude = amplitude
  )

  t_forward = predict(
    fit_forward,
    newdata = new_forward
  )

  t_return = predict(
    fit_return,
    newdata = new_return
  )

  as.numeric(t_forward + 4 + t_return + duration_overhead)
}

#-------------------------------------------------------------------------------
# IN-SAMPLE PREDICTION PERFORMANCE
#-------------------------------------------------------------------------------

performance = function(obs, pred){

  residual = obs - pred

  tibble(
    n = length(obs),
    bias = mean(residual),
    MAE = mean(abs(residual)),
    RMSE = sqrt(mean(residual^2)),
    max_abs_error = max(abs(residual)),
    R2 = cor(obs, pred)^2
  )
}

forward_pred = forward_dat %>%
  mutate(
    pred_duration = predict(
      fit_forward,
      newdata = .
    )
  )

return_pred = return_dat %>%
  mutate(
    pred_duration = predict(
      fit_return,
      newdata = .
    )
  )

in_sample_performance = bind_rows(
  forward = performance(
    forward_pred$duration,
    forward_pred$pred_duration
  ),
  return = performance(
    return_pred$duration,
    return_pred$pred_duration
  ),
  .id = "rotation"
)
print(in_sample_performance)

#-------------------------------------------------------------------------------
# LEAVE-ONE-OUT CROSS-VALIDATION
#-------------------------------------------------------------------------------

loo_forward = map_dfr(
  seq_len(nrow(forward_dat)),
  function(i){

    train = forward_dat[-i, , drop = FALSE]
    test = forward_dat[i, , drop = FALSE]

    fit = lm(
      duration ~ I(amplitude / speed) + speed,
      data = train
    )

    test %>%
      mutate(
        pred_duration = as.numeric(
          predict(fit, newdata = test)
        )
      )
  }
)

loo_return = map_dfr(
  seq_len(nrow(return_dat)),
  function(i){

    train = return_dat[-i, , drop = FALSE]
    test = return_dat[i, , drop = FALSE]

    fit = lm(
      duration ~ amplitude,
      data = train
    )

    test %>%
      mutate(
        pred_duration = as.numeric(
          predict(fit, newdata = test)
        )
      )
  }
)

loo_performance = bind_rows(
  forward = performance(
    loo_forward$duration,
    loo_forward$pred_duration
  ),
  return = performance(
    loo_return$duration,
    loo_return$pred_duration
  ),
  .id = "rotation"
)

print(loo_performance)

dat_processed = dat_processed %>%
  mutate(
    pred_total = calculate_interference_time(amplitude, speed),
    total_error = duration_empirical - pred_total
  )

performance(
  dat_processed$duration_empirical,
  dat_processed$pred_total
)

#-------------------------------------------------------------------------------
# PLOTS
#-------------------------------------------------------------------------------

pred_plot_dat = bind_rows(
  forward_pred %>% mutate(rotation = "forward"),
  return_pred %>% mutate(rotation = "return")
) %>%
  mutate(
    residual = duration - pred_duration
  )

print(
  ggplot(
    pred_plot_dat,
    aes(
      x = pred_duration,
      y = duration,
      colour = rotation
    )
  ) +
    geom_point() +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2
    ) +
    coord_equal() +
    labs(
      x = "Predicted duration (s)",
      y = "Observed duration (s)",
      title = "Observed vs predicted rotation duration"
    )
)

print(
  ggplot(
    pred_plot_dat,
    aes(
      x = pred_duration,
      y = residual,
      colour = rotation
    )
  ) +
    geom_point() +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    labs(
      x = "Predicted duration (s)",
      y = "Observed - predicted (s)",
      title = "Calibration-model residuals"
    )
)

#-------------------------------------------------------------------------------
# SAVE CALIBRATION OBJECTS
#-------------------------------------------------------------------------------

calibration = list(
  fit_forward = fit_forward,
  fit_return = fit_return,
  duration_overhead = audio_length1,
  forward_coefficients = coef(fit_forward),
  return_coefficients = coef(fit_return),
  in_sample_performance = in_sample_performance,
  loo_performance = loo_performance
)

saveRDS(
  calibration,
  "duration_calibration.rds"
)

print(coef(fit_forward))
print(coef(fit_return))

