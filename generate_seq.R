#-------------------------------------------------------------------------------
# GENERATE SEQ/XML FILES BASED ON DESIGN
#-------------------------------------------------------------------------------
AUDIO_DIR = "C:\\audio"

audio_file = function(filename) {
  path = paste0(AUDIO_DIR,"\\", filename)
  
  return(path)
}

trial_silence_file = function(movement_trial_id){
  sprintf(
    "trial_%02d.wav",
    movement_trial_id
  )
}

source("prereg_code/util.R")


#-------------------------------------------------------------------------------
# DEPENDENCIES
#-------------------------------------------------------------------------------
library(xml2)
#-------------------------------------------------------------------------------
# EVENT FUNCTIONS
#-------------------------------------------------------------------------------


new_step = function(id, step, type, ...){
  x = xml_new_root("scriptStep")
  
  defaults = list(
    ID = id,
    Step = step,
    StepType = type,
    LeftArmAngle = "NaN",
    RightArmAngle = "NaN",
    LeftLegAngle = "NaN",
    RightLegAngle = "NaN",
    Speed = 0,
    WaitTime = 0,
    WaitforButton = "false",
    AudioFileName = "",
    LeftArmFreerunState = "false",
    RightArmFreerunState = "false",
    LeftLegFreerunState = "false",
    RightLegFreerunState = "false"
  )
  
  values = modifyList(defaults, list(...))
  
  for (name in names(values)) {
    xml_add_child(x, name, values[[name]])
  }
  
  return(x)
}

set_to_angle = function(id, step, angle, speed){
  new_step(
    id = id,
    step = step,
    type = "SetToAngle",
    LeftArmAngle = angle,
    Speed = speed
  )
}

play_audio = function(id, step, filename){
  new_step(
    id = id,
    step = step,
    type = "PlayAudio",
    AudioFileName = audio_file(filename)
  )
}

#DEPRECIETED
bounded_waits = function(id, start_step, time, max_wait = 10){
  
  # Proprioceptor accepts integer seconds only.
  # Nearest integer minimizes timing discrepancy.
  total_time = as.integer(round(time))
  
  if(total_time < 1){
    stop("Wait time must be at least 1 second.")
  }
  
  # e.g. 23 sec -> c(10, 10, 3)
  n_full = total_time %/% max_wait
  remainder = total_time %% max_wait
  
  wait_times = c(
    rep(max_wait, n_full),
    if(remainder > 0) remainder
  )
  
  waits = lapply(
    seq_along(wait_times),
    function(i){
      wait(
        id = id,
        step = start_step + i - 1,
        time = wait_times[i]
      )
    }
  )
  
  return(waits)
}

wait = function(id, step, time){
  new_step(
    id = id,
    step = step,
    type = "Wait",
    WaitTime = time
  )
}

wait_for_button = function(id, step){
  new_step(
    id = id,
    step = step,
    type = "Wait",
    WaitforButton = "true"
  )
}

set_freerun = function(id, step, state){
  new_step(
    id = id,
    step = step,
    type = "SetFreerunState",
    LeftArmFreerunState = state
  )
}

stop_with_button = function(id, step, angle, speed){
  new_step(
    id = id,
    step = step,
    type = "StopWithButton",
    LeftArmAngle = angle,
    Speed = speed
  )
}

#-------------------------------------------------------------------------------
# TRIAL BUILDER
#-------------------------------------------------------------------------------

build_trial = function(events){
  trial = xml_new_root("ArrayOfScriptStep")
  
  for (event in events) {
    xml_add_child(trial, event)
  }
  
  return(trial)
}

#-------------------------------------------------------------------------------
# TIMING AND MOVEMENT AWARE TRIALS
#-------------------------------------------------------------------------------
encoding_active_trial = function(
    id,
    start,
    target,
    speed,
    interference,
    movement_trial_id
){
  
  silence_file = trial_silence_file(movement_trial_id)

  
  events = c(
    list(
      set_to_angle(id, 1, start, 30),
      play_audio(id, 2, "start_position.wav"),
      play_audio(id, 3, paste0(interference, ".wav")),
      set_to_angle(id, 4, target, speed),
      play_audio(id, 5, "target_position.wav"),
      wait(id, 6, 4),
      set_to_angle(id, 7, start, 30),
      play_audio(id, 8, "interference_end.wav"),
      play_audio(id, 9, silence_file),
      play_audio(id, 10, "reproduction.wav"),
      set_freerun(id, 11, "true"),
      wait_for_button(id, 12),
      set_freerun(id, 13, "false")
    )
  )
  
  build_trial(events)
}

encoding_passive_trial = function(
    id,
    start,
    target,
    speed,
    interference,
    movement_trial_id
){
  
  reproduction_angle = ifelse(target > start, 150, 30)
  
  silence_file = trial_silence_file(movement_trial_id)
  
  events = c(
    list(
      set_to_angle(id, 1, start, 30),
      play_audio(id, 2, "start_position.wav"),
      play_audio(id, 3, paste0(interference, ".wav")),
      set_to_angle(id, 4, target, speed),
      play_audio(id, 5, "target_position.wav"),
      wait(id, 6, 4),
      set_to_angle(id, 7, start, 30),
      play_audio(id, 8, "interference_end.wav"),
      play_audio(id, 9, silence_file),
      play_audio(id, 10, "reproduction.wav"),
      stop_with_button(id,11,reproduction_angle,speed)
    )
  )
  
  build_trial(events)
}

retention_active_trial = function(
    id,
    start,
    target,
    speed,
    interference,
    movement_trial_id
){
  
  silence_file = trial_silence_file(movement_trial_id)
  events = c(
    list(
      set_to_angle(id, 1, start, 30),
      play_audio(id, 2, "start_position.wav"),
      set_to_angle(id, 3, target, speed),
      play_audio(id, 4, "target_position.wav"),
      wait(id, 5, 4),
      set_to_angle(id, 6, start, 30),
      play_audio(id, 7, paste0(interference, ".wav")),
      play_audio(id, 8, silence_file),
      play_audio(id, 9, "interference_end.wav"),
      play_audio(id, 10, "reproduction.wav"),
      set_freerun(id, 11, "true"),
      wait_for_button(id,12),
      set_freerun(id, 13, "false")
    )
  )
  
  build_trial(events)
}

retention_passive_trial = function(
    id,
    start,
    target,
    speed,
    interference,
    movement_trial_id
){
  
  reproduction_angle = ifelse(
    target > start,
    150,
    30
  )
  
  silence_file = trial_silence_file(movement_trial_id)
  
  
  events = c(
    list(
      set_to_angle(id, 1, start, 30),
      play_audio(id, 2, "start_position.wav"),
      set_to_angle(id, 3, target, speed),
      play_audio(id, 4, "target_position.wav"),
      wait(id, 5, 4),
      set_to_angle(id, 6, start, 30),
      play_audio(id, 7, paste0(interference, ".wav")),
      play_audio(id, 8, silence_file),
      play_audio(id, 9, "interference_end.wav"),
      play_audio(id, 10, "reproduction.wav"),
      stop_with_button(id,11,reproduction_angle,speed)
    )
  )
  
  build_trial(events)
}



generate_xml_full = function(randomisation, trials, group) {
  
  #---------------------------------------------------------------------------
  # JOIN RANDOMISATION WITH TRIAL PARAMETERS
  #---------------------------------------------------------------------------
  
  design = merge(
    randomisation,
    trials,
    by = "movement_trial_id",
    sort = FALSE
  )
  design = design[order(design$run_order), ]
  
  
  #---------------------------------------------------------------------------
  # CREATE SEQUENCE ROOT
  #---------------------------------------------------------------------------
  
  sequence = xml_new_root("ArrayOfArrayOfScriptStep")
  
  xml_set_attr(
    sequence,
    "xmlns:xsd",
    "http://www.w3.org/2001/XMLSchema"
  )
  
  xml_set_attr(
    sequence,
    "xmlns:xsi",
    "http://www.w3.org/2001/XMLSchema-instance"
  )
  
  
  #---------------------------------------------------------------------------
  # GENERATE TRIALS IN RANDOMISED ORDER
  #---------------------------------------------------------------------------
  
  for(i in seq_len(nrow(design))){
    
    trial = design[i, ]
    
    id = trial$run_order
    timing = trial$timing
    
    start = trial$start
    target = trial$target
    speed = trial$speed
    interference = trial$interference
    
    movement_trial_id = trial$movement_trial_id
    
    
    #-----------------------------------------------------------------------
    # SELECT APPROPRIATE TRIAL TYPE
    #-----------------------------------------------------------------------
    
    if(timing == "encoding" && group == "active"){
      
      trial_xml = encoding_active_trial(
        id = id,
        start = start,
        target = target,
        speed = speed,
        interference = interference,
        movement_trial_id = movement_trial_id
      )
      
    } else if(timing == "encoding" && group == "passive"){
      
      trial_xml = encoding_passive_trial(
        id = id,
        start = start,
        target = target,
        speed = speed,
        interference = interference,
        movement_trial_id = movement_trial_id
      )
      
    } else if(timing == "retention" && group == "active"){
      
      trial_xml = retention_active_trial(
        id = id,
        start = start,
        target = target,
        speed = speed,
        interference = interference,
        movement_trial_id = movement_trial_id
      )
      
    } else if(timing == "retention" && group == "passive"){
      
      trial_xml = retention_passive_trial(
        id = id,
        start = start,
        target = target,
        speed = speed,
        interference = interference,
        movement_trial_id = movement_trial_id
      )
      
    } 
    
    #-----------------------------------------------------------------------
    # ADD BREAK SIGNAL AFTER FIRST 25 TRIALS
    #-----------------------------------------------------------------------
    if(trial$run_order == 25){
      
      next_step = length(xml_children(trial_xml)) + 1
      
      xml_add_child(
        trial_xml,
        play_audio(
          id = id,
          step = next_step,
          filename = "break.wav"
        )
      )
    }
    #-----------------------------------------------------------------------
    # ADD TRIAL TO SEQUENCE
    #-----------------------------------------------------------------------
    
    xml_add_child(sequence, trial_xml)
  }
  
  
  #---------------------------------------------------------------------------
  # RETURN COMPLETE XML
  #---------------------------------------------------------------------------
  
  return(sequence)
}


save_seq = function(xml, filename, version = "2.3.2.0"){
  xml_text = as.character(xml)
  xml_text = sub(
    '^<\\?xml version="1\\.0" encoding="UTF-8"\\?>',
    '<?xml version="1.0"?>',
    xml_text
  )
  seq_text = paste0(version, xml_text)
  writeLines(seq_text, filename, useBytes = TRUE)
}

#-------------------------------------------------------------------------------
# BATCH GENERATE .seq FILES FOR EVERY DESIGN x GROUP (from randomisation.R output)
#-------------------------------------------------------------------------------
generate_all_seq_files = function(design_trials, trial_params,
                                   out_dir = "proprio_seq"){
  for(id in sort(unique(design_trials$design_id))){

      randomisation = design_trials %>%
        filter(design_id == id) %>%
        select(
          run_order,
          movement_trial_id,
          timing,
          interference
        ) %>%
        arrange(run_order)

      trials = trial_params %>%
        select(
          movement_trial_id,
          start,
          target,
          amplitude,
          speed
        )
      
    for(group in c("active", "passive")){
      seq_xml = generate_xml_full(randomisation = randomisation, trials = trials, group = group)
      save_seq(seq_xml, filename = file.path(out_dir, paste0("design_", id, "_", group, ".seq")))
    }

    print(paste0("design ", id, " done"))
  }
}


library(readr)
library(dplyr)

design_trials = read_csv("material/design_trials.csv")
trial_params = read_csv("material/trial_params.csv")
generate_all_seq_files(design_trials, trial_params)

