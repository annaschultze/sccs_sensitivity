# Program Information  ----------------------------------------------------
#
# Program:      04b_scri_post.R
# Description:  run post-exposure SCRI
# Dependencies: 00_sensitivity_functions 
#               01_clean_data.R
#               02_select_population.R 
# 
# 0. HOUSEKEEPING ---------------------------------------------------------

# specify the study design
design <- "scri_post"

# 1. READ IN DATA ------------------------------------------------------------
load(paste0(dirtemp, "sccs_sensitivity/", "sccs_population", ".RData"))
dir.create(file.path(paste0(diroutput, "sccs_sensitivity/", design)),         showWarnings = FALSE, recursive = TRUE)

# STUDY INPUTS ------------------------------------------------------------
# specify study inputs and create generic variables 

  ## study design options
  control_duration  <- 59
  preexp_start  <- 0
  vax1_end          <- 28
  vax2_end          <- 28
  
  ## outcome 
  sccs_population$outcome_date <- sccs_population$myopericarditis_date
  sccs_population$outcome_days <- round(difftime(sccs_population$outcome_date, as.Date("2020-09-01"), units = "days"),0)
  sccs_population$outcome_days <- as.numeric(sccs_population$outcome_days)
  sccs_population$outcome_binary <- as.numeric(!is.na(sccs_population$outcome_date))
  
# CREATE VARIABLES --------------------------------------------------------
# create variables (in days) needed to fit the chosen study design  
  
  # pre-exposure period start, end and duration 
  sccs_population$preexp <- as.numeric(sccs_population$days_vax1 + preexp_start)
  sccs_population$preexp_end <- as.numeric(sccs_population$days_vax1)
  # start minus end + 1 as the model fit is inclusive (so if start and end is the same, modelled length is 1 day)
  sccs_population$preexp_length <- as.numeric(sccs_population$preexp_end - sccs_population$preexp) + 1
  
  # first risk window start, end and duration. Censor risk period at day before 2nd dose if occurs, as day of 2nd dose should be considered separately 
  sccs_population$risk_d1 <- as.numeric(sccs_population$days_vax1 + 1)
  sccs_population$risk_d1_end <- as.numeric(pmin((sccs_population$days_vax1 + vax1_end), sccs_population$study_exit_days, sccs_population$days_vax2-1, na.rm = T))
  # plus 1 to account for inclusive model fit
  sccs_population$risk_d1_length <- as.numeric((sccs_population$risk_d1_end - sccs_population$risk_d1)) + 1

  # second risk window start, end and duration. Censor risk period at 3rd dose if occurs. 
  sccs_population$risk_d2 <- as.numeric(sccs_population$days_vax2 + 1)
  # 2nd dose end min of dose 2 risk period, study exit and dose 3. Involved call because want to ignore missing dose 3 but return missing if missing dose 2 
  sccs_population$risk_d2_end <- ifelse(!is.na(sccs_population$days_vax2), as.numeric(pmin((sccs_population$days_vax2 + vax2_end), sccs_population$study_exit_days, sccs_population$days_vax3, na.rm = T)), NA)
  # plus 1 to account for inclusive model fit
  sccs_population$risk_d2_length <- as.numeric((sccs_population$risk_d2_end - sccs_population$risk_d2)) + 1

  # time between doses, starts day 1 after dose one risk end or dose 2 whichever is first, ends at dose2 
  # na.rm is F because we want this to be missing if days_vax2 is ever missing 
  sccs_population$between_start <- as.numeric(pmin(sccs_population$risk_d1_end + 1), sccs_population$days_vax2)
  sccs_population$between_start <- ifelse(!is.na(sccs_population$days_vax2), as.numeric(sccs_population$between_start), NA)
  sccs_population$between_end <- as.numeric(pmin(sccs_population$days_vax2, sccs_population$study_exit_days))
  # plus 1 to account for inclusive model fit
  sccs_population$between_length <- as.numeric((sccs_population$between_end - (sccs_population$between_start))) + 1
  
  # control time 
  sccs_population$c_start <- as.numeric(sccs_population$days_last_vax + vax2_end) + 1 
  sccs_population$c_end <- as.numeric(sccs_population$c_start + control_duration)
  sccs_population$c_end <- as.numeric(pmin(sccs_population$study_exit_days, sccs_population$c_end, sccs_population$days_vax3, na.rm = T))
    
  sccs_population$c_length <- as.numeric(sccs_population$c_end - sccs_population$c_start) + 1
  sccs_population$study_length <- as.numeric(sccs_population$c_end)
  
  # start at the first of control window or risk window (depending on design)
  sccs_population$ref_start <- as.numeric(sccs_population$c_start)
  
  # start at the first of control window or risk window (depending on design)
  sccs_population$ref_start <- as.numeric(sccs_population$risk_d1)
  
  # end at the last of the end of the control window, first and second risk windows (depending on design)
  sccs_population$ref_end <- as.numeric(sccs_population$c_end)
  sccs_population$ref_end <- as.numeric(sccs_population$ref_end)

  # NUMBER OF EVENTS --------------------------------------------------------
  # function to summarise the number and duration of each period 
  
  table3a.interim <- lapply(split(sccs_population, sccs_population$type_vax1), describe_scri_periods)
  table3a <- do.call(rbind, table3a.interim )
  
  # add headings
  column_headings <- c("Time period", "N events", "Median Duration (days)", "Q1", "Q3")
  table3a <- rbind(column_headings, table3a)
  
  table3a <- cbind(Names = rownames(table3a), table3a)
  table3a$Names <- gsub("\\..*","",table3a$Names)
  
  write.csv(table3a, paste0(diroutput,"sccs_sensitivity/", design, "/table_3a_describe_period_by_vaccine.csv"), row.names = FALSE)
  
  # FIT MODEL  --------------------------------------------------------------
  
  table3b.interim <- lapply(split(sccs_population, sccs_population$type_vax1), fit_scri)
  table3b <- do.call(rbind, table3b.interim)
  
  # add row headings
  table3b <- cbind(Names = rownames(table3b), table3b)
  table3b$Names <- gsub("\\..*","",table3b$Names)
  
  write.csv(table3b, paste0(diroutput,"sccs_sensitivity/", design, "/table_3b_model_results.csv"), row.names = FALSE)
  
  
