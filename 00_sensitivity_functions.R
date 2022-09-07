# Program Information  ----------------------------------------------------
#
# Program:      00_sensitivity_functions
# Description:  functions to support the sccs sensitivity analyses: 
#               *sccs_data_management
#               *sccs_flowchart
#               *sccs_selection
#               *sccs_table
#               *sccs_analysis 
#
# Dependencies: Function 5 sccs_analysis requires and loads packages tidyverse, SCCS
#
# 0. HOUSEKEEPING ---------------------------------------------------------

# if file is opened through project directory then wd doesn't need to be automatically set 
# current script assumes wd is the project directory 
if(!any(ls()=="thisdir"))   thisdir   <- paste0(getwd(),"/") 
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"g_output/")

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp, "sccs_sensitivity")),           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,"log_files/sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "sccs_sensitivity")),         showWarnings = FALSE, recursive = TRUE)

# 1. FUNCTIONS ---------------------------------------------------------------
# note; these are designed to work for each individual DAP 
# consideration for different BIFAP data sources is not yet applied 

#' Perform SCCS or SCRI data management 
#' 
#' @description 
#' Takes clean data for self-controlled designs and does design/outcome specific data management.  
#' 
#' @param data = Input dataframe, should be loaded.
#' @param outcome = Outcome variable name as a string, referencing a _date variable. 
#' @param reduce_dimensions = T drop unnecessary rows to increase speed. Default is "F", option "T" .  
#' @param design = Type of design - sccs or scri? Default "sccs", option "scri". 
#' @param control_start = Date used to start the control period (days_vax1 or days_last_vax). 
#' @param control_dur = Duration of the control period in days, NEGATIVE for prevaccine, POSITIVE for postvaccine. 
#' @param preexp = Duration of the pre-exposure period, relative to first vaccine date. Default is zero (no pre-exposure period).
#' @param risk1 = Day the risk period after the first vaccine dose ends, relative to the first vaccine date.
#' @param risk2 = Day the risk period after the second dose ends, relative to the second vaccine date.
#' 
#' @returns A dataframe with columns required for sccs or scri added. 
#' 

sccs_data_management <- function(data, 
                                 outcome, 
                                 reduce_dimensions = "F", 
                                 design = "sccs",
                                 control_start = NULL, 
                                 control_dur = NULL, 
                                 preexp = 0, 
                                 risk1, 
                                 risk2) {
  
  # print an update 
  message("starting sccs_data_management")
  
  # load tidyverse
  if(require("tidyverse")) {
    message("tidyverse is loaded")
  } else {
    message("tidyverse is not installed, attempting to install") 
    install.packages("tidyverse") 
    if(require("tidyverse")) {
      message("tidyverse is installed and loaded")
    } else {
      stop("could not install tidyverse, troubleshoot")
    }
  } 
  
  # CHECKS ON THE INPUTS 
  if (!is.character(outcome)) {
    stop("outcome input is not character")
  }
  
  if (!is.character(reduce_dimensions)) {
    stop("reduce_dimensions input is not character")
  }
  
  if (!is.character(design)) {
    stop("design input is not character")
  }
  
  if (design == "sccs") {
    if (!is.null(control_start)) {
      warning("design is sccs, start of control period will be ignored")
    }
    if (!is.null(control_dur)) {
      warning("design is sccs, length of control period will be ignored")
    }
  } 
  
  if (design == "scri") {
    if (is.null(control_start)) {
      stop("design is scri, you need to speficy control_start and control_dur")
    }
    if (is.null(control_dur)) {
      stop("design is scri, you need to speficy control_start and control_dur")
    }
  }
  
  # INPUT VARIABLES 
  control_duration  <- control_dur
  preexp_per_start  <- preexp
  vax1_end          <- risk1
  vax2_end          <- risk2
  
  # slightly hacky solution to create "first of two dates" input, not very generalisable 
  if(outcome == "myopericarditis_date"){
    
    data$outcome_date <- pmin(data$myocarditis_date,data$pericarditis_date,na.rm=T)
    data$outcome_days <- round(difftime(data$outcome_date, as.Date("2020-09-01"), units = "days"),0)
    
  } else {
    
    data$outcome_date <- data[[outcome]]
    data$outcome_days <- round(difftime(data$outcome_date, as.Date("2020-09-01"), units = "days"),0)
    
  }
  
  # DATA MANAGEMENT FOR CREATING RISK AND CONTROL TIME VARIABLES  

  # pre-exposure period start and duration 
  data$preexp_days <- data$days_vax1 + preexp_per_start
  data$preexp_dur <- (data$days_vax1 - data$preexp_days)
  
  # first risk window start, end and duration
  data$risk_d1 <- as.numeric(data$days_vax1 + 1)
  data$riskend_d1 <- pmin((data$days_vax1 + vax1_end), data$study_exit_days, data$days_vax2, na.rm = T)
  data$risk_d1_length <- (data$riskend_d1 - data$days_vax1)
  # start of between period needs to be dose 1 + risk window + 1 (as last day should go to the risk window)
  data$risk_between <- data$risk_d1_length + 1 
  
  # second risk window start, end and duration
  data$risk_d2 <- as.numeric(data$days_vax2 + 1)
  data$riskend_d2 <- pmin((data$days_vax2 + vax2_end), data$study_exit_days, na.rm = T)
  data$risk_d2_length <- (data$riskend_d2 - data$days_vax2)

  ## create a numeric ID variable 
  data$numeric_id <- match(data$person_id, unique(data$person_id))
  
  ## throw a warning if duplicate PIDs 
  
  if(any(duplicated(data$numeric_id))) { 
    warning("there are unexpected duplicate PIDs in the data, check with DAP")
    message("duplicate PIDs dropped to enable function to run")
    data <- data %>% 
      distinct(numeric_id, .keep_all = TRUE) # should really be dropped by DAP 
  }
  
  # DESIGN-SPECIFIC VARIABLES RISK AND CONTROL OPTIONS
  
  if(design == "sccs") {
    
    # length of study and control period 
    data$control_length <- (data$study_exit_days - data$risk_d1_length - data$risk_d2_length - abs(preexp_per_start))
    data$study_length <- data$study_exit_days
    
    # start at the start of study, end at censoring 
    data$ref_start <- as.numeric(data$study_entry_days)
    data$ref_end <- as.numeric(data$study_exit_days)
    
    # when should the time between doses end?
    # for sccs, end at the start (if zero, ignored in model estimation) 
    data$between_end <- data$riskend_d1 + 1 
    
  } else if(design == "scri") {
    
    # control start and end periods 
    # for a prevaccine SCRI, start will be control start plus a negative duration, for a post vacc it'll be control start (preexp will be 0 as does not apply in that)
    data$c_start <- pmin((data[[control_start]] + preexp_per_start), (data[[control_start]] + preexp_per_start + control_duration), na.rm = T)
    data$c_end <- pmax((data[[control_start]] + preexp_per_start), (data[[control_start]] + preexp_per_start + control_duration), na.rm = T)
    data$c_end <- pmin(data$study_exit_days, data$c_end, na.rm = T)
    
    # length of study and control period 
    data$control_length <- (data$c_end - data$c_start)
    data$study_length <- pmin(data$study_exit_days, data$c_end, (data$days_vax1 + vax1_end), (data$days_vax2 + vax2_end), na.rm = T)
    
    # start at the first of control window or risk window (whichever is first)
    data$ref_start <- pmin(data$c_start, data$risk_d1, na.rm = T) 
    data$ref_start <- as.numeric(data$ref_start)
    
    # end at the last of the end of the control window, first and second risk windows (whichever is first)
    data$ref_end <- pmax(data$c_end, data$riskend_d1, data$riskend_d2,  na.rm = T)
    data$ref_end <- as.numeric(data$ref_end)
    
    # time between doses should be taken out of the model 
    data$between_end <- pmin(data$risk_d2, data$study_exit_days)
    
  }
  
  # create a dummy variable for the outcome 
  data$outcome_binary <- as.numeric(!is.na(data$outcome_date))
  
  ## OUTCOME RELATED QC CHECKS 
  # check that event occurs after the start of the observation period 
  if(any((as.numeric(data$outcome_days))<0))
    warning(paste("There are ", sum((as.numeric(data$outcome_days)<0), na.rm=T), " persons with event date before 'study_entry_date' ."))
  
  # if set to true, restrict to only exposed cases (makes data management quicker)
  if(reduce_dimensions == "T"){  
    
    message('Reduce dimensions is set to true, unvaccinated non-cases will be dropped before the flowchart is applied')    
    message(paste(dim(data)[1], 'rows in input data'))
    data <- data[!is.na(data$outcome_date), ]
    message(paste(dim(data)[1], 'rows with non-missing outcome retained'))
    data <- data[!is.na(data$date_vax1), ]
    message(paste(dim(data)[1], 'rows with non-missing vaccination status retained'))
    
  }
  
  # APPLY INCLUSION AND EXCLUSION CRITERIA
  # this section creates a series of logical conditions, and then applies them in turn 
  # defining it like this as opposed to filtering allows you to print out the reduction in the data at each stage later on
  # note that if dimensions are already reduced, this will not reflect the initial selection
  
  # conditions applied to any SCCS
  data$cond_gender        <-  !is.na(data$sex)                  # non-missing gender
  data$cond_age           <-  !is.na(data$age_at_study_entry)   # non-missing age
  data$cond_vax1          <-  !is.na(data$date_vax1)            # non-missing vax1 
  data$cond_type_vax1     <-  !is.na(data$type_vax1)            # non-missing type_vax1
  data$cond_outcome       <-  !is.na(data$outcome_date)           # non-missing myoperi
  data$cond_vax1_outcome  <-  data$cond_vax1 & data$cond_outcome  # non=missing vax1 and myoperi
  data$cond_vax2          <-  !is.na(data$date_vax2)            # non-missing vax2 (not always applied)
  
  # time related conditions 
  # condition for outcome after 01-09-2020:
  data$cond_outcome_from_sept2020 <- data$cond_outcome & (data$outcome_date >= as.Date("2020-09-01"))
  # condition for outcome before censoring 
  data$cond_outcome_before_censor  <-  data$cond_outcome & (data$outcome_days < data$study_exit_days)
  
  if(design == "scri"){  
    
    # three conditions for scri: within control, within risk window 1 or within risk window 2 
    # condition for outcome during the control interval 
    
    data$cond_control_window <- data$cond_vax1_outcome & 
      (data$outcome_days >= (data$c_start)) & 
      (data$outcome_days <= (data$c_end))
    
    # condition for outcome during the first risk interval 
    data$cond_risk_window1 <- data$cond_vax1_outcome  &
      (data$outcome_days > data$days_vax1) & 
      (data$outcome_days <= (data$days_vax1 + risk1))
    
    # condition for outcome during the second risk interval 
    data$cond_risk_window2 <- (data$cond_vax1_outcome) & data$cond_vax2 & 
      (data$outcome_days > data$days_vax2) & 
      (data$outcome_days <= (data$days_vax2 + risk2))
    
    # outcome during either of the three above 
    data$cond_scri_events <- data$cond_control_window | data$cond_risk_window1 | data$cond_risk_window2
    
    # binary indicator for the chosen design 
    data$study <- "scri"
    
  } else {
    
    data$study <- "sccs"
    
  }
  
  # print an update 
  message("completed sccs_data_management")

  # return the edited data frame 
  return(data)

} 

#' Run SCCS/SCRI flowchart
#' 
#' @description 
#' Generate a flowchart for the specified design. 
#' Note not type sensitive, so need to invoke by vaccine type. 
#' 
#' @param data = Input dataframe, should be loaded.
#' @returns Dataframe containing results. 
#' 

sccs_flowchart <- function(data) {
  
  # print an update 
  message("starting sccs_flowchart")
  
  if("sccs" %in% data$study) {
    
    # apply selection, save each step and print as flowchart 
    
    df <- NULL 
    
    df <- cbind("total rows", nrow(data))
    
    data <- data[data$cond_gender,  ]
    df <- rbind(df,(cbind("rows with non-missing gender", nrow(data))))  
    
    data <- data[data$cond_age,  ]
    df <- rbind(df,(cbind("rows with non-missing age", nrow(data))))  
    
    data <- data[data$cond_vax1,  ]
    df <- rbind(df,(cbind("rows with non-missing first dose", nrow(data))))  
    
    data <- data[data$cond_type_vax1,  ]
    df <- rbind(df,(cbind("rows with non-missing type of first dose", nrow(data))))  
    
    data <- data[data$cond_outcome,  ]
    df <- rbind(df,(cbind("rows with non-missing outcome", nrow(data))))  
    
    data <- data[data$cond_outcome_from_sept2020,  ]
    df <- rbind(df,(cbind("rows with outcome after the study start date", nrow(data))))  
    
    data <- data[data$cond_outcome_before_censor,  ]
    df <- rbind(df,(cbind("rows with outcome before the censor date", nrow(data))))  
    
    # name columns and return  
    df <- as.data.frame(df)
    names(df) <- c("Description", "Number")
    
  } else if("scri" %in% data$study){
    
    # apply selection, save each step and print as flowchart 
    
    df <- NULL 
    
    df <- cbind("total rows", nrow(data))
    
    data <- data[data$cond_gender,  ]
    df <- rbind(df,(cbind("rows with non-missing gender", nrow(data))))  
    
    data <- data[data$cond_age,  ]
    df <- rbind(df,(cbind("rows with non-missing age", nrow(data))))  
    
    data <- data[data$cond_vax1,  ]
    df <- rbind(df,(cbind("rows with non-missing first dose", nrow(data))))  
    
    data <- data[data$cond_type_vax1,  ]
    df <- rbind(df,(cbind("rows with non-missing type of first dose", nrow(data))))  
    
    data <- data[data$cond_outcome,  ]
    df <- rbind(df,(cbind("rows with non-missing outcome", nrow(data))))  
    
    data <- data[data$cond_outcome_from_sept2020,  ]
    df <- rbind(df,(cbind("rows with outcome after the study start date", nrow(data))))  
    
    data <- data[data$cond_outcome_before_censor,  ]
    df <- rbind(df,(cbind("rows with outcome before the censor date", nrow(data))))  
    
    data <- data[data$cond_scri_events,  ]
    df <- rbind(df,(cbind("rows with outcome in control, risk 1 or risk 2 intervals", nrow(data))))  
    
    # name columns and return  
    df <- as.data.frame(df)
    names(df) <- c("Description", "Number")
    
  } 

  # print an update 
  message("completed sccs_flowchart")
  
  return(df)
  
}   

#' Apply SCCS/SCRI selection criteria 
#' 
#' @description
#' Small function to apply the selection criteria per specified design. 
#' Note, not type sensitive so will generate a dataset with all vaccine types (less NA). 
#' 
#' @param data = Input dataframe, should be loaded. 
#' @returns Dataframe with inclusion/exclusion applied. 
#' 

sccs_selection <- function(data) {
  
  # print an update 
 message("starting sccs_selection")
  
 if("sccs" %in% data$study){
    
    # apply selection
    data <- data[data$cond_gender,  ]
    data <- data[data$cond_age,  ]
    data <- data[data$cond_vax1,  ]
    data <- data[data$cond_type_vax1,  ]
    data <- data[data$cond_outcome,  ]
    data <- data[data$cond_outcome_from_sept2020,  ]
    data <- data[data$cond_outcome_before_censor,  ]
    
 } else if("scri" %in% data$study){
    
    # apply selection
    data <- data[data$cond_gender,  ]
    data <- data[data$cond_age,  ]
    data <- data[data$cond_vax1,  ]
    data <- data[data$cond_type_vax1,  ]
    data <- data[data$cond_outcome,  ]
    data <- data[data$cond_outcome_from_sept2020,  ]
    data <- data[data$cond_outcome_before_censor,  ]
    data <- data[data$cond_scri_events,  ]
    
  } 
  
  # print an update 
  message("completed sccs_selection")
  
  return(data)
  
}   

#' SCCS baseline table 
#' 
#' @description
#' Summarises number of events in and duration of control and risk periods. 
#' This is important to check, because risk periods and control periods may be truncated by later events. 
#'
#'@param data Input dataframe, should be loaded.
#'@returns Dataframe with descriptives. 
#' 

sccs_table <- function(data) {
  
  # print an update 
  message("starting sccs_table")
  
  # check on input 
  if(nrow(data) <= 0) {
    stop("no observations in input data, no output will generate")
  } 
  
  # Number of Events
  # create variable for end of the "inbetween" period 
  data$between_dose_end <- pmin(data$days_vax2, data$study_exit_days, na.rm = T)
  
  # number of events among those with at least one vaccine dose 
  n_total <- sum(as.numeric(!is.na(data$outcome_days) & (!is.na(data$days_vax1))))
  
  # number of events in risk windows 
  n_risk1 <- sum(as.numeric((data$days_vax1 < data$outcome_days) & (data$outcome_days <= data$riskend_d1)), na.rm = T) 
  n_risk2 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$riskend_d2)), na.rm = T) 
  
  # number of events in "control" period 
 if("sccs" %in% data$study){
    
    # study days starts at zero 
    # double check how this handles missing second doses 
    n_control1 <- sum(as.numeric(data$outcome_days < data$preexp_days), na.rm = T) 
    n_control2 <- sum(as.numeric((data$riskend_d1 < data$outcome_days) & (data$outcome_days < data$between_dose_end)), na.rm = T) 
    n_control3 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <=data$study_exit_days)), na.rm = T) 
    
    n_control <- sum(n_control1, n_control2, n_control3)
    
 } else if("scri" %in% data$study){
    
    n_control <- sum(as.numeric((data$c_start < data$outcome_days) & (data$outcome_days < data$c_end)), na.rm = T) 
    
  }
  
  # Duration
  
  med_control <- as.numeric(median(data$control_length, na.rm = T))
  med_risk1 <- as.numeric(median(data$risk_d1_length, na.rm = T))
  med_risk2 <- as.numeric(median(data$risk_d2_length, na.rm = T))
  med_total <- as.numeric(median(data$study_length, na.rm = T))
  
  min_control <- as.numeric(min(data$control_length, na.rm = T))
  min_risk1 <- as.numeric(min(data$risk_d1_length, na.rm = T))
  min_risk2 <- as.numeric(min(data$risk_d2_length, na.rm = T))
  min_total <- as.numeric(min(data$study_length, na.rm = T))
  
  max_control <- as.numeric(max(data$control_length, na.rm = T))
  max_risk1 <- as.numeric(max(data$risk_d1_length, na.rm = T))
  max_risk2 <- as.numeric(max(data$risk_d2_length, na.rm = T))
  max_total <- as.numeric(max(data$study_length, na.rm = T))
  
  # Create Table 
  df_control <- cbind("control", 
                      n_control, 
                      paste0(med_control, " (", min_control, " -", max_control, ")")
  )
  
  df_risk1 <- cbind("risk 1", 
                    n_risk1, 
                    paste0(med_risk1, " (", min_risk1, " -", max_risk1, ")")
  )
  
  df_risk2 <- cbind("risk 2", 
                    n_risk2, 
                    paste0(med_risk2, " (", min_risk2, " -", max_risk2, ")")
  )
  
  df_total<- cbind("total", 
                   n_total, 
                   paste0(med_total, " (", min_total, " -", max_total, ")")
  )
  
  df_all <- as.data.frame(rbind(df_control, df_risk1, df_risk2, df_total))
  names(df_all) <- c("Time Period", "N Events", "Duration, Days (Q2, min - max")
  
  # print an update 
  message("completed sccs_table")
  
  return(df_all)
  
}   

#' Run SCCS/SCRI analyses
#' 
#' @description
#' Function to run sccs or scri using the sccs package.
#' Requires SCCS and tidyverse packages (and will load).
#' 
#'@param data Input dataset. 
#'@returns Dataframe with formatted results from specified design. 
#' 

sccs_analysis <- function(data) {
  
  # print an update 
  message("starting sccs_analysis")
  
  # require packages 
  if(require("SCCS")) {
    message("SCCS is loaded")
  } else {
    message("SCCS is not installed, attempting to install") 
    install.packages("SCCS") 
    if(require("SCCS")) {
      message("SCCS is installed and loaded")
    } else {
      stop("could not install SCCS, troubleshoot")
    }
  } 
  
  if(require("tidyverse")) {
    message("tidyverse is loaded")
  } else {
    message("tidyverse is not installed, attempting to install") 
    install.packages("tidyverse") 
    if(require("tidyverse")) {
      message("tidyverse is installed and loaded")
    } else {
      stop("could not install tidyverse, troubleshoot")
    }
  } 
  
  # check on data
  if(nrow(data) <= 0) {
    stop("no observations in input data, no output will generate")
  } 
  
  ## extract design 
  design <- data$study[1]
  
  ## create a calendar time variable for calendar time adjustment (30-day interval between start and end)
  max_day = max(data$ref_end, na.rm = T)
  month_cutoffs <- seq(from = 1, to = max_day-1, by = 30)
  
  ## FIT SCCS OR SCRI 
  ### the sccs package does not take strings as input arguments
  ### the variables going into the exposure variable for sccs vs. scri is different, so unfortunately need separate function calls 
  
  if("sccs" %in% data$study) {
    
    # UNADJUSTED 
    
    tryCatch({
      
      message(paste("Fitting unadjusted", design)) 
      
      output.1 <- standardsccs(
        formula = outcome_days ~ risk_d1, 
        indiv  = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind((risk_d1 - preexp_dur), risk_d1, risk_d2), 
        aedrug = cbind(preexp_end, risk_d1_end, risk_d2_end),
        aevent = outcome_days, 
        dataformat = "multi",
        sameexpopar = F, 
        data = data
      )
      
      message("unadjusted model output")
      print(output.1)
      
    }, 
    
      error = function(e) {
        message("Unadjusted model did not fit, this is the error")
        print(e)
      }, 
    
      warning = function(w) {
        message("Unadjusted model did not fit, this is the warning")
        print(w)
      },
    
      finally = {
        
        # check if output exists and create an empty row otherwise 
        if(!exists("output.1")) {
          summary.1 <- data.frame(study_design = design, 
                                  analysis = "unadjusted", 
                                  var = "did not fit", 
                                  irr = "-", 
                                  lci = "-", 
                                  uci = "-")
          
          message("empty summary output table for unadjusted analyses created")
        }
      }
    ) 
  
    # ADJUSTED
  
    tryCatch({
      
      message(paste("fitting adjusted", design)) 
      
      output.2 <- standardsccs(
        formula = outcome_days ~ risk_d1 + age, 
        indiv = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind((risk_d1 - preexp_dur), risk_d1, risk_d2), 
        aedrug = cbind(preexp_end, risk_d1_end, risk_d2_end),
        aevent = outcome_days, 
        agegrp = month_cutoffs, 
        dataformat = "multi",
        sameexpopar = F,  
        data = data
      )
      
      message("adjusted model output")
      print(output.2)

    }, 
    
    error = function(e) {
      message("Adjusted model did not fit, this is the error")
      print(e)
    }, 
    
    warning = function(w) {
      message("Adjusted model did not fit, this is the warning")
      print(w)
    }, 
    
    finally = {
      
      # check if output exists and create an empty row otherwise 
      if(!exists("output.2")) {
        summary.2 <- data.frame(study_design = design, 
                                analysis = "adjusted", 
                                var = "did not fit", 
                                irr = "-", 
                                lci = "-", 
                                uci = "-")
        
        message("empty summary output table for adjusted analyses created")
      }
    
    })
    
  } else if("scri" %in% data$study) {
    
    # UNADJUSTED 
    
    tryCatch({
      
      message(paste("Fitting unadjusted", design)) 
      
      output.1 <- standardsccs(
        formula = outcome_days ~ risk_d1, 
        indiv  = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind((risk_d1 - preexp_dur), risk_d1, (risk_d1_end + 1), risk_d2), 
        aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
        aevent = outcome_days, 
        dataformat = "multi",
        sameexpopar = F, 
        data = data
      )
      
      message("unadjusted model output")
      print(output.1)
      
    }, 
    
    error = function(e) {
      message("Unadjusted model did not fit, this is the error")
      print(e)
    }, 
    
    warning = function(w) {
      message("Unadjusted model did not fit, this is the warning")
      print(w)
    },
    
    finally = {
      
      # check if output exists and create an empty row otherwise 
      if(!exists("output.1")) {
        summary.1 <- data.frame(study_design = design, 
                                analysis = "unadjusted", 
                                var = "did not fit", 
                                irr = "-", 
                                lci = "-", 
                                uci = "-")
        
        message("empty summary output table for unadjusted analyses created")
      }
    }
    ) 
    
    # ADJUSTED
    
    tryCatch({
      
      message(paste("fitting adjusted", design)) 
      
      output.2 <- standardsccs(
        formula = outcome_days ~ risk_d1 + age, 
        indiv = numeric_id, 
        astart = ref_start,  
        aend = ref_end, 
        adrug = cbind((risk_d1 - preexp_dur), risk_d1, (risk_d1_end + 1), risk_d2), 
        aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
        aevent = outcome_days, 
        agegrp = month_cutoffs, 
        dataformat = "multi",
        sameexpopar = F,  
        data = data
      )
      
      message("adjusted model output")
      print(output.2)
      
    }, 
    
    error = function(e) {
      message("Adjusted model did not fit, this is the error")
      print(e)
    }, 
    
    warning = function(w) {
      message("Adjusted model did not fit, this is the warning")
      print(w)
    }, 
    
    finally = {
      
      # check if output exists and create an empty row otherwise 
      if(!exists("output.2")) {
        summary.2 <- data.frame(study_design = design, 
                                analysis = "adjusted", 
                                var = "did not fit", 
                                irr = "-", 
                                lci = "-", 
                                uci = "-")
        
        message("empty summary output table for adjusted analyses created")
      }
      
    })
    
    
  }
  
  # OUTPUT FORMATTING
  
  # unadjusted output formatting
  
  if(!exists("summary.1")) {
    
    summary.1 <- as.data.frame(cbind(output.1$conf.int, output.1$coefficients[,c(1,3)]))
    summary.1$var <- rownames(summary.1)
    summary.1$analysis <- "unadjusted"
    
    rownames(summary.1) <- NULL
    
    summary.1 <- summary.1 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, var, irr, lci, uci)  
    
    message("summary output table for unadjusted analyses created")
    
  } 
  
  # adjusted output 
  
  if(!exists("summary.2")) {
  
    summary.2 <- as.data.frame(cbind(output.2$conf.int, output.2$coefficients[,c(1,3)]))
    summary.2$var <- rownames(summary.2)
    summary.2$analysis <- "adjusted"
    rownames(summary.2) <- NULL
    
    summary.2 <- summary.2 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, var, irr, lci, uci) 
  
    message("summary output table for adjusted analyses created")
  
  }
  
  # bind the two outputs 
  output.total <- rbind(summary.1, summary.2)
    
  # add design specific labels 
  if("sccs" %in% data$study) {
    
    output.total <- output.total %>%  
      mutate(var = case_when(var == "risk_d11" ~ "pre-exposure pindow", 
                             var == "risk_d12" ~ "dose 1 risk window", 
                             var == "risk_d13" ~ "dose 2 risk window", 
                             TRUE ~ var)) %>% 
      slice_head(n = 3)
  
  } else if("scri" %in% data$study) {
    
    output.total <- output.total %>%  
      mutate(var = case_when(var == "risk_d11" ~ "pre-exposure pindow", 
                         var == "risk_d12" ~ "dose 1 risk window", 
                         var == "risk_d13" ~ "between doses", 
                         var == "risk_d14" ~ "dose 2 risk window", 
                         TRUE ~ var)) %>% 
    slice_head(n = 4)
  
  } 
  
  # print an update 
  message("completed sccs_analysis")
  
  # return object 
  return(output.total)

} 

