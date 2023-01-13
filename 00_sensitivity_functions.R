# Program Information  ----------------------------------------------------
#
# Program:      00_sensitivity_functions
# Description:  load packages required to support pipeline 
#               run custom functions required 
#               *utility functions (to plot basic summaries)
#               *event_dependency_table
#               *event_dependency_plot
#               *event_information_table 
#               *sccs_analysis 
#
# Dependencies: Function 5 sccs_analysis requires and loads packages tidyverse, SCCS
#
# 0. HOUSEKEEPING ------------------------------------------------------------

# load packages required by the pipeline 
if(require("tidyverse")) {
  message("tidyverse is loaded")
} else {
  message("tidyverse is not installed, attempting to install") 
  install.packages("tidyverse") 
  if(require("tidyverse")) {
    message("tidyverse is installed and loaded")
  } else {
    stop("could not install tidyverse, troubleshooting needed")
  }
} 

if(require("SCCS")) {
  message("SCCS is loaded")
} else {
  message("SCCS is not installed, attempting to install") 
  install.packages("SCCS") 
  if(require("SCCS")) {
    message("SCCS is installed and loaded")
  } else {
    stop("could not install SCCS, troubleshooting needed")
  }
} 

# 1. FUNCTIONS ---------------------------------------------------------------
# read functions required by the pipeline 

# UTILITY FUNCTIONS

# 1. Function to summarise a vector of binary variables
# Returns a tibble of counts when the variable == 1, and percentages 

indicator_summary <- function(data, ...) {
  
  {{data}} %>% 
    select(all_of(...)) %>% 
    summarise(
      across(
        .cols  = all_of(...),
        list(count = sum, percent = mean), 
        na.rm = TRUE, 
        .names = "{.fn}-{.col}")) %>% 
    mutate(across(contains("percent"), ~round(.x*100, digits = 2))) %>% 
    pivot_longer(everything(), 
                 names_to = c(".value", "variable"), 
                 names_sep = "-") 
  
}

# 2. Function to summarise a vector of continuous variables 
# Returns a tibble of median and IQR (called counts and percentages to enable binding)

cont_summary <- function(data, ...) {

  {{data}} %>% 
    select(all_of(...)) %>% 
    summarise(
      across(
        .cols  = all_of(...),
        .fns = list(median = median,
                    q1 = ~quantile(.,probs =0.25, na.rm = T), 
                    q3 = ~quantile(.,probs =0.75, na.rm = T), 
                    min = min, 
                    max = max), 
        na.rm = TRUE, 
        .names = "{.fn}-{.col}")) %>% 
    mutate(across(everything(), ~as.numeric(.x))) %>% 
    mutate(across(everything(), ~round(.x, digits = 0))) %>% 
    pivot_longer(everything(), 
                 names_to = c(".value", "variable"), 
                 names_sep = "-") 
  
}

# PROJECT SPECIFIC FUNCTIONS 

## TABULATE EVENT-DEPENDENCY INFORMATION 
#' @description Takes clean data and tabulates basic information of interest for assessing event-dependent exposures 
#' Should be run on the SCCS dataset (i.e., don't restrict control time)
#' @param data = Input dataframe, should be loaded into R before invoking, as function does not do this. 
#' @returns A dataframe with basic summaries of the vaccine/event relationship  

event_dependency_table <- function(data) {
  
  # list the variables of interest
  var_list <- c("vax1", "vax2", "vax3", "dose1_afterevent", "dose2_afterevent", "dose3_afterevent")
  # use a utility function from 00_functions to save the N and % of each dose 
  vaccine_summary <- indicator_summary(data, var_list)

  # list continous variables of interest
  cont_list <- c("diff_d1d2", "diff_d2d3", "diff_event_d1", "diff_event_d2", "diff_event_d3")
  # use a utility function from 00_functions to save the median and IQR for each continous variable 
  vaccine_duration_summary <- cont_summary(data, cont_list) 
  
  # Make this into a table 
  vaccine_table <- as.data.frame(bind_rows(vaccine_summary, vaccine_duration_summary))
  names(vaccine_table) <- c("Variable", "N", "%", "Median", "Q1", "Q3", "Min", "Max")
  
  vaccine_table <- vaccine_table %>% 
    mutate(across(everything(), replace_na, replace = "")) %>% 
    mutate(across(everything(), ~ifelse(.x == "Inf", "", .x))) %>%   
    mutate(across(everything(), ~ifelse(.x == "-Inf", "", .x)))
  
  return(vaccine_table) 
  
} 

## PLOT EVENT-DEPENDENCY INFORMATION 
#' @description generate exposure-centered interval plot
#' @param data = Input dataframe, should be loaded into R before invoking, as function does not do this. 
#' @param dose = Dose Number, 1, 2 etc. 
#' @returns ggplot object

event_dependency_plot <- function(data, dose) {

  xvar <- paste0("diff_event_d", dose)
  vaccine_type <- unique(data$type_vax1)
  message(paste("plotting histogram for", xvar, vaccine_type))


  plot <- ggplot(data, aes(x=eval(parse(text = xvar)))) + 
    geom_histogram(colour = "lightslategrey", fill = "lightsteelblue") + 
    theme_minimal() + 
    labs(y = "Number of people", 
         x = "Difference (days)", 
         title = paste("Exposure-Centered Interval Plot", "Dose", dose, "-", vaccine_type),  
         subtitle = "Distribution of time between the event and vaccine dose") + 
    geom_vline(xintercept=0, colour="darkslategrey") 
  
  return(plot)
  
} 

#' SCRI period description
#' @description tabulate number of events and duration of each period
#' @param data = input dataframe
#' @returns dataframe 

describe_scri_periods <- function(data) {

  # describe total number of events in each period 
  n_total <- nrow(data)
  
  # number of events in risk windows 
  n_risk1 <- sum(as.numeric((data$days_vax1 < data$outcome_days) & (data$outcome_days <= data$risk_d1_end)), na.rm = T) 
  n_risk2 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  
  # number of events in pre-exposure period 
  n_preexp <- sum(as.numeric((data$preexp_days < data$outcome_days) & (data$outcome_days <= data$preexp_end)), na.rm = T) 
  
  # number of events in "control" period 
  n_control <- sum(as.numeric((data$c_start < data$outcome_days) & (data$outcome_days < data$c_end)), na.rm = T) 
  
  # number of events in "between" period 
  n_between <- sum(as.numeric((data$between_start < data$outcome_days) & (data$outcome_days < data$between_end)), na.rm = T) 
  
  # Duration
  duration_vars <- c("study_length", "c_length", "preexp_length", "risk_d1_length", "between_length", "risk_d2_length")
  # note, this will throw errors if applied to an empty vector even though it is expected behavior. Warning can be ignored. 
  duration_summary <- cont_summary(data, duration_vars)
  duration_summary <- duration_summary[,-1]
  
  count_summary <- rbind(n_total, n_control, n_preexp, n_risk1, n_between, n_risk2) 
  row_headings <- rbind("Overall", "Control", "Pre-Exposure", "Dose 1 Risk Period", "In between doses", "Dose 2 Risk Period") 
  
  event_summary <- cbind(row_headings,count_summary,duration_summary)
  
  # return the edited data frame 
  return(event_summary)

}   

#' SCCS period description
#' @description tabulate number of events and duration of each period
#' @param data = input dataframe
#' @returns dataframe 

describe_sccs_periods <- function(data) {
  
  # describe total number of events in each period 
  n_total <- nrow(data)
  
  # number of events in risk windows 
  n_risk1 <- sum(as.numeric((data$days_vax1 < data$outcome_days) & (data$outcome_days <= data$risk_d1_end)), na.rm = T) 
  n_risk2 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  n_risk3 <- sum(as.numeric((data$days_vax2 < data$outcome_days) & (data$outcome_days <= data$risk_d2_end)), na.rm = T) 
  
  # number of events in pre-exposure period 
  n_preexp1 <- sum(as.numeric((data$preexp1_days < data$outcome_days) & (data$outcome_days <= data$preexp1_end)), na.rm = T) 
  n_preexp2 <- sum(as.numeric((data$preexp2_days < data$outcome_days) & (data$outcome_days <= data$preexp2_end)), na.rm = T) 
  n_preexp3 <- sum(as.numeric((data$preexp3_days < data$outcome_days) & (data$outcome_days <= data$preexp3_end)), na.rm = T) 
  
  # Duration
  duration_vars <- c("study_length", "preexp1_length", "risk_d1_length", "preexp2_length", "risk_d2_length", "preexp3_length", "risk_d3_length")
  duration_summary <- cont_summary(data, duration_vars)
  duration_summary <- duration_summary[,-1]
  
  count_summary <- rbind(n_total, n_preexp1, n_risk1, n_preexp2, n_risk2, n_preexp3, n_risk3) 
  row_headings <- rbind("Overall", "Pre-Exposure 1", "Dose 1 Risk Period", "Pre-Exposure 2", "Dose 2 Risk Period", "Pre-Exposure 3", "Dose 3 Risk Period") 
  
  event_summary <- cbind(row_headings,count_summary,duration_summary)
  
  # return the edited data frame 
  return(event_summary)
  
}   


#' Run SCRI
#' @description fit an SCRI using the SCCS package and format the output 
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_scri <- function(data) {
  
  ## create a calendar time variable for calendar time adjustment (30-day interval between start and end in the data where model is fit)
  max_day = max(data$ref_end, na.rm = T)
  min_day = min(data$ref_start, na.rm = T)
  month_cutoffs <- seq(from = min_day-1, to = max_day-1, by = 30)
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting SCRI for", vaccine_type))
  
  
  ## wrap the model call to capture errors and warnings
  ## this is not very advanced because neither trycatch or suppresshandlers did exactly what I wanted
  
  try({
    
      message("fitting unadjusted model")
       
        output.1 <- standardsccs(
           formula = outcome_days ~ risk_d1, 
           indiv  = numeric_id, 
           astart = ref_start,  
           aend = ref_end, 
           adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
           aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
           aevent = outcome_days, 
           dataformat = "multi",
           sameexpopar = F, 
           data = data 
       )
      
      message("unadjusted model output")
      print(output.1)
              
    })  

  
  if(!exists("output.1")) {
    
    summary.1 <- data.frame(study_design = design, 
                          analysis = "unadjusted", 
                          n = "n/a", 
                          var = "did not fit", 
                          irr = "-", 
                          lci = "-", 
                          uci = "-")
  
    message("empty summary output table for unadjusted analyses created")  
  
  } 
  
  
  # ADJUSTED
  
  try({
    
    message(paste("fitting adjusted model")) 
    
    output.2 <- standardsccs(
      formula = outcome_days ~ risk_d1 + age, 
      indiv = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind((risk_d1 - preexp_length), risk_d1, between_start, risk_d2), 
      aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
      aevent = outcome_days, 
      agegrp = month_cutoffs, 
      dataformat = "multi",
      sameexpopar = F,  
      data = data
    )
    
    message("adjusted model output")
    print(output.2)

  }) 
  
  if(!exists("output.2")) {
    
    summary.2 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not fit", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for adjusted analyses created")  
    
  } 
  
  if(!exists("summary.1")) {
    
    summary.1 <- as.data.frame(cbind(output.1$n, output.1$nevent, output.1$conf.int, output.1$coefficients[,c(1,3)]))
    summary.1$var <- rownames(summary.1)
    summary.1$analysis <- "unadjusted"
    
    rownames(summary.1) <- NULL
    
    summary.1 <- summary.1 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci)  

  } 
  
  # adjusted output 
  
  if(!exists("summary.2")) {
    
    summary.2 <- as.data.frame(cbind(output.2$n, output.2$nevent, output.2$conf.int, output.2$coefficients[,c(1,3)]))
    summary.2$var <- rownames(summary.2)
    summary.2$analysis <- "adjusted"
    rownames(summary.2) <- NULL
    
    summary.2 <- summary.2 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci) 
    
  }
  
  # bind the two outputs 
  output.total <- rbind(summary.1, summary.2)
  
  # add labels 
  output.total <- output.total %>%  
    mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                           var == "risk_d12" ~ "dose 1 risk window", 
                           var == "risk_d13" ~ "between doses", 
                           var == "risk_d14" ~ "dose 2 risk window", 
                           TRUE ~ var)) 

  # return object 
  return(output.total)
  
}


#' Run SCCS
#' @description fit a standard SCCS
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_standard_sccs <- function(data) {
  
  ## create a calendar time variable for calendar time adjustment (30-day interval between start and end in the data where model is fit)
  max_day = max(data$ref_end, na.rm = T)
  min_day = min(data$ref_start, na.rm = T)
  month_cutoffs <- seq(from = min_day-1, to = max_day-1, by = 30)
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting standard SCCS for", vaccine_type))
  
  ## wrap the model call to capture errors and warnings
  ## this is not very advanced because neither trycatch or suppresshandlers did exactly what I wanted
  
  try({
    
    message("fitting unadjusted model")
    
    output.1 <- standardsccs(
      formula = outcome_days ~ risk_d1, 
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = data 
    )
    
    message("unadjusted model output")
    print(output.1)
    
  })  
  
  
  if(!exists("output.1")) {
    
    summary.1 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not fit", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for unadjusted analyses created")  
    
  } 
  
  
  # ADJUSTED
  
  try({
    
    message(paste("fitting adjusted model")) 
    
    output.2 <- standardsccs(
      formula = outcome_days ~ risk_d1 + age, 
      indiv = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      agegrp = month_cutoffs, 
      dataformat = "multi",
      sameexpopar = F,  
      data = data
    )
    
    message("adjusted model output")
    print(output.2)
    
  }) 
  
  if(!exists("output.2")) {
    
    summary.2 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not fit", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for adjusted analyses created")  
    
  } 
  
  if(!exists("summary.1")) {
    
    summary.1 <- as.data.frame(cbind(output.1$n, output.1$nevent, output.1$conf.int, output.1$coefficients[,c(1,3)]))
    summary.1$var <- rownames(summary.1)
    summary.1$analysis <- "unadjusted"
    
    rownames(summary.1) <- NULL
    
    summary.1 <- summary.1 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci)  
    
  } 
  
  # adjusted output 
  
  if(!exists("summary.2")) {
    
    summary.2 <- as.data.frame(cbind(output.2$n, output.2$nevent, output.2$conf.int, output.2$coefficients[,c(1,3)]))
    summary.2$var <- rownames(summary.2)
    summary.2$analysis <- "adjusted"
    rownames(summary.2) <- NULL
    
    summary.2 <- summary.2 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci) 
    
  }
  
  # bind the two outputs 
  output.total <- rbind(summary.1, summary.2)
  
  # add labels 
  output.total <- output.total %>%  
    mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                           var == "risk_d12" ~ "dose 1 risk window", 
                           var == "risk_d13" ~ "dose 2 pre-exposure", 
                           var == "risk_d14" ~ "dose 2 risk window", 
                           var == "risk_d15" ~ "dose 3 pre-exposure", 
                           var == "risk_d16" ~ "dose 3 risk window", 
                           TRUE ~ var)) 
  
  # return object 
  return(output.total)
  
}



#' Run extended SCCS
#' @description fit an extended SCCS
#' requires input parameters defining start and end of study period
#' these are created outside the function for clarity 
#' @param data = input dataframe
#' @returns dataframe 

fit_extended_sccs <- function(data) {
  
  ## create a calendar time variable for calendar time adjustment (30-day interval between start and end in the data where model is fit)
  max_day = max(data$ref_end, na.rm = T)
  min_day = min(data$ref_start, na.rm = T)
  month_cutoffs <- seq(from = min_day-1, to = max_day-1, by = 30)
  
  ## print info about the call
  vaccine_type <- unique(data$type_vax1)
  message(paste("fitting standard SCCS for", vaccine_type))
  
  ## wrap the model call to capture errors and warnings
  ## this is not very advanced because neither trycatch or suppresshandlers did exactly what I wanted
  
  try({
    
    message("fitting unadjusted model")
    
    output.1 <- eventdepenexp(
      indiv  = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      dataformat = "multi",
      sameexpopar = F, 
      data = data 
    )
    
    message("unadjusted model output")
    print(output.1)
    
  })  
  
  
  if(!exists("output.1")) {
    
    summary.1 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not generate output", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for unadjusted analyses created")  
    
  } 
  
  
  # ADJUSTED
  
  try({
    
    message(paste("fitting adjusted model")) 
    
    output.2 <- eventdepenexp(
      indiv = numeric_id, 
      astart = ref_start,  
      aend = ref_end, 
      adrug = cbind(risk_d1 - 30, risk_d1, preexp2, risk_d2, preexp3, risk_d3), 
      aedrug = cbind(preexp1_end, risk_d1_end, preexp2_end, risk_d2_end, preexp3_end, risk_d3_end),
      aevent = outcome_days, 
      agegrp = month_cutoffs, 
      dataformat = "multi",
      sameexpopar = F,  
      data = data
    )
    
    message("adjusted model output")
    print(output.2)
    
  }) 
  
  if(!exists("output.2")) {
    
    summary.2 <- data.frame(study_design = design, 
                            analysis = "unadjusted", 
                            n = "n/a", 
                            var = "did not generate output", 
                            irr = "-", 
                            lci = "-", 
                            uci = "-")
    
    message("empty summary output table for adjusted analyses created")  
    
  } 
  
  if(!exists("summary.1")) {
    
    summary.1 <- as.data.frame(cbind(output.1$n, output.1$nevent, output.1$conf.int, output.1$coefficients[,c(1,3)]))
    summary.1$var <- rownames(summary.1)
    summary.1$analysis <- "unadjusted"
    
    rownames(summary.1) <- NULL
    
    summary.1 <- summary.1 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci)  
    
  } 
  
  # adjusted output 
  
  if(!exists("summary.2")) {
    
    summary.2 <- as.data.frame(cbind(output.2$n, output.2$nevent, output.2$conf.int, output.2$coefficients[,c(1,3)]))
    summary.2$var <- rownames(summary.2)
    summary.2$analysis <- "adjusted"
    rownames(summary.2) <- NULL
    
    summary.2 <- summary.2 %>%
      rename(irr = `exp(coef)`,
             lci = `lower .95`,
             uci = `upper .95`, 
             n = `V2`) %>%
      mutate(study_design = design) %>% 
      select(study_design, analysis, n, var, irr, lci, uci) 
    
  }
  
  # bind the two outputs 
  output.total <- rbind(summary.1, summary.2)
  
  # add labels 
  output.total <- output.total %>%  
    mutate(var = case_when(var == "risk_d11" ~ "dose 1 pre-exposure", 
                           var == "risk_d12" ~ "dose 1 risk window", 
                           var == "risk_d13" ~ "dose 2 pre-exposure", 
                           var == "risk_d14" ~ "dose 2 risk window", 
                           var == "risk_d15" ~ "dose 3 pre-exposure", 
                           var == "risk_d16" ~ "dose 3 risk window", 
                           TRUE ~ var)) 
  
  # return object 
  return(output.total)
  
}
