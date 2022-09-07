# Program Information  ----------------------------------------------------
#
# Program:      step_12_6_sccs_sens_run.R
# Description:  run sccs under different options for control periods 
# Dependencies: 00_sensitivity_functions 
#               SCCS package 
# 
# 0. HOUSEKEEPING ---------------------------------------------------------

# if file is opened through project directory then wd doesn't need to be automatically set 
# current script assumes wd is the project directory 
if(!any(ls()=="thisdir"))   thisdir   <- paste0(getwd(),"/") 
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"g_output/")
 
# ensure required folders are created  
dir.create(file.path(paste0(dirtemp, "sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,"log_files/sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)

# load custom functions into working memory
source(paste0(thisdir,"00_sensitivity_functions.R"))

# 1. READ IN DATA ------------------------------------------------------------

load(paste0(dirtemp, "sccs_sensitivity/", "sccs_data_extract", ".RData"))

# 2.ANALYSIS -------------------------------------------------------------
# these calls functions for each design option on each strata (i.e vaccine) of interest
# could likely be further simplified as quite a few lines of code for each analysis at the moment 

# 2.1. BASE CASE: SCRI with 90 day pre-vaccination period ----------------

  design <- "scri_pre"
  dir.create(file.path(paste0(diroutput, "sccs_sensitivity/", design)), showWarnings = FALSE, recursive = TRUE)
  
  # STEP 1. APPLY DATA MANAGEMENT
  # this is where all the inputs are set 
  scri_pre <- sccs_data_management(data = sccs_data_extract, 
                                            outcome = "myopericarditis_date", 
                                            reduce_dimensions = "F", 
                                            design = "scri", 
                                            control_start = "days_vax1", 
                                            control_dur = -60,
                                            preexp = -30, 
                                            risk1 = 28,
                                            risk2 = 28)
  
  # STEP 2. PRINT A FLOWCHART
  # by vaccine type 
  # note, don't specify a data input, lapply automatically supplies the first argument as the split)
  table1.interim <- lapply(split(scri_pre, scri_pre$type_vax1), sccs_flowchart)
  table1.split <- do.call(rbind,table1.interim)
  
  # minor tidying of the df before printing out 
  table1.split <- cbind(Names = rownames(table1.split), table1.split)
  table1.split$Names <- substring(table1.split$Names,1,nchar(table1.split$Names)-2)
  
  # apply to the overall dataset for completion 
  table1.all <- sccs_flowchart(data = scri_pre)
  
  # export these two flowcharts
  write.csv(table1.split, paste0(diroutput, "sccs_sensitivity/", design, "/table_1a.csv"), row.names = FALSE)
  write.csv(table1.all, paste0(diroutput, "sccs_sensitivity/", design, "/table_1b.csv"), row.names = FALSE)
  
  # STEP 3. APPLY SELECTION 
  # possibly doesn't need to be run because analysis will automatically drop ineligible cases 
  # should be noted that retaining the entire population will affect age adjustment and table 1 summary stats (not all people will contribute to estimation)

  sccs_analytical_file <- sccs_selection(data = scri_pre)

  # STEP 4. BASELINE TABLE 
  # note, if you haven't run the selection the totals in this table will not sum exactly 
  
  if(nrow(sccs_analytical_file) <= 0) {
    message(paste0("There are NO observations in input data using ", design, " EMPTY OUTPUT CREATED")) 
  } 
  
  table2.interim <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_table)
  table2 <- do.call(rbind, table2.interim )
  
  table2 <- cbind(Names = rownames(table2), table2)
  table2$Names <- substring(table2$Names,1,nchar(table2$Names)-2)

  write.csv(table2, paste0(diroutput,"sccs_sensitivity/", design, "/table_2.csv"), row.names = FALSE)

  # STEP 5. ANALYSIS 
  
  table3.interim  <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_analysis)
  table3.split <- do.call(rbind,table3.interim)
  
  table3.split <- cbind(Names = rownames(table3.split), table3.split)
  table3.split$Names <- substring(table3.split$Names,1,nchar(table3.split$Names)-2)
  
  write.csv(table3.split, paste0(diroutput,"sccs_sensitivity/", design, "/table_3.csv"), row.names = FALSE)
  
# 2.2. ALT 1: SCRI with 90 day post-vaccination period ----------------
  
  design <- "scri_post"
  dir.create(file.path(paste0(diroutput, "sccs_sensitivity/", design)), showWarnings = FALSE, recursive = TRUE)
  
  # STEP 1. APPLY DATA MANAGEMENT
  scri_post <- sccs_data_management(data = sccs_data_extract, 
                                    outcome = "myopericarditis_date", 
                                    reduce_dimensions = "F", 
                                    design = "scri", 
                                    control_start = "days_last_vax", 
                                    control_dur = 60,
                                    risk1 = 28,
                                    risk2 = 28) 
  
  
  # STEP 2. PRINT A FLOWCHART
  table1.interim <- lapply(split(scri_pre, scri_pre$type_vax1), sccs_flowchart)
  table1.split <- do.call(rbind,table1.interim)
  
  # minor tidying of the df before printing out 
  table1.split <- cbind(Names = rownames(table1.split), table1.split)
  table1.split$Names <- substring(table1.split$Names,1,nchar(table1.split$Names)-2)
  
  # apply to the overall dataset for completion 
  table1.all <- sccs_flowchart(data = scri_pre)
  
  # export these two flowcharts
  write.csv(table1.split, paste0(diroutput, "sccs_sensitivity/", design, "/table_1a.csv"), row.names = FALSE)
  write.csv(table1.all, paste0(diroutput, "sccs_sensitivity/", design, "/table_1b.csv"), row.names = FALSE)
  
  # STEP 3. APPLY SELECTION 
  # possibly doesn't need to be run because analysis will automatically drop ineligible cases 
  # should be noted that retaining the entire population will affect age adjustment and table 1 summary stats (not all people will contribute to estimation)

  sccs_analytical_file <- sccs_selection(data = scri_pre)
  
  # STEP 4. BASELINE TABLE 
  # note, if you haven't run the selection the totals in this table will not sum exactly 
  
  if(nrow(sccs_analytical_file) <= 0) {
    message(paste0("There are NO observations in input data using ", design, " EMPTY OUTPUT CREATED")) 
  } 
  
  table2.interim <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_table)
  table2 <- do.call(rbind, table2.interim )
  
  table2 <- cbind(Names = rownames(table2), table2)
  table2$Names <- substring(table2$Names,1,nchar(table2$Names)-2)
  
  write.csv(table2, paste0(diroutput,"sccs_sensitivity/", design, "/table_2.csv"), row.names = FALSE)
  
  # STEP 5. ANALYSIS 
  
  table3.interim  <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_analysis)
  table3.split <- do.call(rbind,table3.interim)
  
  table3.split <- cbind(Names = rownames(table3.split), table3.split)
  table3.split$Names <- substring(table3.split$Names,1,nchar(table3.split$Names)-2)
  
  write.csv(table3.split, paste0(diroutput,"sccs_sensitivity/", design, "/table_3.csv"), row.names = FALSE)
  

# 2.3. ALT 2: SCCS ----------------
  design <- "sccs"
  dir.create(file.path(paste0(diroutput, "sccs_sensitivity/", design)), showWarnings = FALSE, recursive = TRUE)
  
  # STEP 1. APPLY DATA MANAGEMENT
  sccs <- sccs_data_management(data = sccs_data_extract, 
                               outcome = "myopericarditis_date", 
                               reduce_dimensions = "F", 
                               design = "sccs", 
                               risk1 = 28,
                               risk2 = 28)
  
  # STEP 2. PRINT A FLOWCHART
  table1.interim <- lapply(split(scri_pre, scri_pre$type_vax1), sccs_flowchart)
  table1.split <- do.call(rbind,table1.interim)
  
  # minor tidying of the df before printing out 
  table1.split <- cbind(Names = rownames(table1.split), table1.split)
  table1.split$Names <- substring(table1.split$Names,1,nchar(table1.split$Names)-2)
  
  # apply to the overall dataset for completion 
  table1.all <- sccs_flowchart(data = scri_pre)
  
  # export these two flowcharts
  write.csv(table1.split, paste0(diroutput, "sccs_sensitivity/", design, "/table_1a.csv"), row.names = FALSE)
  write.csv(table1.all, paste0(diroutput, "sccs_sensitivity/", design, "/table_1b.csv"), row.names = FALSE)
  
  # STEP 3. APPLY SELECTION 
  # possibly doesn't need to be run because analysis will automatically drop ineligible cases 
  # should be noted that retaining the entire population will affect age adjustment and table 1 summary stats (not all people will contribute to estimation)

  sccs_analytical_file <- sccs_selection(data = scri_pre)

  # STEP 4. BASELINE TABLE 
  # note, if you haven't run the selection the totals in this table will not sum exactly 
  
  if(nrow(sccs_analytical_file) <= 0) {
    message(paste0("There are NO observations in input data using ", design, " EMPTY OUTPUT CREATED")) 
  } 
  
  table2.interim <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_table)
  table2 <- do.call(rbind, table2.interim )
  
  table2 <- cbind(Names = rownames(table2), table2)
  table2$Names <- substring(table2$Names,1,nchar(table2$Names)-2)
  
  write.csv(table2, paste0(diroutput,"sccs_sensitivity/", design, "/table_2.csv"), row.names = FALSE)
  
  # STEP 5. ANALYSIS 
  
  table3.interim  <- lapply(split(sccs_analytical_file, sccs_analytical_file$type_vax1), sccs_analysis)
  table3.split <- do.call(rbind,table3.interim)
  
  table3.split <- cbind(Names = rownames(table3.split), table3.split)
  table3.split$Names <- substring(table3.split$Names,1,nchar(table3.split$Names)-2)
  
  write.csv(table3.split, paste0(diroutput,"sccs_sensitivity/", design, "/table_3.csv"), row.names = FALSE)
