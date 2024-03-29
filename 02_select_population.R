# Program Information  ----------------------------------------------------
#
# Program:      02_select_population.R 
# Description:  restrict study population to vaccinated cases and create a flowchart
# Dependencies: 01_clean_data.R 
#               00_sensitivity_functions 
#               tidyverse, ggplot2
# 
# 0. HOUSEKEEPING ---------------------------------------------------------

# 1. READ IN DATA ------------------------------------------------------------
load(paste0(dirtemp, "sccs_sensitivity/", "sccs_data_extract", ".RData"))

# 2. CREATE CRITERIA ------------------------------------------------------
sccs_data_extract$cond_gender        <-  !is.na(sccs_data_extract$sex)                  # non-missing gender
sccs_data_extract$cond_age           <-  !is.na(sccs_data_extract$age_at_study_entry)   # non-missing age
sccs_data_extract$cond_vax1          <-  !is.na(sccs_data_extract$date_vax1)            # non-missing vax1 
sccs_data_extract$cond_type_vax1     <-  !is.na(sccs_data_extract$type_vax1)            # non-missing type_vax1
sccs_data_extract$cond_type_vax2     <-  ((sccs_data_extract$type_vax1 == sccs_data_extract$type_vax2) & sccs_data_extract$cond_vax1) | is.na(sccs_data_extract$date_vax2) # dose 2 same as dose 1, if exists 
sccs_data_extract$cond_outcome       <-  !is.na(sccs_data_extract$myopericarditis_date)           # non-missing myoperi
sccs_data_extract$cond_vax1_outcome  <-  sccs_data_extract$cond_vax1 & sccs_data_extract$cond_outcome  # non=missing vax1 and myoperi

# time related conditions 
# condition for outcome after 01-09-2020:
sccs_data_extract$cond_outcome_from_sept2020 <- sccs_data_extract$cond_outcome & (sccs_data_extract$myopericarditis_date >= as.Date("2020-09-01"))
# condition for outcome before censoring 
sccs_data_extract$cond_outcome_before_censor  <-  sccs_data_extract$cond_outcome & (sccs_data_extract$myopericarditis_date < sccs_data_extract$study_exit_date)

# 3. APPLY CRITERIA -------------------------------------------------------
# apply selection criteria and create a flowchart 
message("applying SCCS selection criteria")
message(paste("the original data consists of", nrow(sccs_data_extract), "individuals"))

df <- NULL 

df <- cbind("total rows", nrow(sccs_data_extract))

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_gender,  ]
df <- rbind(df,(cbind("rows with non-missing gender", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_age,  ]
df <- rbind(df,(cbind("rows with non-missing age", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_vax1,  ]
df <- rbind(df,(cbind("rows with non-missing first dose", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_type_vax1,  ]
df <- rbind(df,(cbind("rows with non-missing type of first dose", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_type_vax2,  ]
df <- rbind(df,(cbind("if second dose exists, it's the same as the first", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_outcome,  ]
df <- rbind(df,(cbind("rows with non-missing outcome", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_outcome_from_sept2020,  ]
df <- rbind(df,(cbind("rows with outcome after the study start date", nrow(sccs_data_extract))))  

sccs_data_extract <- sccs_data_extract[sccs_data_extract$cond_outcome_before_censor,  ]
df <- rbind(df,(cbind("rows with outcome before the censor date", nrow(sccs_data_extract))))  

message(paste("the selected population consists of", nrow(sccs_data_extract), "individuals"))

# name columns and return  
df <- as.data.frame(df)
names(df) <- c("Description", "Number")

# 4. SAVE DATA ------------------------------------------------------------
# export flowchart 
write.csv(df, paste0(diroutput, "sccs_sensitivity/", "/table_1_flowchart.csv"), row.names = FALSE)

# saving data with new name to prevent overwriting any prior files
sccs_population <- sccs_data_extract 
save(sccs_population,file = paste0(dirtemp,"sccs_sensitivity/","sccs_population.RData"))

