# Program Information  ----------------------------------------------------
#
# Program:      step_12_5_sccs_sens_clean.R 
# Description:  this script does some data cleaning, and outputs a dataset for 
#               SCRI sensitivity analyses. It uses part of the code from step 12_1 
#               but has put this into a function. it removes some output which is not needed (baseline chars), 
# Requirements: dataset in RData format called D3_study_variables_for_SCRI
#               table1.R (userwritten function to tabulate variables)
# Output:       this outputs a single dataset called [placeholder]

# 0. HOUSEKEEPING ---------------------------------------------------------

# if file is opened through project directory then wd doesn't need to be automatically set 
# current script assumes wd is the project directory 
if(!any(ls()=="thisdir"))   thisdir   <- paste0(getwd(),"/") 
if(!any(ls()=="dirtemp"))   dirtemp   <- paste0(thisdir,"g_intermediate/")
if(!any(ls()=="diroutput")) diroutput <- paste0(thisdir,"g_output/")

# specify the input data name 
raw_data_name <- "D3_study_variables_for_SCRI"

# ensure required folders are created  
dir.create(file.path(paste0(dirtemp, "sccs_sensitivity")),           showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(thisdir,"log_files/sccs_sensitivity")), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(paste0(diroutput, "sccs_sensitivity")),         showWarnings = FALSE, recursive = TRUE)

# 1. IMPORT DATA ----------------------------------------------------------

# Import data 
# CONFIRM: IS THIS RDATA OR CSV FOR DAPS??
load(paste0(dirtemp, raw_data_name, ".RData"))
scri_data_extract <- eval(parse(text = raw_data_name))

# 2. DATA CLEANING ------------------------------------------------------
## data cleaning which should apply to all case series, irrespective of DAP and outcome 
## QC and sense checks 

# Check key variables exist and are expected format, otherwise throw warning 
key_vars <- c("date_vax1", "date_vax2", "study_entry_date", "study_exit_date")

# check the vaccine, entry, and exit dates exist 
invisible(lapply(key_vars, function(x) if (!(x %in% colnames(scri_data_extract))) {
  stop(paste(x," does not exist"))
}))

# check the vaccine, entry and exit date are not character
invisible(lapply(scri_data_extract[key_vars], function(x) if (is.character(x)) {
  stop(paste(x,"is character, not date or numeric, and will not be converted by this program"))
}))

# Harmonise variable names; remove underscores
oldnames <- c(names(scri_data_extract[grepl("vax", names(scri_data_extract))]),
              names(scri_data_extract[grepl("covid", names(scri_data_extract))]))
newnames <- gsub("vax_", "vax", oldnames)
newnames <- gsub("covid_19|covid_date", "covid19", newnames)
names(scri_data_extract)[match(oldnames, names(scri_data_extract))] <- newnames

# create list of date variable names for later use (except the study start date)
date_vars <- colnames(scri_data_extract)[grepl("date", colnames(scri_data_extract))==T]
date_vars <- date_vars[date_vars != "start_study_date"]
date_vars_days <- gsub("date","days",date_vars)

# Anna: This is the BIFAP data management option, I moved this up and rewrote to a (conditional) lapply 
# it should be all the date vars as defined by you above, doesn't seem like it applies to start_study_date?
scri_data_extract[date_vars] <- lapply(scri_data_extract[date_vars],
                                       function(x) if (is.numeric(x)) as.Date(x, origin="1970-01-01") else x)

# warning to double check that formats are as expected for one of the date variables
if(is.numeric(scri_data_extract$study_exit_date)) {
  warning("Expected date, but variable is numeric. Check conversion code")
}

# create a variable for study start date 
# only relevant where someone enters late, not applicable to this study
start_scri_date <- as.Date("2020-09-01")
scri_data_extract$start_study_date <- pmax(scri_data_extract$study_entry_date, start_scri_date, na.rm = T)

# get max number of vaccine doses
nvax <- names(scri_data_extract)[ substring(names(scri_data_extract),1,8)=="date_vax" ]
nvax <- max(as.numeric(substring(nvax,9)))

# change the date variables to days from the calendar start date of our observation period 
scri_data_extract[date_vars_days] <- lapply(scri_data_extract[date_vars],
                                            function(x){round(difftime(x, start_scri_date, units = "days"), 0)})

# Check the time difference (in days) between vaccine doses
# It should be at least 14 days; print warning if this is not the case
scri_data_extract$dose_diff = as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1,units="days"))
print("Distribution of time between vaccine doses")
summary(scri_data_extract$dose_diff)
if(any(min(scri_data_extract$dose_diff[!is.na(scri_data_extract$dose_diff)])<14)) 
  warning(paste("There are ", sum(scri_data_extract$dose_diff<14, na.rm=T), " persons with (date_vax2 - date_vax1) < 14 days."))

# Confirm that the first vaccine dose occurs after the study entry date, remove rows if this does not happen
if(any(scri_data_extract$study_entry_days > scri_data_extract$days_vax1 & !is.na(scri_data_extract$days_vax1) & !is.na(scri_data_extract$study_entry_days), na.rm=T )){
  warning(paste("'study_entry_days' after the vax1 for ", sum( scri_data_extract$start > scri_data_extract$days_vax1, na.rm=T), "rows. They are deleted!"))
  nrow0<-nrow(scri_data_extract)  
  scri_data_extract <- scri_data_extract[ scri_data_extract$study_entry_days <= scri_data_extract$days_vax1 & scri_data_extract$vax1==1 & !is.na(scri_data_extract$study_entry_days),]  # if study_entry_days > vax1 ==>  delete rows
  print(c(new_nrow=nrow(scri_data_extract), old_nrow=nrow0)) 
}

# Check whether someone's study exit date occurs before the vaccine, remove rows if this does happen
if(any(scri_data_extract$study_exit_days < scri_data_extract$days_vax1 & !is.na(scri_data_extract$days_vax1), na.rm = T)) {
  warning(paste0("There are ",sum(scri_data_extract$study_exit_days < scri_data_extract$days_vax1,na.rm=T)," persons with 'study_exit_time' before vax1. They are deleted!"))
  nrow0<-nrow(scri_data_extract)  
  scri_data_extract <- scri_data_extract[scri_data_extract$days_vax1 <= scri_data_extract$study_exit_days & !is.na(scri_data_extract$study_entry_days),]  # if study_exit_days < vax1 ==>  delete rows
  print(c(new_nrow=nrow(scri_data_extract), old_nrow=nrow0)) 
}  

if(any( scri_data_extract$study_exit_days < scri_data_extract$days_vax1, na.rm=T )){
  warning(paste0("There are ",sum( scri_data_extract$study_exit_days < scri_data_extract$days_vax1,na.rm=T)," persons with 'study_exit_time' before vax1. They are deleted!"))
  nrow0<-nrow(scri_data_extract)  
  # scri_data_extract <- scri_data_extract[ scri_data_extract$days_vax1<=scri_data_extract$study_exit_days & !is.na(scri_data_extract$study_entry_days),]  # if study_exit_days < vax1 ==>  delete rows
  print(c(new_nrow=nrow(scri_data_extract), old_nrow=nrow0)) 
}  

# number of people with second vaccine dose before end of study (cond is just used in the check for clarity)
if(any((cond <- !is.na(scri_data_extract$days_vax2) & scri_data_extract$study_exit_days < scri_data_extract$days_vax2) )){
  warning(paste0("There are ",sum( cond, na.rm=T)," persons with 'study_exit_time' before vax2!"))
  scri_data_extract$study_exit_days[cond]  <- scri_data_extract$days_vax2[cond] 
}

# cleaning of vaccine doses (only applies if more than 2)
# follows conventions set for myo/perio analyses, so have retained 
if(nvax > 2){
  while(any(cond <- !is.na(scri_data_extract$days_vax2) & (scri_data_extract$days_vax2 - scri_data_extract$days_vax1) < 5)){
    warning(paste(sum(cond),"'dose2' were replace with next dose because dose2 is less than 5 days after dose1."))
    for(iv in 3:nvax){
      scri_data_extract[cond, paste0("date_vax",iv-1) ] <- scri_data_extract[cond, paste0("date_vax",iv) ]
      scri_data_extract[cond, paste0("days_vax",iv-1) ] <- scri_data_extract[cond, paste0("days_vax",iv) ]
      scri_data_extract[cond, paste0("type_vax",iv-1) ] <- scri_data_extract[cond, paste0("type_vax",iv) ]
      scri_data_extract[cond, paste0("vax"     ,iv-1) ] <- scri_data_extract[cond, paste0("vax"     ,iv) ]
    } 
  }  
}

# update the summary of time between doses after the above cleaning
scri_data_extract$dose_diff  <-  as.numeric(difftime(scri_data_extract$date_vax2, scri_data_extract$date_vax1 ,units="days"))

print("Distribution of time between vaccine doses after data cleaning")
summary(scri_data_extract$dose_diff)
if(any(min(scri_data_extract$dose_diff[!is.na(scri_data_extract$dose_diff)])<14)) 
  warning(paste("There are ", sum(scri_data_extract$dose_diff<14, na.rm=T), " persons with (date_vax2 - date_vax1) < 14 days after data cleaning"))

# create a variable for the last possible vaccine date 
scri_data_extract$days_last_vax <- pmax(scri_data_extract$days_vax1, scri_data_extract$days_vax2, na.rm = T)

# 3. SAVE DATA ------------------------------------------------------------
# saving data with new name to prevent overwriting any prior files
sccs_data_extract <- scri_data_extract 
save(sccs_data_extract,file = paste0(dirtemp,"sccs_sensitivity/","sccs_data_extract.RData"))
