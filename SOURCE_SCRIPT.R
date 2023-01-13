# Program Information  ----------------------------------------------------
#
# Program:      SOURCE_SCRIPT
# Description:  calls all required scripts of the ROC20 WP4 SCCS analysis pipeline 

# SCRIPTS --- ------------------------------------------------------------

source("00_sensitivity_functions.R", print.eval = TRUE)
source("01_clean_data.R", print.eval = TRUE)
source("02_select_population.R", print.eval = TRUE)
source("03_describe_event_dependency.R", print.eval = TRUE)
source("04a_scri_pre.R", print.eval = TRUE)
source("04b_scri_post.R", print.eval = TRUE)
source("04c_standard_sccs.R", print.eval = TRUE)
source("04d_standard_sccs.R", print.eval = TRUE)
