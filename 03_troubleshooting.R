# SCCS troubleshooting 

reduced <- scri_pre %>% 
  filter(type_vax1 == "Pfizer")

max_day = max(reduced$ref_end, na.rm = T)
month_cutoffs <- seq(from = -1, to = max_day-1, by = 30)

standardsccs(
  formula = outcome_days ~ risk_d1, 
  indiv  = numeric_id, 
  astart = ref_start,  
  aend = ref_end, 
  adrug = cbind((risk_d1 - preexp_dur), risk_d1, (risk_d1_end + 1), risk_d2), 
  aedrug = cbind(preexp_end, risk_d1_end, between_end, risk_d2_end),
  aevent = outcome_days, 
  dataformat = "multi",
  sameexpopar = F, 
  data = reduced
)

# this should throw error about the non-numeric operator 
standardsccs(
  formula = outcome_days ~ risk_d1 + relevel(age, ref = 2), 
  indiv = numeric_id, 
  astart = ref_start,  
  aend = ref_end, 
  adrug = cbind((risk_d1 - preexp_dur), risk_d1, risk_d2), 
  aedrug = cbind(preexp_end, risk_d1_end, risk_d2_end),
  aevent = outcome_days, 
  agegrp = month_cutoffs, 
  dataformat = "multi",
  sameexpopar = F,  
  data = reduced
) 

# should run 
test <- formatdata(
  indiv = numeric_id, 
  astart = ref_start,  
  aend = ref_end, 
  adrug = cbind((risk_d1 - preexp_dur), risk_d1, risk_d2), 
  aedrug = cbind(preexp_end, risk_d1_end, risk_d2_end),
  aevent = outcome_days, 
  agegrp = month_cutoffs, 
  dataformat = "multi",
  sameexpopar = F,  
  data = reduced_more
) 

# should show if there are categories with zeroes that might explain this... 
tabyl(test, event, age)
