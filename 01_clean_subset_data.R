#...............................................................................................
# Paper: Systolic BP and 6-year mortality in SA
# Purpose: Subset to eligible sample and clean vars for analysis
#...............................................................................................

#........................................
# Load packages ####
#........................................

  library(tidyverse)  # Data management
  library(sjlabelled) # Read variable labels and value lables
  
  # Functions to help clean common missing codes
  non_response2NA <- function(x){
    x[x == -2 | x == -3 | x== -5 | x == -8 | x == -9] <- NA
    x
  }
  
  # Many variables have the "no" as 2, this is to easily make the 0/1
  give_2_zero <- function(x) {
    x[x == 2] <- 0
    x
  }
  
  # This function replicates the SAS proc freq function with the '/ list missing statement'. 
  # Source: https://github.com/asnr/sas-to-r
   pfreq = function(...) {
      dplyr::group_by(...) %>%
        dplyr::summarize(n = n()) %>%
        dplyr::mutate(pct=paste0("(", format(round(100 * n/sum(n), 1), nsmall =1), ")")) 
    }
  
#........................................................................................
# Read NIDS data (just needed vars)
#........................................................................................
  
  sm <- remove_all_labels(readRDS("data/nids-small.rds")) #removing labels to avoid errors.
  
#................................................................................................
# Subset to age-eligible cohort
#................................................................................................
  
  # Age eligible cohort members are those that are 30 years old interviewed in wave 1
  sm$age_eligible.w1 <- ifelse(sm$csm == 1 & sm$w1_ind_outcome == 1 & sm$w1_best_age_yrs >= 30, 1, 0) # N = 10,246
  
  
  # Sensitivity analysis age eligible cohort members are those that are 30 years old interviewed in wave 1.
  sm$age_eligible.w2 <- ifelse(sm$csm == 1 & sm$w1_ind_outcome == 1 & sm$w2_ind_outcome == 1 &
                                 sm$w2_best_age_yrs >= 30, 1, 0) 
    # Among the 8,338 there are 819 participants who were 26-29 years old in wave 1
  
  # For the sake of imputation, we are NOT going to drop people with missing BP
  nids <- sm %>% 
    filter(age_eligible.w2 == 1)
  
#................................................................................................ 
# Clean/create variables for analysis ####
#................................................................................................

  # Create outcome variable: 6-year mortality-----------------------------------------------
  nids <- nids %>%
    mutate(died_6yrs = case_when(wave_died == 3 | wave_died == 4 | wave_died == 5 ~ 1,
                                 w5_ind_outcome == 1 ~ 0,
                                 TRUE ~ as.numeric(NA)))
  
  # Exposure (Systolic BP)-------------------------------------------------------------------------------------
  # For the sake of imputation, allow people to have this be missing. 
  nids <- nids %>% 
    mutate(
           # Correctly code missingness
           w1_a_bpsys_1 = non_response2NA(w1_a_bpsys_1), w1_a_bpsys_2 = non_response2NA(w1_a_bpsys_2), 
           w2_a_bpsys_1 = non_response2NA(w2_a_bpsys_1), w2_a_bpsys_2 = non_response2NA(w2_a_bpsys_2), 
               
           # Average wave 1 and wave 2 Systolic BP 
           avg_sbp = case_when(is.na(w1_a_bpsys_1) == F & is.na(w1_a_bpsys_2) == F &
                                 is.na(w2_a_bpsys_1) == F & is.na(w2_a_bpsys_2) == F ~
                                 (w1_a_bpsys_1 + w1_a_bpsys_2 + w2_a_bpsys_1 + w2_a_bpsys_2) / 4,
                               TRUE ~ as.numeric(NA)),
           avg_sbp_bins = as_factor(case_when(avg_sbp < 130 ~ "<130",
                                              avg_sbp >= 130 & avg_sbp < 140 ~ "130-139",
                                              avg_sbp >= 140 & avg_sbp < 150 ~ "140-149",
                                              avg_sbp >= 150 & avg_sbp < 160 ~ "150-159",
                                              avg_sbp >= 160 & avg_sbp < 170 ~ "160-169",
                                              avg_sbp >= 170 ~ "170+")),
           avg_sbp_bins = fct_relevel(avg_sbp_bins,"<130", "130-139", "140-149", "150-159", "160-169", "170+"),
           avg_sbp_bins120 = as_factor(case_when(avg_sbp < 120 ~ "<120",
                                              avg_sbp >= 120 & avg_sbp < 130 ~ "120-129",
                                              avg_sbp >= 130 & avg_sbp < 140 ~ "130-139",
                                              avg_sbp >= 140 & avg_sbp < 150 ~ "140-149",
                                              avg_sbp >= 150 & avg_sbp < 160 ~ "150-159",
                                              avg_sbp >= 160 & avg_sbp < 170 ~ "160-169",
                                              avg_sbp >= 170 ~ "170+")),
           avg_sbp_bins120 = fct_relevel(avg_sbp_bins120,"<120", "120-129", "130-139", "140-149", "150-159", "160-169", "170+"),
           w1_a_hlbp = give_2_zero(non_response2NA(w1_a_hlbp)),
           w2_a_hlbp = give_2_zero(non_response2NA(w2_a_hlbp)),
           w1_a_hlbp_med = give_2_zero(non_response2NA(w1_a_hlbp_med)),
           w2_a_hlbp_med = give_2_zero(non_response2NA(w2_a_hlbp_med)))
   
  # Covariates: Keep values observed in wave 1, 2, and 3 ----------------------------------------------------
  # Note: Take value observed in wave 2 to serve as baseline variables but keep those observed in wave 1 and 3
  #       for impuations (cary time-invariant vars from Wave 1 forward) and robustness check

  # Covariate list: age, sex, race, marital status, schooling, urbanicty, province,
  #                 self-reported heatlh, current smoking status, exercise, and BMI  
  nids <- nids %>% 
    mutate(age.w1 = w1_best_age_yrs,
           age.w1_group = as_factor(case_when(age.w1 >= 30 & age.w1 < 35 ~ "30-34",
                                              age.w1 >= 35 & age.w1 < 40 ~ "35-39",
                                              age.w1 >= 40 & age.w1 < 45 ~ "40-44",
                                              age.w1 >= 45 & age.w1 < 50 ~ "45-49",
                                              age.w1 >= 50 & age.w1 < 55 ~ "50-54",
                                              age.w1 >= 55 & age.w1 < 60 ~ "55-59",
                                              age.w1 >= 60 & age.w1 < 65 ~ "60-64",
                                              age.w1 >= 65 & age.w1 < 70 ~ "65-69",
                                              age.w1 >= 70 & age.w1 < 75 ~ "70-74",
                                              age.w1 >= 75 & age.w1 < 80 ~ "75-79",
                                              age.w1 >= 80 ~ "80+")),
           age.w2 = w2_best_age_yrs,
           age.w2_group = as_factor(case_when(age.w2 >= 30 & age.w2 < 35 ~ "30-34",
                                              age.w2 >= 35 & age.w2 < 40 ~ "35-39",
                                              age.w2 >= 40 & age.w2 < 45 ~ "40-44",
                                              age.w2 >= 45 & age.w2 < 50 ~ "45-49",
                                              age.w2 >= 50 & age.w2 < 55 ~ "50-54",
                                              age.w2 >= 55 & age.w2 < 60 ~ "55-59",
                                              age.w2 >= 60 & age.w2 < 65 ~ "60-64",
                                              age.w2 >= 65 & age.w2 < 70 ~ "65-69",
                                              age.w2 >= 70 & age.w2 < 75 ~ "70-74",
                                              age.w2 >= 75 & age.w2 < 80 ~ "75-79",
                                              age.w2 >= 80 ~ "80+")))
  
  # sex
  nids <- nids %>% 
    mutate(male.w1 = factor(give_2_zero(nids$w1_best_gen),levels = c(0,1), labels = c("Female", "Male")),
           male.w2 = factor(give_2_zero(nids$w2_best_gen),levels = c(0,1), labels = c("Female", "Male")))
           
  # race
  nids <- nids %>% 
    mutate(race.w1 = factor(non_response2NA(w1_best_race),
                                        levels = c(1,2,3, 4),
                                        labels=c("African", "Coloured","Asian/Indian", "White")),
           race.w2 = factor(non_response2NA(w2_best_race),
                                        levels = c(1,2,3, 4),
                                        labels=c("African", "Coloured","Asian/Indian", "White")))
  
  # Marital status
  nids <- nids %>%
    mutate(married.w1 = as_factor(case_when(w1_a_marstt %in% c(1,2) ~ "Married/partnered",
                                            w1_a_marstt %in% c(3,4) ~ "Widowed/divorced",
                                            w1_a_marstt == 5 ~ "Never married",
                                            TRUE ~ as.character(NA))),
           married.w1 = fct_relevel(married.w1,"Married/partnered","Widowed/divorced","Never married"),
           
           married.w2 = as_factor(case_when(w2_a_marstt %in% c(1,2) ~ "Married/partnered",
                                            w2_a_marstt %in% c(3,4) ~ "Widowed/divorced",
                                            w2_a_marstt == 5 ~ "Never married",
                                            TRUE ~ as.character(NA))),
           married.w2 = fct_relevel(married.w2,"Married/partnered","Widowed/divorced","Never married"))
  
  # Education: Kept only Wave 2 values because of variability from Wave 1 to Wave 2; only 1 missing value
  nids <- nids %>% 
    mutate(schooling.w1 = as_factor(case_when(w1_best_edu %in% c(24,25) ~ "No schooling",
                                              w1_best_edu %in% c(0:12) ~ "Primary/secondary",
                                              w1_best_edu %in% c(13:23) ~ "Tertiary",
                                              TRUE ~ as.character(NA))),
           schooling.w1 = fct_relevel(schooling.w1,"No schooling","Primary/secondary","Tertiary"),
           
           schooling.w2 = as_factor(case_when(w2_best_edu %in% c(24,25) ~ "No schooling",
                                              w2_best_edu %in% c(0:12) ~ "Primary/secondary",
                                              w2_best_edu %in% c(13:23) ~ "Tertiary",
                                              TRUE ~ as.character(NA))),
           schooling.w2 = fct_relevel(schooling.w2,"No schooling","Primary/secondary","Tertiary"))
  
  # Urbanicity
  nids <- nids %>% 
    mutate(urban.w1 = as_factor(case_when(w1_geo2011 == 2 ~ "Urban",
                                          w1_geo2011 == 1 ~ "Traditional",
                                          w1_geo2011 == 3 ~ "Farms")),
           urban.w1 = fct_relevel(urban.w1,"Urban","Traditional","Farms"),
           urban.w2 = as_factor(case_when(w2_geo2011 == 2 ~ "Urban",
                                          w2_geo2011 == 1 ~ "Traditional",
                                          w2_geo2011 == 3 ~ "Farms")),
           urban.w2 = fct_relevel(urban.w2,"Urban","Traditional","Farms"))
  
  # Province
  nids <- nids %>% 
    mutate(province.w1 = factor(nids$w1_prov2011 , levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                            labels=c("Western Cape", "Eastern Cape", "Northern Cape", "Free State",
                                     "KwaZulu-Natal", "North West", "Gauteng", "Mpumalanga", "Limpopo")),
           province.w2 = factor(nids$w2_prov2011 , levels=c(1, 2, 3, 4, 5, 6, 7, 8, 9),
                            labels=c("Western Cape", "Eastern Cape", "Northern Cape", "Free State",
                                     "KwaZulu-Natal", "North West", "Gauteng", "Mpumalanga", "Limpopo")))
  
  # Self-reported health
  nids <- nids %>% 
    mutate(srh.w1 = as_factor(case_when(non_response2NA(w1_a_hldes) == 1 ~ "Excellent",
                                        non_response2NA(w1_a_hldes) == 2 ~ "Very good",
                                        non_response2NA(w1_a_hldes) == 3 ~ "Good",
                                        non_response2NA(w1_a_hldes) == 4 ~ "Fair",
                                        non_response2NA(w1_a_hldes) == 5 ~ "Poor",
                                        TRUE ~ as.character(NA)))) %>% 
    mutate(srh.w2 = as_factor(case_when(non_response2NA(w2_a_hldes) == 1 ~ "Excellent",
                                        non_response2NA(w2_a_hldes) == 2 ~ "Very good",
                                        non_response2NA(w2_a_hldes) == 3 ~ "Good",
                                        non_response2NA(w2_a_hldes) == 4 ~ "Fair",
                                        non_response2NA(w2_a_hldes) == 5 ~ "Poor",
                                        TRUE ~ as.character(NA)))) %>% 
    mutate(srh.w1 = fct_relevel(srh.w1,"Excellent","Very good","Good", "Fair", "Poor")) %>% 
    mutate(srh.w2 = fct_relevel(srh.w2,"Excellent","Very good","Good", "Fair", "Poor"))
  
  # Smoker
  nids <- nids %>% 
    mutate(smoke.w1 = as_factor(give_2_zero(non_response2NA(nids$w1_a_hllfsmk))),
           smoke.w2 = as_factor(give_2_zero(non_response2NA(nids$w2_a_hllfsmk))))
                                
  # Alcohol
  nids <- nids %>% 
    mutate(alcohol.w1 = fct_relevel(as_factor(case_when(non_response2NA(w1_a_hllfalc) %in% c(1:4) ~ "<1",  # <1 day/week",
                                                       non_response2NA(w1_a_hllfalc) %in% c(5) ~   "1-2", # 1-2 days/week",
                                                       non_response2NA(w1_a_hllfalc) %in% c(6:8) ~ "3+",  # 3+ days/week",
                                                       TRUE ~ as.character(NA))),
                                   "<1", "1-2", "3+"),
           alcohol.w2 = fct_relevel(as_factor(case_when(non_response2NA(w2_a_hllfalc) %in% c(1:4) ~ "<1",  # <1 day/week",
                                                       non_response2NA(w2_a_hllfalc) %in% c(5) ~   "1-2", # 1-2 days/week",
                                                       non_response2NA(w2_a_hllfalc) %in% c(6:8) ~ "3+",  # 3+ days/week",
                                                       TRUE ~ as.character(NA))),
                                   "<1", "1-2", "3+"))
  
  # Exercise
  nids <- nids %>% 
    mutate(exercise.w1 = as_factor(case_when(non_response2NA(nids$w1_a_hllfexer) %in% c(1:2) ~ "<1",  # <1 time/wk
                                             non_response2NA(nids$w1_a_hllfexer) %in% c(3:4) ~ "1-2", # 1-2 times/wk
                                             non_response2NA(nids$w1_a_hllfexer) %in% c(5) ~ "3+",    # >=3 times/wk
                                             TRUE ~ as.character(NA)))) %>% 
    mutate(exercise.w2 = as_factor(case_when(non_response2NA(nids$w2_a_hllfexer) %in% c(1:2) ~ "<1",  # <1 time/wk
                                             non_response2NA(nids$w2_a_hllfexer) %in% c(3:4) ~ "1-2", # 1-2 times/wk
                                             non_response2NA(nids$w2_a_hllfexer) %in% c(5)   ~ "3+",  # >=3 times/wk
                                             TRUE ~ as.character(NA)))) %>% 
    mutate(exercise.w1 = fct_relevel(exercise.w1,"<1", "1-2", "3+")) %>% 
    mutate(exercise.w2 = fct_relevel(exercise.w2,"<1", "1-2", "3+"))
  
  # BMI: Referred to Cois & Ehrlich 2018 (https://doi.org/10.1371/journal.pone.0200606) for cleaning 
  # Cois and Ehrlich note in the methods section that "Duplicate measures of of weight and height were recorded, 
  # with a third measurement taken if their difference was greater than 0.5 kg or 0.5 cm, respectively."
  # When I check this for weight measurements for wave 2, there were 169 obs with weight 1 and weight 2 differences
  # greater than 0.5 kg but the third weight measurement was missing for the majority. The difference criteria appears to be
  # a difference of >1 kg. There were 24 obs where the difference was >1 kg and 22 of these observations had a third measurement.
  # As expected there were observations with wave 1 & 2, Wave 1 & 3 or Wave 2 & 3 had closer measurements. However, I have
  # figured out an efficient/valid way to taking the sum of two closest measurements. 
  # This issue is similar for the height measurements where 13 obs had a difference of > 1 cm. 
  # For now I am taking the average of the first two measurements. 
  nids <- nids %>% 
    mutate(weight.w1 = as.numeric(non_response2NA(w1_a_weight_1) + non_response2NA(w1_a_weight_2), na.rm = T) / 2, 
           height.w1 = as.numeric(non_response2NA(w1_a_height_1) + non_response2NA(w1_a_height_2), na.rm = T) / 2,
           bmi.w1 = weight.w1/((height.w1/100)^2)) %>% # height was measured in cm so converted to meters
    mutate(weight.w2 = as.numeric(non_response2NA(w2_a_weight_1) + non_response2NA(w2_a_weight_2), na.rm = T) / 2, 
           height.w2 = as.numeric(non_response2NA(w2_a_height_1) + non_response2NA(w2_a_height_2), na.rm = T) / 2,
           bmi.w2 = weight.w2/((height.w2/100)^2)) # height was measured in cm so converted to meters

  # Antihypertensive drug use: 
  nids <- nids %>% 
    mutate(w1_a_hlbp = give_2_zero(non_response2NA(nids$w1_a_hlbp)),
           w2_a_hlbp = give_2_zero(non_response2NA(nids$w2_a_hlbp)),
           w1_a_hlbp_med = give_2_zero(non_response2NA(nids$w1_a_hlbp_med)),
           w2_a_hlbp_med = give_2_zero(non_response2NA(nids$w2_a_hlbp_med)),
           bp_drugs.w2 = as_factor(case_when((w2_a_hlbp == 0 & is.na(w2_a_hlbp_med)) |
                                               w2_a_hlbp == 1 & w2_a_hlbp_med == 0 ~ 0,
                                             w2_a_hlbp == 1 & w2_a_hlbp_med == 1 ~ 1,
                                             TRUE ~ as.numeric(NA))))
                         
#................................................................................................   
# Export analysis data ####
#................................................................................................ 
 
  # Create dataset with analysis variables to use for assessing missingness and 
  # generate final sample after imputation. 
  analysis <- nids %>% 
    select(pid, w1_ind_outcome, w2_ind_outcome,
          
          # Outcome
          died_6yrs,
          
          # Exposure
          w1_a_bpsys_1, w1_a_bpsys_2, w2_a_bpsys_1, w2_a_bpsys_2,
          avg_sbp, avg_sbp_bins, avg_sbp_bins120, 
   
          #Included for response to reviewers
          w1_a_hlbp, w2_a_hlbp, w1_a_hlbp_med, w2_a_hlbp_med, bp_drugs.w2,
          
          # Covariates
          age.w2, age.w2_group, male.w2,  race.w2,   # age, sex, race
          married.w2, schooling.w2,                 # marital status, schooling
          urban.w2, province.w2,                   # urban, province
          srh.w2, smoke.w2, alcohol.w2,            # current smoker, alcohol, self-reported health status
          exercise.w2, bmi.w2,                   # exercise, BMI
          
          # Design weights
          w2_wgt, w2_dwgt)                    
 
  # Save analysis dataset
  saveRDS(analysis, "data/analysis_data.rds")
  