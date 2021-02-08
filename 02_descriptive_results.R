#...............................................................................................
# Paper: Systolic BP and 6-year mortality in SA
# Purpose: Describe the sample including missingness
#...............................................................................................

#.........................................................................................
# Load packages
#.........................................................................................
  
  library(tidyverse) # Data management & Graphics
  library(visdat)    # Visualization of missing data
  library(naniar)    # Visualization of missing data
  library(ggsci)     # Graphics
  library(tableone)  # Table 1 summary
  library(ggeffects) # Estimate marginal means from regression models
  library(srvyr)     # calculate mean after applying survey weights
  library(survey)    # calculate mean after applying survey weights
  library(splines)   # Modeling splines
 
  #This function replicates the SAS proc freq function with the '/ list missing statement'. 
  #Source: https://github.com/asnr/sas-to-r
  pfreq = function(...) {
    group_by(...) %>%
      summarise(n = n()) %>%
      mutate(pct = round(100 * n/sum(n), 1)) 
  }
  # the fuction is used as: pfreq(mydata, myvariables)
  
#........................................................................................
# Read data
#........................................................................................

  nids <-  readRDS("data/analysis_data.rds") %>% 
    mutate(ltfu = as_factor(case_when(is.na(died_6yrs) == TRUE ~ "Yes",
                                      TRUE ~ as.character("No"))),
           bp_3 = ifelse(avg_sbp_bins120 %in% c("<120", "120-129"), 1,
                         ifelse(avg_sbp_bins120 %in% c("130-139"),2,3)),
           bp_120 = ifelse(avg_sbp_bins120 %in% c("<120", "120-129"), 1, 0),
           bp_130 = ifelse(avg_sbp_bins120 %in% c("130-139"), 1, 0),
           bp_140 = ifelse(avg_sbp_bins120 %in% c("<120", "120-129", "130-139"), 0, 1),
           bp_150 = ifelse(avg_sbp_bins120 %in% c("<120", "120-129", "130-139", "140-149" ), 0, 1),
           bp_160 = ifelse(avg_sbp_bins120 %in% c("<120", "120-129", "130-139", "140-149", "150-159"), 0, 1),
           age_30_49 = as_factor(case_when(age.w2_group %in% c("30-34", "35-39", "40-44", "45-49") ~ "30-49",
                                           TRUE ~ as.character(">=50"))))
  analysis <- na.omit(nids)
  
#........................................................................................
# Missingness Figure Elements & Missing data Exploration
#........................................................................................

  # Number lost to follow-up (missing outcome information) in subsequent waves (3-5)
  pfreq(nids, died_6yrs) # 927 (11.1%) leaving 7,411 with vital information in 2017
  pfreq(nids, ltfu)
  
  # Number missing exposure:
  followed <- nids %>% filter(ltfu == "No") %>% 
    select(-w2_dwgt, -age.w2_group)
  n_miss(followed$avg_sbp) # 2,236 (30.2%), thus Number of age eligible respondets
  # in wave 2 with measured BP in wave 1 & 2: 5,1754
  
  # number of missing at least one confounder: 
  conf <- nids %>% filter(ltfu == "No") %>% 
    select(age.w2, male.w2, race.w2, schooling.w2, province.w2, urban.w2, married.w2,
            srh.w2,  smoke.w2, exercise.w2, alcohol.w2, bmi.w2)
   
  7411 - n_case_complete(conf) #1,355 (18.2%)
  
  # Variables with missing values = 5 / 9; Overall there were 54 Obs missing values for at least one of these counfounders
  gg_miss_upset((conf %>% 
                   select(age.w2, male.w2, race.w2, schooling.w2, province.w2, urban.w2, married.w2,
                          srh.w2,  smoke.w2, exercise.w2, alcohol.w2, bmi.w2)),
                    nsets = 12,
                    nintersects = NA)
  
  # Age-eligible respondents missing BP or a confounder value
  7411 - n_case_complete(followed) #2,418 (32.6%)
  
  # Total non-missing
  n_case_complete(nids) #4,993 (32.6%)
  
#.......................................................................................
# Table 1 Descriptive characteristics of included participants
#.......................................................................................
  
  # tableone package was not recognizing age.w2 as numeric making sure that age is numeric
  analysis$age.w2 <- as.numeric(analysis$age.w2)
  
  # List of variables to summaries
  summary_vars <- c("age.w2", "male.w2", "race.w2", "married.w2", "schooling.w2",
                    "urban.w2", "province.w2", "srh.w2","smoke.w2", "alcohol.w2",
                    "exercise.w2", "bmi.w2")
  
  # List of categorical variables
  catvars <-  c("male.w2", "race.w2", "married.w2", "schooling.w2",
                "urban.w2","province.w2","srh.w2","smoke.w2", 
                "alcohol.w2", "exercise.w2")
  
  # Create table 1 object: continuous vars = mean (SD) by sex
  tbl1 <- CreateTableOne(vars = summary_vars, factorVars = catvars, data = analysis, test = FALSE)

  # Print table 1 to verify calculations
  tbl1_print <- print(tbl1, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

#.......................................................................................
# Descriptive characteristics of participants LTFU
#.......................................................................................
  
  # List of variables to summaries
  summary_ltfu <- c("ltfu", "age.w2", "male.w2", "race.w2", "married.w2", "schooling.w2",
                    "urban.w2", "province.w2", "srh.w2","smoke.w2", "alcohol.w2",
                    "exercise.w2", "bmi.w2")
  
  # List of categorical variables
  catvars_ltfu <-  c("ltfu", "male.w2", "race.w2", "married.w2", "schooling.w2",
                     "urban.w2","province.w2","srh.w2","smoke.w2", 
                     "alcohol.w2", "exercise.w2")
  
  # Create table 1 object: continuous vars = mean (SD) by sex
  tbl1_ltfu <- CreateTableOne(vars = summary_ltfu, factorVars = catvars_ltfu,  strata = "ltfu", data = nids, test = FALSE)

  # Print table 1 to verify calculations
  tbl1_print_ltfu <- print(tbl1_ltfu, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

#........................................................................................
# Figure 1: Description of systolic BP
#........................................................................................ 
  
  # Summarize BP bins
  sbp.bins.summary <- analysis %>%
    filter(is.na(avg_sbp) == FALSE) %>% # drop participants missing BP (N = 2,074)
    select(age.w2_group, avg_sbp_bins120, w2_dwgt) %>% 
    group_by(age.w2_group, avg_sbp_bins120) %>% 
    summarise_each(list(~ sum(., na.rm = TRUE))) %>% 
    mutate(pct_weighted = round(w2_dwgt/sum(w2_dwgt)*100,1))

  # Graph distribution SBP avg
  ggplot(sbp.bins.summary, aes(fill=fct_rev(avg_sbp_bins120), y=pct_weighted, x=age.w2_group)) + 
    geom_bar(position="fill", stat="identity")+
    labs(#title = "Figure 2a: Age-systolic BP relationship" ,
         x = "Age (years)",
         y = "Proportion in SBP categories") +
    guides(fill=guide_legend(title="SBP categories")) +
    scale_fill_nejm() +
    theme_bw() +
    theme(axis.title = element_text(color="black",size=12),
          axis.text = element_text(color = "black", size =  12))
  ggsave("figure1.eps",dpi = 400,units = "in",width = 8,height = 8)

  # Summaries with means and CI for manuscript
  
    # Create survey object using srvyr package
    srvyr_obj <- analysis %>%
      filter(is.na(w2_dwgt) == FALSE) %>% # remove participants missing post-stratification weights (N = 12)
      as_survey_design(ids = pid, weight = w2_dwgt)
    
    age_sbp_summary <- srvyr_obj %>%
      group_by(age_30_49) %>%
      summarize(bp_120_summary = survey_mean(bp_120),
                bp_130_summary = survey_mean(bp_130),
                bp_140_summary = survey_mean(bp_140),
                bp_150_summary = survey_mean(bp_150),
                bp_160_summary = survey_mean(bp_160)) %>% 
      mutate(bp_120_lci = bp_120_summary -  (1.96*bp_120_summary_se),
             bp_120_uci = bp_120_summary +  (1.96*bp_120_summary_se),
             bp_130_lci = bp_130_summary -  (1.96*bp_130_summary_se),
             bp_130_uci = bp_130_summary +  (1.96*bp_130_summary_se),
             bp_140_lci = bp_140_summary -  (1.96*bp_140_summary_se),
             bp_140_uci = bp_140_summary +  (1.96*bp_140_summary_se),
             bp_150_lci = bp_150_summary -  (1.96*bp_150_summary_se),
             bp_150_uci = bp_150_summary +  (1.96*bp_150_summary_se),
             bp_160_lci = bp_160_summary -  (1.96*bp_160_summary_se),
             bp_160_uci = bp_160_summary +  (1.96*bp_160_summary_se)) %>% 
      select(-ends_with("_se"), bp_120_summary, bp_120_lci, bp_120_uci,
             bp_130_summary, bp_130_lci, bp_130_uci,
             bp_140_summary, bp_140_lci, bp_140_uci,
             bp_150_summary, bp_150_lci, bp_150_uci,
             bp_160_summary, bp_160_lci, bp_160_uci)
    
  age_sbp_summary %>% write.table("clipboard", sep="\t", row.names = F)
  

  