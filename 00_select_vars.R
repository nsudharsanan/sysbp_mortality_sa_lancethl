#...............................................................................................
# Paper: Systolic BP and 6-year mortality in SA
# Purpose: Select needed vars to reduce the size of the data in memory for analyses
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
  
  # Many variables have the "no" as 2, this is to easily make them 0/1
  give_2_zero <- function(x) {
    x[x == 2] <- 0
    x
  }
  
#........................................................................................
# Read RDS version of NIDS_1_2_3_4_5.dta ####
#........................................................................................
  
  # Input dataset is generated using the NIDS-provided scripts
  sa = readRDS("data/nids.rds")

#..........................................................................................
# Identify relavent variables ####
#..........................................................................................

  # Variables of interest: 
    # Survey design variables
  		# Person id: pid
  		# 2008 or 2017 refresher sample: sample
  		# Continuing sample member or temporary sample member: csm (tsm's are only in for one wave should be dropped)
  		# Sample weights
  		# Died: wave_died
    # Age: w*_best_age_yrs
  	# Sex: w*_best_gen
    # Marital status: w*_a_marstt (Survey Waves 1-3), W*_best_marstt (Survey Waves 4-5)
  	# Urban/rural (geography type): w*_geo2011 (this is based on the 2011 census there is one for 2001 as well)
  	# Province: w*_prov2011 (based on 2011 census)
  	# Schooling: w*_best_edu
  	# Race: w*_best_race
  	# Self-rated health: w*_a_hldes
  	# Height: w*_a_height_1/2/3 (third measure is only used if discrepancy between first and second measurements)
  	# Weight: w*_a_weight_1/2/3 (third measure only for discrepancies larger than 1 kg)
  	# Physical activity: w*_a_hllfexer
  	# Smoking: w*a_hllfsmk
  	# BP variables
  		# Ever diagnosed: w*_a_hlbp
  		# Taking meds: w*_a_hlbp_med
  		# Sys bp: w*_a_bpsys_1/2
  		# Dia bp: w*_a_bpdia_1/2

  # Create a smaller dataset using key variables because dataset is enormous
  sm <- sa %>%
    dplyr::select(pid, csm, wave_died, sample, ends_with("_best_age_yrs"), ends_with("_best_gen"),
                  ends_with("_a_marstt"), ends_with("_best_marstt"), ends_with("_geo2011"), ends_with("_prov2011"),
                  ends_with("_best_race"), ends_with("_best_edu"), ends_with("_expenditure"),
                  ends_with("_a_height_1"), ends_with("_a_height_2"), ends_with("_a_height_3"),
                  ends_with("_a_weight_1"), ends_with("_a_weight_2"), ends_with("_a_weight_3"),
                  ends_with("_a_hldes"), ends_with("_a_hllfexer"), ends_with("a_hllfalc"), ends_with("a_hllfsmk"),
                  ends_with("_a_hlbp"), ends_with("_a_hlbp_med"), ends_with("_a_bpsys_1"), ends_with("_a_bpsys_2"),
                  ends_with("_a_bpdia_1"), ends_with("_a_bpdia_2"), ends_with("_ind_outcome"),
                  ends_with("_a_refint"), ends_with("_a_refexpl"), ends_with("_wgt"), ends_with("_dwgt")) 

  # Save dataset
  saveRDS(sm, "data/nids_small.rds")
