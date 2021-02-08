#...............................................................................................
# Paper: Systolic BP and 6-year mortality in SA
# Purpose: Do main analysis
#...............................................................................................

#.........................................................................................
# Load packages
#.........................................................................................
  
  library(tidyverse)  # Data management
  library(ggpubr) #Combine ggplots
  library(cowplot) #Combine ggplots

#........................................................................................
# Read cleaned data
#........................................................................................
  
  nids = readRDS("data/analysis_data.rds")
  
#........................................................................................
# Pre-cleaning
#........................................................................................
  
  # Keep only necessary variables
  analysis <- nids %>% 
    select(c("died_6yrs", "avg_sbp", "age.w2", "male.w2", "race.w2", "married.w2", "schooling.w2",
             "urban.w2", "province.w2", "srh.w2","smoke.w2", "alcohol.w2",
             "exercise.w2", "bmi.w2", "age.w2_group","w2_dwgt")) 
  
  # Subset to complete cases
  analysis <- analysis[complete.cases(analysis),]
  
#.......................................................................................
# Bootstrap parameters
#.......................................................................................
  
  # Number of bootstrap loops
  bsize <- 500
  
  # Bootstrap seed
  set.seed(15062020)
  
  # Within each bootstrap loop we need to save
    # 1. Mortality at 120, 130, 140, 150, 160
    # 2. Mean SBP shift
    # 3. Total share under each scenario
    # 4. Population effect
    
  # Levels of stratification variable
  strata <- 10 # Overall, Sex (Men, Women), Race (African, Coloured, Asian/Indian, White), Urbanicity (Urban, Traditional, Farms)
  
  # Number of scenarioss
  nscenarios <- 4
  
  # Strata * 4 scenarios * 500 loops
  bs_pop <- matrix(ncol = strata*nscenarios, nrow = bsize)
  colnames(bs_pop) <- c("Overall-120",  "Overall-130", "Overall-140", "Overall-150",
                        "Men-120",  "Men-130", "Men-140", "Men-150",
                        "Women-120",  "Women-130", "Women-140", "Women-150",
                        "African-120",  "African-130", "African-140", "African-150",
                        "Coloured-120",  "Coloured-130", "Coloured-140", "Coloured-150",
                        "Asian-120",  "Asian-130", "Asian-140", "Asian-150",
                        "White-120",  "White-130", "White-140", "White-150",
                        "Urban-120",  "Urban-130", "Urban-140", "Urban-150",
                        "Trad-120",  "Trad-130", "Trad-140", "Trad-150",
                        "Farms-120",  "Farms-130", "Farms-140", "Farms-150") 
  
  #Vectors to save output
  bs_160 <- bs_150 <- bs_140 <- bs_130 <- bs_120 <- matrix(nrow = bsize, ncol = 1)
  bs_share <- bs_shift <- bs_pop
  
#.......................................................................................
# Start the bootstrap
#.......................................................................................

  # G-formula estimates: Start the bootstrap
  for(b in 1:bsize) {
    print(b)
     
    # Make separate male and female data so that we can do a stratified bootstrap
    male <- filter(analysis, male.w2 == "Male")
    female <- filter(analysis, male.w2 == "Female")
    
    # Weighted draw
    index_male <- sample(length(male$died_6yrs), replace = T, prob=male$w2_dwgt)
    male_bs <- male[index_male,]
    index_female <- sample(length(female$died_6yrs), replace = T, prob=female$w2_dwgt)
    female_bs <- female[index_female,]
    bs_data <- rbind(male_bs,female_bs)
    
#.......................................................................................
# Fit the outcome model
#.......................................................................................

  fit_died_6yrs <- glm(died_6yrs ~ avg_sbp + I(avg_sbp*avg_sbp) + I(avg_sbp*avg_sbp*avg_sbp) +
                         age.w2 + I(age.w2^2) + male.w2 + race.w2 +
                         married.w2 + schooling.w2 + urban.w2 + province.w2 + srh.w2 +
                         smoke.w2 + alcohol.w2 + exercise.w2 + bmi.w2 + I(bmi.w2 ^2),
                       data = bs_data,
                       family = binomial)
   
#.......................................................................................
# Predict mortality at the different BP levels to generate incremental mortality reductions figure
#.......................................................................................
    
  #Dataset to hold all predictions
  predictions <- bs_data
    
    #160
    predictions$avg_sbp <- 160
    bs_160[b,1] <- mean(predict(fit_died_6yrs,newdata = predictions, type = "response"))
    
    #150
    predictions$avg_sbp<- 150
    bs_150[b,1] <- mean(predict(fit_died_6yrs,newdata = predictions, type = "response"))
    
    #140
    predictions$avg_sbp <- 140
    bs_140[b,1] <- mean(predict(fit_died_6yrs,newdata = predictions, type = "response"))
    
    #130
    predictions$avg_sbp <- 130
    bs_130[b,1] <- mean(predict(fit_died_6yrs,newdata = predictions, type = "response"))
    
    #120
    predictions$avg_sbp <- 120
    bs_120[b,1] <- mean(predict(fit_died_6yrs,newdata = predictions, type = "response"))

#.......................................................................................
# Predict the natural course
#.......................................................................................    

  # Predict mortality - can be done without rbinom since we are aggregating anyway
  bs_data$died <- predict(fit_died_6yrs, newdata = bs_data, type = "response")
  
#.......................................................................................
# Flag people who are eligible to receive each BP intervention
#.......................................................................................    
  
  bs_data <- bs_data %>% 
    
    # Intervention 1: 120
    mutate(eligible1 = case_when(avg_sbp > 120 ~ 1, 
                                 avg_sbp <=120 ~ 0)) %>% 
    
    # Intervention 2: 130
    mutate(eligible2 = case_when(avg_sbp > 130 ~ 1, 
                                 avg_sbp <=130 ~ 0)) %>% 
    
    # Intervention 3: 140
    mutate(eligible3 = case_when(avg_sbp > 140 ~ 1, 
                                 avg_sbp <=140 ~ 0)) %>%  
    
    # Intervention 4: 150
    mutate(eligible4 = case_when(avg_sbp > 150 ~ 1, 
                                 avg_sbp <=150 ~ 0))
    
#.......................................................................................
# Intervention 1: Target of 120
#.......................................................................................    
  
  int1_data <- bs_data
  
    # Do the intervention
    int1_data$avg_sbp <- ifelse(int1_data$avg_sbp > 120, 120,int1_data$avg_sbp)
    
    # Predict mortality - can be done without rbinom since we are aggregating anyway
    int1_data$died <- predict(fit_died_6yrs, newdata = int1_data, type = "response")
  
#.......................................................................................
# Intervention 2: Target of 130
#.......................................................................................    

  int2_data <- bs_data

    # Do the intervention
    int2_data$avg_sbp <- ifelse(int2_data$avg_sbp > 130, 130,int2_data$avg_sbp)
    
    # Predict mortality - can be done without rbinom since we are aggregating anyway
    int2_data$died <- predict(fit_died_6yrs, newdata = int2_data, type = "response")
    
#.......................................................................................
# Intervention 3: Target of 140
#.......................................................................................

  int3_data <- bs_data
  
    # Do the intervention
    int3_data$avg_sbp <- ifelse(int3_data$avg_sbp > 140, 140,int3_data$avg_sbp)
    
    # Predict mortality - can be done without rbinom since we are aggregating anyway
    int3_data$died <- predict(fit_died_6yrs, newdata = int3_data, type = "response")
    
#.......................................................................................
# Intervention 4: Target of 150
#.......................................................................................

  int4_data <- bs_data
  
    # Do the intervention
    int4_data$avg_sbp <- ifelse(int4_data$avg_sbp > 150, 150,int4_data$avg_sbp)
    
    # Predict mortality - can be done without rbinom since we are aggregating anyway
    int4_data$died <- predict(fit_died_6yrs, newdata = int4_data, type = "response")    
        
#.......................................................................................
# Share who would need to be on treatment
#.......................................................................................
  
  # Overall
  bs_share[b,1] <- mean(bs_data$eligible1)
  bs_share[b,2] <- mean(bs_data$eligible2)
  bs_share[b,3] <- mean(bs_data$eligible3)
  bs_share[b,4] <- mean(bs_data$eligible4)
  
  # Men
  bs_share[b,5] <- mean(bs_data$eligible1[bs_data$male.w2=="Male"])
  bs_share[b,6] <- mean(bs_data$eligible2[bs_data$male.w2=="Male"])
  bs_share[b,7] <- mean(bs_data$eligible3[bs_data$male.w2=="Male"])
  bs_share[b,8] <- mean(bs_data$eligible4[bs_data$male.w2=="Male"])
  
  # Female
  bs_share[b,9]  <- mean(bs_data$eligible1[bs_data$male.w2=="Female"])
  bs_share[b,10] <- mean(bs_data$eligible2[bs_data$male.w2=="Female"])
  bs_share[b,11] <- mean(bs_data$eligible3[bs_data$male.w2=="Female"])
  bs_share[b,12] <- mean(bs_data$eligible4[bs_data$male.w2=="Female"])
  
  # African
  bs_share[b,13] <- mean(bs_data$eligible1[bs_data$race.w2=="African"])
  bs_share[b,14] <- mean(bs_data$eligible2[bs_data$race.w2=="African"])
  bs_share[b,15] <- mean(bs_data$eligible3[bs_data$race.w2=="African"])
  bs_share[b,16] <- mean(bs_data$eligible4[bs_data$race.w2=="African"])
  
  # Coloured
  bs_share[b,17] <- mean(bs_data$eligible1[bs_data$race.w2=="Coloured"])
  bs_share[b,18] <- mean(bs_data$eligible2[bs_data$race.w2=="Coloured"])
  bs_share[b,19] <- mean(bs_data$eligible3[bs_data$race.w2=="Coloured"])
  bs_share[b,20] <- mean(bs_data$eligible4[bs_data$race.w2=="Coloured"])
  
  # Asian/Indian
  bs_share[b,21] <- mean(bs_data$eligible1[bs_data$race.w2=="Asian/Indian"])
  bs_share[b,22] <- mean(bs_data$eligible2[bs_data$race.w2=="Asian/Indian"])
  bs_share[b,23] <- mean(bs_data$eligible3[bs_data$race.w2=="Asian/Indian"])
  bs_share[b,24] <- mean(bs_data$eligible4[bs_data$race.w2=="Asian/Indian"])
  
  # White
  bs_share[b,25] <- mean(bs_data$eligible1[bs_data$race.w2=="White"])
  bs_share[b,26] <- mean(bs_data$eligible2[bs_data$race.w2=="White"])
  bs_share[b,27] <- mean(bs_data$eligible3[bs_data$race.w2=="White"])
  bs_share[b,28] <- mean(bs_data$eligible4[bs_data$race.w2=="White"])
  
  # Urban
  bs_share[b,29] <- mean(bs_data$eligible1[bs_data$urban.w2=="Urban"])
  bs_share[b,30] <- mean(bs_data$eligible2[bs_data$urban.w2=="Urban"])
  bs_share[b,31] <- mean(bs_data$eligible3[bs_data$urban.w2=="Urban"])
  bs_share[b,32] <- mean(bs_data$eligible4[bs_data$urban.w2=="Urban"])
  
  # Traditional
  bs_share[b,33] <- mean(bs_data$eligible1[bs_data$urban.w2=="Traditional"])
  bs_share[b,34] <- mean(bs_data$eligible2[bs_data$urban.w2=="Traditional"])
  bs_share[b,35] <- mean(bs_data$eligible3[bs_data$urban.w2=="Traditional"])
  bs_share[b,36] <- mean(bs_data$eligible4[bs_data$urban.w2=="Traditional"])
  
  # Farms
  bs_share[b,37] <- mean(bs_data$eligible1[bs_data$urban.w2=="Farms"])
  bs_share[b,38] <- mean(bs_data$eligible2[bs_data$urban.w2=="Farms"])
  bs_share[b,39] <- mean(bs_data$eligible3[bs_data$urban.w2=="Farms"])
  bs_share[b,40] <- mean(bs_data$eligible4[bs_data$urban.w2=="Farms"])
  
#.......................................................................................
# Mean shift in systolic BP
#.......................................................................................
  
  # Overall
  bs_shift[b,1] <- mean(int1_data$avg_sbp) - mean(bs_data$avg_sbp)
  bs_shift[b,2] <- mean(int2_data$avg_sbp) - mean(bs_data$avg_sbp)
  bs_shift[b,3] <- mean(int3_data$avg_sbp) - mean(bs_data$avg_sbp)
  bs_shift[b,4] <- mean(int4_data$avg_sbp) - mean(bs_data$avg_sbp)
  
  # Men
  bs_shift[b,5] <- mean(int1_data$avg_sbp[bs_data$male.w2=="Male"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Male"])
  bs_shift[b,6] <- mean(int2_data$avg_sbp[bs_data$male.w2=="Male"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Male"])
  bs_shift[b,7] <- mean(int3_data$avg_sbp[bs_data$male.w2=="Male"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Male"])
  bs_shift[b,8] <- mean(int4_data$avg_sbp[bs_data$male.w2=="Male"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Male"])
  
  # Women
  bs_shift[b,9] <- mean(int1_data$avg_sbp[bs_data$male.w2=="Female"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Female"])
  bs_shift[b,10] <- mean(int2_data$avg_sbp[bs_data$male.w2=="Female"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Female"])
  bs_shift[b,11] <- mean(int3_data$avg_sbp[bs_data$male.w2=="Female"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Female"])
  bs_shift[b,12] <- mean(int4_data$avg_sbp[bs_data$male.w2=="Female"]) - mean(bs_data$avg_sbp[bs_data$male.w2=="Female"])
 
  # African
  bs_shift[b,13] <- mean(int1_data$avg_sbp[bs_data$race.w2=="African"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="African"])
  bs_shift[b,14] <- mean(int2_data$avg_sbp[bs_data$race.w2=="African"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="African"])
  bs_shift[b,15] <- mean(int3_data$avg_sbp[bs_data$race.w2=="African"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="African"])
  bs_shift[b,16] <- mean(int4_data$avg_sbp[bs_data$race.w2=="African"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="African"])
  
  # Coloured
  bs_shift[b,17] <- mean(int1_data$avg_sbp[bs_data$race.w2=="Coloured"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Coloured"])
  bs_shift[b,18] <- mean(int2_data$avg_sbp[bs_data$race.w2=="Coloured"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Coloured"])
  bs_shift[b,19] <- mean(int3_data$avg_sbp[bs_data$race.w2=="Coloured"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Coloured"])
  bs_shift[b,20] <- mean(int4_data$avg_sbp[bs_data$race.w2=="Coloured"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Coloured"])
  
  # Asian/Indian
  bs_shift[b,21] <- mean(int1_data$avg_sbp[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Asian/Indian"])
  bs_shift[b,22] <- mean(int2_data$avg_sbp[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Asian/Indian"])
  bs_shift[b,23] <- mean(int3_data$avg_sbp[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Asian/Indian"])
  bs_shift[b,24] <- mean(int4_data$avg_sbp[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="Asian/Indian"])
  
  # White
  bs_shift[b,25] <- mean(int1_data$avg_sbp[bs_data$race.w2=="White"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="White"])
  bs_shift[b,26] <- mean(int2_data$avg_sbp[bs_data$race.w2=="White"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="White"])
  bs_shift[b,27] <- mean(int3_data$avg_sbp[bs_data$race.w2=="White"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="White"])
  bs_shift[b,28] <- mean(int4_data$avg_sbp[bs_data$race.w2=="White"]) - mean(bs_data$avg_sbp[bs_data$race.w2=="White"])
  
  # Urban
  bs_shift[b,29] <- mean(int1_data$avg_sbp[bs_data$urban.w2=="Urban"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Urban"])
  bs_shift[b,30] <- mean(int2_data$avg_sbp[bs_data$urban.w2=="Urban"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Urban"])
  bs_shift[b,31] <- mean(int3_data$avg_sbp[bs_data$urban.w2=="Urban"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Urban"])
  bs_shift[b,32] <- mean(int4_data$avg_sbp[bs_data$urban.w2=="Urban"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Urban"])
  
  # Traditional
  bs_shift[b,33] <- mean(int1_data$avg_sbp[bs_data$urban.w2=="Traditional"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Traditional"])
  bs_shift[b,34] <- mean(int2_data$avg_sbp[bs_data$urban.w2=="Traditional"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Traditional"])
  bs_shift[b,35] <- mean(int3_data$avg_sbp[bs_data$urban.w2=="Traditional"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Traditional"])
  bs_shift[b,36] <- mean(int4_data$avg_sbp[bs_data$urban.w2=="Traditional"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Traditional"])
  
  # Farms
  bs_shift[b,37] <- mean(int1_data$avg_sbp[bs_data$urban.w2=="Farms"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Farms"])
  bs_shift[b,38] <- mean(int2_data$avg_sbp[bs_data$urban.w2=="Farms"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Farms"])
  bs_shift[b,39] <- mean(int3_data$avg_sbp[bs_data$urban.w2=="Farms"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Farms"])
  bs_shift[b,40] <- mean(int4_data$avg_sbp[bs_data$urban.w2=="Farms"]) - mean(bs_data$avg_sbp[bs_data$urban.w2=="Farms"])

#.......................................................................................
# Population reduction in mortality
#.......................................................................................

  # Overall
  bs_pop[b,1] <- mean(int1_data$died) - mean(bs_data$died)
  bs_pop[b,2] <- mean(int2_data$died) - mean(bs_data$died)
  bs_pop[b,3] <- mean(int3_data$died) - mean(bs_data$died)
  bs_pop[b,4] <- mean(int4_data$died) - mean(bs_data$died)
  
  # Men
  bs_pop[b,5] <- mean(int1_data$died[bs_data$male.w2=="Male"]) - mean(bs_data$died[bs_data$male.w2=="Male"])
  bs_pop[b,6] <- mean(int2_data$died[bs_data$male.w2=="Male"]) - mean(bs_data$died[bs_data$male.w2=="Male"])
  bs_pop[b,7] <- mean(int3_data$died[bs_data$male.w2=="Male"]) - mean(bs_data$died[bs_data$male.w2=="Male"])
  bs_pop[b,8] <- mean(int4_data$died[bs_data$male.w2=="Male"]) - mean(bs_data$died[bs_data$male.w2=="Male"])
  
  # Women
  bs_pop[b,9] <- mean(int1_data$died[bs_data$male.w2=="Female"]) - mean(bs_data$died[bs_data$male.w2=="Female"])
  bs_pop[b,10] <- mean(int2_data$died[bs_data$male.w2=="Female"]) - mean(bs_data$died[bs_data$male.w2=="Female"])
  bs_pop[b,11] <- mean(int3_data$died[bs_data$male.w2=="Female"]) - mean(bs_data$died[bs_data$male.w2=="Female"])
  bs_pop[b,12] <- mean(int4_data$died[bs_data$male.w2=="Female"]) - mean(bs_data$died[bs_data$male.w2=="Female"])
 
  # African
  bs_pop[b,13] <- mean(int1_data$died[bs_data$race.w2=="African"]) - mean(bs_data$died[bs_data$race.w2=="African"])
  bs_pop[b,14] <- mean(int2_data$died[bs_data$race.w2=="African"]) - mean(bs_data$died[bs_data$race.w2=="African"])
  bs_pop[b,15] <- mean(int3_data$died[bs_data$race.w2=="African"]) - mean(bs_data$died[bs_data$race.w2=="African"])
  bs_pop[b,16] <- mean(int4_data$died[bs_data$race.w2=="African"]) - mean(bs_data$died[bs_data$race.w2=="African"])
  
  # Coloured
  bs_pop[b,17] <- mean(int1_data$died[bs_data$race.w2=="Coloured"]) - mean(bs_data$died[bs_data$race.w2=="Coloured"])
  bs_pop[b,18] <- mean(int2_data$died[bs_data$race.w2=="Coloured"]) - mean(bs_data$died[bs_data$race.w2=="Coloured"])
  bs_pop[b,19] <- mean(int3_data$died[bs_data$race.w2=="Coloured"]) - mean(bs_data$died[bs_data$race.w2=="Coloured"])
  bs_pop[b,20] <- mean(int4_data$died[bs_data$race.w2=="Coloured"]) - mean(bs_data$died[bs_data$race.w2=="Coloured"])
  
  # Asian/Indian
  bs_pop[b,21] <- mean(int1_data$died[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$died[bs_data$race.w2=="Asian/Indian"])
  bs_pop[b,22] <- mean(int2_data$died[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$died[bs_data$race.w2=="Asian/Indian"])
  bs_pop[b,23] <- mean(int3_data$died[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$died[bs_data$race.w2=="Asian/Indian"])
  bs_pop[b,24] <- mean(int4_data$died[bs_data$race.w2=="Asian/Indian"]) - mean(bs_data$died[bs_data$race.w2=="Asian/Indian"])
  
  # White
  bs_pop[b,25] <- mean(int1_data$died[bs_data$race.w2=="White"]) - mean(bs_data$died[bs_data$race.w2=="White"])
  bs_pop[b,26] <- mean(int2_data$died[bs_data$race.w2=="White"]) - mean(bs_data$died[bs_data$race.w2=="White"])
  bs_pop[b,27] <- mean(int3_data$died[bs_data$race.w2=="White"]) - mean(bs_data$died[bs_data$race.w2=="White"])
  bs_pop[b,28] <- mean(int4_data$died[bs_data$race.w2=="White"]) - mean(bs_data$died[bs_data$race.w2=="White"])
  
  # Urban
  bs_pop[b,29] <- mean(int1_data$died[bs_data$urban.w2=="Urban"]) - mean(bs_data$died[bs_data$urban.w2=="Urban"])
  bs_pop[b,30] <- mean(int2_data$died[bs_data$urban.w2=="Urban"]) - mean(bs_data$died[bs_data$urban.w2=="Urban"])
  bs_pop[b,31] <- mean(int3_data$died[bs_data$urban.w2=="Urban"]) - mean(bs_data$died[bs_data$urban.w2=="Urban"])
  bs_pop[b,32] <- mean(int4_data$died[bs_data$urban.w2=="Urban"]) - mean(bs_data$died[bs_data$urban.w2=="Urban"])
  
  # Traditional
  bs_pop[b,33] <- mean(int1_data$died[bs_data$urban.w2=="Traditional"]) - mean(bs_data$died[bs_data$urban.w2=="Traditional"])
  bs_pop[b,34] <- mean(int2_data$died[bs_data$urban.w2=="Traditional"]) - mean(bs_data$died[bs_data$urban.w2=="Traditional"])
  bs_pop[b,35] <- mean(int3_data$died[bs_data$urban.w2=="Traditional"]) - mean(bs_data$died[bs_data$urban.w2=="Traditional"])
  bs_pop[b,36] <- mean(int4_data$died[bs_data$urban.w2=="Traditional"]) - mean(bs_data$died[bs_data$urban.w2=="Traditional"])
  
  # Urban
  bs_pop[b,37] <- mean(int1_data$died[bs_data$urban.w2=="Farms"]) - mean(bs_data$died[bs_data$urban.w2=="Farms"])
  bs_pop[b,38] <- mean(int2_data$died[bs_data$urban.w2=="Farms"]) - mean(bs_data$died[bs_data$urban.w2=="Farms"])
  bs_pop[b,39] <- mean(int3_data$died[bs_data$urban.w2=="Farms"]) - mean(bs_data$died[bs_data$urban.w2=="Farms"])
  bs_pop[b,40] <- mean(int4_data$died[bs_data$urban.w2=="Farms"]) - mean(bs_data$died[bs_data$urban.w2=="Farms"])
  
#.......................................................................................
# END BOOTSTRAP
#.......................................................................................    

}
    
#.......................................................................................
# P-value estimates
#.......................................................................................
  
  # P-value: Mortality ####
  mort_pvalue <-  tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                         strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                                                  "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
                         mort_pvalue = as.double(NA))
  
    # Mortality 120
    mort_pvalue[1,3] <- (2 * pnorm(abs(mean(bs_pop[,1])/sd(bs_pop[,1])), lower.tail = FALSE))#Overall
    mort_pvalue[2,3] <- (2 * pnorm(abs(mean(bs_pop[,5])/sd(bs_pop[,5])), lower.tail = FALSE))#Men
    mort_pvalue[3,3] <- (2 * pnorm(abs(mean(bs_pop[,9])/sd(bs_pop[,9])), lower.tail = FALSE))#Women
    mort_pvalue[4,3] <- (2 * pnorm(abs(mean(bs_pop[,13])/sd(bs_pop[,13])), lower.tail = FALSE))#African
    mort_pvalue[5,3] <- (2 * pnorm(abs(mean(bs_pop[,17])/sd(bs_pop[,17])), lower.tail = FALSE))#Coloured
    mort_pvalue[6,3] <- (2 * pnorm(abs(mean(bs_pop[,21])/sd(bs_pop[,21])), lower.tail = FALSE))#Asian/Indian
    mort_pvalue[7,3] <- (2 * pnorm(abs(mean(bs_pop[,25])/sd(bs_pop[,25])), lower.tail = FALSE))#White
    mort_pvalue[8,3] <- (2 * pnorm(abs(mean(bs_pop[,29])/sd(bs_pop[,29])), lower.tail = FALSE))#Urban
    mort_pvalue[9,3] <- (2 * pnorm(abs(mean(bs_pop[,33])/sd(bs_pop[,33])), lower.tail = FALSE))#traditional
    mort_pvalue[10,3] <- (2 * pnorm(abs(mean(bs_pop[,37])/sd(bs_pop[,37])), lower.tail = FALSE))#farms
  
    # Mortality 130
    mort_pvalue[11,3] <- (2 * pnorm(abs(mean(bs_pop[,2])/sd(bs_pop[,2])), lower.tail = FALSE))#Overall
    mort_pvalue[12,3] <- (2 * pnorm(abs(mean(bs_pop[,6])/sd(bs_pop[,6])), lower.tail = FALSE))#Men
    mort_pvalue[13,3] <- (2 * pnorm(abs(mean(bs_pop[,10])/sd(bs_pop[,10])), lower.tail = FALSE))#Women
    mort_pvalue[14,3] <- (2 * pnorm(abs(mean(bs_pop[,14])/sd(bs_pop[,14])), lower.tail = FALSE))#African
    mort_pvalue[15,3] <- (2 * pnorm(abs(mean(bs_pop[,18])/sd(bs_pop[,18])), lower.tail = FALSE))#Coloured
    mort_pvalue[16,3] <- (2 * pnorm(abs(mean(bs_pop[,22])/sd(bs_pop[,22])), lower.tail = FALSE))#Asian/Indian
    mort_pvalue[17,3] <- (2 * pnorm(abs(mean(bs_pop[,26])/sd(bs_pop[,26])), lower.tail = FALSE))#White
    mort_pvalue[18,3] <- (2 * pnorm(abs(mean(bs_pop[,30])/sd(bs_pop[,30])), lower.tail = FALSE))#Urban
    mort_pvalue[19,3] <- (2 * pnorm(abs(mean(bs_pop[,34])/sd(bs_pop[,34])), lower.tail = FALSE))#traditional
    mort_pvalue[20,3] <- (2 * pnorm(abs(mean(bs_pop[,38])/sd(bs_pop[,38])), lower.tail = FALSE))#farms
  
    # Mortality 140
    mort_pvalue[21,3] <- (2 * pnorm(abs(mean(bs_pop[,3])/sd(bs_pop[,3])), lower.tail = FALSE))#Overall
    mort_pvalue[22,3] <- (2 * pnorm(abs(mean(bs_pop[,7])/sd(bs_pop[,7])), lower.tail = FALSE))#Men
    mort_pvalue[23,3] <- (2 * pnorm(abs(mean(bs_pop[,11])/sd(bs_pop[,11])), lower.tail = FALSE))#Women
    mort_pvalue[24,3] <- (2 * pnorm(abs(mean(bs_pop[,15])/sd(bs_pop[,15])), lower.tail = FALSE))#African
    mort_pvalue[25,3] <- (2 * pnorm(abs(mean(bs_pop[,19])/sd(bs_pop[,19])), lower.tail = FALSE))#Coloured
    mort_pvalue[26,3] <- (2 * pnorm(abs(mean(bs_pop[,23])/sd(bs_pop[,23])), lower.tail = FALSE))#Asian/Indian
    mort_pvalue[27,3] <- (2 * pnorm(abs(mean(bs_pop[,27])/sd(bs_pop[,27])), lower.tail = FALSE))#White
    mort_pvalue[28,3] <- (2 * pnorm(abs(mean(bs_pop[,31])/sd(bs_pop[,31])), lower.tail = FALSE))#Urban
    mort_pvalue[29,3] <- (2 * pnorm(abs(mean(bs_pop[,35])/sd(bs_pop[,35])), lower.tail = FALSE))#traditional
    mort_pvalue[30,3] <- (2 * pnorm(abs(mean(bs_pop[,39])/sd(bs_pop[,39])), lower.tail = FALSE))#farms
  
    # Mortality 150
    mort_pvalue[31,3] <- (2 * pnorm(abs(mean(bs_pop[,4])/sd(bs_pop[,4])), lower.tail = FALSE))#Overall
    mort_pvalue[32,3] <- (2 * pnorm(abs(mean(bs_pop[,8])/sd(bs_pop[,8])), lower.tail = FALSE))#Men
    mort_pvalue[33,3] <- (2 * pnorm(abs(mean(bs_pop[,12])/sd(bs_pop[,12])), lower.tail = FALSE))#Women
    mort_pvalue[34,3] <- (2 * pnorm(abs(mean(bs_pop[,16])/sd(bs_pop[,16])), lower.tail = FALSE))#African
    mort_pvalue[35,3] <- (2 * pnorm(abs(mean(bs_pop[,20])/sd(bs_pop[,20])), lower.tail = FALSE))#Coloured
    mort_pvalue[36,3] <- (2 * pnorm(abs(mean(bs_pop[,24])/sd(bs_pop[,24])), lower.tail = FALSE))#Asian/Indian
    mort_pvalue[37,3] <- (2 * pnorm(abs(mean(bs_pop[,28])/sd(bs_pop[,28])), lower.tail = FALSE))#White
    mort_pvalue[38,3] <- (2 * pnorm(abs(mean(bs_pop[,32])/sd(bs_pop[,32])), lower.tail = FALSE))#Urban
    mort_pvalue[39,3] <- (2 * pnorm(abs(mean(bs_pop[,36])/sd(bs_pop[,36])), lower.tail = FALSE))#traditional
    mort_pvalue[40,3] <- (2 * pnorm(abs(mean(bs_pop[,40])/sd(bs_pop[,40])), lower.tail = FALSE))#farms
    
  # P-value: Share ####
  share_pvalue <-  tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                             "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
              share_pvalue = as.double(NA))
  
    # Share 120
    share_pvalue[1,3] <- (2 * pnorm(abs(mean(bs_share[,1])/sd(bs_share[,1])), lower.tail = FALSE))#Overall
    share_pvalue[2,3] <- (2 * pnorm(abs(mean(bs_share[,5])/sd(bs_share[,5])), lower.tail = FALSE))#Men
    share_pvalue[3,3] <- (2 * pnorm(abs(mean(bs_share[,9])/sd(bs_share[,9])), lower.tail = FALSE))#Women
    share_pvalue[4,3] <- (2 * pnorm(abs(mean(bs_share[,13])/sd(bs_share[,13])), lower.tail = FALSE))#African
    share_pvalue[5,3] <- (2 * pnorm(abs(mean(bs_share[,17])/sd(bs_share[,17])), lower.tail = FALSE))#Coloured
    share_pvalue[6,3] <- (2 * pnorm(abs(mean(bs_share[,21])/sd(bs_share[,21])), lower.tail = FALSE))#Asian/Indian
    share_pvalue[7,3] <- (2 * pnorm(abs(mean(bs_share[,25])/sd(bs_share[,25])), lower.tail = FALSE))#White
    share_pvalue[8,3] <- (2 * pnorm(abs(mean(bs_share[,29])/sd(bs_share[,29])), lower.tail = FALSE))#Urban
    share_pvalue[9,3] <- (2 * pnorm(abs(mean(bs_share[,33])/sd(bs_share[,33])), lower.tail = FALSE))#traditional
    share_pvalue[10,3] <- (2 * pnorm(abs(mean(bs_share[,37])/sd(bs_share[,37])), lower.tail = FALSE))#farms
    
    # Share 130
    share_pvalue[11,3] <- (2 * pnorm(abs(mean(bs_share[,2])/sd(bs_share[,2])), lower.tail = FALSE))#Overall
    share_pvalue[12,3] <- (2 * pnorm(abs(mean(bs_share[,6])/sd(bs_share[,6])), lower.tail = FALSE))#Men
    share_pvalue[13,3] <- (2 * pnorm(abs(mean(bs_share[,10])/sd(bs_share[,10])), lower.tail = FALSE))#Women
    share_pvalue[14,3] <- (2 * pnorm(abs(mean(bs_share[,14])/sd(bs_share[,14])), lower.tail = FALSE))#African
    share_pvalue[15,3] <- (2 * pnorm(abs(mean(bs_share[,18])/sd(bs_share[,18])), lower.tail = FALSE))#Coloured
    share_pvalue[16,3] <- (2 * pnorm(abs(mean(bs_share[,22])/sd(bs_share[,22])), lower.tail = FALSE))#Asian/Indian
    share_pvalue[17,3] <- (2 * pnorm(abs(mean(bs_share[,26])/sd(bs_share[,26])), lower.tail = FALSE))#White
    share_pvalue[18,3] <- (2 * pnorm(abs(mean(bs_share[,30])/sd(bs_share[,30])), lower.tail = FALSE))#Urban
    share_pvalue[19,3] <- (2 * pnorm(abs(mean(bs_share[,34])/sd(bs_share[,34])), lower.tail = FALSE))#traditional
    share_pvalue[20,3] <- (2 * pnorm(abs(mean(bs_share[,38])/sd(bs_share[,38])), lower.tail = FALSE))#farms
    
    # Share 140
    share_pvalue[21,3] <- (2 * pnorm(abs(mean(bs_share[,3])/sd(bs_share[,3])), lower.tail = FALSE))#Overall
    share_pvalue[22,3] <- (2 * pnorm(abs(mean(bs_share[,7])/sd(bs_share[,7])), lower.tail = FALSE))#Men
    share_pvalue[23,3] <- (2 * pnorm(abs(mean(bs_share[,11])/sd(bs_share[,11])), lower.tail = FALSE))#Women
    share_pvalue[24,3] <- (2 * pnorm(abs(mean(bs_share[,15])/sd(bs_share[,15])), lower.tail = FALSE))#African
    share_pvalue[25,3] <- (2 * pnorm(abs(mean(bs_share[,19])/sd(bs_share[,19])), lower.tail = FALSE))#Coloured
    share_pvalue[26,3] <- (2 * pnorm(abs(mean(bs_share[,23])/sd(bs_share[,23])), lower.tail = FALSE))#Asian/Indian
    share_pvalue[27,3] <- (2 * pnorm(abs(mean(bs_share[,27])/sd(bs_share[,27])), lower.tail = FALSE))#White
    share_pvalue[28,3] <- (2 * pnorm(abs(mean(bs_share[,31])/sd(bs_share[,31])), lower.tail = FALSE))#Urban
    share_pvalue[29,3] <- (2 * pnorm(abs(mean(bs_share[,35])/sd(bs_share[,35])), lower.tail = FALSE))#traditional
    share_pvalue[30,3] <- (2 * pnorm(abs(mean(bs_share[,39])/sd(bs_share[,39])), lower.tail = FALSE))#farms
    
    # Share 150
    share_pvalue[31,3] <- (2 * pnorm(abs(mean(bs_share[,4])/sd(bs_share[,4])), lower.tail = FALSE))#Overall
    share_pvalue[32,3] <- (2 * pnorm(abs(mean(bs_share[,8])/sd(bs_share[,8])), lower.tail = FALSE))#Men
    share_pvalue[33,3] <- (2 * pnorm(abs(mean(bs_share[,12])/sd(bs_share[,12])), lower.tail = FALSE))#Women
    share_pvalue[34,3] <- (2 * pnorm(abs(mean(bs_share[,16])/sd(bs_share[,16])), lower.tail = FALSE))#African
    share_pvalue[35,3] <- (2 * pnorm(abs(mean(bs_share[,20])/sd(bs_share[,20])), lower.tail = FALSE))#Coloured
    share_pvalue[36,3] <- (2 * pnorm(abs(mean(bs_share[,24])/sd(bs_share[,24])), lower.tail = FALSE))#Asian/Indian
    share_pvalue[37,3] <- (2 * pnorm(abs(mean(bs_share[,28])/sd(bs_share[,28])), lower.tail = FALSE))#White
    share_pvalue[38,3] <- (2 * pnorm(abs(mean(bs_share[,32])/sd(bs_share[,32])), lower.tail = FALSE))#Urban
    share_pvalue[39,3] <- (2 * pnorm(abs(mean(bs_share[,36])/sd(bs_share[,36])), lower.tail = FALSE))#traditional
    share_pvalue[40,3] <- (2 * pnorm(abs(mean(bs_share[,40])/sd(bs_share[,40])), lower.tail = FALSE))#farms
  
  # P-value: shift ####
  shift_pvalue <-  tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                             "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
              shift_pvalue = as.double(NA))
  
    # shift 120
    shift_pvalue[1,3] <- (2 * pnorm(abs(mean(bs_shift[,1])/sd(bs_shift[,1])), lower.tail = FALSE))#Overall
    shift_pvalue[2,3] <- (2 * pnorm(abs(mean(bs_shift[,5])/sd(bs_shift[,5])), lower.tail = FALSE))#Men
    shift_pvalue[3,3] <- (2 * pnorm(abs(mean(bs_shift[,9])/sd(bs_shift[,9])), lower.tail = FALSE))#Women
    shift_pvalue[4,3] <- (2 * pnorm(abs(mean(bs_shift[,13])/sd(bs_shift[,13])), lower.tail = FALSE))#African
    shift_pvalue[5,3] <- (2 * pnorm(abs(mean(bs_shift[,17])/sd(bs_shift[,17])), lower.tail = FALSE))#Coloured
    shift_pvalue[6,3] <- (2 * pnorm(abs(mean(bs_shift[,21])/sd(bs_shift[,21])), lower.tail = FALSE))#Asian/Indian
    shift_pvalue[7,3] <- (2 * pnorm(abs(mean(bs_shift[,25])/sd(bs_shift[,25])), lower.tail = FALSE))#White
    shift_pvalue[8,3] <- (2 * pnorm(abs(mean(bs_shift[,29])/sd(bs_shift[,29])), lower.tail = FALSE))#Urban
    shift_pvalue[9,3] <- (2 * pnorm(abs(mean(bs_shift[,33])/sd(bs_shift[,33])), lower.tail = FALSE))#traditional
    shift_pvalue[10,3] <- (2 * pnorm(abs(mean(bs_shift[,37])/sd(bs_shift[,37])), lower.tail = FALSE))#farms
    
    # shift 130
    shift_pvalue[11,3] <- (2 * pnorm(abs(mean(bs_shift[,2])/sd(bs_shift[,2])), lower.tail = FALSE))#Overall
    shift_pvalue[12,3] <- (2 * pnorm(abs(mean(bs_shift[,6])/sd(bs_shift[,6])), lower.tail = FALSE))#Men
    shift_pvalue[13,3] <- (2 * pnorm(abs(mean(bs_shift[,10])/sd(bs_shift[,10])), lower.tail = FALSE))#Women
    shift_pvalue[14,3] <- (2 * pnorm(abs(mean(bs_shift[,14])/sd(bs_shift[,14])), lower.tail = FALSE))#African
    shift_pvalue[15,3] <- (2 * pnorm(abs(mean(bs_shift[,18])/sd(bs_shift[,18])), lower.tail = FALSE))#Coloured
    shift_pvalue[16,3] <- (2 * pnorm(abs(mean(bs_shift[,22])/sd(bs_shift[,22])), lower.tail = FALSE))#Asian/Indian
    shift_pvalue[17,3] <- (2 * pnorm(abs(mean(bs_shift[,26])/sd(bs_shift[,26])), lower.tail = FALSE))#White
    shift_pvalue[18,3] <- (2 * pnorm(abs(mean(bs_shift[,30])/sd(bs_shift[,30])), lower.tail = FALSE))#Urban
    shift_pvalue[19,3] <- (2 * pnorm(abs(mean(bs_shift[,34])/sd(bs_shift[,34])), lower.tail = FALSE))#traditional
    shift_pvalue[20,3] <- (2 * pnorm(abs(mean(bs_shift[,38])/sd(bs_shift[,38])), lower.tail = FALSE))#farms
    
    # shift 140
    shift_pvalue[21,3] <- (2 * pnorm(abs(mean(bs_shift[,3])/sd(bs_shift[,3])), lower.tail = FALSE))#Overall
    shift_pvalue[22,3] <- (2 * pnorm(abs(mean(bs_shift[,7])/sd(bs_shift[,7])), lower.tail = FALSE))#Men
    shift_pvalue[23,3] <- (2 * pnorm(abs(mean(bs_shift[,11])/sd(bs_shift[,11])), lower.tail = FALSE))#Women
    shift_pvalue[24,3] <- (2 * pnorm(abs(mean(bs_shift[,15])/sd(bs_shift[,15])), lower.tail = FALSE))#African
    shift_pvalue[25,3] <- (2 * pnorm(abs(mean(bs_shift[,19])/sd(bs_shift[,19])), lower.tail = FALSE))#Coloured
    shift_pvalue[26,3] <- (2 * pnorm(abs(mean(bs_shift[,23])/sd(bs_shift[,23])), lower.tail = FALSE))#Asian/Indian
    shift_pvalue[27,3] <- (2 * pnorm(abs(mean(bs_shift[,27])/sd(bs_shift[,27])), lower.tail = FALSE))#White
    shift_pvalue[28,3] <- (2 * pnorm(abs(mean(bs_shift[,31])/sd(bs_shift[,31])), lower.tail = FALSE))#Urban
    shift_pvalue[29,3] <- (2 * pnorm(abs(mean(bs_shift[,35])/sd(bs_shift[,35])), lower.tail = FALSE))#traditional
    shift_pvalue[30,3] <- (2 * pnorm(abs(mean(bs_shift[,39])/sd(bs_shift[,39])), lower.tail = FALSE))#farms
    
    # shift 150
    shift_pvalue[31,3] <- (2 * pnorm(abs(mean(bs_shift[,4])/sd(bs_shift[,4])), lower.tail = FALSE))#Overall
    shift_pvalue[32,3] <- (2 * pnorm(abs(mean(bs_shift[,8])/sd(bs_shift[,8])), lower.tail = FALSE))#Men
    shift_pvalue[33,3] <- (2 * pnorm(abs(mean(bs_shift[,12])/sd(bs_shift[,12])), lower.tail = FALSE))#Women
    shift_pvalue[34,3] <- (2 * pnorm(abs(mean(bs_shift[,16])/sd(bs_shift[,16])), lower.tail = FALSE))#African
    shift_pvalue[35,3] <- (2 * pnorm(abs(mean(bs_shift[,20])/sd(bs_shift[,20])), lower.tail = FALSE))#Coloured
    shift_pvalue[36,3] <- (2 * pnorm(abs(mean(bs_shift[,24])/sd(bs_shift[,24])), lower.tail = FALSE))#Asian/Indian
    shift_pvalue[37,3] <- (2 * pnorm(abs(mean(bs_shift[,28])/sd(bs_shift[,28])), lower.tail = FALSE))#White
    shift_pvalue[38,3] <- (2 * pnorm(abs(mean(bs_shift[,32])/sd(bs_shift[,32])), lower.tail = FALSE))#Urban
    shift_pvalue[39,3] <- (2 * pnorm(abs(mean(bs_shift[,36])/sd(bs_shift[,36])), lower.tail = FALSE))#traditional
    shift_pvalue[40,3] <- (2 * pnorm(abs(mean(bs_shift[,40])/sd(bs_shift[,40])), lower.tail = FALSE))#farms  
  
  mort_pvalue1 <- mort_pvalue %>% 
      select(target, strata, mort_pvalue) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "mort",
                  names_sep =  "_",
                  values_from = mort_pvalue) 
  
  share_pvalue1 <- share_pvalue %>% 
      select(target, strata, share_pvalue) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "share",
                  names_sep =  "_",
                  values_from = share_pvalue) 
  
  shift_pvalue1 <- shift_pvalue %>% 
      select(target, strata, shift_pvalue) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "shift",
                  names_sep =  "_",
                  values_from = shift_pvalue) 
          
    pvalue_table <- mort_pvalue1 %>% 
    left_join(share_pvalue1, by= c("strata")) %>%
    left_join(shift_pvalue1, by= c("strata")) 
    write.csv(pvalue_table, "pop_pvalue.csv")   
      
  # Ratio estimates pvalues
  est_pvalue <- tibble(start_bp = c(160,160,160,160,160,150,150,150,150,140,140,140,130,130),
              end_bp = c(160,150,140,130,120,150,140,130,120,140,130,120,130,120),
              pvalue = as.double(NA))
  est_pvalue[1,3] <- 1
  est_pvalue[2,3] <- (2 * pnorm(abs((log(mean(bs_150/bs_160)))/(sd(bs_150/bs_160))), lower.tail = FALSE))
  est_pvalue[3,3] <- (2 * pnorm(abs((log(mean(bs_140/bs_160)))/(sd(bs_140/bs_160))), lower.tail = FALSE))
  est_pvalue[4,3] <- (2 * pnorm(abs((log(mean(bs_130/bs_160)))/(sd(bs_130/bs_160))), lower.tail = FALSE))
  est_pvalue[5,3] <- (2 * pnorm(abs((log(mean(bs_120/bs_160)))/(sd(bs_120/bs_160))), lower.tail = FALSE))
  est_pvalue[6,3] <- 1
  est_pvalue[7,3] <- (2 * pnorm(abs((log(mean(bs_140/bs_150)))/(sd(bs_140/bs_150))), lower.tail = FALSE))
  est_pvalue[8,3] <- (2 * pnorm(abs((log(mean(bs_130/bs_150)))/(sd(bs_130/bs_150))), lower.tail = FALSE))
  est_pvalue[9,3] <- (2 * pnorm(abs((log(mean(bs_120/bs_150)))/(sd(bs_120/bs_150))), lower.tail = FALSE))
  est_pvalue[10,3] <- 1
  est_pvalue[11,3] <- (2 * pnorm(abs((log(mean(bs_130/bs_140)))/(sd(bs_130/bs_140))), lower.tail = FALSE))
  est_pvalue[12,3] <- (2 * pnorm(abs((log(mean(bs_120/bs_140)))/(sd(bs_120/bs_140))), lower.tail = FALSE))
  est_pvalue[13,3] <- 1
  est_pvalue[12,3] <- (2 * pnorm(abs((log(mean(bs_120/bs_130)))/(sd(bs_120/bs_130))), lower.tail = FALSE))
  est_pvalue[14,3] <- 1
  
  write.csv(est_pvalue, "clinical_pvalue.csv")  
          
#.......................................................................................
# Clinical graph
#.......................................................................................    

  clinical <- tibble(start_bp = c(160,160,160,160,160,150,150,150,150,140,140,140,130,130),
                    end_bp = c(160,150,140,130,120,150,140,130,120,140,130,120,130,120),
                    est = as.double(NA),
                    lcl = as.double(NA),
                    ucl = as.double(NA))
 
  
  clinical[1,3] <- 1
  clinical[2,3] <- mean(bs_150/bs_160)
  clinical[2,4] <- mean(bs_150/bs_160) - 1.96*sd(bs_150/bs_160)
  clinical[2,5] <- mean(bs_150/bs_160) + 1.96*sd(bs_150/bs_160)
  
  clinical[3,3] <- mean(bs_140/bs_160)
  clinical[3,4] <- mean(bs_140/bs_160) - 1.96*sd(bs_140/bs_160)
  clinical[3,5] <- mean(bs_140/bs_160) + 1.96*sd(bs_140/bs_160)
  
  clinical[4,3] <- mean(bs_130/bs_160)
  clinical[4,4] <- mean(bs_130/bs_160) - 1.96*sd(bs_130/bs_160)
  clinical[4,5] <- mean(bs_130/bs_160) + 1.96*sd(bs_130/bs_160)
  
  clinical[5,3] <- mean(bs_120/bs_160)
  clinical[5,4] <- mean(bs_120/bs_160) - 1.96*sd(bs_120/bs_160)
  clinical[5,5] <- mean(bs_120/bs_160) + 1.96*sd(bs_120/bs_160)
  
  clinical[6,3] <- 1
  
  clinical[7,3] <- mean(bs_140/bs_150)
  clinical[7,4] <- mean(bs_140/bs_150) - 1.96*sd(bs_140/bs_150)
  clinical[7,5] <- mean(bs_140/bs_150) + 1.96*sd(bs_140/bs_150)
  
  clinical[8,3] <- mean(bs_130/bs_150)
  clinical[8,4] <- mean(bs_130/bs_150) - 1.96*sd(bs_130/bs_150)
  clinical[8,5] <- mean(bs_130/bs_150) + 1.96*sd(bs_130/bs_150)
  
  clinical[9,3] <- mean(bs_120/bs_150)
  clinical[9,4] <- mean(bs_120/bs_150) - 1.96*sd(bs_120/bs_150)
  clinical[9,5] <- mean(bs_120/bs_150) + 1.96*sd(bs_120/bs_150)

  clinical[10,3] <- 1
  
  clinical[11,3] <- mean(bs_130/bs_140)
  clinical[11,4] <- mean(bs_130/bs_140) - 1.96*sd(bs_130/bs_140)
  clinical[11,5] <- mean(bs_130/bs_140) + 1.96*sd(bs_130/bs_140)
  
  clinical[12,3] <- mean(bs_120/bs_140)
  clinical[12,4] <- mean(bs_120/bs_140) - 1.96*sd(bs_120/bs_140)
  clinical[12,5] <- mean(bs_120/bs_140) + 1.96*sd(bs_120/bs_140)
  
  clinical[13,3] <- 1
  
  clinical[14,3] <- mean(bs_120/bs_130)
  clinical[14,4] <- mean(bs_120/bs_130) - 1.96*sd(bs_120/bs_130)
  clinical[14,5] <- mean(bs_120/bs_130) + 1.96*sd(bs_120/bs_130)
  
  clinical %>% write.table("clipboard", sep="\t", row.names = F)
  
  ggplot(data = clinical, aes(x = end_bp,y=est,ymin=lcl,ymax=ucl)) +
    geom_point() +
    geom_errorbar(width = 0.5) +
    geom_hline(yintercept = 1,linetype = "dashed") +
    scale_x_continuous(limits = c(120,160),breaks = seq(120,160,10)) +
    scale_y_continuous(limits = c(0.7,1.1),breaks = seq(0.7,1.1,0.1)) +
    labs(y = "Risk ratio", x = "Ending SBP level (mmHg)") +
    facet_wrap(~fct_rev(as_factor(start_bp))) + theme_classic() + coord_flip() +
    theme(text = element_text(size = 14))
  
  ggsave("figure3.jpg",units = "in",height = 7,width = 5)
  
#.......................................................................................
# Pop-level estimates by strata
#....................................................................................... 
      
  #Mortality reduction####
  mort_results <- tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                          strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                                       "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
                        est = as.double(NA),
                        lcl = as.double(NA),
                        ucl = as.double(NA))
  
    #Fill in: Overall

      #120
      mort_results[1,3] <- mean(bs_pop[,1])
      mort_results[1,4] <- mean(bs_pop[,1]) - 1.96*sd(bs_pop[,1])
      mort_results[1,5] <- mean(bs_pop[,1]) + 1.96*sd(bs_pop[,1])
      
      
      #130
      mort_results[11,3] <- mean(bs_pop[,2])
      mort_results[11,4] <- mean(bs_pop[,2]) - 1.96*sd(bs_pop[,2])
      mort_results[11,5] <- mean(bs_pop[,2]) + 1.96*sd(bs_pop[,2])
      
      #140
      mort_results[21,3] <- mean(bs_pop[,3])
      mort_results[21,4] <- mean(bs_pop[,3]) - 1.96*sd(bs_pop[,3])
      mort_results[21,5] <- mean(bs_pop[,3]) + 1.96*sd(bs_pop[,3])
      
      #150
      mort_results[31,3] <- mean(bs_pop[,4])
      mort_results[31,4] <- mean(bs_pop[,4]) - 1.96*sd(bs_pop[,4])
      mort_results[31,5] <- mean(bs_pop[,4]) + 1.96*sd(bs_pop[,4])
      
    #Fill in: Men
    
      #120
      mort_results[2,3] <- mean(bs_pop[,5])
      mort_results[2,4] <- mean(bs_pop[,5]) - 1.96*sd(bs_pop[,5])
      mort_results[2,5] <- mean(bs_pop[,5]) + 1.96*sd(bs_pop[,5])
      
      #130
      mort_results[12,3] <- mean(bs_pop[,6])
      mort_results[12,4] <- mean(bs_pop[,6]) - 1.96*sd(bs_pop[,6])
      mort_results[12,5] <- mean(bs_pop[,6]) + 1.96*sd(bs_pop[,6])
      
      #140
      mort_results[22,3] <- mean(bs_pop[,7])
      mort_results[22,4] <- mean(bs_pop[,7]) - 1.96*sd(bs_pop[,7])
      mort_results[22,5] <- mean(bs_pop[,7]) + 1.96*sd(bs_pop[,7])
      
      #150
      mort_results[32,3] <- mean(bs_pop[,8])
      mort_results[32,4] <- mean(bs_pop[,8]) - 1.96*sd(bs_pop[,8])
      mort_results[32,5] <- mean(bs_pop[,8]) + 1.96*sd(bs_pop[,8])
      
    #Fill in: Women
     
      #120
      mort_results[3,3] <- mean(bs_pop[,9])
      mort_results[3,4] <- mean(bs_pop[,9]) - 1.96*sd(bs_pop[,9])
      mort_results[3,5] <- mean(bs_pop[,9]) + 1.96*sd(bs_pop[,9])
      
      #130
      mort_results[13,3] <- mean(bs_pop[,10])
      mort_results[13,4] <- mean(bs_pop[,10]) - 1.96*sd(bs_pop[,10])
      mort_results[13,5] <- mean(bs_pop[,10]) + 1.96*sd(bs_pop[,10])
      
      #140
      mort_results[23,3] <- mean(bs_pop[,11])
      mort_results[23,4] <- mean(bs_pop[,11]) - 1.96*sd(bs_pop[,11])
      mort_results[23,5] <- mean(bs_pop[,11]) + 1.96*sd(bs_pop[,11])
      
      #150
      mort_results[33,3] <- mean(bs_pop[,12])
      mort_results[33,4] <- mean(bs_pop[,12]) - 1.96*sd(bs_pop[,12])
      mort_results[33,5] <- mean(bs_pop[,12]) + 1.96*sd(bs_pop[,12])
      
    #Fill in: African
     
      #120
      mort_results[4,3] <- mean(bs_pop[,13])
      mort_results[4,4] <- mean(bs_pop[,13]) - 1.96*sd(bs_pop[,13])
      mort_results[4,5] <- mean(bs_pop[,13]) + 1.96*sd(bs_pop[,13])
      
      #130
      mort_results[14,3] <- mean(bs_pop[,14])
      mort_results[14,4] <- mean(bs_pop[,14]) - 1.96*sd(bs_pop[,14])
      mort_results[14,5] <- mean(bs_pop[,14]) + 1.96*sd(bs_pop[,14])
      
      #140
      mort_results[24,3] <- mean(bs_pop[,15])
      mort_results[24,4] <- mean(bs_pop[,15]) - 1.96*sd(bs_pop[,15])
      mort_results[24,5] <- mean(bs_pop[,15]) + 1.96*sd(bs_pop[,15])
      
      #150
      mort_results[34,3] <- mean(bs_pop[,16])
      mort_results[34,4] <- mean(bs_pop[,16]) - 1.96*sd(bs_pop[,16])
      mort_results[34,5] <- mean(bs_pop[,16]) + 1.96*sd(bs_pop[,16])
      
    #Fill in: Coloured
     
      #120
      mort_results[5,3] <- mean(bs_pop[,17])
      mort_results[5,4] <- mean(bs_pop[,17]) - 1.96*sd(bs_pop[,17])
      mort_results[5,5] <- mean(bs_pop[,17]) + 1.96*sd(bs_pop[,17])
      
      #130
      mort_results[15,3] <- mean(bs_pop[,18])
      mort_results[15,4] <- mean(bs_pop[,18]) - 1.96*sd(bs_pop[,18])
      mort_results[15,5] <- mean(bs_pop[,18]) + 1.96*sd(bs_pop[,18])
      
      #140
      mort_results[25,3] <- mean(bs_pop[,19])
      mort_results[25,4] <- mean(bs_pop[,19]) - 1.96*sd(bs_pop[,19])
      mort_results[25,5] <- mean(bs_pop[,19]) + 1.96*sd(bs_pop[,19])
      
      #150
      mort_results[35,3] <- mean(bs_pop[,20])
      mort_results[35,4] <- mean(bs_pop[,20]) - 1.96*sd(bs_pop[,20])
      mort_results[35,5] <- mean(bs_pop[,20]) + 1.96*sd(bs_pop[,20])
      
    #Fill in: Asian/Indian
     
      #120
      mort_results[6,3] <- mean(bs_pop[,21])
      mort_results[6,4] <- mean(bs_pop[,21]) - 1.96*sd(bs_pop[,21])
      mort_results[6,5] <- mean(bs_pop[,21]) + 1.96*sd(bs_pop[,21])
      
      #130
      mort_results[16,3] <- mean(bs_pop[,22])
      mort_results[16,4] <- mean(bs_pop[,22]) - 1.96*sd(bs_pop[,22])
      mort_results[16,5] <- mean(bs_pop[,22]) + 1.96*sd(bs_pop[,22])
      
      #140
      mort_results[26,3] <- mean(bs_pop[,23])
      mort_results[26,4] <- mean(bs_pop[,23]) - 1.96*sd(bs_pop[,23])
      mort_results[26,5] <- mean(bs_pop[,23]) + 1.96*sd(bs_pop[,23])
      
      #150
      mort_results[36,3] <- mean(bs_pop[,24])
      mort_results[36,4] <- mean(bs_pop[,24]) - 1.96*sd(bs_pop[,24])
      mort_results[36,5] <- mean(bs_pop[,24]) + 1.96*sd(bs_pop[,24])
      
    #Fill in: White
      
      #120
      mort_results[7,3] <- mean(bs_pop[,25])
      mort_results[7,4] <- mean(bs_pop[,25]) - 1.96*sd(bs_pop[,25])
      mort_results[7,5] <- mean(bs_pop[,25]) + 1.96*sd(bs_pop[,25])
      
      #130
      mort_results[17,3] <- mean(bs_pop[,26])
      mort_results[17,4] <- mean(bs_pop[,26]) - 1.96*sd(bs_pop[,26])
      mort_results[17,5] <- mean(bs_pop[,26]) + 1.96*sd(bs_pop[,26])
      
      #140
      mort_results[27,3] <- mean(bs_pop[,27])
      mort_results[27,4] <- mean(bs_pop[,27]) - 1.96*sd(bs_pop[,27])
      mort_results[27,5] <- mean(bs_pop[,27]) + 1.96*sd(bs_pop[,27])
      
      #150
      mort_results[37,3] <- mean(bs_pop[,28])
      mort_results[37,4] <- mean(bs_pop[,28]) - 1.96*sd(bs_pop[,28])
      mort_results[37,5] <- mean(bs_pop[,28]) + 1.96*sd(bs_pop[,28])
      
    #Fill in: Urban
     
      #120
      mort_results[8,3] <- mean(bs_pop[,29])
      mort_results[8,4] <- mean(bs_pop[,29]) - 1.96*sd(bs_pop[,29])
      mort_results[8,5] <- mean(bs_pop[,29]) + 1.96*sd(bs_pop[,29])
      
      #130
      mort_results[18,3] <- mean(bs_pop[,30])
      mort_results[18,4] <- mean(bs_pop[,30]) - 1.96*sd(bs_pop[,30])
      mort_results[18,5] <- mean(bs_pop[,30]) + 1.96*sd(bs_pop[,30])
      
      #140
      mort_results[28,3] <- mean(bs_pop[,31])
      mort_results[28,4] <- mean(bs_pop[,31]) - 1.96*sd(bs_pop[,31])
      mort_results[28,5] <- mean(bs_pop[,31]) + 1.96*sd(bs_pop[,31])
      
      #150
      mort_results[38,3] <- mean(bs_pop[,32])
      mort_results[38,4] <- mean(bs_pop[,32]) - 1.96*sd(bs_pop[,32])
      mort_results[38,5] <- mean(bs_pop[,32]) + 1.96*sd(bs_pop[,32])
      
    #Fill in: Traditional
     
      #120
      mort_results[9,3] <- mean(bs_pop[,33])
      mort_results[9,4] <- mean(bs_pop[,33]) - 1.96*sd(bs_pop[,33])
      mort_results[9,5] <- mean(bs_pop[,33]) + 1.96*sd(bs_pop[,33])
      
      #130
      mort_results[19,3] <- mean(bs_pop[,34])
      mort_results[19,4] <- mean(bs_pop[,34]) - 1.96*sd(bs_pop[,34])
      mort_results[19,5] <- mean(bs_pop[,34]) + 1.96*sd(bs_pop[,34])
      
      #140
      mort_results[29,3] <- mean(bs_pop[,35])
      mort_results[29,4] <- mean(bs_pop[,35]) - 1.96*sd(bs_pop[,35])
      mort_results[29,5] <- mean(bs_pop[,35]) + 1.96*sd(bs_pop[,35])
      
      #150
      mort_results[39,3] <- mean(bs_pop[,36])
      mort_results[39,4] <- mean(bs_pop[,36]) - 1.96*sd(bs_pop[,36])
      mort_results[39,5] <- mean(bs_pop[,36]) + 1.96*sd(bs_pop[,36]) 
      
    #Fill in: Farms
     
      #120
      mort_results[10,3] <- mean(bs_pop[,37])
      mort_results[10,4] <- mean(bs_pop[,37]) - 1.96*sd(bs_pop[,37])
      mort_results[10,5] <- mean(bs_pop[,37]) + 1.96*sd(bs_pop[,37])
      
      #130
      mort_results[20,3] <- mean(bs_pop[,38])
      mort_results[20,4] <- mean(bs_pop[,38]) - 1.96*sd(bs_pop[,38])
      mort_results[20,5] <- mean(bs_pop[,38]) + 1.96*sd(bs_pop[,38])
      
      #140
      mort_results[30,3] <- mean(bs_pop[,39])
      mort_results[30,4] <- mean(bs_pop[,39]) - 1.96*sd(bs_pop[,39])
      mort_results[30,5] <- mean(bs_pop[,39]) + 1.96*sd(bs_pop[,39])
      
      #150
      mort_results[40,3] <- mean(bs_pop[,40])
      mort_results[40,4] <- mean(bs_pop[,40]) - 1.96*sd(bs_pop[,40])
      mort_results[40,5] <- mean(bs_pop[,40]) + 1.96*sd(bs_pop[,40]) 
   
    mort_summary <- mort_results %>% 
      mutate(est = format(round(est*1000,1), nsmall = 1),
             lcl = format(round(lcl*1000,1), nsmall = 1),
             ucl = format(round(ucl*1000, 1), nsmall = 1),
             stest1 = str_c(est, lcl, sep = " ("),
             stest2 = str_c(stest1, ucl, sep = ", "),
             estimate = paste0(stest2, ")")) %>% 
      select(target, strata, estimate) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "mort",
                  names_sep =  "_",
                  values_from = estimate)
    
  #Share treated
  share_results <- tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                          strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                                       "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
                        est = as.double(NA),
                        lcl = as.double(NA),
                        ucl = as.double(NA))
  
    #Fill in: Overall
      
      #120
      share_results[1,3] <- mean(bs_share[,1])
      share_results[1,4] <- mean(bs_share[,1]) - 1.96*sd(bs_share[,1])
      share_results[1,5] <- mean(bs_share[,1]) + 1.96*sd(bs_share[,1])
      
      #130
      share_results[11,3] <- mean(bs_share[,2])
      share_results[11,4] <- mean(bs_share[,2]) - 1.96*sd(bs_share[,2])
      share_results[11,5] <- mean(bs_share[,2]) + 1.96*sd(bs_share[,2])
      
      #140
      share_results[21,3] <- mean(bs_share[,3])
      share_results[21,4] <- mean(bs_share[,3]) - 1.96*sd(bs_share[,3])
      share_results[21,5] <- mean(bs_share[,3]) + 1.96*sd(bs_share[,3])
      
      #150
      share_results[31,3] <- mean(bs_share[,4])
      share_results[31,4] <- mean(bs_share[,4]) - 1.96*sd(bs_share[,4])
      share_results[31,5] <- mean(bs_share[,4]) + 1.96*sd(bs_share[,4])
      
    #Fill in: Men
     
      #120
      share_results[2,3] <- mean(bs_share[,5])
      share_results[2,4] <- mean(bs_share[,5]) - 1.96*sd(bs_share[,5])
      share_results[2,5] <- mean(bs_share[,5]) + 1.96*sd(bs_share[,5])
      
      #130
      share_results[12,3] <- mean(bs_share[,6])
      share_results[12,4] <- mean(bs_share[,6]) - 1.96*sd(bs_share[,6])
      share_results[12,5] <- mean(bs_share[,6]) + 1.96*sd(bs_share[,6])
      
      #140
      share_results[22,3] <- mean(bs_share[,7])
      share_results[22,4] <- mean(bs_share[,7]) - 1.96*sd(bs_share[,7])
      share_results[22,5] <- mean(bs_share[,7]) + 1.96*sd(bs_share[,7])
      
      #150
      share_results[32,3] <- mean(bs_share[,8])
      share_results[32,4] <- mean(bs_share[,8]) - 1.96*sd(bs_share[,8])
      share_results[32,5] <- mean(bs_share[,8]) + 1.96*sd(bs_share[,8])
      
    #Fill in: Women
     
      #120
      share_results[3,3] <- mean(bs_share[,9])
      share_results[3,4] <- mean(bs_share[,9]) - 1.96*sd(bs_share[,9])
      share_results[3,5] <- mean(bs_share[,9]) + 1.96*sd(bs_share[,9])
      
      #130
      share_results[13,3] <- mean(bs_share[,10])
      share_results[13,4] <- mean(bs_share[,10]) - 1.96*sd(bs_share[,10])
      share_results[13,5] <- mean(bs_share[,10]) + 1.96*sd(bs_share[,10])
      
      #140
      share_results[23,3] <- mean(bs_share[,11])
      share_results[23,4] <- mean(bs_share[,11]) - 1.96*sd(bs_share[,11])
      share_results[23,5] <- mean(bs_share[,11]) + 1.96*sd(bs_share[,11])
      
      #150
      share_results[33,3] <- mean(bs_share[,12])
      share_results[33,4] <- mean(bs_share[,12]) - 1.96*sd(bs_share[,12])
      share_results[33,5] <- mean(bs_share[,12]) + 1.96*sd(bs_share[,12])
      
    #Fill in: African
     
      #120
      share_results[4,3] <- mean(bs_share[,13])
      share_results[4,4] <- mean(bs_share[,13]) - 1.96*sd(bs_share[,13])
      share_results[4,5] <- mean(bs_share[,13]) + 1.96*sd(bs_share[,13])
      
      #130
      share_results[14,3] <- mean(bs_share[,14])
      share_results[14,4] <- mean(bs_share[,14]) - 1.96*sd(bs_share[,14])
      share_results[14,5] <- mean(bs_share[,14]) + 1.96*sd(bs_share[,14])
      
      #140
      share_results[24,3] <- mean(bs_share[,15])
      share_results[24,4] <- mean(bs_share[,15]) - 1.96*sd(bs_share[,15])
      share_results[24,5] <- mean(bs_share[,15]) + 1.96*sd(bs_share[,15])
      
      #150
      share_results[34,3] <- mean(bs_share[,16])
      share_results[34,4] <- mean(bs_share[,16]) - 1.96*sd(bs_share[,16])
      share_results[34,5] <- mean(bs_share[,16]) + 1.96*sd(bs_share[,16])
      
    #Fill in: Coloured 
     
      #120
      share_results[5,3] <- mean(bs_share[,17])
      share_results[5,4] <- mean(bs_share[,17]) - 1.96*sd(bs_share[,17])
      share_results[5,5] <- mean(bs_share[,17]) + 1.96*sd(bs_share[,17])
      
      #130
      share_results[15,3] <- mean(bs_share[,18])
      share_results[15,4] <- mean(bs_share[,18]) - 1.96*sd(bs_share[,18])
      share_results[15,5] <- mean(bs_share[,18]) + 1.96*sd(bs_share[,18])
      
      #140
      share_results[25,3] <- mean(bs_share[,19])
      share_results[25,4] <- mean(bs_share[,19]) - 1.96*sd(bs_share[,19])
      share_results[25,5] <- mean(bs_share[,19]) + 1.96*sd(bs_share[,19])
      
      #150
      share_results[35,3] <- mean(bs_share[,20])
      share_results[35,4] <- mean(bs_share[,20]) - 1.96*sd(bs_share[,20])
      share_results[35,5] <- mean(bs_share[,20]) + 1.96*sd(bs_share[,20])
      
    #Fill in: Asian/Indian
     
      #120
      share_results[6,3] <- mean(bs_share[,21])
      share_results[6,4] <- mean(bs_share[,21]) - 1.96*sd(bs_share[,21])
      share_results[6,5] <- mean(bs_share[,21]) + 1.96*sd(bs_share[,21])
      
      #130
      share_results[16,3] <- mean(bs_share[,22])
      share_results[16,4] <- mean(bs_share[,22]) - 1.96*sd(bs_share[,22])
      share_results[16,5] <- mean(bs_share[,22]) + 1.96*sd(bs_share[,22])
      
      #140
      share_results[26,3] <- mean(bs_share[,23])
      share_results[26,4] <- mean(bs_share[,23]) - 1.96*sd(bs_share[,23])
      share_results[26,5] <- mean(bs_share[,23]) + 1.96*sd(bs_share[,23])
      
      #150
      share_results[36,3] <- mean(bs_share[,24])
      share_results[36,4] <- mean(bs_share[,24]) - 1.96*sd(bs_share[,24])
      share_results[36,5] <- mean(bs_share[,24]) + 1.96*sd(bs_share[,24])
      
    #Fill in: White
      
      #120
      share_results[7,3] <- mean(bs_share[,25])
      share_results[7,4] <- mean(bs_share[,25]) - 1.96*sd(bs_share[,25])
      share_results[7,5] <- mean(bs_share[,25]) + 1.96*sd(bs_share[,25])
      
      #130
      share_results[17,3] <- mean(bs_share[,26])
      share_results[17,4] <- mean(bs_share[,26]) - 1.96*sd(bs_share[,26])
      share_results[17,5] <- mean(bs_share[,26]) + 1.96*sd(bs_share[,26])
      
      #140
      share_results[27,3] <- mean(bs_share[,27])
      share_results[27,4] <- mean(bs_share[,27]) - 1.96*sd(bs_share[,27])
      share_results[27,5] <- mean(bs_share[,27]) + 1.96*sd(bs_share[,27])
      
      #150
      share_results[37,3] <- mean(bs_share[,28])
      share_results[37,4] <- mean(bs_share[,28]) - 1.96*sd(bs_share[,28])
      share_results[37,5] <- mean(bs_share[,28]) + 1.96*sd(bs_share[,28])
      
    #Fill in: Urban
      
      #120
      share_results[8,3] <- mean(bs_share[,29])
      share_results[8,4] <- mean(bs_share[,29]) - 1.96*sd(bs_share[,29])
      share_results[8,5] <- mean(bs_share[,29]) + 1.96*sd(bs_share[,29])
      
      #130
      share_results[18,3] <- mean(bs_share[,30])
      share_results[18,4] <- mean(bs_share[,30]) - 1.96*sd(bs_share[,30])
      share_results[18,5] <- mean(bs_share[,30]) + 1.96*sd(bs_share[,30])
      
      #140
      share_results[28,3] <- mean(bs_share[,31])
      share_results[28,4] <- mean(bs_share[,31]) - 1.96*sd(bs_share[,31])
      share_results[28,5] <- mean(bs_share[,31]) + 1.96*sd(bs_share[,31])
      
      #150
      share_results[38,3] <- mean(bs_share[,32])
      share_results[38,4] <- mean(bs_share[,32]) - 1.96*sd(bs_share[,32])
      share_results[38,5] <- mean(bs_share[,32]) + 1.96*sd(bs_share[,32])
      
    #Fill in: Traditional 
     
      #120
      share_results[9,3] <- mean(bs_share[,33])
      share_results[9,4] <- mean(bs_share[,33]) - 1.96*sd(bs_share[,33])
      share_results[9,5] <- mean(bs_share[,33]) + 1.96*sd(bs_share[,33])
      
      #130
      share_results[19,3] <- mean(bs_share[,34])
      share_results[19,4] <- mean(bs_share[,34]) - 1.96*sd(bs_share[,34])
      share_results[19,5] <- mean(bs_share[,34]) + 1.96*sd(bs_share[,34])
      
      #140
      share_results[29,3] <- mean(bs_share[,35])
      share_results[29,4] <- mean(bs_share[,35]) - 1.96*sd(bs_share[,35])
      share_results[29,5] <- mean(bs_share[,35]) + 1.96*sd(bs_share[,35])
      
      #150
      share_results[39,3] <- mean(bs_share[,36])
      share_results[39,4] <- mean(bs_share[,36]) - 1.96*sd(bs_share[,36])
      share_results[39,5] <- mean(bs_share[,36]) + 1.96*sd(bs_share[,36]) 
      
    #Fill in: Farms
     
      #120
      share_results[10,3] <- mean(bs_share[,37])
      share_results[10,4] <- mean(bs_share[,37]) - 1.96*sd(bs_share[,37])
      share_results[10,5] <- mean(bs_share[,37]) + 1.96*sd(bs_share[,37])
      
      #130
      share_results[20,3] <- mean(bs_share[,38])
      share_results[20,4] <- mean(bs_share[,38]) - 1.96*sd(bs_share[,38])
      share_results[20,5] <- mean(bs_share[,38]) + 1.96*sd(bs_share[,38])
      
      #140
      share_results[30,3] <- mean(bs_share[,39])
      share_results[30,4] <- mean(bs_share[,39]) - 1.96*sd(bs_share[,39])
      share_results[30,5] <- mean(bs_share[,39]) + 1.96*sd(bs_share[,39])
      
      #150
      share_results[40,3] <- mean(bs_share[,40])
      share_results[40,4] <- mean(bs_share[,40]) - 1.96*sd(bs_share[,40])
      share_results[40,5] <- mean(bs_share[,40]) + 1.96*sd(bs_share[,40]) 
   
    share_summary <- share_results %>% 
      mutate(est = format(round(est*100,1), nsmall = 1),
             lcl = format(round(lcl*100,1), nsmall = 1),
             ucl = format(round(ucl*100, 1), nsmall = 1),
             stest1 = str_c(est, lcl, sep = " ("),
             stest2 = str_c(stest1, ucl, sep = ", "),
             estimate = paste0(stest2, ")")) %>% 
      select(target, strata, estimate) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "share",
                  names_sep =  "_",
                  values_from = estimate)
  #BP Shift
  shift_results <- tibble(target = as_factor(rep(c("120","130","140","150"),each=10)),
                          strata = as_factor(rep(c("Overall","Men","Women","African","Coloured",
                                       "Asian/Indian","White", "Urban", "Traditional", "Farms"),4)),
                        est = as.double(NA),
                        lcl = as.double(NA),
                        ucl = as.double(NA))
    
    #Fill in: Overall
      
      #120
      shift_results[1,3] <- mean(bs_shift[,1])
      shift_results[1,4] <- mean(bs_shift[,1]) - 1.96*sd(bs_shift[,1])
      shift_results[1,5] <- mean(bs_shift[,1]) + 1.96*sd(bs_shift[,1])
      
      #130
      shift_results[11,3] <- mean(bs_shift[,2])
      shift_results[11,4] <- mean(bs_shift[,2]) - 1.96*sd(bs_shift[,2])
      shift_results[11,5] <- mean(bs_shift[,2]) + 1.96*sd(bs_shift[,2])
      
      #140
      shift_results[21,3] <- mean(bs_shift[,3])
      shift_results[21,4] <- mean(bs_shift[,3]) - 1.96*sd(bs_shift[,3])
      shift_results[21,5] <- mean(bs_shift[,3]) + 1.96*sd(bs_shift[,3])
      
      #150
      shift_results[31,3] <- mean(bs_shift[,4])
      shift_results[31,4] <- mean(bs_shift[,4]) - 1.96*sd(bs_shift[,4])
      shift_results[31,5] <- mean(bs_shift[,4]) + 1.96*sd(bs_shift[,4])
      
    #Fill in: Men
     
      #120
      shift_results[2,3] <- mean(bs_shift[,5])
      shift_results[2,4] <- mean(bs_shift[,5]) - 1.96*sd(bs_shift[,5])
      shift_results[2,5] <- mean(bs_shift[,5]) + 1.96*sd(bs_shift[,5])
      
      #130
      shift_results[12,3] <- mean(bs_shift[,6])
      shift_results[12,4] <- mean(bs_shift[,6]) - 1.96*sd(bs_shift[,6])
      shift_results[12,5] <- mean(bs_shift[,6]) + 1.96*sd(bs_shift[,6])
      
      #140
      shift_results[22,3] <- mean(bs_shift[,7])
      shift_results[22,4] <- mean(bs_shift[,7]) - 1.96*sd(bs_shift[,7])
      shift_results[22,5] <- mean(bs_shift[,7]) + 1.96*sd(bs_shift[,7])
      
      #150
      shift_results[32,3] <- mean(bs_shift[,8])
      shift_results[32,4] <- mean(bs_shift[,8]) - 1.96*sd(bs_shift[,8])
      shift_results[32,5] <- mean(bs_shift[,8]) + 1.96*sd(bs_shift[,8])
      
    #Fill in: Women
      
      #120
      shift_results[3,3] <- mean(bs_shift[,9])
      shift_results[3,4] <- mean(bs_shift[,9]) - 1.96*sd(bs_shift[,9])
      shift_results[3,5] <- mean(bs_shift[,9]) + 1.96*sd(bs_shift[,9])
      
      #130
      shift_results[13,3] <- mean(bs_shift[,10])
      shift_results[13,4] <- mean(bs_shift[,10]) - 1.96*sd(bs_shift[,10])
      shift_results[13,5] <- mean(bs_shift[,10]) + 1.96*sd(bs_shift[,10])
      
      #140
      shift_results[23,3] <- mean(bs_shift[,11])
      shift_results[23,4] <- mean(bs_shift[,11]) - 1.96*sd(bs_shift[,11])
      shift_results[23,5] <- mean(bs_shift[,11]) + 1.96*sd(bs_shift[,11])
      
      #150
      shift_results[33,3] <- mean(bs_shift[,12])
      shift_results[33,4] <- mean(bs_shift[,12]) - 1.96*sd(bs_shift[,12])
      shift_results[33,5] <- mean(bs_shift[,12]) + 1.96*sd(bs_shift[,12])
      
    #Fill in: African
      
      #120
      shift_results[4,3] <- mean(bs_shift[,13])
      shift_results[4,4] <- mean(bs_shift[,13]) - 1.96*sd(bs_shift[,13])
      shift_results[4,5] <- mean(bs_shift[,13]) + 1.96*sd(bs_shift[,13])
      
      #130
      shift_results[14,3] <- mean(bs_shift[,14])
      shift_results[14,4] <- mean(bs_shift[,14]) - 1.96*sd(bs_shift[,14])
      shift_results[14,5] <- mean(bs_shift[,14]) + 1.96*sd(bs_shift[,14])
      
      #140
      shift_results[24,3] <- mean(bs_shift[,15])
      shift_results[24,4] <- mean(bs_shift[,15]) - 1.96*sd(bs_shift[,15])
      shift_results[24,5] <- mean(bs_shift[,15]) + 1.96*sd(bs_shift[,15])
      
      #150
      shift_results[34,3] <- mean(bs_shift[,16])
      shift_results[34,4] <- mean(bs_shift[,16]) - 1.96*sd(bs_shift[,16])
      shift_results[34,5] <- mean(bs_shift[,16]) + 1.96*sd(bs_shift[,16])
      
    #Fill in: Coloured
      
      #120
      shift_results[5,3] <- mean(bs_shift[,17])
      shift_results[5,4] <- mean(bs_shift[,17]) - 1.96*sd(bs_shift[,17])
      shift_results[5,5] <- mean(bs_shift[,17]) + 1.96*sd(bs_shift[,17])
      
      #130
      shift_results[15,3] <- mean(bs_shift[,18])
      shift_results[15,4] <- mean(bs_shift[,18]) - 1.96*sd(bs_shift[,18])
      shift_results[15,5] <- mean(bs_shift[,18]) + 1.96*sd(bs_shift[,18])
      
      #140
      shift_results[25,3] <- mean(bs_shift[,19])
      shift_results[25,4] <- mean(bs_shift[,19]) - 1.96*sd(bs_shift[,19])
      shift_results[25,5] <- mean(bs_shift[,19]) + 1.96*sd(bs_shift[,19])
      
      #150
      shift_results[35,3] <- mean(bs_shift[,20])
      shift_results[35,4] <- mean(bs_shift[,20]) - 1.96*sd(bs_shift[,20])
      shift_results[35,5] <- mean(bs_shift[,20]) + 1.96*sd(bs_shift[,20])
      
    #Fill in: Asian/Indian
      
      #120
      shift_results[6,3] <- mean(bs_shift[,21])
      shift_results[6,4] <- mean(bs_shift[,21]) - 1.96*sd(bs_shift[,21])
      shift_results[6,5] <- mean(bs_shift[,21]) + 1.96*sd(bs_shift[,21])
      
      #130
      shift_results[16,3] <- mean(bs_shift[,22])
      shift_results[16,4] <- mean(bs_shift[,22]) - 1.96*sd(bs_shift[,22])
      shift_results[16,5] <- mean(bs_shift[,22]) + 1.96*sd(bs_shift[,22])
      
      #140
      shift_results[26,3] <- mean(bs_shift[,23])
      shift_results[26,4] <- mean(bs_shift[,23]) - 1.96*sd(bs_shift[,23])
      shift_results[26,5] <- mean(bs_shift[,23]) + 1.96*sd(bs_shift[,23])
      
      #150
      shift_results[36,3] <- mean(bs_shift[,24])
      shift_results[36,4] <- mean(bs_shift[,24]) - 1.96*sd(bs_shift[,24])
      shift_results[36,5] <- mean(bs_shift[,24]) + 1.96*sd(bs_shift[,24])
      
    #Fill in: White
      
      #120
      shift_results[7,3] <- mean(bs_shift[,25])
      shift_results[7,4] <- mean(bs_shift[,25]) - 1.96*sd(bs_shift[,25])
      shift_results[7,5] <- mean(bs_shift[,25]) + 1.96*sd(bs_shift[,25])
      
      #130
      shift_results[17,3] <- mean(bs_shift[,26])
      shift_results[17,4] <- mean(bs_shift[,26]) - 1.96*sd(bs_shift[,26])
      shift_results[17,5] <- mean(bs_shift[,26]) + 1.96*sd(bs_shift[,26])
      
      #140
      shift_results[27,3] <- mean(bs_shift[,27])
      shift_results[27,4] <- mean(bs_shift[,27]) - 1.96*sd(bs_shift[,27])
      shift_results[27,5] <- mean(bs_shift[,27]) + 1.96*sd(bs_shift[,27])
      
      #150
      shift_results[37,3] <- mean(bs_shift[,28])
      shift_results[37,4] <- mean(bs_shift[,28]) - 1.96*sd(bs_shift[,28])
      shift_results[37,5] <- mean(bs_shift[,28]) + 1.96*sd(bs_shift[,28])
      
    #Fill in: Urban
      
      #120
      shift_results[8,3] <- mean(bs_shift[,29])
      shift_results[8,4] <- mean(bs_shift[,29]) - 1.96*sd(bs_shift[,29])
      shift_results[8,5] <- mean(bs_shift[,29]) + 1.96*sd(bs_shift[,29])
      
      #130
      shift_results[18,3] <- mean(bs_shift[,30])
      shift_results[18,4] <- mean(bs_shift[,30]) - 1.96*sd(bs_shift[,30])
      shift_results[18,5] <- mean(bs_shift[,30]) + 1.96*sd(bs_shift[,30])
      
      #140
      shift_results[28,3] <- mean(bs_shift[,31])
      shift_results[28,4] <- mean(bs_shift[,31]) - 1.96*sd(bs_shift[,31])
      shift_results[28,5] <- mean(bs_shift[,31]) + 1.96*sd(bs_shift[,31])
      
      #150
      shift_results[38,3] <- mean(bs_shift[,32])
      shift_results[38,4] <- mean(bs_shift[,32]) - 1.96*sd(bs_shift[,32])
      shift_results[38,5] <- mean(bs_shift[,32]) + 1.96*sd(bs_shift[,32])
      
    #Fill in: Traditional
      
      #120
      shift_results[9,3] <- mean(bs_shift[,33])
      shift_results[9,4] <- mean(bs_shift[,33]) - 1.96*sd(bs_shift[,33])
      shift_results[9,5] <- mean(bs_shift[,33]) + 1.96*sd(bs_shift[,33])
      
      #130
      shift_results[19,3] <- mean(bs_shift[,34])
      shift_results[19,4] <- mean(bs_shift[,34]) - 1.96*sd(bs_shift[,34])
      shift_results[19,5] <- mean(bs_shift[,34]) + 1.96*sd(bs_shift[,34])
      
      #140
      shift_results[29,3] <- mean(bs_shift[,35])
      shift_results[29,4] <- mean(bs_shift[,35]) - 1.96*sd(bs_shift[,35])
      shift_results[29,5] <- mean(bs_shift[,35]) + 1.96*sd(bs_shift[,35])
      
      #150
      shift_results[39,3] <- mean(bs_shift[,36])
      shift_results[39,4] <- mean(bs_shift[,36]) - 1.96*sd(bs_shift[,36])
      shift_results[39,5] <- mean(bs_shift[,36]) + 1.96*sd(bs_shift[,36]) 
      
    #Fill in: Farms
     
      #120
      shift_results[10,3] <- mean(bs_shift[,37])
      shift_results[10,4] <- mean(bs_shift[,37]) - 1.96*sd(bs_shift[,37])
      shift_results[10,5] <- mean(bs_shift[,37]) + 1.96*sd(bs_shift[,37])
      
      #130
      shift_results[20,3] <- mean(bs_shift[,38])
      shift_results[20,4] <- mean(bs_shift[,38]) - 1.96*sd(bs_shift[,38])
      shift_results[20,5] <- mean(bs_shift[,38]) + 1.96*sd(bs_shift[,38])
      
      #140
      shift_results[30,3] <- mean(bs_shift[,39])
      shift_results[30,4] <- mean(bs_shift[,39]) - 1.96*sd(bs_shift[,39])
      shift_results[30,5] <- mean(bs_shift[,39]) + 1.96*sd(bs_shift[,39])
      
      #150
      shift_results[40,3] <- mean(bs_shift[,40])
      shift_results[40,4] <- mean(bs_shift[,40]) - 1.96*sd(bs_shift[,40])
      shift_results[40,5] <- mean(bs_shift[,40]) + 1.96*sd(bs_shift[,40])    

  shift_summary <- shift_results %>% 
      mutate(est = format(round(est,1), nsmall = 1),
             lcl = format(round(lcl,1), nsmall = 1),
             ucl = format(round(ucl, 1), nsmall = 1),
             stest1 = str_c(est, lcl, sep = " ("),
             stest2 = str_c(stest1, ucl, sep = ", "),
             estimate = paste0(stest2, ")")) %>% 
      select(target, strata, estimate) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "shift",
                  names_sep =  "_",
                  values_from = estimate)
  
  #NNT
  nnt_results <- mort_results %>% left_join(share_results, by = c('target', "strata")) %>% 
    mutate(nnt =  ceiling((est.y / est.x)/-1))
  
  nnt_summary <- nnt_results %>% 
      select(target, strata, nnt) %>% 
      pivot_wider(id_cols = strata,
                  names_from = "target",
                  names_prefix = "nnt",
                  names_sep =  "_",
                  values_from = nnt)
  
#.......................................................................................
# Population results table summary
#.......................................................................................   
  
  pop_summary_table <- mort_summary %>% 
    left_join(share_summary) %>%
    left_join(shift_summary) 
  write.csv(pop_summary_table, "pop_summary_table.csv")
  
#.......................................................................................
# Population Figures
#.......................................................................................       
  
  mort <- ggplot(data = filter(mort_results, strata == "Overall"), aes(x = fct_rev(as_factor(target)), y = est*1000, ymin = lcl*1000, ymax = ucl*1000, label = format(round(est*1000,digits = 1),nsmall =1))) +
    geom_errorbar(width = 0.3, size = 1) +
    geom_label() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "SBP target (mmHg)", y = "Deaths per 1,000", title = "Deaths averted") +
    coord_flip() + 
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14),
          axis.text = element_text(size = 12, color = "black"))

  share <- ggplot(data = filter(share_results, strata == "Overall"), aes(x = fct_rev(as_factor(target)), y = est*100, ymin = lcl*100, ymax = ucl*100, label = paste(round(est*100,digits = 0),"%",sep = ""))) +
    geom_errorbar(width = 0.5, size = 1) +
    geom_label() +
    labs(y = "Percentage", title = "Share in need of care") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14))

  nnt <- ggplot(data = filter(nnt_results, strata == "Overall"), aes(x = fct_rev(as_factor(target)), y = nnt, ymin = nnt, ymax = nnt, label = paste(round(nnt,digits = 0)))) +
    geom_errorbar(width = 0.0, size = 0) +
    geom_label() +
    labs(y = "Number", title = "Number Needed to Treat") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14));

  shift <- ggplot(data = filter(shift_results, strata == "Overall") , aes(x = fct_rev(as_factor(target)), y = est, ymin = lcl, ymax = ucl, label = format(round(est,digits = 1),nsmall =1))) +
    geom_errorbar(width = 0.5, size = 1) +
    geom_label() +
    labs(y = "mmHg", title = "Mean SBP reduction") +
    coord_flip() +
    theme_classic() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5),
          text = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_text(size = 12, color = "black"))

  #Combine
  ggarrange(mort,share,shift,nrow = 1,widths = c(1.2,1.1,1.1))
  ggsave("figure4.pdf",dpi = 320,units = "in",width = 11.5,height = 5)
