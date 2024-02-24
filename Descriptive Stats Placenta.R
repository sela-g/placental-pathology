#########################################
####  BC CANCOVID PLACENTAL PATH     ####
####   Sela Grays - Feb 23 2024      ####
#########################################
rm(list = ls())

## SET WORKING DIRECTORY
#setwd("/Users/selagrays/dev/cancov data cleaning/")
getwd()
`%notin%` <- Negate(`%in%`)

library(Gmisc)
library(broman)
library(dplyr)
library(DescTools)
library(tidyverse)
library(sjPlot)
library(here)
library(janitor)
library(zoo)
library(lubridate)
library(haven)
library(ghibli)
library(ggpubr)

#### LOAD BC QIQA DATA ####
data_bc <- read.csv("CanadianSurveillance_DATA_2024-02-23_1315.csv", header = TRUE)

data_bc[data_bc == ""] <- NA
data_bc_1 <- data_bc %>% filter(redcap_event_name == "demographics_arm_1")
data_bc_2 <- data_bc %>% filter(redcap_event_name == "sarscov2_arm_1")
data_bc_3 <- data_bc %>% filter(redcap_event_name == "antepartum_arm_1")
data_bc_4 <- data_bc %>% filter(redcap_event_name == "intrapartum_arm_1")
data_bc_5 <- data_bc %>% filter(redcap_event_name == "other_forms_arm_1")
data_bc_6 <- data_bc %>% filter(redcap_event_name == "postpartum_arm_1")

data_bc_1 <- remove_empty(data_bc_1, which = "cols")
data_bc_2 <- remove_empty(data_bc_2, which = "cols")
data_bc_3 <- remove_empty(data_bc_3, which = "cols")
data_bc_4 <- remove_empty(data_bc_4, which = "cols")
data_bc_5 <- remove_empty(data_bc_5, which = "cols")
data_bc_6 <- remove_empty(data_bc_6, which = "cols")

data_bc <- left_join(data_bc_1,data_bc_2, by = "a_record_id")
data_bc <- left_join(data_bc,data_bc_3, by = "a_record_id")
data_bc <- left_join(data_bc,data_bc_4, by = "a_record_id")
data_bc <- left_join(data_bc,data_bc_5, by = "a_record_id")
data_bc <- left_join(data_bc,data_bc_6, by = "a_record_id")

dim(data_bc)


twin_placenta_cases <- read.csv("CanadianSurveillance-NonSingletons_DATA_2024-02-23_1300.csv", header = TRUE)

twin_placenta_cases[twin_placenta_cases == ""] <- NA

twin_placenta_cases_3 <- twin_placenta_cases %>% filter(redcap_event_name == "antepartum_arm_1")
twin_placenta_cases_4 <- twin_placenta_cases %>% filter(redcap_event_name == "intrapartum_arm_1")

twin_placenta_cases_3 <- remove_empty(twin_placenta_cases_3, which = "cols")
twin_placenta_cases_4 <- remove_empty(twin_placenta_cases_4, which = "cols")

twin_placenta_cases <- left_join(twin_placenta_cases_4,twin_placenta_cases_3, by = "a_record_id")

dim(twin_placenta_cases)

twin_placenta_cases <- twin_placenta_cases$a_record_id

placentas <- read.csv("CanadianCOVID19InPre_DATA_2024-02-23_1254.csv", header = TRUE)

placentas <- placentas$cancovid_study_id

twin_placenta <- placentas[which(placentas %in% twin_placenta_cases)]

placentas_no_twins <- data_bc[which(data_bc$a_record_id %in% placentas
                              & data_bc$a_record_id %notin% twin_placenta_cases),]

### DEFINE CATAGORIES ###

## fix some diagnosis dates
placentas_no_twins$e_diagnosis <- as.Date(placentas_no_twins$e_diagnosis)
placentas_no_twins$r_dob <- as.Date(placentas_no_twins$r_dob)

placentas_no_twins$e_diagnosis[which(placentas_no_twins$e_diagnosis == "9999-09-09")] <- NA
placentas_no_twins$r_dob[which(placentas_no_twins$r_dob == "9999-09-09")] <- NA
placentas_no_twins$time_del <- as.numeric(as.Date(placentas_no_twins$r_dob) - placentas_no_twins$e_diagnosis)
summary(placentas_no_twins$time_del)
placentas_no_twins$i_lmp <- ifelse(placentas_no_twins$i_lmp == "9999-09-09", NA, as.character(placentas_no_twins$i_lmp))

placentas_no_twins$i_deliverydate_est <- as.Date(placentas_no_twins$i_deliverydate_est)

placentas_no_twins$i_deliverydate_est[which(is.na(placentas_no_twins$i_deliverydate_est) & !is.na(placentas_no_twins$i_lmp))] <- as.Date(placentas_no_twins$i_lmp[which(is.na(placentas_no_twins$i_deliverydate_est) & !is.na(placentas_no_twins$i_lmp))]) + 280

placentas_no_twins$time_del <- as.numeric(as.Date(placentas_no_twins$r_dob) - placentas_no_twins$e_diagnosis)
summary(placentas_no_twins$time_del)


## Ethnicity ##

placentas_no_twins <- placentas_no_twins %>% 
  rowwise() %>% 
  mutate(eth = case_when(
    all(c(b_ethnicity___1, b_ethnicity___2, b_ethnicity___3, b_ethnicity___4, b_ethnicity___5, b_ethnicity___6, b_ethnicity___7, b_ethnicity___8, b_ethnicity___998, b_ethnicity___999) == 0) ~ "Missing",
    b_ethnicity___999 == 1 ~ "Unknown",
    b_ethnicity___998 == 1 ~ "Other",
    b_ethnicity___1 == 1 ~ "White",
    b_ethnicity___2 == 1 ~ "African/Carribean/Black",
    b_ethnicity___3 == 1 ~ "Hispanic/Latino",
    b_ethnicity___4 == 1 ~ "East Asian",
    b_ethnicity___5 == 1 ~ "South Asian",
    b_ethnicity___6 == 1 ~ "South East Asian",
    b_ethnicity___7 == 1 ~ "Middle East",
    b_ethnicity___8 == 1 ~ "Indigenous"
  ))

describeFactors(placentas_no_twins$eth)
placentas_no_twins$eth <- factor(placentas_no_twins$eth)

levels(placentas_no_twins$eth)[which(levels(placentas_no_twins$eth) == "Missing" | levels(placentas_no_twins$eth) == "Unknown")] <- NA

placentas_no_twins$eth_cat <- placentas_no_twins$eth
levels(placentas_no_twins$eth_cat)[which(levels(placentas_no_twins$eth) == "East Asian" | levels(placentas_no_twins$eth) == "South East Asian")] <- "East or SE Asian"
placentas_no_twins$eth_cat <- relevel(placentas_no_twins$eth_cat, ref = "White")
describeFactors(placentas_no_twins$eth_cat)



## AGE
summary(placentas_no_twins$b_age) # check that Missing are coded as NA and not as 999 or something like that. Should be NA
placentas_no_twins$b_age <- ifelse(placentas_no_twins$b_age == 999, NA, placentas_no_twins$b_age)


### calculate gestational age at testing
# EDD is only available on the placentas_no_twins data set, therefore, some of the people tested will be missing GA at testing/diagnosis
# placentas_no_twins$i_deliverydate_est # EDD
placentas_no_twins$i_deliverydate_est <- as.Date(as.character(placentas_no_twins$i_deliverydate_est), format = "%Y-%m-%d")
summary(placentas_no_twins$i_deliverydate_est)

# fix
placentas_no_twins <- merge(placentas_no_twins, placentas_no_twins[, c("a_record_id", "e_diagnosis")], by = "a_record_id")

placentas_no_twins$e_diagnosis <- as.Date(placentas_no_twins$e_diagnosis.x)
# replace 9999-09-09 with NA
placentas_no_twins$i_deliverydate_est <- replace(placentas_no_twins$i_deliverydate_est, which(placentas_no_twins$i_deliverydate_est == "9999-09-09"), NA)

# split maternal age  
placentas_no_twins$age_cat <- factor(case_when(
  placentas_no_twins$b_age <25 ~ "<25 years",
  placentas_no_twins$b_age <30 ~ "25-29 years",
  placentas_no_twins$b_age <36 ~ "30-35 years",
  placentas_no_twins$b_age <40 ~ "36-39 years",
  placentas_no_twins$b_age >=40 ~ "≥40 years"
),
levels = c("<25 years", "25-29 years", "30-35 years", "36-39 years", "≥40 years"))

describeFactors(placentas_no_twins$age_cat)

summary(placentas_no_twins$b_age)

# GA at diagnosis into categories


# maternal BMI
placentas_no_twins$i_weight <- replace(placentas_no_twins$i_weight, which(placentas_no_twins$i_weight == 999 | placentas_no_twins$i_weight == 666), NA)
placentas_no_twins$i_height <- replace(placentas_no_twins$i_height, which(placentas_no_twins$i_height == 999 | placentas_no_twins$i_height == 666), NA)

placentas_no_twins$BMI <- placentas_no_twins$i_weight/(placentas_no_twins$i_height/100)^2
describeMedian(placentas_no_twins$BMI, iqr = FALSE) # check for out of range values - and remove or correct

# BMI >= 30 variable
placentas_no_twins$BMI_cat <- factor(case_when(
  placentas_no_twins$BMI < 18.5 ~ "<18.5",
  placentas_no_twins$BMI < 25 ~ "18.5-24",
  placentas_no_twins$BMI < 30 ~ "25-29",
  placentas_no_twins$BMI >= 30 ~ "≥30"
))
describeFactors(placentas_no_twins$BMI_cat)
placentas_no_twins$BMI_cat <- factor(placentas_no_twins$BMI_cat, levels = c("<18.5", "18.5-24", "25-29", "≥30"))


## Comorbs ##

# determine the number missing
placentas_no_twins <- placentas_no_twins %>%
  rowwise() %>%
  mutate(cvs = case_when(
    all(is.na(c(j_cns, j_cvs, j_resp, j_gi, j_gu, j_repro, j_endo, j_ms, j_hem, j_mh, j_aai))) & a_record_id %in% placentas_no_twins$a_record_id ~ "No",
    all(is.na(c(j_cns, j_cvs, j_resp, j_gi, j_gu, j_repro, j_endo, j_ms, j_hem, j_mh, j_aai))) & j_none___1 == 0 ~ "No entry",
    j_none___1 == 1 ~ "No",
    j_cvs == 1 ~ "Yes",
    l_htn == 1 ~ "Yes",
    j_cvs == 0 ~ "No",
    is.na(j_cvs) ~ "No"
  ))


# #replace the No entry with NA
# cp.ant <- cp.ant %>%
#   mutate(across(.cols = c(cvs, cns, resp, eentm, gi, gu, repro, endo, ms, hem, mh, aai, other_comor), ~ifelse(.x == "No entry", NA, .x)))
placentas_no_twins <- placentas_no_twins %>%
  mutate(across(.cols = c(cvs), ~ifelse(.x == "No entry", NA, .x)))

# check to make sure
describeFactors(placentas_no_twins$cvs)

# hypertension pre-existing
placentas_no_twins$htn <- case_when(
  placentas_no_twins$j_cvs_htn___1 == 1 ~ "Yes",
  is.na(placentas_no_twins$cvs) ~ NA_character_,
  !is.na(placentas_no_twins$cvs) ~ "No"
)
describeFactors(placentas_no_twins$htn)

# diabetes 1 or 2
placentas_no_twins$diabetes <- case_when(
  placentas_no_twins$j_endo_diabt1___1 == 1 | placentas_no_twins$j_endo_diabt2___1 == 1 ~ "Yes",
  is.na(placentas_no_twins$cvs) ~ NA_character_,
  !is.na(placentas_no_twins$cvs) ~ "No"
)
describeFactors(placentas_no_twins$diabetes)

## infant outcomes ##

# Gravida
placentas_no_twins$gravida <- placentas_no_twins$h_gravida
placentas_no_twins$gravida <- factor(case_when(
  placentas_no_twins$gravida == 0 ~ "0",
  placentas_no_twins$gravida == 1 ~ "1",
  placentas_no_twins$gravida >= 2 ~ "2+"
))

## losses/stillbirth
placentas_no_twins$p_outcome <- case_when(
  placentas_no_twins$p_outcome == 2 ~ "Stillbirth",
  placentas_no_twins$p_outcome == 3 ~ "Livebirth"
)

factor(placentas_no_twins$p_outcome)

describeFactors(placentas_no_twins$p_outcome)

# mode of delivery
# Vag/CS
placentas_no_twins$mode_del <- case_when(
  placentas_no_twins$p_mode == 1 ~ "Vaginal",
  placentas_no_twins$p_mode == 2 ~ "CS"
)

describeFactors(placentas_no_twins$mode_del)

# 5 minute apgar
placentas_no_twins$apgar5 <- factor(ifelse(placentas_no_twins$s_apgar_5 < 7, "<7", "≥7"))

describeFactors(placentas_no_twins$apgar5)

# birth weight
placentas_no_twins$bw_cat <- case_when(
  placentas_no_twins$s_bw_gm < 2500 ~ "<2500",
  placentas_no_twins$s_bw_gm >= 2500 & placentas_no_twins$s_bw_gm <=4000 ~ "2500-4000",
  placentas_no_twins$s_bw_gm > 4000 ~ ">4000"
)

describeFactors(placentas_no_twins$bw_cat)

summary(placentas_no_twins$s_bw_gm)

# NICU admission
# assuming if we have birth weight then we should know if baby admitted to NICU or not
placentas_no_twins$NICU <- case_when(
  placentas_no_twins$t_nicu == 1 ~ "Yes",
  placentas_no_twins$t_nicu == 0 ~ "No",
  is.na(placentas_no_twins$s_bw_gm) == FALSE ~ "No")

describeFactors(placentas_no_twins$NICU)

# duration of admission
placentas_no_twins$nicu_dur <- as.numeric(as.Date(placentas_no_twins$t_nicu_discharge) - as.Date(placentas_no_twins$t_nicu_admission))

diag_date <- case_when(
       !is.na(placentas_no_twins$d_naso1_collect) ~ as.Date(placentas_no_twins$d_naso1_collect),
       .default = as.Date(placentas_no_twins$d_other_collect)
       )
placentas_no_twins$ga_at_diag <- as.numeric((280 - (as.Date(placentas_no_twins$i_deliverydate_est) - diag_date))/7) # in weeks
# describeMedian(placentas_no_twins$ga_at_diag, iqr = FALSE) # check the range and fix any that are not possible - or remove. <0 and >43 weeks
placentas_no_twins$ga_at_diag <- replace(placentas_no_twins$ga_at_diag, which(placentas_no_twins$ga_at_diag >43 | placentas_no_twins$ga_at_diag < 0), NA)

placentas_no_twins$ga_dx_cat <- factor(case_when(
  placentas_no_twins$ga_at_diag <= 14 ~ "<=14 weeks",
  placentas_no_twins$ga_at_diag < 27 ~ "15-27 weeks",
  placentas_no_twins$ga_at_diag < 38 ~ "28-37 weeks",
  placentas_no_twins$ga_at_diag >=38 ~ "≥38 weeks"
), levels = c("<=14 weeks", "15-27 weeks", "28-37 weeks", "≥38 weeks"))

any(which(placentas_no_twins$ga_at_diag < 37 & placentas_no_twins$ga_at_diag > 27 ))

describeFactors(placentas_no_twins$ga_dx_cat)

summary(placentas_no_twins$ga_at_diag)


placentas_no_twins <- placentas_no_twins %>% mutate(covid_period = case_when(
  e_diagnosis < as.Date("2021-04-04") ~ "pre-Delta",
  e_diagnosis >= as.Date("2021-04-04") & e_diagnosis < as.Date("2021-12-19") ~ "Delta",
  e_diagnosis >= as.Date("2021-12-19") ~ "Omicron"
))

describeFactors(placentas_no_twins$covid_period)

## ga at delivery
placentas_no_twins$ga_at_del <- as.numeric((280 - (as.Date(placentas_no_twins$i_deliverydate_est) - as.Date(placentas_no_twins$r_dob)))/7)


placentas_no_twins$ga_del_cat <- case_when(
  placentas_no_twins$ga_at_del < 28 ~ "extremely preterm",
  placentas_no_twins$ga_at_del < 32 ~ "very preterm",
  placentas_no_twins$ga_at_del < 34 ~ "moderate preterm",
  placentas_no_twins$ga_at_del < 37 ~ "late preterm",
  placentas_no_twins$ga_at_del >= 37 ~ "term",
)


describeFactors(placentas_no_twins$ga_del_cat) # missings

# COVID outcomes

## hospitalizations
describeFactors(placentas_no_twins$e_hosp)

#describeFactors(placentas_no_twins$abx_pne)
describeFactors(placentas_no_twins$e_coag)
describeFactors(placentas_no_twins$g_icu)
describeFactors(placentas_no_twins$e_oxygen___1)
describeFactors(placentas_no_twins$e_inv___1)
describeFactors(placentas_no_twins$e_sepsis)


describeFactors(placentas_no_twins$l_htn)

describeFactors(placentas_no_twins$l_diab)

describeFactors(placentas_no_twins$n_covid)

placentas_no_twins$vacc_count <- case_when(
  (!is.na(placentas_no_twins$n_covid_date1) & !is.na(placentas_no_twins$n_covid_date2)) ~ "two doses",
  (!is.na(placentas_no_twins$n_covid_date1) & is.na(placentas_no_twins$n_covid_date2)) ~ "one doses"
)

describeFactors(placentas_no_twins$vacc_count)
