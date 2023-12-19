#***************************** 5_model_subgroup_adm_cat ***********************#
#*
#******************************************************************************#
#*
#***************************** Parameters *************************************#
# load in the neccesary packages
source("1_packages.R")
# set the wished output theme
theme_set(theme_bw())
# get the helper functions
source("2_function_helpers.R")
source("3_dir_create.R")
#******************************************************************************#
#*
#*************************** loading in data **********************************#

patient <- readRDS("data/patient.Rds")
daily <- readRDS("data/daily.Rds")
mergedData <- readRDS("data/mergedAndCleanedData.Rds")

summary(patient)
summary(daily)

patient <- patient %>%
  filter(DaysMechVent != 0)

daily <- daily[daily$CombinedID %in% patient$CombinedID, ]

ll_fun <- function(t, tz) { t >= tz + 4 & t <= tz * 3 + 10 }
ll <- get_laglead(
  0:60,
  tz = 1:11,
  ll_fun = ll_fun)
ll_dynamic <- gg_laglead(ll) + scale_x_discrete(guide = guide_axis(angle = 90))

ggsave(ll_dynamic, filename = paste0(main_subfolder,"/results/figures/ll.png"))
saveRDS(ll_fun, paste0(main_subfolder, "/data/ll.Rds"))

#******************************************************************************#
#*
#**************************** Creation of new Variable ************************#
# Creation of the variables mv_status and event_time related to our study characteristic: extubation

# mv_status:
# 2 died before extubation
# 1 extubation
# 0 censored

patient <- patient %>%
  mutate(mv_status = case_when(DaysMechVent >= 60 | DaysMechVent > Disc0To60 ~ 0,
                               DaysMechVent <= Disc0To60 | DaysMechVent < Surv0To60 ~ 1,
                               DaysMechVent >= Surv0To60 ~ 2))

# event_time:
patient <- patient %>%
  mutate(event_time = case_when(mv_status == 1 ~ DaysMechVent,
                                mv_status == 2 ~ Surv0To60,
                                mv_status == 0 ~ 61))


# delta_1:
# 1 Death before extubation
# 0 else

# delta_2:
# 1 Extubation
# 0 else

patient <- patient %>%
  mutate(delta_1 = 1L * (mv_status == 2),
         delta_2 = 1L * (mv_status == 1))

#******************************************************************************#
#*
#***************************** ADM cat  model *************************************#
# AdmCatID removed from predictors
formula = as.formula('ped_status ~ s(int_mid, bs = "ps") +
  Year +
  ApacheIIScore + ApacheIIScore:int_mid +
  s(Age, by = int_mid, bs = "ps") +
  s(BMI, by = int_mid, bs = "ps") + DiagID2+ Gender +
  inMV2_4 + Propofol2_4 + OralIntake2_4 + PN2_4 +
  s(CombinedicuID, bs = "re", by = icuByDummy) +
 te(t, tz, by = I(LL * proteinCat3),
  bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal") +
  te(t, tz, by = I(LL * proteinCat2),
   bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal")
  ')

############## Model 2 ###############ein
patient_medical <- filter(patient, AdmCatID == "Medical")
### Model A: Death before extubation
ped_2_A_medical <- as_ped(
  data    = patient_medical,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_medical <- ped_2_A_medical %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_medical$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_medical$int_mid <- 0.5 * (ped_2_A_medical$tstart + ped_2_A_medical$tend)
ped_2_A_medical <- ped_2_A_medical[ped_2_A_medical$tend >= 5, ]
saveRDS(ped_2_A_medical, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_medical.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_medical <- bam(
  formula = formula,
  data    = ped_2_A_medical,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_medical,  paste0(main_subfolder,"/models/m_2_A_medical.Rds"))

### Model B: Extubation
ped_2_B_medical <- as_ped(
  data    = patient_medical,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_medical <- ped_2_B_medical %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_medical$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_medical$int_mid <- 0.5 * (ped_2_B_medical$tstart + ped_2_B_medical$tend)
ped_2_B_medical <- ped_2_B_medical[ped_2_B_medical$tend >= 5, ]
saveRDS(ped_2_B_medical, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_medical.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_medical <- bam(
  formula = formula,
  data    = ped_2_B_medical,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_medical, paste0(main_subfolder,"/models/m_2_B_medical.Rds"))


############## subgroup surgical emergency #####################################
patient_emergency <- filter(patient, AdmCatID == "Surgical/Emeregency")
### Model A: Death before extubation
ped_2_A_emergency <- as_ped(
  data    = patient_emergency,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_emergency <- ped_2_A_emergency %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_emergency$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_emergency$int_mid <- 0.5 * (ped_2_A_emergency$tstart + ped_2_A_emergency$tend)
ped_2_A_emergency <- ped_2_A_emergency[ped_2_A_emergency$tend >= 5, ]
saveRDS(ped_2_A_emergency, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_emergency.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_emergency <- bam(
  formula = formula,
  data    = ped_2_A_emergency,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_emergency,  paste0(main_subfolder,"/models/m_2_A_emergency.Rds"))

### Model B: Extubation
ped_2_B_emergency <- as_ped(
  data    = patient_emergency,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_emergency <- ped_2_B_emergency %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_emergency$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_emergency$int_mid <- 0.5 * (ped_2_B_emergency$tstart + ped_2_B_emergency$tend)
ped_2_B_emergency <- ped_2_B_emergency[ped_2_B_emergency$tend >= 5, ]
saveRDS(ped_2_B_emergency, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_emergency.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_emergency <- bam(
  formula = formula,
  data    = ped_2_B_emergency,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_emergency, paste0(main_subfolder,"/models/m_2_B_emergency.Rds"))





############## subgroup surgical elective #####################################
patient_elective <- filter(patient, AdmCatID == "Surgical/Elective")
### Model A: Death before extubation
ped_2_A_elective <- as_ped(
  data    = patient_elective,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_elective <- ped_2_A_elective %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_elective$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_elective$int_mid <- 0.5 * (ped_2_A_elective$tstart + ped_2_A_elective$tend)
ped_2_A_elective <- ped_2_A_elective[ped_2_A_elective$tend >= 5, ]
saveRDS(ped_2_A_elective, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_elective.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_elective <- bam(
  formula = formula,
  data    = ped_2_A_elective,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_elective,  paste0(main_subfolder,"/models/m_2_A_elective.Rds"))

### Model B: Extubation
ped_2_B_elective <- as_ped(
  data    = patient_elective,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_elective <- ped_2_B_elective %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_elective$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_elective$int_mid <- 0.5 * (ped_2_B_elective$tstart + ped_2_B_elective$tend)
ped_2_B_elective <- ped_2_B_elective[ped_2_B_elective$tend >= 5, ]
saveRDS(ped_2_B_elective, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_elective.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_elective <- bam(
  formula = formula,
  data    = ped_2_B_elective,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_elective, paste0(main_subfolder,"/models/m_2_B_elective.Rds"))
