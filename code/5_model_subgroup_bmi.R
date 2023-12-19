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

patient <- patient %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5 ~ "bmi<18.5",
    BMI > 30   ~ "bmi>30",
    TRUE       ~ "bmi18.5-30"
  ))

#******************************************************************************#
#*
#***************************** ADM cat  model *************************************#
# BMI removed from predictors
formula = as.formula('ped_status ~ s(int_mid, bs = "ps") +
  Year +
  ApacheIIScore + ApacheIIScore:int_mid +
  s(Age, by = int_mid, bs = "ps") + DiagID2 + Gender +
  inMV2_4 + Propofol2_4 + OralIntake2_4 + PN2_4 +
  s(CombinedicuID, bs = "re", by = icuByDummy) +
 te(t, tz, by = I(LL * proteinCat3),
  bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal") +
  te(t, tz, by = I(LL * proteinCat2),
   bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal")
  ')

############## Model 2 ###############ein
patient_bmi1 <- filter(patient, bmi_cat == "bmi<18.5")
### Model A: Death before extubation
ped_2_A_bmi1 <- as_ped(
  data    = patient_bmi1,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_bmi1 <- ped_2_A_bmi1 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi1$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_bmi1$int_mid <- 0.5 * (ped_2_A_bmi1$tstart + ped_2_A_bmi1$tend)
ped_2_A_bmi1 <- ped_2_A_bmi1[ped_2_A_bmi1$tend >= 5, ]
saveRDS(ped_2_A_bmi1, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_bmi1.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_bmi1 <- bam(
  formula = formula,
  data    = ped_2_A_bmi1,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_bmi1,  paste0(main_subfolder,"/models/m_2_A_bmi1.Rds"))

### Model B: Extubation
ped_2_B_bmi1 <- as_ped(
  data    = patient_bmi1,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_bmi1 <- ped_2_B_bmi1 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi1$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_bmi1$int_mid <- 0.5 * (ped_2_B_bmi1$tstart + ped_2_B_bmi1$tend)
ped_2_B_bmi1 <- ped_2_B_bmi1[ped_2_B_bmi1$tend >= 5, ]
saveRDS(ped_2_B_bmi1, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_bmi1.Rds"))

model_2_B_bmi1 <- bam(
  formula = formula,
  data    = ped_2_B_bmi1,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_bmi1, paste0(main_subfolder,"/models/m_2_B_bmi1.Rds"))


############## subgroup bmi > 30 #####################################
patient_bmi3 <- filter(patient, bmi_cat == "bmi>30")
### Model A: Death before extubation
ped_2_A_bmi3 <- as_ped(
  data    = patient_bmi3,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_bmi3 <- ped_2_A_bmi3 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi3$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_bmi3$int_mid <- 0.5 * (ped_2_A_bmi3$tstart + ped_2_A_bmi3$tend)
ped_2_A_bmi3 <- ped_2_A_bmi3[ped_2_A_bmi3$tend >= 5, ]
saveRDS(ped_2_A_bmi3, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_bmi3.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_bmi3 <- bam(
  formula = formula,
  data    = ped_2_A_bmi3,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_bmi3,  paste0(main_subfolder,"/models/m_2_A_bmi3.Rds"))

### Model B: Extubation
ped_2_B_bmi3 <- as_ped(
  data    = patient_bmi3,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_bmi3 <- ped_2_B_bmi3 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi3$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_bmi3$int_mid <- 0.5 * (ped_2_B_bmi3$tstart + ped_2_B_bmi3$tend)
ped_2_B_bmi3 <- ped_2_B_bmi3[ped_2_B_bmi3$tend >= 5, ]
saveRDS(ped_2_B_bmi3, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_bmi3.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_bmi3 <- bam(
  formula = formula,
  data    = ped_2_B_bmi3,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_bmi3, paste0(main_subfolder,"/models/m_2_B_bmi3.Rds"))





############## subgroup bmi 18.5 - 30 #####################################
patient_bmi2 <- filter(patient, bmi_cat == "bmi18.5-30")
### Model A: Death before extubation
ped_2_A_bmi2 <- as_ped(
  data    = patient_bmi2,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_bmi2 <- ped_2_A_bmi2 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi2$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_bmi2$int_mid <- 0.5 * (ped_2_A_bmi2$tstart + ped_2_A_bmi2$tend)
ped_2_A_bmi2 <- ped_2_A_bmi2[ped_2_A_bmi2$tend >= 5, ]
saveRDS(ped_2_A_bmi2, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_bmi2.Rds"))

model_2_A_bmi2 <- bam(
  formula = formula,
  data    = ped_2_A_bmi2,
  family  = "poisson",
  offset  = offset
  )

saveRDS(model_2_A_bmi2,  paste0(main_subfolder,"/models/m_2_A_bmi2.Rds"))

### Model B: Extubation
ped_2_B_bmi2 <- as_ped(
  data    = patient_bmi2,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_bmi2 <- ped_2_B_bmi2 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_bmi2$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_bmi2$int_mid <- 0.5 * (ped_2_B_bmi2$tstart + ped_2_B_bmi2$tend)
ped_2_B_bmi2 <- ped_2_B_bmi2[ped_2_B_bmi2$tend >= 5, ]
saveRDS(ped_2_B_bmi2, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_bmi2.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_bmi2 <- bam(
  formula = formula,
  data    = ped_2_B_bmi2,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_bmi2, paste0(main_subfolder,"/models/m_2_B_bmi2.Rds"))
