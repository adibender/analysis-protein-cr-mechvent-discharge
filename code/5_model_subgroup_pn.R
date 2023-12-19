#***************************** 5_model_subgroup_pn ****************************#
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
  mutate(pn_bin = case_when(
    PN2_4 %in% c(0, 1) ~ "pn01",
    TRUE               ~ "pn23"
  ))
table(patient$pn_bin)
#******************************************************************************#
#*
#****************************** PN  model *************************************#
# PN2_4 removed from predictors
formula = as.formula('ped_status ~ s(int_mid, bs = "ps") +
  Year +
  ApacheIIScore + ApacheIIScore:int_mid +
  s(Age, by = int_mid, bs = "ps") +
  s(BMI, by = int_mid, bs = "ps") + DiagID2+ Gender +
  inMV2_4 + Propofol2_4 + OralIntake2_4 +
  s(CombinedicuID, bs = "re", by = icuByDummy) +
 te(t, tz, by = I(LL * proteinCat3),
  bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal") +
  te(t, tz, by = I(LL * proteinCat2),
   bs = "ps", m = list(c(2, 1), c(2, 1)), id = "cal")
  ')

############## Model 2 ###############
patient_pn01 <- filter(patient, pn_bin == "pn01")
### Model A: Death before extubation
ped_2_A_pn01 <- as_ped(
  data    = patient_pn01,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_pn01 <- ped_2_A_pn01 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_pn01$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_pn01$int_mid <- 0.5 * (ped_2_A_pn01$tstart + ped_2_A_pn01$tend)
ped_2_A_pn01 <- ped_2_A_pn01[ped_2_A_pn01$tend >= 5, ]
saveRDS(ped_2_A_pn01, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_pn01.Rds"))


model_2_A_pn01 <- bam(
  formula = formula,
  data    = ped_2_A_pn01,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_pn01,  paste0(main_subfolder,"/models/m_2_A_pn01.Rds"))

### Model B: Extubation
ped_2_B_pn01 <- as_ped(
  data    = patient_pn01,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_pn01 <- ped_2_B_pn01 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_pn01$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_pn01$int_mid <- 0.5 * (ped_2_B_pn01$tstart + ped_2_B_pn01$tend)
ped_2_B_pn01 <- ped_2_B_pn01[ped_2_B_pn01$tend >= 5, ]
saveRDS(ped_2_B_pn01, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_pn01.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_pn01 <- bam(
  formula = formula,
  data    = ped_2_B_pn01,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_pn01, paste0(main_subfolder,"/models/m_2_B_pn01.Rds"))


############## subgroup pn23 #####################################
patient_pn23 <- filter(patient, pn_bin == "pn23")
### Model A: Death before extubation
ped_2_A_pn23 <- as_ped(
  data    = patient_pn23,
  formula = Surv(event_time, delta_1) ~ Year + DiagID2 + AdmCatID + Gender +
    ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 + PN2_4 + Age +
    CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_A_pn23 <- ped_2_A_pn23 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_pn23$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_A_pn23$int_mid <- 0.5 * (ped_2_A_pn23$tstart + ped_2_A_pn23$tend)
ped_2_A_pn23 <- ped_2_A_pn23[ped_2_A_pn23$tend >= 5, ]
saveRDS(ped_2_A_pn23, paste0(main_subfolder, "/data/ped-data-death-hosp-sub5_pn23.Rds"))

# Ucan can use this code as a basis for the parallelization if needed
# parallezation for the bam function
# numCores <- detectCores()-6
# cluster_1 <- makeCluster(numCores)
# for this you also need to set the cluster argument of the BAM function to
# cluster = cluster_1

model_2_A_pn23 <- bam(
  formula = formula,
  data    = ped_2_A_pn23,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_A_pn23,  paste0(main_subfolder,"/models/m_2_A_pn23.Rds"))

### Model B: Extubation
ped_2_B_pn23 <- as_ped(
  data    = patient_pn23,
  formula = Surv(event_time, delta_2) ~ Year + DiagID2 + AdmCatID +
    Gender + ApacheIIScore + BMI + Propofol2_4 + inMV2_4 + OralIntake2_4 +
    PN2_4 + Age + CombinedicuID + icuByDummy,
  cut     = 0:60, id = "CombinedID")
ped_2_B_pn23 <- ped_2_B_pn23 %>%
  add_cumulative_eff_vec(
    filter(daily, CombinedID %in% patient_pn23$CombinedID),
    "proteinCat2", "proteinCat3", LL = ll)
ped_2_B_pn23$int_mid <- 0.5 * (ped_2_B_pn23$tstart + ped_2_B_pn23$tend)
ped_2_B_pn23 <- ped_2_B_pn23[ped_2_B_pn23$tend >= 5, ]
saveRDS(ped_2_B_pn23, paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5_pn23.Rds"))

# parallelization
# numCores <- detectCores()-6
#
# cluster_1 <- makeCluster(numCores)

model_2_B_pn23 <- bam(
  formula = formula,
  data    = ped_2_B_pn23,
  family  = "poisson",
  offset  = offset
  )

# stopCluster()
saveRDS(model_2_B_pn23, paste0(main_subfolder,"/models/m_2_B_pn23.Rds"))
