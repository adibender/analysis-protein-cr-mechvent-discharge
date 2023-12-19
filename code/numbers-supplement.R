source("1_packages.R")
source("2_function_helpers.R")
source("3_dir_create.R")

# Data import
ped_death <- readRDS(paste0(main_subfolder, "/data/ped-data-death-hosp-sub5.Rds"))
ped_extubation <- readRDS(paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5.Rds"))

patient <- readRDS(paste0(main_subfolder, "/data/patient.Rds"))
daily <- readRDS(paste0(main_subfolder, "/data/daily.Rds"))

### Table S1 of patient at risk/dying/discharged per interval ###
smry_death <- ped_death %>%
  group_by(interval) %>%
  summarize(at_risk = n(), dying = sum(ped_status))
smry_extubation <- ped_extubation %>%
  group_by(interval) %>%
  summarize(at_risk = n(), extubation = sum(ped_status)) %>%
  select(-at_risk)

smry_tab <- smry_death %>% left_join(smry_extubation, by = c("interval"))

readr::write_excel_csv(smry_tab, paste0(main_subfolder, "/results/supplement-table-death-extubation-per-interval.csv"))

### Table S2 of baseline characteristics ###
# Year
table(patient$Year)
round((table(patient$Year) / nrow(patient)) * 100, 1)

# Table 60-day mortality under MV (%)
ped_death_last <- ped_death %>%
  group_by(CombinedID) %>%
  slice(n()) %>%
  left_join(select(patient, EN2_4, CombinedID))
ped_extubation_last <- ped_extubation %>%
  group_by(CombinedID) %>%
  slice(n()) %>%
  left_join(select(patient, EN2_4, CombinedID))

ped_death_last %>%
  group_by(Year) %>%
  summarize(mdeath = mean(ped_status) * 100)
ped_death_last %>%
  group_by(Gender) %>%
  summarize(mdeath = mean(ped_status) * 100)
ped_death_last %>%
  group_by(inMV2_4) %>%
  summarize(mdeath = mean(ped_status) * 100)
ped_death_last %>%
  group_by(OralIntake2_4) %>%
  summarize(mdeath = mean(ped_status) * 100)
ped_death_last %>%
  group_by(Propofol2_4) %>%
  summarize(mdeath = mean(ped_status) * 100)
ped_death_last %>%
  group_by(PN2_4) %>%
  summarize(mdeath = mean(ped_status) * 100)

ped_death_last %>%
  group_by(EN2_4) %>%
  summarize(mdeath = mean(ped_status) * 100)

ped_death_last %>%
  group_by(AdmCatID) %>%
  summarize(mdeath = mean(ped_status) * 100)

ped_death_last %>%
  group_by(DiagID2) %>%
  summarize(mdeath = mean(ped_status) * 100)


res <- rep(0, length(unique(patient$Year)))
names(res) <- unique(patient$Year)
i <- 0
for (y in unique(patient$Year)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Year == y])
}
round(res * 100, 1)

# Gender
table(patient$Gender)
round((table(patient$Gender) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$Gender)))
names(res) <- unique(patient$Gender)
i <- 0
for (y in unique(patient$Gender)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Gender == y])
}
round(res * 100, 1)

# MV
MV <- patient$DaysMechVent * 24
MV <- cut(MV , c(0, 24, 48, max(MV)), right = TRUE, include.lowest = TRUE)
table(MV)
round((table(MV) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(MV)))
names(res) <- unique(MV)
i <- 0
for (y in unique(MV)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[MV == y])
}
round(res * 100, 1)

#OI
table(patient$OralIntake2_4)
round((table(patient$OralIntake2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$OralIntake2_4)))
names(res) <- sort(unique(patient$OralIntake2_4))
i <- 0
for (y in sort(unique(patient$OralIntake2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$OralIntake2_4 == y])
}
round(res * 100, 1)

#MV 2_4
table(patient$inMV2_4)
round((table(patient$inMV2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$inMV2_4)))
names(res) <- sort(unique(patient$inMV2_4))
i <- 0
for (y in sort(unique(patient$inMV2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$inMV2_4 == y])
}
round(res * 100, 1)

#PF
table(patient$Propofol2_4)
round((table(patient$Propofol2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$Propofol2_4)))
names(res) <- sort(unique(patient$Propofol2_4))
i <- 0
for (y in sort(unique(patient$Propofol2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$Propofol2_4 == y])
}
round(res * 100, 1)

#PN
table(patient$PN2_4)
round((table(patient$PN2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$PN2_4)))
names(res) <- sort(unique(patient$PN2_4))
i <- 0
for (y in sort(unique(patient$PN2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$PN2_4 == y])
}
round(res * 100, 1)

#EN # empty
table(patient$EN2_4)
round((table(patient$EN2_4) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$EN2_4)))
names(res) <- sort(unique(patient$EN2_4))
i <- 0
for (y in sort(unique(patient$EN2_4))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$EN2_4 == y])
}
round(res * 100, 1)

#Admin cat
table(patient$AdmCatID)
round((table(patient$AdmCatID) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$AdmCatID)))
names(res) <- names(table(patient$AdmCatID))
i <- 0
for (y in unique(patient$AdmCatID)) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$AdmCatID == y])
}
round(res * 100, 1)

# Diag
table(patient$DiagID2)
round((table(patient$DiagID2) / nrow(patient)) * 100, 1)

res <- rep(0, length(unique(patient$DiagID2)))
names(res) <- names(table(patient$DiagID2))
i <- 0
for (y in names(table(patient$DiagID2))) {
  i <- i + 1
  res[i] <- mean(patient$PatientDied[patient$DiagID2 == y])
}
round(res * 100, 1)

### Table S3: summaries ###
summary(patient[, c("Age", "ApacheIIScore", "BMI")])

### Table S4: Low versus standard protein intake ###
# HR per interval (death and discharge) for all comparisons
frames2A <- readRDS(paste0(main_subfolder, "/models/m_2_A_6frames.Rds"))
frames2B <- readRDS(paste0(main_subfolder, "/models/m_2_B_6frames.Rds"))

frames2A <- map(
  .x = frames2A,
  .f = ~{
    .x %>%
      mutate_at(
        c("fit", "lo", "hi"),
        ~round(.x, 2)) %>%
      mutate(CI = paste0("[", lo, ", ", hi, "]"))
  }
)

frames2B <- map(
  .x = frames2B,
  .f = ~{
    .x %>%
      mutate_at(
        c("fit", "lo", "hi"),
        ~round(.x, 2)) %>%
      mutate(CI = paste0("[", lo, ", ", hi, "]"))
  }
)

iwalk(frames2A, ~ readr::write_excel_csv(.x, paste0(main_subfolder, "/results/main/tab-death-", .y, ".csv")))
iwalk(frames2B, ~ readr::write_excel_csv(.x, paste0(main_subfolder, "/results/main/tab-extubation-", .y, ".csv")))



## Tables according to PN

ped_death_last <- ped_death_last %>%
  mutate(pn_bin = case_when(
    PN2_4 %in% c(0, 1) ~ "pn01",
    TRUE               ~ "pn23"
  ))

ped_death_last <- ped_death_last %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5 ~ "bmi<18.5",
    BMI > 30   ~ "bmi>30",
    TRUE       ~ "bmi18.5-30"
  ))
ped_extubation_last <- ped_extubation_last %>%
  mutate(pn_bin = case_when(
    PN2_4 %in% c(0, 1) ~ "pn01",
    TRUE               ~ "pn23"
  ))

ped_extubation_last <- ped_extubation_last %>%
  mutate(bmi_cat = case_when(
    BMI < 18.5 ~ "bmi<18.5",
    BMI > 30   ~ "bmi>30",
    TRUE       ~ "bmi18.5-30"
  ))

tab_pn <- ped_death_last %>%
  group_by(pn_bin) %>%
  summarize(
    m_age = mean(Age),
    IQR_age = paste0("[", quantile(Age, .25),", ", quantile(Age, .75), "]"),
    m_apache = mean(ApacheIIScore),
    IQR_apache = paste0("[", quantile(ApacheIIScore, .25), ", ", quantile(ApacheIIScore, .75), "]"),
    m_bmi = round(mean(BMI), 1),
    IQR_bmi = paste0("[", round(quantile(BMI, .25),1),", ", round(quantile(BMI, .75), 1), "]"),
    sex_n = paste0("male = ", sum(Gender == "Male"), ", female = ", sum(Gender == "Female")),
    sex_percentage = paste0("male = ", round(mean(Gender == "Male"), 2), ", female = ", round(mean(Gender == "Female"), 2)),
    admcat_n = paste0("medical = ", sum(AdmCatID == "Medical"), ", elective = ", sum(AdmCatID == "Surgical/Elective"), ", emergency = ", sum(AdmCatID == "Surgical/Emeregency")),
    admcat_percentage = paste0("medical = ", round(mean(AdmCatID == "Medical"), 2), ", elective = ", round(mean(AdmCatID == "Surgical/Elective"), 2), ", emergency = ", round(mean(AdmCatID == "Surgical/Emeregency"),2))
  )
tab_pn
readr::write_csv(tab_pn, "tab_pn.csv")


## cases subgroups 
ped_death_last %>% 
  group_by(AdmCatID) %>% 
  summarize(n = n(), cases = sum(ped_status))

ped_death_last %>% 
  group_by(pn_bin) %>% 
  summarize(n = n(), cases = sum(ped_status))

ped_death_last %>% 
  group_by(bmi_cat) %>% 
  summarize(n = n(), cases = sum(ped_status))


# time-to-death under mv 
ped_death_last %>% 
  filter(ped_status == 1) %>% 
  pull(tend) %>% 
  summary()