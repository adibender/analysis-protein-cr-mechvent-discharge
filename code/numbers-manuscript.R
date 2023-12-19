source("1_packages.R")
source("2_function_helpers.R")
source("3_dir_create.R")

patient_orig      <- readRDS("data/patient.Rds")
patient           <- readRDS(paste0(main_subfolder, "/data/patient.Rds"))
daily             <- readRDS(paste0(main_subfolder, "/data/daily.Rds"))

ped_death         <- readRDS(paste0(main_subfolder, "/data/ped-data-death-hosp-sub5.Rds"))
ped_extubation    <- readRDS(paste0(main_subfolder, "/data/ped-data-extubation-hosp-sub5.Rds"))

####### Results ######
### number of subjects in data base
nrow(patient_orig)


### Study participants ###
uid <- unique(ped_death$CombinedID)

# number of subjects in final model
uid_n <- length(unique(uid))
uid_n

# number extubation & deaths before day 60
ped_extubation %>% pull(ped_status) %>% sum()
ped_extubation %>% pull(ped_status) %>% sum() / uid_n

ped_death %>% pull(ped_status) %>% sum()
ped_death %>% pull(ped_status) %>% sum() / uid_n

uid_n - (ped_extubation %>% pull(ped_status) %>% sum()) - (ped_death %>% pull(ped_status) %>% sum())
(uid_n - (ped_extubation %>% pull(ped_status) %>% sum()) - (ped_death %>% pull(ped_status) %>% sum())) / uid_n

# number of subject died within 48 after weaning
patient %>%
  filter(PatientDied == 1,
         delta_2 == 1,
         Surv0To60 - DaysMechVent <= 2) %>%
  nrow()

# number of subject died within 48 - 72 after weaning
patient %>%
  filter(PatientDied == 1,
         delta_2 == 1,
         Surv0To60 - DaysMechVent > 2,
         Surv0To60 - DaysMechVent <= 3) %>%
  nrow()

# days on ventilation, in ICU, in hospital
patient %>%
  filter(delta_2 == 1) %>%
  summarise(summary(DaysMechVent), summary(DaysInICU), summary(event))

### Protein intake ###
# number of days on diet
dietdays <- table(daily$incomplete_day)[1]
dietdays
dietdays / uid_n

# days with parenteral amino acids intake
table(daily$PN)
table(daily$PN) / dietdays

# Percentage EN of total protein intake
daily %>%
  summarise(sum(EN)/sum(EN, PN))

# number of days with protein intake
proteindays <- daily %>%
  filter(incomplete_day == FALSE,
         proteinAdjustedPercentage != 0) %>%
  summarise(n())
proteindays
# number patients with at least one day with low amount of protein
daily %>%
  filter(incomplete_day == FALSE,
         proteinCat2 == 0,
         proteinCat3 == 0) %>%
  summarise(length(unique(CombinedID)), length(unique(CombinedID))/uid_n)

# days of protein & average intake & avg. calorie intake: level I, II, III
daily %>%
  filter(incomplete_day == FALSE,
         proteinGproKG < 0.8 & proteinGproKG != 0) %>%
  summarise("days protein I:" = n()/proteindays,
            "daily protein intake" = summary(proteinGproKG),
            "daily total calorie intake" = summary(calproKg))

daily %>%
  filter(incomplete_day == FALSE,
         proteinGproKG >= 0.8 & proteinGproKG <= 1.2) %>%
  summarise("days protein II:" = n()/proteindays,
            "daily protein intake" = summary(proteinGproKG),
            "daily total calorie intake" = summary(calproKg))

daily %>%
  filter(incomplete_day == FALSE,
         proteinGproKG > 1.2) %>%
  summarise("days protein III:" = n()/proteindays,
            "daily protein intake" = summary(proteinGproKG),
            "daily total calorie intake" = summary(calproKg))

# comparison of hazard ratios
f6$`late standard vs. early standard` %>%
  filter(fit == min(fit))

f6$`exclusively low vs. late standard` %>%
  filter(fit == min(fit))

f6$`exclusively low vs. early standard` %>%
  filter(fit == min(fit))



f6$`early standard vs. early high` %>%
  filter(fit == max(fit))
