source("1_packages.R")
source("2_function_helpers.R")
source("3_dir_create.R")

# set maximal follow-up
max.follow <- 60
# for preprocessing, set survival time of patients with survival > max.follow to max.time
max.time <- 60.1
# interval break points (for piece-wise exponential model format)
brks <- c(0:floor(max.time))

## for pattern_label function
# define number of days time-dependent covariate was provided
maxdays.nutri <- 11
# define protocols for comparisons

protocol1 <- rep("< 0.8 g/kg", maxdays.nutri)
protocol2 <- c(rep("< 0.8 g/kg", 4), rep("0.8 - 1.2 g/kg", 7))
protocol3 <- rep("0.8 - 1.2 g/kg", maxdays.nutri)
protocol4 <- c(rep("0.8 - 1.2 g/kg", 4), rep("> 1.2 g/kg", 7))
protocol5 <- rep("> 1.2 g/kg", maxdays.nutri)

comparisons_df <-
  data.frame(
    comparison = LETTERS[1:6],
    p1 = c("protocol1", "protocol2", "protocol1", "protocol3", "protocol4", "protocol3"),
    p2 = c("protocol2", "protocol3", "protocol3", "protocol4", "protocol5", "protocol5"),
    p1.name = c("exclusively low", "late standard", "exclusively low", "early standard", "late high", "early standard"),
    p2.name = c("late standard", "early standard", "early standard", "late high", "early high", "early high"),
    p1.count = c(1, 2, 1, 3, 4, 3),
    p2.count = c(2, 3, 3, 4, 5, 5)
  )

list_p_scheme <- purrr::map(
  .x = seq_len(nrow(comparisons_df)),
  .f = ~ggcomparison_scheme(
    protocol1 = get(comparisons_df[.x,"p1"]),
    protocol2 = get(comparisons_df[.x, "p2"]),
    protocol1.name = comparisons_df[.x, "p1.name"],
    protocol2.name = comparisons_df[.x, "p2.name"],
    legend = TRUE
  )
)

table(patient$AdmCatID)
# exclude Surgical/Elective -> too few cases

# model objects
# subgroup: medical
m_2_A_medical <- readRDS(paste0(main_subfolder, "/models/m_2_A_medical.Rds")) # cause specific: death
m_2_B_medical <- readRDS(paste0(main_subfolder, "/models/m_2_B_medical.Rds")) # cause specific: extubation
# subgroup: emergency
m_2_A_emergency <- readRDS(paste0(main_subfolder, "/models/m_2_A_emergency.Rds")) # cause specific: death
m_2_B_emergency <- readRDS(paste0(main_subfolder, "/models/m_2_B_emergency.Rds")) # cause specific: extubation

daily   <- readRDS(paste0(main_subfolder, "/data/daily.Rds"))
patient <- readRDS(paste0(main_subfolder, "/data/patient.Rds"))
ll      <- readRDS(paste0(main_subfolder, "/data/ll.Rds"))

# subgroup_data
patient_medical <- filter(patient, AdmCatID == "Medical")
patient_emergency <- filter(patient, AdmCatID == "Surgical/Emeregency")
daily_medical <- daily %>% filter(CombinedID %in% patient_medical$CombinedID)
daily_emergency <- daily %>% filter(CombinedID %in% patient_emergency$CombinedID)

## subgroup_medical
f6 <- make_six_frames(m_2_A_medical, patient_medical, daily_medical, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_extubation <- make_six_frames(m_2_B_medical, patient_medical, daily_medical, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)

list_p_effect <- purrr::map(seq_along(f6), ~single_plot(f6[[.x]], NULL))
list_p_discharge <- purrr::map(
  seq_along(f6_extubation),
  ~single_plot(f6_extubation[[.x]], NULL))

names(f6) <- names(f6_extubation) <-
  paste0(comparisons_df$p1.name, " vs. ", comparisons_df$p2.name)


p_1 <- list_p_scheme[[1]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[1]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_discharge[[1]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[2]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[2]] + list_p_discharge[[2]] +
  list_p_scheme[[3]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[3]] + list_p_discharge[[3]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_1
ggsave(p_1,
  file = paste0(main_subfolder,
    "/results/figures/supplement/sub_medical_protein_effect_comparisons_1.png"),
  height = 9,
  width = 9)

p_2 <- list_p_scheme[[4]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[4]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_discharge[[4]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[5]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[5]] + list_p_discharge[[5]] +
  list_p_scheme[[6]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[6]] + list_p_discharge[[6]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_2
ggsave(p_2,
  file = paste0(main_subfolder,
    "/results/figures/supplement/sub_medical_protein_effect_comparisons_2.png"),
  height = 9,
  width = 9)

## subgroup emergency
f6 <- make_six_frames(m_2_A_emergency, patient_emergency, daily_emergency, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)
f6_extubation <- make_six_frames(m_2_B_emergency, patient_emergency, daily_emergency, ll,
  var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
  surv_max = max.follow)

list_p_effect <- purrr::map(seq_along(f6), ~single_plot(f6[[.x]], NULL))
list_p_discharge <- purrr::map(
  seq_along(f6_extubation),
  ~single_plot(f6_extubation[[.x]], NULL))

names(f6) <- names(f6_extubation) <-
  paste0(comparisons_df$p1.name, " vs. ", comparisons_df$p2.name)


p_1 <- list_p_scheme[[1]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[1]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_discharge[[1]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[2]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[2]] + list_p_discharge[[2]] +
  list_p_scheme[[3]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[3]] + list_p_discharge[[3]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_1
ggsave(p_1,
  file = paste0(main_subfolder,
    "/results/figures/supplement/sub_emergency_protein_effect_comparisons_1.png"),
  height = 9,
  width = 9)

p_2 <- list_p_scheme[[4]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[4]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_discharge[[4]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[5]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[5]] + list_p_discharge[[5]] +
  list_p_scheme[[6]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[6]] + list_p_discharge[[6]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_2
ggsave(p_2,
  file = paste0(main_subfolder,
    "/results/figures/supplement/sub_emergency_protein_effect_comparisons_2.png"),
  height = 9,
  width = 9)
