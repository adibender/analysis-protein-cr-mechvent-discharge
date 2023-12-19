#************************ cif_hazard*******************************************#
# STEP 5 for sensitivity analysis
# This Script takes care of data preparation for generating hazard ratio plots
# and CIF
#******************************************************************************#
source("1_packages.R")
source("2_function_helpers.R")
#********************************* Parameters *********************************# 
# set maximal follow-up time
max.follow <- 60
# for preprocessing, set extubaton time of patients with extubation time > max.follow to max.time
max.time <- 60.1
# define number of days time-dependent covariate was provided
maxdays.nutri <- 11
# set max extubation time for cif
surv_max <- 60
#******************************************************************************#
#*
#*#*************************** loading in data ********************************#

daily   <- readRDS(paste0(main_subfolder, "/data/daily.Rds"))
patient <- readRDS(paste0(main_subfolder, "/data/patient.Rds"))
ll      <- readRDS(paste0(main_subfolder, "/data/ll.Rds"))
sens_m_2_A   <- readRDS(paste0(main_subfolder, "/models/sens_m_2_A.Rds"))
sens_m_2_B   <- readRDS(paste0(main_subfolder, "/models/sens_m_2_B.Rds"))
sens_ped_death <- readRDS(paste0(main_subfolder, "/data/sens-ped-data-death-hosp-sub5.Rds"))
sens_ped_extubation <- readRDS(paste0(main_subfolder, "/data/sens-ped-data-extubation-hosp-sub5.Rds"))

#******************* Data Preperation for Hazard Ratio Plots ******************#

sens_f6 <- make_six_frames(sens_m_2_A, patient, daily, ll,
                           var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
                           surv_max = max.follow)
sens_f6_extubation <- make_six_frames(sens_m_2_B, patient, daily, ll,
                                      var1 = "proteinCat2", var2 = "proteinCat3", type = "1", effect = "proteinCat",
                                      surv_max = max.follow)


protocol1 <- rep("< 0.8 g/kg", maxdays.nutri)
protocol2 <- c(rep("< 0.8 g/kg", 4), rep("0.8 - 1.2 g/kg", 7))
protocol3 <- rep("0.8 - 1.2 g/kg", maxdays.nutri)
protocol4 <- c(rep("0.8 - 1.2 g/kg", 4), rep("> 1.2 g/kg", 7))
protocol5 <- rep("> 1.2 g/kg", maxdays.nutri)


list_p_effect <- purrr::map(seq_along(sens_f6), ~single_plot(sens_f6[[.x]], NULL)) # list of plots Competing Risk
list_p_extubation <- purrr::map(seq_along(sens_f6_extubation),
                                ~single_plot(sens_f6_extubation[[.x]], NULL)) # list of plots Extubation

comparisons_df <-
  data.frame(
    comparison = LETTERS[1:6],
    p1 = c("protocol1", "protocol2", "protocol1", "protocol3", "protocol4", "protocol3"),
    p2 = c("protocol2", "protocol3", "protocol3", "protocol4", "protocol5", "protocol5"),
    p1.name = c("exclusively low", "late standard", "exclusively low", "early standard", "late high", "early standard"),
    p2.name = c("late standard", "early standard", "early standard", "late high", "early high", "early high"),
    p1.count = c(1, 2, 1, 3, 4, 3),
    p2.count = c(2, 3, 3, 4, 5, 5)
  ) %>%
  mutate(comparison_text = paste0(comparison, ": ", p1.name, " vs. ", p2.name)) # dataframe with interested comparisons

names(sens_f6) <- names(sens_f6_extubation) <- paste0(comparisons_df$p1.name, " vs. ", comparisons_df$p2.name)
#******************************************************************************#
#*
#************************ Data Preperation for CIF ****************************#
only_low_ped <- prepare_predict_frame(
  patient,
  daily,
  ll_fun   = ll,
  scheme   = list(C2 = 0, C3 = 0),
  surv_max = surv_max,
  type     = 1,
  var1     = "proteinCat2",
  var2     = "proteinCat3")

# data for profile: early standard
only_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = 1, C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

# data for profile: late standard
low_to_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(0, 4), rep(1, 7)), C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

# data for profile: late high
med_to_high_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(1, 4), rep(0, 7)), C3 = c(rep(0, 4), rep(1, 7))),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

# data for profile: early high
only_high_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = 0, C3 = 1),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")


# create list of comparisons
comp1 <- list(only_low_ped, low_to_med_ped)
names(comp1) <- c("low", "low_to_med")
comp2 <- list(low_to_med_ped, only_med_ped)
names(comp2)<- c("low_to_med", "med")
comp3 <-  list(only_low_ped, only_med_ped)
names(comp3)<- c("low", "med")
comp4 <-  list(only_med_ped, med_to_high_ped)
names(comp4)<- c("med", "med_to_high")
comp5 <-  list(med_to_high_ped, only_high_ped)
names(comp5)<- c("med_to_high", "high")
comp6 <-  list(only_med_ped, only_high_ped)
names(comp6)<- c("med", "high")

comparisons_list <- list(comp1, comp2, comp3, comp4, comp5, comp6)
names(comparisons_list) <- purrr::map_chr(comparisons_list, ~paste0(names(.x), collapse = " vs. "))

res_df_cif <- purrr::imap_dfr(
  .x = comparisons_list,
  .f = ~{
    get_cif_df(.x, sens_m_2_A, sens_m_2_B, comparison_name = .y)
  })

res_df_cif <- res_df_cif %>%
  mutate(comparison = factor(comparison,
                             levels = c("low vs. low_to_med", "low_to_med vs. med", "low vs. med",
                                        "med vs. med_to_high", "med_to_high vs. high", "med vs. high"),
                             labels = comparisons_df$comparison_text))
res_df_cif <- res_df_cif %>%
  mutate(
    protocol = factor(
      protocol,
      levels = c("low", "low_to_med", "med", "med_to_high", "high"),
      labels = c("exclusively low", "late standard", "early standard", "late high", "early high")))

#******************** Plot for manuscript: HR protein effect ******************#
#*
#* code for main model: figures-manuscript.R
#* 
#******************************************************************************#


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
saveRDS(list_p_scheme, paste0(main_subfolder, "/data/list_p_scheme_sens.Rds"))

p_4 <- list_p_scheme[[1]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[1]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_extubation[[1]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[2]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[2]] + list_p_extubation[[2]] +
  list_p_scheme[[3]] + xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[3]] + list_p_extubation[[3]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_4
ggsave(p_4, file = paste0(main_subfolder, "/results/figures/sens_protein_effect_comparisons_1.png"), height = 9, width = 9)

p_5 <- list_p_scheme[[4]] + xlab("Days on diet") + ylab("g protein / kg day") +
  ggtitle("comparison of\nhypothetical protein diets") + theme(plot.title = element_text(face = "bold")) +
  list_p_effect[[4]] + ggtitle("death under ventilation") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_extubation[[4]] + ggtitle("successful weaning") + theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  list_p_scheme[[5]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[5]] + list_p_extubation[[5]] +
  list_p_scheme[[6]] +  xlab("Days on diet") + ylab("g protein / kg day") +
  list_p_effect[[6]] + list_p_extubation[[6]] +
  plot_layout(nrow = 3, byrow = TRUE)
p_5
ggsave(p_5, file = paste0(main_subfolder, "/results/figures/sens_protein_effect_comparisons_2.png"), height = 9, width = 9)


#*************************** Plot for additional file *************************#
#*
#* code for main model: figures-supplement.R
#* 
#******************************************************************************#
### Figure S9: Fixed coefficients categorial ###
sens_p_coef_tab_2_A_B_cat <- coefPlotGAM_cat(
  models = list(sens_m_2_A, sens_m_2_B),
  modelnames = c("Risk of death under mechanical ventilation", "Chance of successful weaning"))

ggsave(sens_p_coef_tab_2_A_B_cat, file=paste0(main_subfolder, "/results/figures/supplement/sens_fixed_coefs_main_categorial.png"), width =  7, height = 8)

### Figure S10: Fixed coefficients metric ###
sens_p_coef_tab_2_A_B_met <- coefPlotGAM_met(
  models = list(sens_m_2_A, sens_m_2_B),
  modelnames = c("Risk of death under mechanical ventilation", "Chance of successful weaning"))

ggsave(sens_p_coef_tab_2_A_B_met, file=paste0(main_subfolder, "/results/figures/supplement/sens_fixed_coefs_main_metric.png"), width =  7, height = 8)

### Figure: smooth effects baseline, age, bmi ###
sens_p_4_panel_2_A <- gg_4_panel(sens_m_2_A)
sens_p_4_panel_2_B <- gg_4_panel(sens_m_2_B)

ggsave(sens_p_4_panel_2_A, file = paste0(main_subfolder, "/results/figures/supplement/sens_p_4_panel_2_A.png"), height = 9, width = 9)
ggsave(sens_p_4_panel_2_B, file = paste0(main_subfolder, "/results/figures/supplement/sens_p_4_panel_2_B.png"), height = 9, width = 9)

### Figure: Cumulative Incidence Function ###

sens_p_cs_death_extubation <- ggplot(
  data = res_df_cif %>% filter(cause == "death"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause-specific (death under ventilation)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")

sens_p_cs_extubation <- ggplot(
  data = res_df_cif %>% filter(cause == "discharge"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause-specific (successful weaning)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")


ggsave(sens_p_cs_death_extubation, filename = paste0(main_subfolder, "/results/figures/supplement/sens_p-cs-cif-death-extubation.png"),
       height = 6, width = 9)
ggsave(sens_p_cs_extubation, filename = paste0(main_subfolder, "/results/figures/supplement/sens_p-cs-cif-extubation.png"),
       height = 6, width = 9)


#*************************** Numbers for additional file **********************#
#*
#* code for main model: numbers-supplement.R
#* 
#******************************************************************************#
### Table of patient at risk/dying/discharged per interval ###
smry_death <- sens_ped_death %>%
  group_by(interval) %>%
  summarize(at_risk = n(), dying = sum(ped_status))
smry_extubation <- sens_ped_extubation %>%
  group_by(interval) %>%
  summarize(at_risk = n(), extubation = sum(ped_status)) %>%
  select(-at_risk)

smry_tab <- smry_death %>% left_join(smry_extubation, by = c("interval"))

readr::write_excel_csv(smry_tab, paste0(main_subfolder, "/results/sens-supplement-table-death-extubation-per-interval.csv"))

### Table Low versus standard protein intake ###
# HR per interval (death and discharge) for all comparisons
sens_f6 <- map(
  .x = sens_f6,
  .f = ~{
    .x %>%
      mutate_at(
        c("fit", "lo", "hi"),
        ~round(.x, 2)) %>%
      mutate(CI = paste0("[", lo, ", ", hi, "]"))
  }
)

sens_f6_extubation <- map(
  .x = sens_f6_extubation,
  .f = ~{
    .x %>%
      mutate_at(
        c("fit", "lo", "hi"),
        ~round(.x, 2)) %>%
      mutate(CI = paste0("[", lo, ", ", hi, "]"))
  }
)

iwalk(sens_f6, ~ readr::write_excel_csv(.x, paste0(main_subfolder, "/results/sens/sens-tab-death-", .y, ".csv")))
iwalk(sens_f6_extubation, ~ readr::write_excel_csv(.x, paste0(main_subfolder, "/results/sens/sens-tab-extubation-", .y, ".csv")))

