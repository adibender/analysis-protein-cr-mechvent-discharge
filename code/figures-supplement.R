# global variables

surv_max <- 60

## objects
# data frames that contain hazard ratios for 6 comparisons of different nutrition
# protocols for the different models
f6                 <- readRDS(paste0(main_subfolder, "/models/m_2_A_6frames.Rds"))
f6_extubation      <- readRDS(paste0(main_subfolder, "/models/m_2_B_6frames.Rds"))

# data frame that summarizes the 6 comparisons
comparisons_df <- readRDS(paste0(main_subfolder, "/results/comparisons_df.Rds")) %>%
  mutate(comparison_text = paste0(p1.name, " vs. ", p2.name))

# model objects
mod_2_A <- readRDS(paste0(main_subfolder, "/models/m_2_A.Rds")) # cause specific: death
mod_2_B <- readRDS(paste0(main_subfolder, "/models/m_2_B.Rds")) # cause specific: extubation

### Figure S2: Fixed coefficients categorial ###
p_coef_tab_2_A_B_cat <- coefPlotGAM_cat(
  models = list(mod_2_A, mod_2_B),
  modelnames = c("Risk of death under mechanical ventilation", "Chance of successful weaning"))

ggsave(p_coef_tab_2_A_B_cat, file=paste0(main_subfolder, "/results/figures/supplement/fixed_coefs_main_categorial.png"), width =  7, height = 8)

### Figure S3: Fixed coefficients metric ###
p_coef_tab_2_A_B_met <- coefPlotGAM_met(
  models = list(mod_2_A, mod_2_B),
  modelnames = c("Risk of death under mechanical ventilation", "Chance of successful weaning"))

ggsave(p_coef_tab_2_A_B_met, file=paste0(main_subfolder, "/results/figures/supplement/fixed_coefs_main_metric.png"), width =  7, height = 8)

### Figure S4: smooth effects baseline, age, bmi ###
p_4_panel_2_A <- gg_4_panel(mod_2_A)
p_4_panel_2_B <- gg_4_panel(mod_2_B)

ggsave(p_4_panel_2_A, file = paste0(main_subfolder, "/results/figures/supplement/p_4_panel_2_A.png"), width =  7, height = 8)
ggsave(p_4_panel_2_B, file = paste0(main_subfolder, "/results/figures/supplement/p_4_panel_2_B.png"), width =  7, height = 8)

### Figure S5 & S6: Cumulative Incidence Function ###
only_low_ped <- prepare_predict_frame(
  patient,
  daily,
  ll_fun   = ll,
  scheme   = list(C2 = 0, C3 = 0),
  surv_max = surv_max,
  type     = 1,
  var1     = "proteinCat2",
  var2     = "proteinCat3")

only_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = 1, C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

low_to_med_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(0, 4), rep(1, 7)), C3 = 0),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

med_to_high_ped <- prepare_predict_frame(
  patient, daily, ll_fun=ll,
  scheme = list(C2 = c(rep(1, 4), rep(0, 7)), C3 = c(rep(0, 4), rep(1, 7))),
  surv_max = surv_max, type=1, var1="proteinCat2", var2 = "proteinCat3")

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
    get_cif_df(.x, mod_2_A, mod_2_B, comparison_name = .y)
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
      labels = c("exclusively low", "late standard", "early standard", "late high", "early high"))
  )

p_cs_death_extubation <- ggplot(
  data = res_df_cif %>% filter(cause == "death"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause-specific (death under ventilation)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")

p_cs_extubation <- ggplot(
  data = res_df_cif %>% filter(cause == "discharge"),
  aes(x = int_mid, y = cif, col = protocol)) +
  geom_line() +
  geom_ribbon(aes(ymin = cif_lower, ymax = cif_upper, fill = protocol), alpha = .3) +
  facet_wrap(~comparison) +
  ylab("CIF") + xlab("Days after ICU admission") +
  ggtitle("Cause-specific (successful weaning)") +
  labs(color = "hypothetical diet", fill = "hypothetical diet")


ggsave(p_cs_death_extubation, filename = paste0(main_subfolder, "/results/figures/supplement/p-cs-cif-death-extubation.png"),
  height = 6, width = 9)
ggsave(p_cs_extubation, filename = paste0(main_subfolder, "/results/figures/supplement/p-cs-cif-extubation.png"),
  height = 6, width = 9)


