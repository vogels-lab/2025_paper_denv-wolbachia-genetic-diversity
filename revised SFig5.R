# ============================================================
# LIBRARIES
# ============================================================
library(tidyverse)
library(lme4)
library(ggeffects)
library(ggpubr)

# ============================================================
# REBUILD DATA CLEANLY
# ============================================================
Dataset_final <- readxl::read_xlsx("DENV1-4_2026.xlsx", sheet = "Sheet1")
#Dataset_final <- readxl::read_xlsx("Source Data File 2026.xlsx", sheet = "Fig2, Fig4 & Supplementary Fig3")


CT_stock <- Dataset_final %>%
  select(Virus, `GE/uL_stock`, Serotype, Genotype) %>%
  distinct()

infection_rate <- Dataset_final %>%
  #mutate(`infection(<38)` = as.numeric(`infection(<38)`)) %>%
  group_by(Virus, Mosquito_strain, Serotype) %>%
  summarise(
    infected    = sum(Body_Infection, na.rm = TRUE),
    sample_size = n(),
    .groups     = "drop"
  )

combined_infection_rate_stock <- infection_rate %>%
  left_join(CT_stock, by = c("Virus", "Serotype"))

# Verify: should be 60 isolates x 3 strains = 180 rows
nrow(combined_infection_rate_stock)
length(unique(combined_infection_rate_stock$Virus))

# ============================================================
# PREPARE VARIABLES
# ============================================================
combined_infection_rate_stock$Mosquito_strain <- factor(
  combined_infection_rate_stock$Mosquito_strain,
  levels = c("WT", "wAlbB", "wMel")
)

combined_infection_rate_stock <- combined_infection_rate_stock %>%
  mutate(Mosquito_strain = recode(Mosquito_strain, "wMel" = "wMelM"))



combined_infection_rate_stock$log10_titre <-
  log10(combined_infection_rate_stock$`GE/uL_stock`)

titre_mean <- mean(combined_infection_rate_stock$log10_titre, na.rm = TRUE)
titre_sd   <- sd(combined_infection_rate_stock$log10_titre,   na.rm = TRUE)

combined_infection_rate_stock$log10_titre_scaled <-
  as.numeric(scale(combined_infection_rate_stock$log10_titre))

combined_infection_rate_stock$Serotype_label <-
  paste0("DENV-", combined_infection_rate_stock$Serotype)

combined_infection_rate_stock$Serotype <- factor(
  combined_infection_rate_stock$Serotype)

# ============================================================
# FIT MODELS
# ============================================================

# --- Random intercept only (used for caterpillar plot) ---
model_ri <- glmer(
  cbind(infected, sample_size - infected) ~
    log10_titre_scaled * Mosquito_strain + Serotype +
    (1 | Virus),
  family  = binomial(link = "logit"),
  data    = combined_infection_rate_stock,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

 


model_rs_strain <- glmer(
  cbind(infected, sample_size - infected) ~
    log10_titre_scaled * Mosquito_strain + Serotype +
    (1 + Mosquito_strain | Virus),
  family  = binomial(link = "logit"),
  data    = combined_infection_rate_stock,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model_rs_strain)

# --- Same model, correlation constrained to zero ---
# Resolves singularity; used for plotting (stable CIs)
# Biologically equivalent: still allows Wolbachia inhibition to vary by isolate
model_rs_strain_nocor <- glmer(
  cbind(infected, sample_size - infected) ~
    log10_titre_scaled * Mosquito_strain + Serotype +
    (1 | Virus) + (0 + Mosquito_strain | Virus),
  family  = binomial(link = "logit"),
  data    = combined_infection_rate_stock,
  control = glmerControl(optimizer  = "nloptwrap",
                         optCtrl    = list(maxfun = 2e5))
)

summary(model_rs_strain_nocor)

isSingular(model_rs_strain_nocor)  # should be FALSE
VarCorr(model_rs_strain_nocor)

# ============================================================
# COLOUR SCHEME
# ============================================================
mosquito_strain_colors <- c("WT"    = "#33b89a",
                            "wAlbB" = "#b8cc00",
                            "wMelM"  = "#3b5323")

serotype_colors <- c("DENV-1" = "#227DCB",
                     "DENV-2" = "#b020e0",
                     "DENV-3" = "#ffce9d",
                     "DENV-4" = "#954226")


# ============================================================
# REVISED FIGURE A
# Predicted infection by titre, strain, serotype
# Uses model_rs_strain_nocor with type="fixed" for stable CIs
# ============================================================

# Generate predictions (fixed effects only)

pred <- ggpredict(
  model_rs_strain_nocor,
  terms = c("log10_titre_scaled [all]", "Mosquito_strain", "Serotype"),
  type  = "fixed"
)

# Back-transform x to original log10 scale
pred$x_orig <- pred$x * titre_sd + titre_mean

# Rename pred facet column to match observed data column name
pred$Serotype_label <- paste0("DENV-", pred$facet)

Fig_revised_A <- ggplot() +
  # Confidence ribbon
  geom_ribbon(
    data = pred,
    aes(x    = x_orig,
        y    = predicted * 100,
        ymin = conf.low  * 100,
        ymax = conf.high * 100,
        fill = group),
    alpha = 0.2
  ) +
  # Fitted lines
  geom_line(
    data = pred,
    aes(x     = x_orig,
        y     = predicted * 100,
        color = group),
    linewidth = 1.4
  ) +
  # Observed data points — matching original figure style exactly
  geom_point(
    data  = combined_infection_rate_stock,
    aes(x    = log10_titre,
        y    = 100 * infected / sample_size,
        fill = Mosquito_strain),
    color = "black",
    size  = 3.5,
    shape = 21
  ) +
  facet_wrap(~ Serotype_label, nrow = 1) +
  scale_fill_manual(
    values = mosquito_strain_colors,
    limits = c("WT", "wAlbB", "wMelM"),
    labels = c("WT", "wAlbB", "wMelM"),
    name   = "Mosquito colony"
  ) +
  scale_color_manual(
    values = mosquito_strain_colors,
    limits = c("WT", "wAlbB", "wMelM"),
    labels = c("WT", "wAlbB", "wMelM"),
    name   = "Mosquito colony"
  ) +
  scale_y_continuous(
    name   = "Infection rate (%)",
    breaks = c(0, 25, 50, 75, 100),
    limits = c(0, 100)
  ) +
  scale_x_continuous(
    name   = expression(paste("Input stock virus GE/µL (log"[10], ")")),
    breaks = seq(2, 5, by = 1)
  ) +
  guides(color = "none") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(color = "black", size = 20, face="bold"),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    strip.text = element_text(color = "black", size = 20, face="bold"),
    axis.text.x = element_text(color = "black", size = 20, face="bold"),
    axis.text.y = element_text(color = "black", size = 20, face="bold"),
    axis.title.x = element_text(color = "black", size = 20, face="bold"),
    axis.title.y = element_text(color = "black", size = 20, face="bold"),
    axis.line.x = element_line(color="black", linewidth = 1),
    axis.line.y = element_line(color="black", linewidth = 1),
    plot.title = element_text(size = 20, face="bold", hjust = 0.5)
  )  



print(Fig_revised_A)

# ============================================================
# REVISED FIGURE B
# Caterpillar plot of isolate-level random intercepts (model_ri)
# Residual between-isolate variability after accounting for
# titre and Wolbachia strain
# ============================================================

# --- Extract fixed effects ---
fe        <- fixef(model_rs_strain)
fe_wAlbB  <- fe["Mosquito_strainwAlbB"]
fe_wMelM  <- fe["Mosquito_strainwMelM"]

# --- Extract random slopes + SE ---
re_obj <- ranef(model_rs_strain, condVar = TRUE)
re     <- re_obj$Virus
pv     <- attr(re_obj$Virus, "postVar")

# pv is a 3x3xN array: dim 1/2 = (intercept, wAlbB, wMelM), dim 3 = isolate
re_df_slopes <- data.frame(
  Virus       = rownames(re),
  # Total effect = population mean + isolate deviation
  eff_wAlbB   = fe_wAlbB + re[, "Mosquito_strainwAlbB"],
  eff_wMelM   = fe_wMelM + re[, "Mosquito_strainwMelM"],
  se_wAlbB    = sqrt(pv[2, 2, ]),
  se_wMelM    = sqrt(pv[3, 3, ])
)
# Create serotype lookup table
serotype_lookup <- unique(
  combined_infection_rate_stock[, c("Virus", "Serotype_label")]
)
# --- Add serotype info ---
re_df_slopes <- merge(re_df_slopes, serotype_lookup, by = "Virus")

# --- Pivot to long format for ggplot ---
re_long <- re_df_slopes %>%
  pivot_longer(
    cols      = c(eff_wAlbB, eff_wMelM),
    names_to  = "Comparison",
    values_to = "log_odds_inhibition"
  ) %>%
  mutate(
    se = ifelse(Comparison == "eff_wAlbB", se_wAlbB, se_wMelM),
    Comparison = recode(Comparison,
                        "eff_wAlbB" = "wAlbB vs WT",
                        "eff_wMelM" = "wMelM vs WT")
  )

# --- Order viruses by mean inhibition effect ---
virus_order <- re_long %>%
  group_by(Virus) %>%
  summarise(mean_eff = mean(log_odds_inhibition)) %>%
  arrange(mean_eff) %>%
  pull(Virus)

re_long$Virus <- factor(re_long$Virus, levels = virus_order)

# --- Plot ---
wolbachia_colors <- c(
  "wAlbB vs WT" = "#b8cc00",
  "wMelM vs WT" = "#3b5323"
)

Fig_revised_B <- ggplot(
  re_long,
  aes(x     = Virus,
      y     = log_odds_inhibition,
      color = Comparison,
      shape = Comparison)
) +
  geom_hline(
    yintercept = 0, linetype = "dashed",
    color = "grey50", linewidth = 0.8
  ) +
  geom_errorbar(
    aes(ymin = log_odds_inhibition - 1.96 * se,
        ymax = log_odds_inhibition + 1.96 * se),
    width    = 0.4,
    linewidth = 0.6,
    position = position_dodge(width = 0.6)
  ) +
  geom_point(
    size     = 3,
    position = position_dodge(width = 0.6)
  ) +
  coord_flip() +
  facet_wrap(~ Serotype_label, scales = "free_y") +
  scale_color_manual(values = wolbachia_colors, name = "Comparison") +
  scale_shape_manual(
    values = c("wAlbB vs WT" = 16, "wMelM vs WT" = 17),
    name   = "Comparison"
  ) +
  labs(
    x = "Virus isolate",
    y = "Wolbachia inhibition effect\n(log-odds relative to WT, titre-adjusted)"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(color = "black", size = 20, face="bold"),
    legend.title = element_text(color = "black", size = 20, face="bold"),
    strip.text = element_text(color = "black", size = 20, face="bold"),
    axis.text.x = element_text(color = "black", size = 20, face="bold"),
    axis.text.y = element_text(color = "black", size = 20, face="bold"),
    axis.title.x = element_text(color = "black", size = 20, face="bold"),
    axis.title.y = element_text(color = "black", size = 20, face="bold"),
    axis.line.x = element_line(color="black", linewidth = 1),
    axis.line.y = element_line(color="black", linewidth = 1),
    plot.title = element_text(size = 20, face="bold", hjust = 0.5)
  )  


# ============================================================
# COMBINE AND SAVE
# ============================================================
Fig_revised <- ggarrange(
  Fig_revised_A,
  Fig_revised_B,
  labels     = c("A", "B"),
  nrow       = 2,
  heights    = c(0.4, 0.6),
  font.label = list(size = 20)
)

print(Fig_revised)


or_table <- data.frame(
  OR    = exp(fixef(model_rs_strain)),
  lower = exp(fixef(model_rs_strain) - 1.96 * sqrt(diag(vcov(model_rs_strain)))),
  upper = exp(fixef(model_rs_strain) + 1.96 * sqrt(diag(vcov(model_rs_strain))))
)
round(or_table, 3)




model_ri_serotype_interaction <- glmer(
  cbind(infected, sample_size - infected) ~
    log10_titre_scaled * Mosquito_strain +
    log10_titre_scaled * Serotype +          # titre x serotype interaction
    (1 | Virus),
  family  = binomial(link = "logit"),
  data    = combined_infection_rate_stock,
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

summary(model_ri_serotype_interaction)

AIC(model_ri_serotype_interaction, 
    model_rs_strain)

#ggsave("Fig_revised.svg", Fig_revised, width = 14, height = 14)
