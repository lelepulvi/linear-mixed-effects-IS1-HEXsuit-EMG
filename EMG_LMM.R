###############################
# LMM for EMG (2×2 design)
# Author: Emanuele Pulvirenti
# Date: 12/09/2025
###############################

# ---- 0) Packages ----
need <- c("tidyverse","lme4","lmerTest","emmeans","effectsize","broom","broom.mixed")
to_install <- need[!need %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)
invisible(lapply(need, library, character.only = TRUE))

# ---- 1) Load long-format data ----
DATA_PATH_PREF  <- "C:/Users/ep15603/OneDrive - University of Bristol/Desktop/VIVO hub/4. ICEE.space project/2. EMG data/Gait_Outcomes_Master_20250911_1103.xlsx"       # <-- EDIT ME: put your .xlsx here (use forward slashes)
SHEET_NAME <- "EMG_Stride_Long"

DATA_PATH <- if (file.exists(DATA_PATH_PREF)) DATA_PATH_PREF else DATA_PATH_ALT
if(!file.exists(DATA_PATH)){
  stop("Excel file not found. Tried: ", DATA_PATH_PREF, " and ", DATA_PATH_ALT)
}

df <- readxl::read_excel(DATA_PATH, sheet = SHEET_NAME) %>% as.data.frame()

# ---- 1a) Harmonise names & basic checks ----
# Expected: Subject, Condition, Muscle, StrideIndex, Activation
required <- c("Subject","Condition","Muscle","StrideIndex","Activation")
missing  <- setdiff(required, names(df))
if(length(missing)){
  stop("Missing required columns: ", paste(missing, collapse=", "))
}

# Coerce columns into the types the model needs:
# - Subject: factor (categorical) → enables random effects (1|Subject)
# - Condition: factor with explicit level order (so results/contrasts are stable/reproducible)
# - Muscle: factor (categorical).
# - Activation: numeric (response variable; per-stride EMG metric)
df <- df %>%
  mutate(
    Subject   = factor(Subject),
    Condition = factor(Condition,
                       levels = c("Baseline","NoIS1_HEXOff","NoIS1_HEXOn","IS1On_NoHEX","IS1On_HEXOn")),
    Muscle    = factor(Muscle),
    Activation = as.numeric(Activation)  # ensure numeric even if Excel read as character
  )
# Quick counts
cat("\n== Data summary ==\n")
print(df %>% count(Condition) %>% arrange(Condition))
print(df %>% count(Muscle) %>% arrange(Muscle))

# ---- 2) Paired t-test: Baseline vs NoIS1_HEXOff -------------------------------
# Goal: decide whether Baseline and NoIS1_HEXOff can be treated as the same OFF/OFF bucket
# Method: avoid pseudoreplication by averaging per Subject × Muscle × Condition first,
#         then run a paired t-test across those paired (Subject×Muscle) means.

pair_df <- df %>%
  filter(Condition %in% c("Baseline","NoIS1_HEXOff")) %>%      # keep only the two OFF/OFF candidates
  group_by(Subject, Muscle, Condition) %>%                      # define the pairing unit
  summarise(emg_mean = mean(Activation, na.rm = TRUE), .groups = "drop")  # average within unit

# Reshape to wide so we have two columns: Baseline and NoIS1_HEXOff, paired by Subject×Muscle
wide_pair <- pair_df %>%
  pivot_wider(names_from = Condition, values_from = emg_mean) %>%
  drop_na(Baseline, NoIS1_HEXOff)  # drop incomplete pairs

# If we have enough pairs, run the paired t-test and compute a paired Hedges' g
if (nrow(wide_pair) < 3) {
  warning("Not enough paired (Subject×Muscle) observations for a stable t-test.")
  merge_offoff <- FALSE
} else {
  tt <- t.test(wide_pair$Baseline, wide_pair$NoIS1_HEXOff, paired = TRUE)   # paired t-test
  es <- effectsize::hedges_g(wide_pair$Baseline, wide_pair$NoIS1_HEXOff, paired = TRUE)
  
  cat("\n== Baseline vs NoIS1_HEXOff (paired t-test on Subject×Muscle means) ==\n")
  print(tt)
  cat("\nEffect size (Hedges' g, paired):\n")
  print(es)
  
  # Decision rule: merge only if both statistically (p≥0.05) AND practically (|g|<0.2) similar
  merge_offoff <- (tt$p.value >= 0.05) && (abs(es$Hedges_g[1]) < 0.20)
  cat("\nMerge decision (Baseline + NoIS1_HEXOff into single OFF/OFF bucket): ",
      ifelse(merge_offoff, "YES", "NO"), "\n")
}

# ---- 3) Build clean 2×2 factors; drop DEACT -----------------------------------
# Map 5 conditions to a spacesuit × exosuit grid:
# - Baseline       → spacesuit OFF, exosuit OFF
# - NoIS1_HEXOff   → spacesuit OFF, exosuit OFF
# - NoIS1_HEXOn    → spacesuit OFF, exosuit ON
# - IS1On_NoHEX    → spacesuit ON,  exosuit OFF
# - IS1On_HEXOn    → spacesuit ON,  exosuit ON
map_2x2 <- tribble(
  ~Condition,        ~spacesuit, ~exosuit,
  "Baseline",        "OFF",      "OFF",
  "NoIS1_HEXOff",    "OFF",      "OFF",
  "NoIS1_HEXOn",     "OFF",      "ON",
  "IS1On_NoHEX",     "ON",       "OFF",
  "IS1On_HEXOn",     "ON",       "ON"
)

# Join the mapping, create factors for spacesuit/exosuit, and drop the DEACT rows
df2 <- df %>%
  left_join(map_2x2, by = "Condition") %>%                 # add spacesuit/exosuit columns
  mutate(
    spacesuit = factor(spacesuit, levels = c("OFF","ON")),# ensure consistent ordering
    exosuit   = factor(exosuit,   levels = c("OFF","ON"))
  ) %>%
  filter(exosuit %in% c("OFF","ON"))                      

# Create a merged condition label for plots/tables (not used in the model fit)
if (isTRUE(merge_offoff)) {
  df2 <- df2 %>%
    mutate(Condition_Merged = forcats::fct_collapse(Condition,
                                                    OFF_OFF = c("Baseline","NoIS1_HEXOff"),
                                                    OFF_ON  = "NoIS1_HEXOn",
                                                    ON_ON   = "IS1On_HEXOn"))
} else {
  df2 <- df2 %>%
    mutate(Condition_Merged = forcats::fct_relevel(Condition,
                                                   "Baseline","NoIS1_HEXOff","NoIS1_HEXOn","IS1On_NoHEX","IS1On_HEXOn"))
}

# Show the 2×2 cell counts so you can confirm design balance
cat("\n== 2×2 counts after mapping ==\n")
print(df2 %>% count(spacesuit, exosuit))

# ---- 4) Fit the LMM -----------------------------------------------------------
# Model structure:
#   Activation ~ exosuit*spacesuit + Muscle + (1|Subject) + (1|Subject:Muscle)
# Rationale:
# - Fixed effects: exosuit, spacesuit, and their interaction (what we care about),
#                  plus Muscle (to adjust for systematic differences between muscles).
# - Random effects: random intercepts for each Subject, and for each Muscle within Subject,
#                   which absorb baseline differences so suit effects aren't confounded.
# Coding choice:
# - Use effect coding for Muscle so suit coefficients represent *average effects across muscles*.

# Apply effect coding to Muscle (sum-to-zero contrasts)
contrasts(df2$Muscle) <- contr.sum(nlevels(df2$Muscle))

# Fit the model (REML = TRUE is standard for LMMs when you care about variance components)
#m2 <- lmer(Activation ~ exosuit*spacesuit + Muscle +
 #            (1|Subject) + (1|Subject:Muscle),
  #         data = df2, REML = TRUE)

m2 <- lmer(Activation ~ exosuit*spacesuit + Muscle +
                         (1|Subject) + (0 + Muscle|Subject),
                      data = df2, REML = TRUE)
# Summaries: fixed effects (coefficients, SEs, t, p) and random effects variances
cat("\n== LMM summary ==\n")
print(summary(m2))

# F-tests / p-values for fixed effects (Type II/III-like from lmerTest)
cat("\n== Tests of fixed effects (Type II/III style) ==\n")
print(anova(m2))

# ---- 5) Marginal means & contrasts via emmeans --------------------------------
# Marginal means for exosuit (averaging over spacesuit and muscles)
emm_exo    <- emmeans(m2, ~ exosuit)
# Marginal means for spacesuit (averaging over exosuit and muscles)
emm_suit   <- emmeans(m2, ~ spacesuit)
# Cell means for the full 2×2 (exosuit × spacesuit)
emm_inter  <- emmeans(m2, ~ exosuit * spacesuit)

# Print the marginal means and their pairwise contrasts
cat("\n== Marginal means: exosuit ==\n");   print(emm_exo)
cat("\nPairwise contrast (exosuit):\n");    print(pairs(emm_exo))

cat("\n== Marginal means: spacesuit ==\n"); print(emm_suit)
cat("\nPairwise contrast (spacesuit):\n");  print(pairs(emm_suit))

# Print the 2×2 cell means and simple effects (aka conditional effects)
cat("\n== Cell means (exosuit × spacesuit) and simple effects ==\n")
print(emm_inter)
cat("\nSimple effects of exosuit within spacesuit levels:\n")
print(contrast(emm_inter, interaction = "pairwise", by = "spacesuit"))
cat("\nSimple effects of spacesuit within exosuit levels:\n")
print(contrast(emm_inter, interaction = "pairwise", by = "exosuit"))


# ---- 6) Quick diagnostics ------------------------------------------------------
# Basic numerical checks on residuals; plots are commented for non-interactive runs
cat("\n== Quick diagnostics ==\n")
res <- resid(m2)                   # model residuals
fitted_vals <- fitted(m2)          # model fitted values
cat("Residual summary:\n"); print(summary(res))
cat("Correlation(resid, fitted): ", suppressWarnings(cor(res, fitted_vals)), "\n")

# If you run interactively, you can uncomment for visuals:
# par(mfrow=c(1,2))
# qqnorm(res); qqline(res)         # should be roughly straight if residuals ~ normal
# plot(fitted_vals, res,
#      main="Residuals vs Fitted", xlab="Fitted", ylab="Residuals"); abline(h=0,lty=2)




# ---- 7) Interpretation crib sheet ---------------------------------------------
cat("\n================ HOW TO READ THE RESULTS ================\n")
cat("
A) Paired t-test (Baseline vs NoIS1_HEXOff):
   - We averaged within Subject×Muscle first to avoid step-level pseudoreplication.
   - If p ≥ 0.05 and |Hedges g| < 0.20 → they behave similarly; narrate as one OFF/OFF bucket.
   - If not, report them separately in descriptions; the 2×2 model is still valid.

B) LMM:
   - exosuitON: average change in Activation with exosuit ON vs OFF,
     averaged across muscles (effect coding) and accounting for Subject/Muscle baselines.
   - spacesuitON: average change with spacesuit ON vs OFF.
   - exosuitON:spacesuitON: whether the exosuit effect is different when the spacesuit is ON (synergy/attenuation).

C) emmeans:
   - ~ exosuit and ~ spacesuit: marginal means + OFF vs ON contrasts (clean main effects).
   - ~ exosuit*spacesuit: the four cell means (OFF/OFF, ON/OFF, OFF/ON, ON/ON).
     'by=' contrasts are simple effects (e.g., exosuit effect when spacesuit is ON only).

D) Diagnostics:
   - Residuals should look roughly normal and show no strong pattern vs fitted values.
   - If skewed/heteroscedastic, consider transforming Activation (e.g., log1p)
     and then request emmeans on the response scale: emmeans(..., type='response').

E) Reporting:
   - Provide estimates, 95% CIs, and p-values for main effects and interaction.
   - Translate into physiology: e.g., a positive, significant exosuit main effect suggests
     higher muscle activation on average with resistive loading.
")

# ============================
# Pretty EMM plots (SVG, white bg)
# ============================
if (!dir.exists("outputs")) dir.create("outputs", recursive = TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(emmeans)
})

# Ensure SVG device is available
if (!"svglite" %in% rownames(installed.packages())) install.packages("svglite")
library(svglite)

# Robust CI extractor (always yields lower.CL / upper.CL)
safe_emm_df <- function(emm_obj, level = 0.95) {
  df <- as.data.frame(confint(emm_obj, level = level))
  nm <- names(df)
  if (all(c("asymp.LCL","asymp.UCL") %in% nm)) {
    df <- dplyr::rename(df, lower.CL = asymp.LCL, upper.CL = asymp.UCL)
  } else if (all(c("LCL","UCL") %in% nm)) {
    df <- dplyr::rename(df, lower.CL = LCL, upper.CL = UCL)
  }
  if ("exosuit"   %in% nm) df$exosuit   <- factor(df$exosuit,   levels = c("OFF","ON"))
  if ("spacesuit" %in% nm) df$spacesuit <- factor(df$spacesuit, levels = c("OFF","ON"))
  df
}

# Theme & palette (explicit white background)
theme_pub <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background  = element_rect(fill = "grey95", colour = NA),
      strip.text        = element_text(face = "bold"),
      plot.title        = element_text(face = "bold"),
      axis.title        = element_text(face = "bold"),
      panel.background  = element_rect(fill = "white", colour = NA),
      plot.background   = element_rect(fill = "white", colour = NA)
    )
}
pal2 <- c("OFF" = "#1b9e77", "ON" = "#d95f02")

# Build plotting data
stopifnot(exists("emm_inter"), exists("emm_exo"), exists("emm_suit"))
df_cells <- safe_emm_df(emm_inter) %>% filter(is.finite(emmean))
df_exo   <- safe_emm_df(emm_exo)    %>% filter(is.finite(emmean))
df_suit  <- safe_emm_df(emm_suit)   %>% filter(is.finite(emmean))

# Helpful to check you actually have rows:
# print(df_cells); print(df_exo); print(df_suit)

# 1) 2×2 cell means
p_cells <- ggplot(df_cells, aes(x = exosuit, y = emmean,
                                group = spacesuit, colour = spacesuit)) +
  geom_line(position = position_dodge(width = 0.25), linewidth = 0.9) +
  geom_point(position = position_dodge(width = 0.25), size = 2.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1, linewidth = 0.6,
                position = position_dodge(width = 0.25)) +
  scale_colour_manual(values = pal2, name = "Spacesuit") +
  # Small values: label nicely
  scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
  labs(x = "Exosuit", y = "EMG (estimated marginal mean)",
       title = "EMG EMMs by Exosuit × Spacesuit (averaged across muscles)") +
  theme_pub()

ggsave("outputs/fig_emm_cells_2x2.svg", p_cells, width = 7, height = 4, bg = "white", device = "svg")
ggsave("outputs/fig_emm_cells_2x2.png", p_cells, width = 7, height = 4, dpi = 300, bg = "white")

# 2) Marginal means: Exosuit
p_exo <- ggplot(df_exo, aes(x = exosuit, y = emmean, colour = exosuit)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.6) +
  scale_colour_manual(values = pal2, guide = "none") +
  scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
  labs(x = "Exosuit", y = "EMG (EMM)", title = "Marginal means: Exosuit") +
  theme_pub()

ggsave("outputs/fig_emm_marginal_exosuit.svg", p_exo, width = 4, height = 3.2, bg = "white", device = "svg")
ggsave("outputs/fig_emm_marginal_exosuit.png",  p_exo, width = 4, height = 3.2, dpi = 300, bg = "white")

# 3) Marginal means: Spacesuit
p_suit <- ggplot(df_suit, aes(x = spacesuit, y = emmean, colour = spacesuit)) +
  geom_point(size = 2.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.08, linewidth = 0.6) +
  scale_colour_manual(values = pal2, guide = "none") +
  scale_y_continuous(labels = scales::label_scientific(digits = 2)) +
  labs(x = "Spacesuit", y = "EMG (EMM)", title = "Marginal means: Spacesuit") +
  theme_pub()

ggsave("outputs/fig_emm_marginal_spacesuit.svg", p_suit, width = 4, height = 3.2, bg = "white", device = "svg")
ggsave("outputs/fig_emm_marginal_spacesuit.png",  p_suit, width = 4, height = 3.2, dpi = 300, bg = "white")

