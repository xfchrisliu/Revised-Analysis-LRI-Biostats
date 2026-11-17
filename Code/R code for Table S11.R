
# ============================================================
# Longitudinal Body Weight Analysis
# Linear Mixed-Effects Models
#
# Outcome:
#   Body Weight
#
# ------------------------------------------------------------
# 1. Combined GF + MPF analysis
#
# Fixed effects:
#   Microbiota * Diet + Time
#
# Model:
#   Weight ~ Microbiota * Diet + Time + (1 | Mouse)
#
# Includes:
#   Microbiota
#   Diet
#   Time
#   Microbiota × Diet
#
# Does NOT include:
#   Microbiota × Time
#   Diet × Time
#   Microbiota × Diet × Time
#
# ------------------------------------------------------------
# 2. GF-only analysis
#
# Model:
#   Weight ~ Diet + Time + (1 | Mouse)
#
# ------------------------------------------------------------
# 3. MPF-only analysis
#
# Model:
#   Weight ~ Diet + Time + (1 | Mouse)
# ============================================================


# ============================================================
# 0. Load packages
# ============================================================

library(readxl)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)


# ============================================================
# 1. GF data
# ============================================================

gf_raw <- read_excel(
  "GF_bodyweight.xlsx",
  col_names = FALSE
)

# Time
gf_time <- as.numeric(
  gf_raw[-1, 1][[1]]
)

# Diet for each mouse
gf_diet <- as.character(
  gf_raw[1, -1]
)

# Weight matrix
gf_weight <- gf_raw[-1, -1]

gf_weight <- gf_weight %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  )

# Mouse names
colnames(gf_weight) <- paste0(
  "GF_mouse_",
  1:ncol(gf_weight)
)


# ============================================================
# 2. GF long format
# ============================================================

gf_long <- gf_weight %>%
  mutate(
    Time = gf_time
  ) %>%
  pivot_longer(
    cols = starts_with("GF_mouse"),
    names_to = "Mouse",
    values_to = "Weight"
  )


# ============================================================
# 3. GF mouse-diet mapping
# ============================================================

gf_mouse_info <- tibble(
  Mouse = colnames(gf_weight),
  Diet = as.character(gf_diet)
)

gf_long <- gf_long %>%
  left_join(
    gf_mouse_info,
    by = "Mouse"
  ) %>%
  mutate(
    Microbiota = "GF"
  )


# ============================================================
# 4. MPF data
# ============================================================

mpf_raw <- read_excel(
  "MPF_bodyweight.xlsx",
  col_names = FALSE
)

# Time
mpf_time <- as.numeric(
  mpf_raw[-1, 1][[1]]
)

# Diet
mpf_diet <- as.character(
  mpf_raw[1, -1]
)

# Weight matrix
mpf_weight <- mpf_raw[-1, -1]

mpf_weight <- mpf_weight %>%
  mutate(
    across(
      everything(),
      as.numeric
    )
  )

# Mouse names
colnames(mpf_weight) <- paste0(
  "MPF_mouse_",
  1:ncol(mpf_weight)
)


# ============================================================
# 5. MPF long format
# ============================================================

mpf_long <- mpf_weight %>%
  mutate(
    Time = mpf_time
  ) %>%
  pivot_longer(
    cols = starts_with("MPF_mouse"),
    names_to = "Mouse",
    values_to = "Weight"
  )


# ============================================================
# 6. MPF mouse-diet mapping
# ============================================================

mpf_mouse_info <- tibble(
  Mouse = colnames(mpf_weight),
  Diet = as.character(mpf_diet)
)

mpf_long <- mpf_long %>%
  left_join(
    mpf_mouse_info,
    by = "Mouse"
  ) %>%
  mutate(
    Microbiota = "MPF"
  )


# ============================================================
# 7. Combine GF + MPF
# ============================================================

bw <- bind_rows(
  gf_long,
  mpf_long
)


# ============================================================
# 8. Clean variables
# ============================================================

bw <- bw %>%
  mutate(

    # Germ status
    Microbiota = factor(
      Microbiota,
      levels = c("GF", "MPF")
    ),

    # Diet
    Diet = factor(
      Diet,
      levels = c(
        "palm",
        "borage",
        "echium",
        "fish",
        "high linoleic safflower",
        "high oleic safflower"
      )
    ),

    # Mouse
    Mouse = factor(
      Mouse
    ),

    # Keep numeric time
    Time_numeric = as.numeric(Time),

    # Treat time as categorical
    Time = factor(
      Time_numeric,
      levels = sort(
        unique(Time_numeric)
      )
    )
  )


# ============================================================
# 9. Check data structure
# ============================================================

cat("\n============================================================\n")
cat("SAMPLE SIZE\n")
cat("============================================================\n")

print(
  table(
    bw$Microbiota,
    bw$Diet
  )
)


cat("\n============================================================\n")
cat("TIME POINTS\n")
cat("============================================================\n")

print(
  table(
    bw$Time
  )
)


cat("\n============================================================\n")
cat("MICROBIOTA × DIET × TIME\n")
cat("============================================================\n")

print(
  table(
    bw$Microbiota,
    bw$Diet,
    bw$Time
  )
)


cat("\n============================================================\n")
cat("NUMBER OF MICE\n")
cat("============================================================\n")

print(
  length(
    unique(bw$Mouse)
  )
)


# ============================================================
# 10. Check repeated measurements per mouse
# ============================================================

mouse_time_check <- bw %>%
  group_by(Mouse) %>%
  summarise(
    N_timepoints = n_distinct(Time),
    .groups = "drop"
  )

cat("\n============================================================\n")
cat("NUMBER OF TIME POINTS PER MOUSE\n")
cat("============================================================\n")

print(
  table(
    mouse_time_check$N_timepoints
  )
)


# ============================================================
# 11. COMBINED ANALYSIS: GF + MPF
#
# Chris requested:
#   Microbiota
#   Diet
#   Time
#   Microbiota × Diet
#
# No interactions involving Time.
# ============================================================

lmm_combined <- lmer(
  Weight ~
    Microbiota * Diet +
    Time +
    (1 | Mouse),
  data = bw,
  REML = FALSE
)


# ============================================================
# 12. Combined model summary
# ============================================================

cat("\n============================================================\n")
cat("COMBINED LMM SUMMARY\n")
cat("============================================================\n")

print(
  summary(lmm_combined)
)


# ============================================================
# 13. Combined model ANOVA
# ============================================================

cat("\n============================================================\n")
cat("COMBINED LMM ANOVA\n")
cat("============================================================\n")

anova_combined <- anova(
  lmm_combined,
  type = 3
)

print(
  anova_combined
)


# ============================================================
# 14. Combined analysis: key p-values
# ============================================================

cat("\n============================================================\n")
cat("COMBINED ANALYSIS: KEY P-VALUES\n")
cat("============================================================\n")

combined_key_results <- anova_combined[
  c(
    "Microbiota",
    "Diet",
    "Time",
    "Microbiota:Diet"
  ),
  ,
  drop = FALSE
]

print(
  combined_key_results
)


# ============================================================
# 15. Combined analysis: estimated marginal means
#
# Overall Microbiota × Diet comparison,
# averaged over Time.
# ============================================================

emm_micro_diet <- emmeans(
  lmm_combined,
  ~ Microbiota * Diet
)

cat("\n============================================================\n")
cat("COMBINED EMM: MICROBIOTA × DIET\n")
cat("============================================================\n")

print(
  emm_micro_diet
)


# ============================================================
# 16. Combined analysis:
#     GF vs MPF within each Diet
# ============================================================

microbiota_within_diet <- emmeans(
  lmm_combined,
  pairwise ~ Microbiota | Diet,
  adjust = "tukey"
)

cat("\n============================================================\n")
cat("GF vs MPF WITHIN EACH DIET\n")
cat("============================================================\n")

print(
  microbiota_within_diet$contrasts
)


# ============================================================
# 17. GF-ONLY ANALYSIS
#
# Model:
#   Weight ~ Diet + Time + (1 | Mouse)
# ============================================================

bw_gf <- bw %>%
  filter(
    Microbiota == "GF"
  )


lmm_gf <- lmer(
  Weight ~
    Diet +
    Time +
    (1 | Mouse),
  data = bw_gf,
  REML = FALSE
)


# ============================================================
# 18. GF-only model summary
# ============================================================

cat("\n============================================================\n")
cat("GF-ONLY LMM SUMMARY\n")
cat("============================================================\n")

print(
  summary(lmm_gf)
)


# ============================================================
# 19. GF-only ANOVA
# ============================================================

cat("\n============================================================\n")
cat("GF-ONLY LMM ANOVA\n")
cat("============================================================\n")

anova_gf <- anova(
  lmm_gf,
  type = 3
)

print(
  anova_gf
)


# ============================================================
# 20. GF-only key p-values
# ============================================================

cat("\n============================================================\n")
cat("GF-ONLY: KEY P-VALUES\n")
cat("============================================================\n")

gf_key_results <- anova_gf[
  c(
    "Diet",
    "Time"
  ),
  ,
  drop = FALSE
]

print(
  gf_key_results
)


# ============================================================
# 21. GF-only estimated marginal means for Diet
#
# Averaged over Time.
# ============================================================

emm_gf_diet <- emmeans(
  lmm_gf,
  ~ Diet
)

cat("\n============================================================\n")
cat("GF-ONLY EMM: DIET\n")
cat("============================================================\n")

print(
  emm_gf_diet
)


# ============================================================
# 22. GF-only pairwise Diet comparisons
# ============================================================

gf_diet_pairs <- pairs(
  emm_gf_diet,
  adjust = "tukey"
)

cat("\n============================================================\n")
cat("GF-ONLY: PAIRWISE DIET COMPARISONS\n")
cat("============================================================\n")

print(
  gf_diet_pairs
)


# ============================================================
# 23. MPF-ONLY ANALYSIS
#
# Model:
#   Weight ~ Diet + Time + (1 | Mouse)
# ============================================================

bw_mpf <- bw %>%
  filter(
    Microbiota == "MPF"
  )


lmm_mpf <- lmer(
  Weight ~
    Diet +
    Time +
    (1 | Mouse),
  data = bw_mpf,
  REML = FALSE
)


# ============================================================
# 24. MPF-only model summary
# ============================================================

cat("\n============================================================\n")
cat("MPF-ONLY LMM SUMMARY\n")
cat("============================================================\n")

print(
  summary(lmm_mpf)
)


# ============================================================
# 25. MPF-only ANOVA
# ============================================================

cat("\n============================================================\n")
cat("MPF-ONLY LMM ANOVA\n")
cat("============================================================\n")

anova_mpf <- anova(
  lmm_mpf,
  type = 3
)

print(
  anova_mpf
)


# ============================================================
# 26. MPF-only key p-values
# ============================================================

cat("\n============================================================\n")
cat("MPF-ONLY: KEY P-VALUES\n")
cat("============================================================\n")

mpf_key_results <- anova_mpf[
  c(
    "Diet",
    "Time"
  ),
  ,
  drop = FALSE
]

print(
  mpf_key_results
)


# ============================================================
# 27. MPF-only estimated marginal means for Diet
#
# Averaged over Time.
# ============================================================

emm_mpf_diet <- emmeans(
  lmm_mpf,
  ~ Diet
)

cat("\n============================================================\n")
cat("MPF-ONLY EMM: DIET\n")
cat("============================================================\n")

print(
  emm_mpf_diet
)


# ============================================================
# 28. MPF-only pairwise Diet comparisons
# ============================================================

mpf_diet_pairs <- pairs(
  emm_mpf_diet,
  adjust = "tukey"
)

cat("\n============================================================\n")
cat("MPF-ONLY: PAIRWISE DIET COMPARISONS\n")
cat("============================================================\n")

print(
  mpf_diet_pairs
)


# ============================================================
# 29. OPTIONAL: Time-specific EMMs
#
# These are descriptive/follow-up comparisons only.
# They are NOT required for Chris's requested model.
#
# Uncomment if you want to inspect the model estimates
# at individual time points.
# ============================================================

# emm_micro_diet_time <- emmeans(
#   lmm_combined,
#   ~ Microbiota * Diet | Time
# )
#
# print(emm_micro_diet_time)


# ============================================================
# 30. OPTIONAL: Longitudinal plot
#
# This is descriptive only and does not change the model.
# ============================================================

ggplot(
  bw,
  aes(
    x = Time_numeric,
    y = Weight,
    group = Mouse,
    color = Microbiota
  )
) +
  geom_line(
    alpha = 0.25
  ) +
  geom_point(
    alpha = 0.4
  ) +
  stat_summary(
    aes(
      group = Microbiota
    ),
    fun = mean,
    geom = "line",
    linewidth = 1.2
  ) +
  facet_wrap(
    ~ Diet
  ) +
  labs(
    x = "Time",
    y = "Body Weight",
    color = "Germ Status"
  ) +
  theme_bw()


# ============================================================
# END OF ANALYSIS
# ============================================================
