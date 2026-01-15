# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT
### PART II: FIGURES

### Libraries ------------------
library(patchwork) #for layout

# Graphics setup -----------------------------

## Colorpalette
simplecomplexcolors <- c("#c7d6d5", "#d4c1e3")

# FIGURE 3 -----------------------------------

# Boxplot of innovation time for each trial
graph_data <- innovationsuccess
ymax <- max(graph_data$time)
ggplot(graph_data, aes(x = factor(trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2")),
                             y = time, fill = env)) +
  geom_boxplot(coef = Inf, position = position_dodge(width = 0.75)) +
  scale_fill_manual(
    values = c("s" = "#c7d6d5", "c" = "#d4c1e3"),
    labels = c("Simple", "Complex"),
    name = "Environment"
  )+
  theme_minimal() +
  labs(x = "Trial", y = "Time to Solve") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  annotate("text", x = 1, y = ymax + 2, label = "A", size = 6, color = "black") +
  annotate("text", x = 2, y = ymax + 2, label = "A", size = 6, color = "black") +
  annotate("text", x = 3, y = ymax + 2, label = "B", size = 6, color = "black") +
  annotate("text", x = 4, y = ymax + 2, label = "B", size = 6, color = "black")


# FIGURE 4 -----------------------------------

#Three panel plot
# Step 1: Get average solving time per bee
bee_avg <- innovatedatalong %>%
  group_by(BeeID) %>%
  summarise(
    avg_time = mean(time, na.rm = TRUE),
    SRI = first(SRI),
    resp = first(resp),
    HB10 = first(HB10)
  )

# Step 2: Reshape for faceted plot
bee_longplot <- bee_avg %>%
  pivot_longer(cols = c(SRI, resp, HB10), names_to = "trait", values_to = "score") %>%
  mutate(trait = recode(trait,
                        SRI = "Routine-ness (SRI)",
                        resp = "Responsiveness",
                        HB10 = "Exploration"))

library(scales)  # for number_format()

# Step 1: Reorder the 'trait' factor
bee_longplot <- bee_longplot %>%
  mutate(trait = factor(trait, 
                        levels = c("Responsiveness", "Routine-ness (SRI)", "Exploration")))

# Step 2: Plot with formatted x-axis and reordered facets
ggplot(bee_longplot, aes(x = score, y = avg_time)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ trait, scales = "free_x") +
  scale_x_continuous(labels = number_format(accuracy = 0.1)) +
  labs(
    x = "Trait score",
    y = "Mean solving time (s)"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  )


####barplot for land and solve####
library(dplyr)
library(ggplot2)

# STEP 1: Summarize landing outcomes
landing_summary <- bee_long %>%
  filter(!is.na(landed)) %>%
  group_by(Env, outcome = ifelse(landed == 1, "Landed", "Did not land")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Landing")

# STEP 2: Summarize solving outcomes (only for bees that landed)
solve_summary <- bee_long %>%
  filter(landed == 1 & !is.na(solve)) %>%
  group_by(Env, outcome = ifelse(solve == 1, "Solved", "Did not solve")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Solving")

# STEP 3: Combine
plot_data <- bind_rows(landing_summary, solve_summary)

# STEP 4: Set outcome as a factor to control legend order
plot_data$outcome <- factor(
  plot_data$outcome,
  levels = c("Did not land", "Landed", "Did not solve", "Solved")
)

# STEP 5: Define muted custom colors
custom_colors <- c(
  "Did not land" = "#c7d6d5",
  "Landed" = "#5a9a8f",
  "Did not solve" = "#d4c1e3",
  "Solved" = "#7b6ea8"
)

# STEP 6: Plot
sample_sizes <- plot_data %>%
  group_by(Env, category) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  mutate(label = paste0("n = ", total))

ggplot(plot_data, aes(x = Env, y = count, fill = outcome)) +
  geom_col(position = "fill", width = 0.6) +
  facet_wrap(~ category, scales = "free_x") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = c("c" = "Complex", "s" = "Simple"))+
  labs(
    x = "Environment",
    y = "Proportion of bees",
    fill = NULL
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18)
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    strip.text = element_text(size = 14)
  ) +
  geom_text(
    aes(label = count),
    position = position_fill(vjust = 0.5),
    color = "black",  # <-- black text labels
    size = 4
  ) +
  geom_text(  # Optional: asterisk annotation
    data = data.frame(Env = "s", category = "Landing", count = 1.05),
    aes(x = Env, y = count, label = "*"),
    inherit.aes = FALSE,
    size = 8
  )


####traits and land/solve ####
library(tidyverse)

# Step 1: Reshape data to long format
trait_long <- trait_summary %>%
  select(BeeID, Env,SRI, HB10, resp, prop_landed, prop_solved) %>%
  pivot_longer(cols = c(SRI, HB10, resp), names_to = "trait", values_to = "score") %>%
  pivot_longer(cols = c(prop_landed, prop_solved), names_to = "outcome", values_to = "prop") %>%
  mutate(
    trait = recode(trait,
                   SRI = "Routine-ness (SRI)",
                   HB10 = "Exploration",
                   resp = "Responsiveness"),
    outcome = recode(outcome,
                     prop_landed = "Proportion landed",
                     prop_solved = "Proportion solved")
  )
trait_long <- trait_long %>%
  mutate(
    trait = factor(trait, levels = c("Responsiveness", "Routine-ness (SRI)", "Exploration")),
    outcome = factor(outcome, levels = c("Proportion landed", "Proportion solved"))
  )
# Step 2: Plot
ggplot(trait_long, aes(x = score, y = prop, shape = Env)) +
  geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +
  facet_grid(rows = vars(outcome), cols = vars(trait), scales = "free_x", switch = "y") +
  scale_x_continuous(labels = number_format(accuracy = 0.1), n.breaks = 5) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  scale_shape_manual(
    values = c("c" = 16, "s" = 17),  # same shapes you want
    labels = c("c" = "Complex", "s" = "Simple"),
    name = "Environment"
  )+
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 14),
    strip.text.x = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

####Search time and innovation ####
ggplot(bumpy_folded_search, aes(x = search_time, y = time, color = trial)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Search time",
    y = "Solving Time",
    color = "Trial"
  ) +
  theme_minimal(base_size = 14)

####env and search####

ggplot(bee_avg_search, aes(x = Env, y = mean_search, fill = Env)) +
  geom_boxplot(width = 0.6, alpha = 0.7, coef = Inf) +
  labs(
    x = "Environment",
    y = "Mean search time per bee"
  ) +
  scale_fill_manual(
    values = c("s" = "#66c2a5", "c" = "#fc8d62"),
    labels = c("s" = "Simple", "c" = "Complex")
  ) +
  scale_x_discrete(labels = c("s" = "Simple", "c" = "Complex")) +
  theme_minimal(base_size = 20) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    legend.position = "none",
  )

#handling time first flower
ggplot(innovatedatalong, aes(x = H_F1_T1, y = time, color = trial)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Handling Time on First Flower",
    y = "Time to Solve Novel Flower"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2")


# Plot: SRI vs. Mean Time to Solve
bee_avg_time <- innovationlong %>%
  group_by(BeeID, SRI, trial) %>%
  summarise(mean_time = mean(time, na.rm = TRUE)) %>%
  ungroup()

ggplot(bee_avg_time, aes(x = SRI, y = mean_time)) +
  geom_point(aes(shape = trial), size = 3, alpha = 0.7) +  # Slightly smaller points
  geom_smooth(method = "lm", se = FALSE, color = "blue", aes(group = 1)) +  # One overall trend line
  scale_y_continuous(limits = c(0, 60)) +
  labs(x = "Routine-ness (SRI)", y = "Mean Time to Solve") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )