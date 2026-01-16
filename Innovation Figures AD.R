# Shannon McWaters, James J. Kearsley, David W. Kikuchi, 
# Timothy J. Polnaszek, Anna Dornhaus
# R script for data analysis and figures for 
# Innovation project: bumble bees, inherent traits, and environment

# Paper reference: XXX (to be updated on acceptance of the paper)
### PAPER VERSION    DRAFT
### PART II: FIGURES

### Libraries ------------------
#library(patchwork) #for layout
library(scales)  # for number_format
library(tidyverse)

# Graphics setup -----------------------------

## Colorpalette
simplecomplexcolors <- c("#c7d6d5", "#d4c1e3")
simcomcolors_dark <- c("#5a9a8f", "#7b6ea8")

# Data format ----------------------------
# Reshape per-bee values for faceted plot
bee_longplot <- beedata %>%
  pivot_longer(cols = c(SRI, resp, HB10), names_to = "trait", values_to = "score") %>%
  pivot_longer(cols = c(prop_landed, prop_solved), names_to = "outcome", values_to = "prop") %>%
  mutate(trait = recode(trait,
                        SRI = "Routine formation (SRI)",
                        resp = "Responsiveness",
                        HB10 = "Exploration")
         , outcome = recode(outcome,
                          prop_landed = "Proportion landed",
                          prop_solved = "Proportion solved")
         ) %>%
  mutate(trait = factor(trait, levels = c("Responsiveness", "Routine formation (SRI)", "Exploration"))
         , outcome = factor(outcome, levels = c("Proportion landed", "Proportion solved"))
  )

# FIGURE 3 -----------------------------------

# Boxplot of innovation time for each trial
graph_data <- innovationsuccess
ymax <- max(graph_data$time)
offset <- 0.2
N_s <- table(subset(graph_data, env=="s")$trial)
N_c <- table(subset(graph_data, env=="c")$trial)

### ggplot version
ggplot(graph_data, aes(x = factor(trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2")),
                             y = time, fill = env)) +
  geom_boxplot(coef = Inf, position = position_dodge(width = 0.75)) +
  scale_fill_manual(
    values = c("s" = simplecomplexcolors[1], "c" = simplecomplexcolors[2]),
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
  annotate("text", x = 1, y = ymax + 2, label = "A", size = 5, color = "black") +
  annotate("text", x = 2, y = ymax + 2, label = "A", size = 5, color = "black") +
  annotate("text", x = 3, y = ymax + 2, label = "B", size = 5, color = "black") +
  annotate("text", x = 4, y = ymax + 2, label = "B", size = 5, color = "black") +
  annotate("text", x = 1+offset, y = - 2, label = N_c[[1]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 2+offset, y = - 2, label = N_c[[2]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 3+offset, y = - 2, label = N_c[[3]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 4+offset, y = - 2, label = N_c[[4]], size = 5, color = simcomcolors_dark[2]) +
  annotate("text", x = 1-offset, y = - 2, label = N_s[[1]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 2-offset, y = - 2, label = N_s[[2]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 3-offset, y = - 2, label = N_s[[3]], size = 5, color = simcomcolors_dark[1]) +
  annotate("text", x = 4-offset, y = - 2, label = N_s[[4]], size = 5, color = simcomcolors_dark[1]) 



# FIGURE 4 -----------------------------------
# Three panel plot generated automatically by ggplot
# Plot with formatted x-axis and reordered facets
ggplot(bee_longplot, aes(x = score, y = avg_time, color = env)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("s" = simcomcolors_dark[1], "c" = simcomcolors_dark[2]),
                     labels = c("Simple", "Complex"),
                     name = "Environment"
  ) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ trait, scales = "free_x") +
  scale_x_continuous(labels = number_format(accuracy = 0.1)) +
  labs(
    x = "Trait score",
    y = "Mean solving time for each bee (s)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# FIGURE 5 -----------------------------------

# STEP 1: Summarize landing outcomes
landing_summary <- innovationlong %>%
  filter(!is.na(landed)) %>%
  group_by(env, outcome = ifelse(landed == 1, "Landed", "Did not land")) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "Landing")

# STEP 2: Summarize solving outcomes (only for bees that landed)
solve_summary <- innovationlong %>%
  filter(landed == 1 & !is.na(solve)) %>%
  group_by(env, outcome = ifelse(solve == 1, "Solved", "Did not solve")) %>%
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
  "Did not land" = simplecomplexcolors[1],
  "Landed" = simcomcolors_dark[1],
  "Did not solve" = simplecomplexcolors[2],
  "Solved" = simcomcolors_dark[2]
)

# STEP 6: Plot
sample_sizes <- plot_data %>%
  group_by(env, category) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  mutate(label = paste0("n = ", total))

ggplot(plot_data, aes(x = env, y = count, fill = outcome)) +
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
    data = data.frame(env = "s", category = "Landing", count = 1.05),
    aes(x = env, y = count, label = "*"),
    inherit.aes = FALSE,
    size = 8
  )


# FIGURE 6 ------------------------------

## This needs work (see manuscript)

ggplot(bee_longplot, aes(x = score, y = prop, color = env)) +
  geom_point(alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black") +
  facet_grid(rows = vars(outcome), cols = vars(trait), scales = "free_x", switch = "y") +
  scale_x_continuous(labels = number_format(accuracy = 0.1), n.breaks = 5) +
  labs(x = "Trait score", y = NULL) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("s" = simcomcolors_dark[1], "c" = simcomcolors_dark[2]),
    labels = c("s" = "Simple", "c" = "Complex"),
    name = "Environment"
  ) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90, size = 14),
    strip.text.x = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.spacing = unit(1, "lines"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

# FIGURE 7 --------------------------------
# Needs work, fig doesn't make sense as-is: solving time can't be negative,
# clearly the max point is missing
graph_data <- subset(bumpy_folded_search, time>0)
ggplot(graph_data, aes(x = log(search_time), y = log(time), color = env)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "Search time",
    y = "Solving Time",
    color = "Environment"
  ) +
  scale_color_manual(values = c("s" = simcomcolors_dark[1], "c" = simcomcolors_dark[2]),
                     labels = c("s" = "Simple", "c" = "Complex"),
                     name = "Environment"
  ) +
  theme_minimal(base_size = 14)

####env and search####

ggplot(bee_avg_search, aes(x = env, y = mean_search, fill = env)) +
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
ggplot(innovationlong, aes(x = H_F1_T1, y = time, color = trial)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Handling Time on First Flower",
    y = "Time to Solve Novel Flower"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2")


# Plot: SRI vs. Mean Time to Solve
# Doesn't work - what is 'trial' supposed to be? Time is averaged over all trials.
# Changed it to env. 
ggplot(beedata, aes(x = SRI, y = avg_time)) +
  geom_point(aes(color = env), size = 3, alpha = 0.7) +  # Slightly smaller points
  geom_smooth(method = "lm", se = FALSE, aes(group = 1)) +  # One overall trend line
  scale_y_continuous(limits = c(0, 60)) +
  labs(x = "Routine formation (SRI)", y = "Mean Time to Solve") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )