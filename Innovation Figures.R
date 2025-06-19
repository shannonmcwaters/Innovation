#Plot: SRI vs. Mean Time to Solve
bee_avg_time <- innovatedatalong %>%
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

#Boxplot of innovation time for each trial
innovatedatalong$Env <- factor(innovatedatalong$Env, levels = c("s", "c"))

ggplot(innovatedatalong, aes(x = factor(trial, levels = c("Bumpy", "Folded", "Cap1", "Cap2")),
                             y = time, fill = Env)) +
  geom_boxplot(coef = Inf, position = position_dodge(width = 0.75)) +
  theme_minimal() +
  labs(x = "Trial", y = "Time to Solve") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  annotate("text", x = 1, y = max(innovatedatalong$time) + 2, label = "A", size = 6, color = "black") +
  annotate("text", x = 2, y = max(innovatedatalong$time) + 2, label = "A", size = 6, color = "black") +
  annotate("text", x = 3, y = max(innovatedatalong$time) + 2, label = "B", size = 6, color = "black") +
  annotate("text", x = 4, y = max(innovatedatalong$time) + 2, label = "B", size = 6, color = "black")


library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)#for layout
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
bee_long <- bee_avg %>%
  pivot_longer(cols = c(SRI, resp, HB10), names_to = "trait", values_to = "score") %>%
  mutate(trait = recode(trait,
                        SRI = "Routine-ness (SRI)",
                        resp = "Responsiveness",
                        HB10 = "Exploration"))

# Step 3: Make the 3-panel plot
ggplot(bee_long, aes(x = score, y = avg_time)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ trait, scales = "free_x") +
  labs(
    x = "Trait score",
    y = "Mean solving time (s)"
  ) +
  theme_minimal(base_size = 14)