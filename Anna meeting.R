innovation_subset <- innovation %>%
  group_by(trial) %>%
  mutate(rank = rank(time, ties.method = "first")) %>%
  ungroup()  # Ungroup after the operation
innovation_subset_all <- innovation_subset %>%
  group_by(BeeID) %>%
  filter(n_distinct(trial) == 4) %>%  # Keep only bees with 4 unique trials
  ungroup() 

innovation_subset_vw<- innovation_subset %>%
  filter(trial %in% c("v", "w"))
innovation_subset_vw <- innovation_subset_vw %>%
  group_by(BeeID) %>%
  filter(n_distinct(trial) == 2) %>%  # Keep only bees with 4 unique trials
  ungroup() 

innovation_subset_xy <- innovation_subset %>%
  filter(trial %in% c("x", "y"))
innovation_subset_xy <- innovation_subset_xy %>%
  group_by(BeeID) %>%
  filter(n_distinct(trial) == 2) %>%  # Keep only bees with 4 unique trials
  ungroup() 

innovation_subset_x <- innovation %>%
  filter(trial %in% "x")
summary(lm(log(time) ~ resp + env, data=innovation_subset_x))


ggplot(innovation_subset_all, aes(x = trial, y = time, group = BeeID, color = as.factor(env))) +
  geom_line() +
  geom_point() +
  scale_y_reverse() +  # Reverse the y-axis if lower ranks are better
  labs(x = "Trial",
       y = "Time to access nectar",
       color = "Environment") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(innovation_subset_vw, aes(x = trial, y = time, group = BeeID, color = as.factor(env))) +
  geom_line() +
  geom_point() +
  scale_y_reverse() +  # Reverse the y-axis if lower ranks are better
  labs(x = "Trial",
       y = "Time to access nectar",
       color = "Environment") +
  theme_minimal() +
  theme(legend.position = "right")

ggplot(innovation_subset_xy, aes(x = trial, y = time, group = BeeID, color = as.factor(env))) +
  geom_line() +
  geom_point() +
  scale_y_reverse() +  # Reverse the y-axis if lower ranks are better
  labs(x = "Trial",
       y = "Time to access nectar",
       color = "Environment") +
  theme_minimal() +
  theme(legend.position = "right")










innovation_reshaped <- innovation_subset_xy %>%
  select(BeeID, trial, time) %>%
  pivot_wider(names_from = trial, values_from = time)

regression_model <- lm(y ~ x, data = innovation_reshaped)
summary(regression_model)
plot(innovation_reshaped$x~innovation_reshaped$y)


solve_env_summary <- env_chi %>%
  group_by(trial, env) %>%
  summarize(
    solved = sum(solve == 1),
    not_solved = sum(solve == 0)
  )

solve_env_summary

solve_env_summary <- env_chi %>%
  group_by(BeeID, env) %>%
  summarize(
    solved_x = sum(trial == "x" & solve == 1),
    solved_y = sum(trial == "y" & solve == 1)
  ) %>%
  mutate(
    solve_both = ifelse(solved_x > 0 & solved_y > 0, 1, 0),
    solve_only_x = ifelse(solved_x > 0 & solved_y == 0, 1, 0),
    solve_only_y = ifelse(solved_y > 0 & solved_x == 0, 1, 0),
    solve_neither = ifelse(solved_x == 0 & solved_y == 0, 1, 0)
  )

# Count the number of bees in each solving category within each environment
solution_counts_by_env <- solve_env_summary %>%
  group_by(env) %>%
  summarize(
    both_solved = sum(solve_both),
    only_x_solved = sum(solve_only_x),
    only_y_solved = sum(solve_only_y),
    neither_solved = sum(solve_neither)
  )

solution_counts_by_env

####

# Calculate the per-bee average time to solve for each environment
avg_time_per_bee <- innovation_subset_xy %>%
  group_by(BeeID, env) %>%
  summarize(avg_time = mean(time), .groups = "drop")

test_result <- wilcox.test(avg_time ~ env, data = avg_time_per_bee)

boxplot(avg_time_per_bee$avg_time~avg_time_per_bee$env)

##########
complex = innovatedata  %>%
  filter(Env %in% "c")
simple = innovatedata  %>%
  filter(Env %in% "s")
wilcox.test(complex$resp,simple$resp)
boxplot(complex$resp,simple$resp)
t.test(complex$resp,simple$resp)
shapiro.test(log(complex$resp+1))
shapiro.test(log(simple$resp + 1))




innovation_subset_x <- innovatedatalong %>%
  filter(trial %in% "x")
summary(lm(log(time) ~ resp + env, data=innovation_subset_x))

coxph(Surv(time, solve) ~ env + resp, data = innovation_subset_x)
