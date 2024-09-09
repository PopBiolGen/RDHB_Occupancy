hive.removals <- select(df, date.time, hive.removed) %>%
  st_drop_geometry() %>%
  mutate(month = month(date.time),
         year = year(date.time)) %>%
  arrange(year, month) %>%
  group_by(year, month) %>%
  summarise(effort = n(),
            catch = sum(hive.removed),
            mid.date = median(date.time)) %>%
  ungroup() %>%
  mutate(cpue = catch/effort, cum.catch = cumsum(catch)) %>%
  filter(effort > 50) # exclude months with low effort

# catch per unit effort plot
p <- ggplot(hive.removals, aes(x = cum.catch, y = cpue)) +
  geom_point() +  # Scatter plot
  geom_smooth(method = "lm", formula = y ~ x, level = 0.95) +  
  labs(x = "Cumulative removals", y = "Catch per unit effort") +
  theme_minimal() 
ggsave(filename = "out/cpue-cumulative-removals.png")

# effort over time
p <- ggplot(hive.removals, aes(x = mid.date, y = effort)) +
  geom_point() +  # Scatter plot
  labs(x = "Date", y = "Effort") +
theme_minimal() 
ggsave(filename = "out/effort-date.png")
