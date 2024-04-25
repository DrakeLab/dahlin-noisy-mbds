# Fixing the confidence interval calculations

# Load library
library(DescTools)

# Load data
jop <- read.csv("results/julia_mean_test.csv")
jop$R0 <- jop$Thv
Thv_list <- unique(jop$Thv)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[1], 0.75)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[2], 0.95)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[3], 1.05)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[4], 1.25)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[5], 2)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[6], 4)
jop$R0 <- replace(jop$R0, jop$R0 == Thv_list[7], 6.5)

# Total number of simulations
num_runs = 1000
num_batches = 100
num_sims = num_runs * num_batches

fns_list <- list(
  mean = ~ BinomCI(.x, num_sims, conf.level = 0.99995, method = "wilsoncc")[1], 
  lwr.ci = ~ BinomCI(.x, num_sims, conf.level = 0.99995, method = "wilsoncc")[2],
  upr.ci = ~ BinomCI(.x, num_sims, conf.level = 0.99995, method = "wilsoncc")[3]
)

jop2 <- select(jop, Thv, sigma, prob_end_mean, prob_out10_mean, prob_out100_mean, R0) %>%
  mutate(across(ends_with("_mean"),  ~ .x * num_sims)) %>% 
  rename_with(~gsub("[_]mean$","",.x)) %>% 
  rowwise() %>% 
  mutate(across(.cols = c(prob_end, prob_out10, prob_out100), 
                .fns = fns_list,
                .names = "{col}_{fn}"))

plot_df <- jop2 %>% select(-Thv) %>% 
  pivot_longer(cols = -c(sigma, R0)) %>% 
  mutate(output = stringr::str_extract(name, "[^_]*_[^_]*"),
         type = stringi::stri_reverse(str_split_i(stringi::stri_reverse(name), '_', i = 1))
  ) %>% 
  mutate(output = case_when(output == "prob_end" ~ "(a) Probability disease is endemic",
                            output == "prob_out10" ~ "(b) Probability of an outbreak > 10 hosts",
                            output == "prob_out100" ~ "(c) Probability of an outbreak > 100 hosts")) %>% 
  select(-name) %>%
  mutate(across(output, ~factor(., levels = c("(a) Probability disease is endemic", "(b) Probability of an outbreak > 10 hosts","(c) Probability of an outbreak > 100 hosts")))) %>% 
  pivot_wider(id_cols = c(sigma, R0, output), names_from = type) %>% 
  select(-c(end, out10, out100)) %>%
  arrange(sigma, R0)

plot <- ggplot(plot_df, aes(x = sigma)) +
  geom_line(aes(y =  mean, color = as.factor(R0), group = as.factor(R0))) +
  geom_point(aes(y =  mean, color = as.factor(R0), group = as.factor(R0))) +
  geom_ribbon(aes(ymin = `lwr.ci`, ymax = `upr.ci`, fill = as.factor(R0)), alpha = 0.3) +
  # scale_color_viridis(discrete=TRUE, guide = "legend")+
  # scale_fill_viridis(discrete=TRUE)+
  scale_color_manual(values = R0_cols, breaks = R0_vals) + 
  scale_fill_manual(values = R0_cols, breaks = R0_vals) +
  facet_wrap(~ output, ncol = 1, scales = "free", strip.position = "top") +
  theme_cowplot() +
  labs(color = unname(TeX("$R_0$")), fill = unname(TeX("$R_0$"))) +
  xlab(unname(TeX("Environmental noise strength $(\\sigma)$"))) +
  ylab("") +
  theme(strip.text.x = element_text(hjust = 0),
        strip.background = element_blank())
plot

