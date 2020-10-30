library(tidyverse)
library(rstanarm)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

df <- read_csv("analysis/doctrinal-required/data/doctrinal-required_raw.csv") %>%
  # convert 2016-D class to 2017 since these students entered the same year
  mutate(class = str_replace(class, "2016-D", "2017"),
         class = str_replace(class, "2015-D", "2016")) %>%
  # calculate percent ranks in doctrinal classes for individual years
  group_by(class) %>%
  # calculate percentile ranks
  mutate_at(vars(contains("cum_gpa_")), list(rank = ~percent_rank(.))) %>%
  mutate(perc_rank = percent_rank(cum_gpa)) %>%
  mutate(cum_gpa_exp = scale(cum_gpa)) %>%
  ungroup() %>%
  #rename(perc_rank = cum_gpa_rank) %>%
  # center all percent rank variables
  mutate_at(vars(contains('_rank')), ~scale(., center = T, scale = F)) 

basic_mod <- stan_glmer(result ~ cum_gpa_exp + (1 | jur) + (1 | class) + state_rate_std, 
                        data = df, adapt_delta = .995,
                        family = binomial("logit"),
                        seed = 12345
                        )

sq_basic_mod <- stan_glmer(result ~ poly(cum_gpa_exp, 2) + (1 | jur) + (1 | class) + state_rate_std, 
                          data = df, adapt_delta = .995,
                          family = binomial("logit"),
                          seed = 12345)

perc_mod <- stan_glmer(result ~ perc_rank + (1 | jur) + (1 | class) + state_rate_std, 
                       data = df, adapt_delta = .995,
                       family = binomial("logit"),
                       seed = 12345)

sq_perc_mod <- stan_glmer(result ~ poly(perc_rank, 2) + (1| jur) + (1 | class) + state_rate_std, 
                         data = df, adapt_delta = .995,
                         family = binomial("logit"),
                         seed = 12345) 

loo_basic <- loo(basic_mod)
loo_sq_basic <- loo(sq_basic_mod)
loo_perc <- loo(perc_mod)
loo_sq_perc <- loo(sq_perc_mod)

loo_compare(loo_basic, loo_sq_basic, loo_perc, loo_sq_perc)

# results: model with scaled GPA does best, so we will use it.

# elpd_diff se_diff
# basic_mod     0.0       0.0   
# sq_basic_mod -1.0       0.3   
# perc_mod     -8.1       2.7   
# sq_perc_mod  -9.5       2.7   

# remove state rates
gpa_no_rate <- stan_glmer(result ~ cum_gpa_exp + (1 | jur) + (1 | class), 
                             data = df, adapt_delta = .995,
                             family = binomial("logit"),
                             seed = 12345)

loo_gpa_no_rate <- loo(gpa_no_rate)

loo_compare(loo_basic, loo_gpa_no_rate)

# model with rates does better
