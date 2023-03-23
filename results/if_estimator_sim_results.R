library(tidyverse)
library(knitr)
results = read_csv('results/if_estimator_sim.csv')
truth = read_csv('results/if_estimator_sim_truth.csv')


results_summary = results %>%
  pivot_longer(-sim:-method, values_to = 'est', names_to=c('name','par'), names_sep = '_') %>%
  pivot_wider(names_from=par, values_from=est) %>%
  left_join(truth %>% select(name, truth=psi_est)) %>%
  group_by(method, name) %>%
  mutate(error = mean-truth,
         bias = mean(error),
         sim_var = var(mean),
         mean_var = mean(var))



# figures -----------------------------------------------------------------



#estimates
results_summary %>%
  ggplot(aes(x=error)) +
  geom_histogram() +
  facet_grid(method ~ name, scales='free')

#variances
results_summary %>%
  ggplot(aes(x=var)) +
  geom_histogram() +
  geom_vline(aes(xintercept = sim_var)) +
  facet_grid(method ~ name, scales='free')



# tables ------------------------------------------------------------------

results_summary %>%
  select(method, bias, sim_var, mean_var) %>%
  unique() %>%
  pivot_wider(names_from=name, values_from = c(bias, sim_var, mean_var)) %>%
  set_names(c('method',
  r"(bias ($\widehat\Psi_1$))",
  r"(bias ($\widehat\Psi_2$))",
  r"($\widehat V_{sim}(\widehat\Psi_1)$)",
  r"($\widehat V_{sim}(\widehat\Psi_2)$)",
  r"($\widehat V(\widehat\Psi_1)$)",
  r"($\widehat V(\widehat\Psi_2)$)")) %>%
  kable(format = 'latex',digits = 4, escape=FALSE)

