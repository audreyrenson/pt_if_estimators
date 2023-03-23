library(tidyverse)
library(simstudy)
library(rsample)
library(rlang)
library(SuperLearner)
library(furrr)
library(didsim)

nsims=3
seed=7147147 #this has to be seven digits for L'Ecuyer-CMRG seed for furrr
cv_folds = 10 #this is used by sl_wrapper()

#plan(cluster)
plan(multisession, workers=2)

# funs -----------------------------------------------------------------
source('if_estimator_sim_helpers.R')
source('superlearner_helpers.R')
source('crossfit_estimator.R')

# global variables --------------------------------------------------------
formulas = list(
  w1 = 'w1{t-1} + a{t-1}', #note no intercepts - these will be added automatically
  w2 = 'w1{t}   + w2{t-1} + a{t-1}',
  w3 = 'w1{t}   + w2{t}   + w3{t-1} + a{t-1}',
  a  = 'w1{t}   + cos(w2{t})   + I(w3{t}^2)  + u',
  y  = 'sin(w1{t})   + w2{t}*w3{t}   + a{t} + u'
)

dists =  c(
  w1 = 'normal',
  w2 = 'normal',
  w3 = 'normal',
  a  = 'binary',
  y  = 'normal'
)

variances = c(
  w1 = .1,
  w2 = .1,
  w3=  .1,
  a = 1, #this has no effect
  y = 0.1
)


#have covariates related to the outcome using e.g. sin()
#example in targeted learning book
#variance of parametric vs. more flexible estimators. How much efficiency are we losing by going nonparametric
#if we could have been unbiased parametrically?
#confirming root-n consistency under stated assumption
#Comparing estimated standard errors to true observed standard errors across simulated datasets. Does it require both models correctly specified?
#Bootstrap variance as a comparison
#Check robustness to unmeasured confounding / deviations from parallel trends? If changing how u enters the outcome model
# doesn't cause the estimator to be unbiased, something is wrong

#Scenarios (gut check - i.e. check double robust, root n consistent, variance consistent):
# 1. neither set of models correctly specified
# 2. g correct, 3. Q correct, 4. both correct

#Scenarios (efficiency)
# how much efficiency is lost when using ML when parametric was fine


#functional forms: use mixture models and nonlinear terms. e.g. covarites affect the outcome using a sin curve
# if w1 breaks some threshold then this is the effect of w2


# simulate observed data ---------------------------------------------------------------------

#generate a def object based on the above specifications
def = pt_def(formulas, dists, links=default_link(dists),
             variances, periods = 3, u_coef = 1)

#add balancing intercepts to control marginal means
def = balance_intercepts(def,
                         new_means = c(a1_conditional =0.3, a2_conditional = 0.3, y0=1, y1=1.5, y2=-1),
                         omit=c('a1','a2'),
                         default = 0)

df = genData(n = 1e4, dtDefs = def) %>% select(-ends_with('conditional'), -id, -u)

#check that balancing intercepts worked
round(colMeans(df),3)


# simulate under the intervention distribution -------------------------------------------------

def0 = intervene(def, varnames = c('a1','a2'), values = c(0,0))
df0 = genData(n=1e5, dtDefs = def0)

psi_true = tibble(name=c('psi1','psi2'), psi_est = c(mean(df0$y1), mean(df0$y2)), psi_var = c(var(df0$y1), var(df0$y2))/nrow(df0))

# setup for estimators ----------------------------------

y_lib_flex = c("SL.mean", "SL.glmnet","SL.ksvm","SL.earth")
a_lib_flex =  c("SL.mean","SL.glmnet","SL.ksvm", "SL.earth")

nuissance_markov = function(y_library, a_library) {
  list(
    Q22t = function(df) sl_wrapper(y2~ .,    df %>% select(w12:y2), y_library),
    Q21t = function(df) sl_wrapper(Q22t ~ ., df %>% select(w11:a1, Q22t), y_library),
    Q20t = function(df) sl_wrapper(Q21t ~ ., df %>% select(w10:w30, Q21t), y_library),
    Q22tmin1 = function(df) sl_wrapper(y1~ .,        df %>% select(y1:a2), y_library),
    Q21tmin1 = function(df) sl_wrapper(Q22tmin1 ~ ., df %>% select(w11:a1, Q22tmin1), y_library),
    Q20tmin1 = function(df) sl_wrapper(Q21tmin1 ~ ., df %>% select(w10:w30, Q21tmin1), y_library),
    Q11t= function(df) sl_wrapper(y1 ~ .,    df %>% select(w11:y1), y_library),
    Q10t = function(df) sl_wrapper(Q11t ~ ., df %>% select(w10:w30, Q11t), y_library),
    Q11tmin1 = function(df) sl_wrapper(y0 ~ .,       df %>% select(y0:a1), y_library),
    Q10tmin1 = function(df) sl_wrapper(Q11tmin1 ~ ., df %>% select(w10:w30, Q11tmin1), y_library),
    g1 = function(df) sl_wrapper(a1 ~ ., df %>% select(w11:a1), a_library, family=binomial()),
    g2 = function(df) sl_wrapper(a2 ~ ., df %>% select(w12:a2, a1), a_library, family=binomial())
  )
}

nuissance_nonmarkov = function(y_library, a_library) {
  list(
    Q22t = function(df) sl_wrapper(y2~ .,    df %>% select_upto(y2), y_library),
    Q21t = function(df) sl_wrapper(Q22t ~ ., df %>% select_upto(a1, Q22t), y_library),
    Q20t = function(df) sl_wrapper(Q21t ~ ., df %>% select_upto(w30, Q21t), y_library),
    Q22tmin1 = function(df) sl_wrapper(y1~ .,        df %>% select_upto(a2, y1), y_library),
    Q21tmin1 = function(df) sl_wrapper(Q22tmin1 ~ ., df %>% select_upto(a1, Q22tmin1), y_library),
    Q20tmin1 = function(df) sl_wrapper(Q21tmin1 ~ ., df %>% select_upto(w30, Q21tmin1), y_library),
    Q11t= function(df) sl_wrapper(y1 ~ .,    df %>% select_upto(y1), y_library),
    Q10t = function(df) sl_wrapper(Q11t ~ ., df %>% select_upto(w30, Q11t), y_library),
    Q11tmin1 = function(df) sl_wrapper(y0 ~ .,       df %>% select_upto(a1, y0), y_library),
    Q10tmin1 = function(df) sl_wrapper(Q11tmin1 ~ ., df %>% select_upto(w30, Q11tmin1), y_library),
    g1 = function(df) sl_wrapper(a1 ~ ., df %>% select_upto(a1), a_library, family=binomial()),
    g2 = function(df) sl_wrapper(a2 ~ ., df %>% select_upto(a2), a_library, family=binomial())
  )
}



# # e.g. for one run:
# nuissance_fits = list(
#   true = nuissance_true_parametric(formulas, dists, periods=3),
#   gfal = nuissance_true_parametric(list(y=formulas$y, a='1'), dists, periods=3),
#   qfal = nuissance_true_parametric(list(y='1', a=formulas$a), dists, periods=3),
#   bfal = nuissance_true_parametric(list(y='1', a='1'),  dists,  periods=3),
#   super = nuissance_nonmarkov(y_lib_flex, a_lib_flex) # could use nuissance_markov here too
# ) %>%
#   map(fit_nuiss, df_cv = vfold_cv(df, v=2)) # currently only using two splits. chernozhukov discusses reasons for generalizing to multiple splits
#
# psi = map(nuissance_fits, psi_from_df_nuiss)
#
#
# psi %>%
#   bind_rows(.id='method') %>%
#   select(method, name, psi_est, psi_var) %>%
#   left_join(psi_true %>% select(name, psi_true = psi_est)) %>%
#   mutate(bias=psi_est - psi_true)



# sim --------------------------------------------------------------




nuissance_specs = list(
  true = nuissance_parametric(formulas, dists, periods=3),
  gfal = nuissance_parametric(list(y=formulas$y, a='1'), dists, periods=3),
  qfal = nuissance_parametric(list(y='1', a=formulas$a), dists, periods=3),
  bfal = nuissance_parametric(list(y='1', a='1'),  dists,  periods=3)#,
  # super = nuissance_nonmarkov(y_library = c("SL.mean", "SL.glmnet",
  #                                           "SL.ksvm","SL.earth"),
  #                             a_library = c("SL.mean", "SL.glmnet",
  #                                           "SL.ksvm","SL.earth")) # why can't I put y_lib_flex here?
)

#put time trackers on each of the components

df_sims = expand_grid(n=c(1e3),
                      sim=1:nsims) %>%
  mutate(data = future_map(n,~genData(n = ., dtDefs = def) %>% select(-ends_with('conditional'), -id, -u),
                           .options = furrr_options(seed=seed)),
         data_split = future_map(data, vfold_cv, v=2, .options = furrr_options(seed=seed))) %>%
  unnest(data_split)



df_nuissance = df_sims %>%
  expand_grid(tibble(method = names(nuissance_specs),
                     nuissance_spec = nuissance_specs)) %>%
  mutate(nuissance_fits = future_map2(splits,
                                      nuissance_spec,
                                      ~get_nuissance_fits(data=analysis(.x),
                                                          nuissance_models = .y,
                                                          first_colname = 'y2',
                                                          intervene = function(df) df %>% mutate(across(starts_with('a'), ~0))),
                                      .options = furrr_options(seed=seed)))

df_psi = df_nuissance %>%
  group_by(n, sim, method) %>%
  nest() %>%
  mutate(psi = map(data, psi_from_df_nuiss)) %>%
           #future_map(data, psi_from_df_nuiss,  .options = furrr_options(seed=seed))) %>% #parallel not worth the overhead here, at least for nsims=5
  unnest(psi) %>%
  select(-data) %>%
  group_by(n, sim, method) %>%
  summarise(across(psi1:psi2, .fns = c(est=mean, var=var))) %>%
  mutate(across(ends_with('var'), ~./n))

write_csv(df_psi, file = 'if_estimator_sim.csv')
write_csv(psi_true, file='if_estimator_sim_truth.csv')



df_psi %>%
  left_join(psi_true %>% select(name, psi_true=psi_est)) %>%
  mutate(error = psi_est - psi_true) %>%
  group_by(method, name) %>%
  mutate(bias = mean(error)) %>%
  ggplot(aes(x=error)) +
  geom_histogram() +
  geom_vline(aes(xintercept=bias)) +
  facet_wrap(~method)


