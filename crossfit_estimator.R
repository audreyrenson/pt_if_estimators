get_nuissance_fits = function(data,
                              nuissance_models,
                              first_colname,
                              intervene,# = function(df) df %>% mutate(across(starts_with('a'), ~0)),
                              return_data = FALSE) {

  predict.sl_wrapper = function(object, newdata, ...) {
    #this function solves two problems.
    #First, sl_wrapper allows a formula, which will lead to potential basis functions being added to the
    #data frame ultimately passed to SuperLearner, but predict.SuperLearner doesn't know this.
    #Second, if there are NAs in the dataset, some SL algorithms (like glmnet) return a vector of length of the complete cases,
    #and others (like glm) return a vector of length nrow(data) with NAs in the missing obs slots. This function
    #makes that behavior uniform, and always returns a vector of length nrow(data) in the glm-style.

    newdata =  model.matrix.lm(object$formula[-2], newdata, na.action = na.pass) #object$formula[-2] drops the outcome from the formula
    newdata = newdata[, object$xnames, drop=FALSE] #keep only columns in model fit, needed by e.g. ksvm
    preds = rep(NA, nrow(newdata))
    preds[complete.cases(newdata)] = predict.SuperLearner(object, newdata[complete.cases(newdata), ], ..., onlySL = TRUE)$pred[,1]
    preds
  } #needed because future_map cannot find predict.sl_wrapper. I think the real solution is for this to be in a package.

  #a bundle is a list containing a model fit and the data with predictions from that model, to be used by the next fit
  get_next_bundle = function(last_bundle, model,  colname) {

    last_data = last_bundle$data
    fit = model(last_data)
    predict.fun = attr(model, 'predict.fun') %>% {if(is.null(.)) predict else .}
    preds = predict.fun(fit, newdata=intervene(last_data), type='response')

    #when get_nuissance_fits is run through map(), predict returns a vector as usual,
    #but when through future_map, it returns a list with $pred and $library.pred (idk why)
    #we want $pred.
    if(is.list(preds)) preds = c(preds[[1]])


    colname = enexpr(colname)
    #could the below prediction make use of add_outsamp_pred() below? Instead of hardcording the intervention (above)
    list(data = last_data %>% mutate(!!colname := preds),
         fit = fit)
  }

  bundle_init = list(data = data,
                     fit = pseudomodel(data[[first_colname]]))

  bundles = accumulate2(.x = nuissance_models,
                        .y = parse_exprs(names(nuissance_models)),
                        .f = get_next_bundle,
                        .init=bundle_init)[-1] #nuissance_models does not have the first element, which is just a pseudomodel

  nuissance_fits = map(bundles, 'fit') %>% set_names(names(nuissance_models))

  if(return_data) attr(nuissance_fits, 'data') <- bundles[[length(bundles)]]$data

  return( nuissance_fits )

}

# A pseudomodel is just a vector packaged as an object whose predict() method returns that vector.
# The purpose is to facilitate accumulate2()ing bundles of data and models, by allowing the initial value
# (.init) to have the same form as the subsequent values.
pseudomodel = function(y, data) {
  res = match.call()
  class(res) = c('pseudomodel', class(res))
  return(res)
}
predict.pseudomodel = function(object, newdata, ...) {
  eval(object$y, envir=if(missing(newdata)) eval(object$data, envir=parent.frame()) else newdata)
}

add_outsamp_pred = function(outsamp_data,
                            nuissance_fit,
                            colname,
                            intervene = function(df) df %>% mutate(across(starts_with('a'), ~0))) {

  predict.sl_wrapper = function(object, newdata, ...) {
    #this function solves two problems.
    #First, sl_wrapper allows a formula, which will lead to potential basis functions being added to the
    #data frame ultimately passed to SuperLearner, but predict.SuperLearner doesn't know this.
    #Second, if there are NAs in the dataset, some SL algorithms (like glmnet) return a vector of length of the complete cases,
    #and others (like glm) return a vector of length nrow(data) with NAs in the missing obs slots. This function
    #makes that behavior uniform, and always returns a vector of length nrow(data) in the glm-style.

    newdata =  model.matrix.lm(object$formula[-2], newdata, na.action = na.pass) #object$formula[-2] drops the outcome from the formula
    newdata = newdata[, object$xnames, drop=FALSE] #keep only columns in model fit, needed by e.g. ksvm
    preds = rep(NA, nrow(newdata))
    preds[complete.cases(newdata)] = predict.SuperLearner(object, newdata[complete.cases(newdata), ], ..., onlySL = TRUE)$pred[,1]
    preds
  } #needed because future_map cannot find predict.sl_wrapper. I think the real solution is for this to be in a package.


  colname = enexpr(colname)
  outsamp_data %>%
    mutate(!!colname := predict(nuissance_fit, newdata=intervene(outsamp_data), type='response'))
}

# repair_g takes predictions from the exposure  models and converts them to appropriate weights
repair_g = function(df) mutate(df, g2 = a2*g1*g2 + (1-a2)*(1-g1)*(1-g2),  g1= a1*g1 + (1-a1)*(1-g1))

#calculating the influence functions
calc_if = function(df) {
  df %>% mutate(
    phi11 = ( 1*(a1==0)/g1 )*(y1 - Q11t)     + Q10t,
    phi10 = ( 1*(a1==0)/g1 )*(y0 - Q11tmin1) + Q10tmin1,
    phi22 = ( 1*(a2==0)/g2 )*(y2 - Q22t)     + (1*(a1==0)/g1)*(Q22t - Q21t)         + Q20t,
    phi21 = ( 1*(a2==0)/g2 )*(y1 - Q22tmin1) + (1*(a1==0)/g1)*(Q22tmin1 - Q21tmin1) + Q20tmin1,
    psi1 =  y0 + phi11 - phi10,
    psi2  = y0 + phi11 - phi10 + phi22 - phi21
  )
}

calc_eif_pooled = function(df, g_trunc = 0) {
  require(dtplyr)
  require(data.table)
  df_eif = df %>%
    as.data.table() %>%
    select(id, starts_with('t'), starts_with('a'), starts_with('g'), starts_with('Q'), event, lag_event) %>%
    group_by(id) %>%
    arrange(id, t) %>%
    mutate(g = cumprod(replace(g, is.na(g), 1))) %>%
    left_join(select(., id, t_lag1 = t, g_lag1 = g)) %>%
    left_join(select(., id, t_lag2 = t, g_lag2 = g)) %>%
    left_join(select(., id, t_lag3 = t, g_lag3 = g)) %>%
    left_join(select(filter(., t==71), id, event_0 = lag_event)) %>%
    mutate(Q2 = coalesce(Q2, Q3),
           Q1 = coalesce(Q1, Q2, Q3),
           Q0 = coalesce(Q0, Q1, Q2, Q3),
           across(starts_with('Q'), ~(.*2 - 1)),
           across(starts_with('g_lag'), ~ifelse(is.na(.), 1, .)),
           across(starts_with('g'), ~ifelse(.<g_trunc, g_trunc, .)),
           across(starts_with('a_lag'), ~ifelse(is.na(.), 0, .)), #this causes all irrelevant Q's to zero out of the influence function
           phi_t = (1*(a==1)/g)*(event - lag_event - Q3)
           + (1*(a_lag1==1)/g_lag1)*(Q3 - Q2)
           + (1*(a_lag2==1)/g_lag2)*(Q2 - Q1)
           + (1*(a_lag3==1)/g_lag3)*(Q1 - Q0) + Q0,
           psi_t = event_0 + cumsum(phi_t)
    ) %>%
    as_tibble()
}

#get point and variance estimates from the influence function
summarise_psi = function(df_with_if, .fns = list(mean=mean, var=var)) {
  df_with_if %>%
    select(psi1, psi2) %>%
    pivot_longer(everything(), values_to='psi') %>%
    group_by(name) %>%
    summarise(across(psi, .fns))
}

#zivich method for estimating the variance
zivich_estimator = function(df_if_fold) {
  #df_if_fold has influence function estimates and a fold id variable called 'fold'
  df_if_fold %>%
    group_by(fold) %>%
    nest() %>%
    mutate(psi_est = map(data, summarise_psi)) %>%
    unnest(psi_est) %>%
    group_by(name) %>%
    summarise(psi_est_ziv = median(psi_mean),
              psi_var_ziv = median(psi_var + (psi_mean - psi_est_ziv)^2))
}


fit_nuiss = function(df_cv, nuissance) {
  df_cv %>%
    mutate(nuissance_fits = map(splits, ~get_nuissance_fits(data=analysis(.),
                                                            nuissance_models = nuissance,
                                                            first_colname = 'y2')))
}

# # 'nuissance' is a named list of functions that take in a data frame and output
# # an object that has a predict method. Names are the column name that will be predicted from that model in the data frame.
# # For example:
#
# nuissance = list(
#   Q22t = function(df) lm(y2~ w12+w22+w32+a2, data=df),
#   Q21t = function(df) lm(Q22t ~ w11+w21+w31+a1, data=df),
#   Q20t = function(df) lm(Q21t ~ w10+w20+w30, data=df),
#   Q22tmin1 = function(df) lm(y1 ~ w12+w22+w32+a2, data=df),
#   Q21tmin1 = function(df) lm(Q22tmin1 ~ w11+w21+w31+a1, data=df),
#   Q20tmin1 = function(df) lm(Q21tmin1 ~ w10+w20+w30, data=df),
#   Q11t= function(df) lm(y1 ~ w11+w21+w31+a1, data=df),
#   Q10t = function(df) lm(Q11t ~ w10+w20+w30, data=df),
#   Q11tmin1 = function(df) lm(y0 ~ w11+w21+w31+a1, data=df),
#   Q10tmin1 = function(df) lm(Q11tmin1 ~ w10+w20+w30, data=df),
#   g1 = function(df) glm(a1 ~ w11 + w21 + w31, binomial, df),
#   g2 = function(df) glm(a2 ~ w21 + w22 + w32, binomial, df, subset = a1==0)
# )
#
# df_crossfit = vfold_cv(df, v=2)  %>%
# mutate(nuissance_fits = map(splits, ~get_nuissance_fits(data = analysis(.), nuissance_models = nuissance, first_colname = 'y2')),
#        outsamp_nuissance_estimates = map2(nuissance_fits, splits,
#                                           ~reduce2(.x = .x,
#                                                    .y = parse_exprs(names(.x)),
#                                                    .f = add_outsamp_pred,
#                                                    .init = assessment(.y)) %>%
#                                             repair_g()),
#        psi_est = map(outsamp_nuissance_estimates, estimate_psi)) %>%
#   unnest(psi_est)
