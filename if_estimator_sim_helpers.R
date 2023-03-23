library(glue)
library(tidyverse)

# select_upto() is intended for nonparametric modeling. If columns are temporally ordered,
# pass df %>% select_upto(varname) when modeling varname to have the whole history in the model.
# With parallel trends estimators we typically want to drop any previous outcomes (drop_previous='y')
select_upto = function(df, varname,  ..., drop_previous='y') {
  if('id' %in% names(df)) df = df %>% select(-id)
  if('u' %in% names(df)) df = df %>% select(-u)
  varname = enexpr(varname)
  df %>% select(
    1:!!varname, #select everything prior to and including `varname`
    -(if(!is.null(drop_previous)) starts_with(drop_previous)),  #drop anything starting with `drop_previous`
    !!varname, #make sure to include `varname` (previous line will drop it if it starts with `drop_previous`)
    ... #any other variables the user would like to select (e.g., Q22t)
  )
}


nuissance_parametric = function(formulas, dists, periods=3) {
  # formulas and dists are those specified for pt_def()
  # this function will return a list of nuissance functions with glms
  # using the same kind of specification as for pt_def().
  # Thus it is convenient for fitting correctly specified parametric models
  # but also for misspecified models by changing 'formulas'

  tt = periods -1

  formulas_no_u = formulas %>%
    map(str_replace_all, ' ', '') %>%
    map(str_replace_all, '\\+u', '')

  formulas_no_u$y = paste0(formulas_no_u$y, '-a0') # a hacky fix to remove a0 - this won't work if there are interactions

  vars_to_model = map(tt:1, function(k) c(map(k:0, function(m) c(glue('Q{k}{m}t'))),
                                         map(k:0, function(m) glue('Q{k}{m}tmin1')))) %>%
    c(map(tt:1, function(k) glue('g{k}'))) %>%
    flatten_chr()

  glm_formulas = list()

  for(k in tt:1) {
    #first adding all the Qkt functions
    glm_formulas = c(glm_formulas,
                     glue('y{t} ~ ', formulas_no_u$y, .envir = list(t=k)),
                     map((k-1):0, function(m) glue('Q{k}{m+1}t ~', formulas_no_u$y, .envir=list(m=m,k=k, t=m))),
                     glue('y{t-1} ~ ', formulas_no_u$y, .envir = list(t=k)), #possibly this should be list(t=s-1). E.g. does y0 given history up to t=1 depend on (a1,w1) or (a0, w0)?
                     map((k-1):0, function(m) glue('Q{k}{m+1}tmin1 ~', formulas_no_u$y, .envir=list(m=m,k=k, t=m))))#ditto
  }
  for(k in tt:1) {
    #then adding all the g functions. This order is required for get_nuissance_fits()
    glm_formulas = c(glm_formulas,
                     glue('a{t} ~ ', formulas_no_u$a, .envir = list(t=k)))
  }

  names(glm_formulas) = vars_to_model

  glm_dists = vars_to_model %>%
    map(startsWith, 'Q') %>%
    ifelse(dists['y'], dists['a']) %>%
    map_lgl(~.=='normal') %>%
    ifelse('gaussian', 'quasibinomial') %>%
    set_names(vars_to_model)

  glm_dfs = suppressWarnings( ifelse(startsWith(vars_to_model, 'g'),
                                     paste0('df[df$a', as.numeric(str_sub(vars_to_model, 2)) - 1, '==0, ]'),
                                     'df') )

  pmap(list(glm_formulas, glm_dists, glm_dfs), function(x,y,z) parse_expr(paste0('function(df) glm(', x, ', family=',y, ', data=', z,')')) %>% eval())
}

