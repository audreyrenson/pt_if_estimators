library(glue)
library(tidyverse)



# sim functions -----------------------------------------------------------

pt_def = function(formulas, dists, links = default_link(dists), variances = rep(1, length(formulas)), periods, u_coef = 1) {

  stopifnot(length(formulas) == length(dists) & length(dists==length(variances)))

  varnames <- names(dists) <- names(links) <- names(variances) <- names(formulas)


  if(!check_varnames(varnames)) stop('variable names in formula argument must be exactly a, y, and w{i} where i are numbers')
  stopifnot(length(u_coef) %in% c(1, periods)) #allowing u_coef to vary for violations of parallel trends

  n_covs =  varnames %>% str_starts('w') %>% sum()

  #t=0 variables
  def = defData(varname = 'u', dist = 'normal', formula = 0, variance = 1)
  for(j in 1:n_covs) {
    def = defData(def,
                  varname = as.character(glue('w{j}{t}', .envir = list(t = 0, j=j))),
                  dist = dists[glue('w{j}')],
                  formula = glue_formula(formulas[[glue('w{j}')]], glue_vars = list(t=0)) %>% add_coefs_to_formula(),
                  variance = variances[glue('w{j}')])
  }

  def = defData(def,
                varname = as.character(glue('a{t}', .envir = list(t = 0))),
                dist = 'nonrandom', formula = 0)
  def = defData(def,
                varname = as.character(glue('y{t}', .envir = list(t=0))),
                dist = dists['y'],
                link = links['y'],
                formula = glue_formula(formulas$y, glue_vars = list(t=0)) %>% add_coefs_to_formula(u_coef = u_coef[1]),
                variance = variances[glue('y')])

  #t = 1:tau variables
  for(k in 1:(periods - 1)) {
    #w first
    for(varname in names(formulas))
    {
      def = defData(def,
                    varname = paste0(varname, k, if(varname=='a') '_conditional'),#first generate what A{k} should be conditional on A{k-1}=0
                    formula = glue_formula(formulas[varname], glue_vars = list(t=k)) %>%
                      add_coefs_to_formula(u_coef = if(varname !='y') NULL else if(length(u_coef)==1) u_coef else u_coef[k + 1]),
                    dist = dists[varname],
                    link = links[varname],
                    variance = variances[varname])

      if(varname == 'a')  def=defData(def,  #then A{k} marginal
                                      varname = paste0(varname, k),
                                      formula = glue_formula('0 + a{t-1} + (1-a{t-1})*a{t}_conditional', glue_vars = list(t=k)),
                                      dist = 'nonrandom')
    }
  }

  return(def)

}

check_varnames = function(varnames) {
  a_varnames = varnames[str_starts(varnames, 'a')]
  y_varnames = varnames[str_starts(varnames, 'y')]
  w_varnames = varnames[str_starts(varnames, 'w')]

  #check that all variable names start with w, a, or y
  starts_with_way = varnames %>% str_starts('w|a|y') %>% all()

  #contains all w, a, and y
  contains_way = 0 < length(a_varnames)*length(y_varnames)*length(w_varnames)

  #check that w variables are of form w{n} where n is a number
  w_form_correct = w_varnames %>% str_match('[^w+[:digit:]]') %>% is.na() %>% all()

  #check that a and y are just a and y
  a_form_correct = a_varnames %>% str_match('[^a]') %>% is.na() %>% all()
  y_form_correct = y_varnames %>% str_match('[^y]') %>% is.na() %>% all()

  return (starts_with_way & contains_way & w_form_correct & a_form_correct & y_form_correct)
}


balance_intercepts = function(def, new_means, omit, default=0, round_digits=3) {

  get_balancing_intercept = function(formula, mean, var_means, link) {

    if(!grepl('\\+', formula) & link=='identity') {
      return(mean) #in this case the formula does not have covariates, just return the mean
    }
    if(!grepl('\\+', formula) & link=='logit') return(plogis(mean))

    mean_replacer = set_names(as.character(var_means), names(var_means))

    formula_nointercept =  str_split(formula, '\\+', n=2, simplify = TRUE)[,2]
    formula_nointercept_with_means = str_replace_all(formula_nointercept, mean_replacer)

    balancer = if(link=='identity') mean else if(link=='logit') -log(1/mean - 1) else stop('link must be "logit" or "identity"')

    new_intercept = balancer - eval(parse(text = formula_nointercept_with_means))

    return( if(is.infinite(round_digits)) new_intercept else round(new_intercept, round_digits) )

  }

  replace_intercept = function(formula, new_intercept) {
    if(!grepl('\\+', formula)) {
      return(new_intercept)
    } else {
      split_formula = str_split(formula, pattern='\\+', n=2, simplify=TRUE)
      split_formula[,1] = new_intercept
      return(paste(split_formula, collapse='+') )
    }
  }

  indices_keep = if(!missing(omit)) seq_along(def$varname)[-match(omit, def$varname)] else seq_along(def$varname)
  varnames = def$varname[indices_keep]
  formulas = def$formula[indices_keep] %>% set_names(varnames)
  links = def$link[indices_keep] %>% set_names(varnames)

  #fill in defaults
  new_mean_vector = rep(default, length(def$varname)) %>% set_names(def$varname)
  new_mean_vector[names(new_means)] = new_means

  #currently this will not get the means exactly right for anything with ak in it
  #it will just treat those as zeros. Probably not far off and good enough for our purposes.
  balancing_intercepts = sapply(varnames, function(v) get_balancing_intercept(formulas[v],
                                                                              mean = new_mean_vector[v],
                                                                              var_means = new_mean_vector,
                                                                              link = links[v])) %>%
    set_names(varnames)



  new_formulas = sapply(varnames, function(v) replace_intercept(formulas[v], balancing_intercepts[v])) %>%
    set_names(varnames)

  def$formula[indices_keep] = new_formulas
  def
}

# aux functions -----------------------------------------------------------



glue_formula = function(string, glue_vars) {

  trimmed_string = stringr::str_remove_all(string, '[:space:]')
  glued_trimmed_string = glue::glue(trimmed_string, .envir = glue_vars)

  remove_negatives(glued_trimmed_string)
}

remove_negatives = function(glued_string) {
  #this is to get rid of references to variables before time 0
  result = glued_string %>%
    stringr::str_remove_all('\\+[:alpha:]+-[:digit:]+') %>%#first remove those starting with +
    stringr::str_remove_all('\\*[:alpha:]+-[:digit:]+') %>%#or a *
    stringr::str_remove_all('\\:[:alpha:]+-[:digit:]+') %>%#etc.
    stringr::str_remove_all('[:alpha:]+-[:digit:]+') %>%
    stringr::str_remove_all('\\+[:alpha:]+[:digit:]+-[:digit:]+') %>%#first remove those starting with +
    stringr::str_remove_all('\\*[:alpha:]+[:digit:]+-[:digit:]+') %>%#or a *
    stringr::str_remove_all('\\:[:alpha:]+[:digit:]+-[:digit:]+') %>%#etc.
    stringr::str_remove_all('[:alpha:]+[:digit:]+-[:digit:]+')

  if(result == '') return(character(0)) else return(result)
}

add_coefs_to_formula = function(chr_formula, coefs = NULL, u_coef = NULL){

  if(length(chr_formula) == 0) {

    return(generate_coefficients(p = 1))

  } else {

    terms = terms.glue_formula(chr_formula)

    if(is.null(coefs)) {
      coefs = generate_coefficients(p = length(terms) + 1) #+1 for an intercept
    }

    if(!is.null(u_coef)) {
      coefs[ which(terms == 'u') + 1] = u_coef
    }

    stopifnot(length(coefs) == length(terms) + 1)
    coefs_times_terms = paste0(coefs, c('', paste0('*', terms)))

    return(paste(coefs_times_terms, collapse=' + '))

  }

}
terms.glue_formula = function(glued_string) {
  glued_string %>%
    paste0('~', .) %>%
    as.formula() %>%
    terms() %>%
    attr('term.labels') %>%
    str_replace_all(':', '*')
}

generate_coefficients = function(p) round(rnorm(p), 2)

default_link = function(dist) case_when(dist == 'normal' ~ 'identity',
                                        dist %in% c('binomial','binary') ~ 'logit',
                                        TRUE ~ 'identity')

intervene = function(def, varnames, values) {
  stopifnot(length(varnames) == length(values))
  for(i in seq_along(varnames)) {
    def = updateDef(def, changevar = varnames[i], newformula = values[i], newdist = 'nonrandom')
  }
  return(def)
}

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

