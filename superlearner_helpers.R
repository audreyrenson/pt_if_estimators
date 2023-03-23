sl_wrapper = function(formula, data, SL.library, family=gaussian(),  parallel=TRUE, subset, weights, ...) {
  #sl_wrapper is designed to have an interface and generally act more like a typical linear model object.
  #e.g., you pass a formula and data (this helps when including parametric regressions), you can subset and weight
  #using quoted expressions, and predict.sl_wrapper responds to NAs in the same way as predict.lm.


  if(!missing(subset)) {
    r = eval( enexpr(subset), data, parent.frame())
    data = data[r, , drop=FALSE]
  }

  y = model.frame(formula, as.data.frame(data))[,1]
  x = as.data.frame(model.matrix(formula, as.data.frame(data)))[,-1, drop=FALSE]

  if(!missing(weights)) {
    keep = rowSums(is.na(model.matrix.lm(formula, as.data.frame(data), na.action = na.pass))) == 0
    obsWeights = eval(enexpr(weights), data[keep, ])
  } else {
    obsWeights = rep(1, nrow(x))
  }

  if(parallel)
    sl = mcSuperLearner(Y = y, X =  x,  SL.library=SL.library, family=family, obsWeights = obsWeights, ...) else
      sl = SuperLearner(Y = y, X =  x,  SL.library=SL.library, family=family, obsWeights = obsWeights, ...)

  sl$formula = formula
  sl$xnames = colnames(x)
  class(sl) = c('sl_wrapper', class(sl))

  sl
}
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
}



SL.glmnet.interaction = function (Y, X, newX, family, obsWeights, id, alpha = 1, nfolds = 10,
          nlambda = 100, useMin = TRUE, loss = "deviance", ...)
{
  require("glmnet")
  if (!is.matrix(X)) {
    X <- model.matrix(~-1 + (.)^2, X)
    newX <- model.matrix(~-1 + (.)^2, newX)
  }
  fitCV <- glmnet::cv.glmnet(x = X, y = Y, weights = obsWeights,
                             lambda = NULL, type.measure = loss, nfolds = nfolds,
                             family = family$family, alpha = alpha, nlambda = nlambda,
                             ...)
  pred <- predict(fitCV, newx = newX, type = "response", s = ifelse(useMin,
                                                                    "lambda.min", "lambda.1se"))
  fit <- list(object = fitCV, useMin = useMin)
  class(fit) <- "SL.glmnet"
  out <- list(pred = pred, fit = fit)
  return(out)
}

SL.earth.morefamilies = function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3,
          nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0,
          ncross = 1, minspan = 0, endspan = 0, ...)
  #slight modification to allow other families besides guassian and binomial
{
  require("earth")
  if (family$family == "gaussian") {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
                              ncross = ncross, minspan = minspan, endspan = endspan)
  } else {
    fit.earth <- earth::earth(x = X, y = Y, degree = degree,
                              nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold,
                              ncross = ncross, minspan = minspan, endspan = endspan,
                              glm = list(family = family$family))
  }
  pred <- predict(fit.earth, newdata = newX, type = "response")
  fit <- list(object = fit.earth)
  out <- list(pred = pred, fit = fit)
  class(out$fit) <- c("SL.earth")
  return(out)
}




