# Fit one model; mod: output of .modify_formula_and_data()
.fit_one_model <- function(mod, weights = NULL, time_var, val_time, family) {
  
  
  # Effectively subset without dropping any observations
  if (is_null(weights)) {
    weights <- rep.int(1, nrow(mod$data))
  }
  
  is.na(weights)[mod$data[[time_var]] >= val_time] <- TRUE
  
  # Model fitting function; note: need do.call() to correctly process `weights`
  fit_fun <- {
    if (identical(family$family, "Negative Binomial"))
      function(.formula, .data, .family, .weights = NULL) {
        do.call("glm.nb", list(.formula, data = .data,
                               weights = .weights,
                               na.action = "na.exclude",
                               x = TRUE, y = TRUE))
      }
    else if (identical(family$family, "gaussian") &&
             identical(family$link, "identity"))
      function(.formula, .data, .family, .weights = NULL) {
        do.call("lm", list(.formula, data = .data,
                           weights = .weights,
                           na.action = "na.exclude",
                           x = TRUE, y = TRUE))
      }
    else 
      function(.formula, .data, .family, .weights = NULL) {
        out <- do.call("glm", list(.formula, data = .data,
                                   family = .family,
                                   weights = .weights,
                                   na.action = "na.exclude",
                                   x = TRUE, y = TRUE))
        
        out$call$family <- str2lang(sprintf('%s("%s")',
                                            .family$family,
                                            .family$link))
        out
      }
  }
  
  out <- fit_fun(.formula = mod$formula,
                 .data = mod$data,
                 .family = family,
                 .weights = weights)
  
  out$call$data <- quote(.data)
  out$call$na.action <- NULL
  out$call$weights <- quote(.weights)
  
  out
}

# Modifies formula and dataset to account for model features
.modify_formula_and_data <- function(model, data, group_var, unit_var, time_var) {
  
  formula <- model$formula
  
  # Create log and lagged variables
  outcome_name <- as.character(formula[[2L]])
  
  if (model$time_trend == 1) {
    formula <- update(formula, sprintf(". ~ . + %s", time_var))
  }
  else if (model$time_trend > 0) {
    formula <- update(formula, sprintf(". ~ . + poly(%s, %s)", time_var, model$time_trend))
  }
  
  # Add interaction with group
  formula <- update(formula, sprintf(". ~ %s * (.)", group_var))
  
  outcome <- NULL
  # Log outcome if requested
  if (model$log) {
    outcome <- model.response(model.frame(update(formula, . ~ 1), data = data))
    if (min(outcome) <= 0) {
      chk::err("`log` cannot be `TRUE` when the outcome takes on values of 0 or lower")
    }
    
    formula <- update(formula, log(.) ~ .)
  }
  
  if (model$diff_k > 0 && model$family$link == "log") {
    if (is_null(outcome)) {
      outcome <- model.response(model.frame(update(formula, . ~ 1), data = data))
    }
    
    if (min(outcome) <= 0) {
      chk::err("no model can have a log link and an outcome lag when the outcome takes on values of 0 or lower")
    }
  }
  
  # Add lagged outcome or offset thereof if requested
  if (model$lag > 0 || model$diff_k > 0) {
    # 
    # if (lag_outcome_name %in% all.vars(formula)) {
    #     chk::wrn(sprintf("the variable named %s will be replaced. Give this variable a different name before running",
    #                      lag_outcome_name))
    # }
    
    for (i in seq_len(max(model$diff_k, model$lag))) {
      if (i > model$lag && i != model$diff_k) next
      
      lag_outcome_name <- sprintf("%s_lag_%s", outcome_name, i)
      
      lag_i <- data[[outcome_name]]
      is.na(lag_i)[] <- TRUE
      
      beginning <- seq_len(i)
      for (u in levels(data[[unit_var]])) {
        #We can lag here because data is ordered by time_var already
        #Note: assumes complete time series for each unit, uses previous value in dataset (ignoring actual
        #      value of time var)
        ending <- sum(data[[unit_var]] == u) + 1 - seq_len(i)
        lag_i[data[[unit_var]] == u][-beginning] <- data[[outcome_name]][data[[unit_var]] == u][-ending]
      }
      
      data[[lag_outcome_name]] <- lag_i
      
      if (i <= model$lag) {
        #Add lag as predictor
        formula <- {
          if (model$log) {
            update(formula, sprintf(". ~ . + log(%s)", lag_outcome_name))
          }
          else {
            update(formula, sprintf(". ~ . + %s", lag_outcome_name))
          }
        }
      }
      else if (i == model$diff_k) {
        #Add lag as offset
        formula <- {
          if (model$log || model$family$link == "log") {
            update(formula, sprintf(". ~ . + offset(log(%s))", lag_outcome_name))
          }
          else {
            update(formula, sprintf(". ~ . + offset(%s)", lag_outcome_name))
          }
        }
      }
    }
  }
  
  # Add unit fixed effects if requested, remove group var main effect
  if (model$fixef) {
    formula <- update(formula, sprintf(". ~ . + %s - %s", unit_var, group_var))
  }
  
  list(formula = formula, data = data)
}

.get_y <- function(models, data) {
  model.response(model.frame(update(models[[1L]]$formula, . ~ 1), data = data))
}

# Prepares validation data to be used in .predict_quick() to compute predictions
.make_predict_prep <- function(fit, newdata) {
  tt <- terms(fit)
  Terms <- delete.response(tt)
  
  m <- model.frame(Terms, newdata,
                   na.action = na.pass, 
                   xlev = fit$xlevels)
  
  cl <- attr(Terms, "dataClasses")
  if (!is_null(cl)) {
    .checkMFClasses(cl, m)
  }
  
  X <- model.matrix(Terms, m, contrasts.arg = fit$contrasts)
  
  offset <- model.offset(m)
  addO <- fit$call$offset
  if (!is_null(addO)) {
    addO <- eval(addO, newdata, environment(tt))
    offset <- if (is_null(offset)) addO else offset + addO
  }
  
  piv <- fit$qr$pivot[seq_len(fit$rank)]
  
  list(X = X[, piv, drop = FALSE],
       offset = offset)
}

# Computes predictions quickly from coefficients (beta), data (output of
# .make_predict_prep()), and inverse link (fit$family$linkinv)
.predict_quick <- function(beta, predict_prep, linkinv = NULL) {
  
  if (is_null(linkinv)) {
    linkinv <- identity
  }
  
  eta <- drop(predict_prep$X %*% beta)
  
  if (length(predict_prep$offset) > 0) {
    eta <- eta + predict_prep$offset
  }
  
  linkinv(eta)
}

.refit_with_weights <- function(fit, weights) {
  if (length(fit$na.action) > 0) {
    weights <- weights[-fit$na.action]
  }
  
  if (class(fit)[1L] == "lm") {
    x <- model.matrix(fit)
    y <- {
      if (is_null(fit[["y"]])) model.response(model.frame(fit))
      else fit[["y"]]
    }
    
    new_fit <- lm.wfit(x = x, y = y, w = weights, offset = fit$offset)
    
    fit <- .list_modify(fit, new_fit)
  }
  else if (class(fit)[1L] == "glm") {
    x <- model.matrix(fit)
    y <- {
      if (is_null(fit[["y"]])) model.response(model.frame(fit))
      else fit[["y"]]
    }
    
    start <- fit$coefficients
    
    new_fit <- do.call(fit$method,
                       list(x = x, y = y, weights = weights,
                            start = start, offset = fit$offset,
                            family = fit$family,
                            control = fit$control))
    
    fit <- .list_modify(fit, new_fit)
  }
  else {
    fit <- do.call("update", list(fit, weights = weights))
  }
  
  fit
}