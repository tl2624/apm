#' Generate models used to fit outcomes
#' 
#' @description `apm_mod()` generates a list of models characterized by a basic model formulas and other options (e.g., lags, families, etc.) that are supplied to [apm_pre()]. These values are completely crossed to create a grid of model specifications, and multiple sets of model specifications can be combined using `c()` (see Examples).
#' 
#' @param formula_list a list of model formulas with the outcome on the left side and predictions (or just an intercept) on the right side.
#' @param family a list of family specifications; see [family()] for allowable options. These will eventually be passed to [glm()] when fitting the models in [apm_pre()]. `"negbin"` can also be supplied to request a negative binomial model with a log link fit using [MASS::glm.nb()]. Default is `"gaussian"` to specify a linear model.
#' @param lag a vector of integers indicating the desired outcome lags to be used as predictors. For example, a `lag` value of 3 means the outcome lagged once, twice, and three times will be included as predictors. Default is 0 for no lags.
#' @param diff_k a vector of integers indicating the desired outcome lag to be used a an offset For example, a `diff_k` value of 1 means the prior time point's outcome will be included as an offset, equivalent to using the outcome minus its corresponding lag as the outcome of the corresponding model. Default is 0 for no lags. Any models with a `diff_k` value less than a `lag` value will be removed automatically. When used with a family with a log link, the lags are automatically log-transformed; an error will be thrown by `apm_pre()` if nonpositive values are present in the outcome.
#' @param log a logical vector indicating whether the outcome should be log-transformed. Default is `FALSE` to use the original outcome. When `lag` or `diff_k` are greater than 0, the outcome lags will also be log-transformed if `TRUE`. When the family has a log link and `diff_k` is greater than zero, the lag in the offset will be log transformed.
#' @param time_trend a vector of integers indicating the desired powers to be included in a time trend. For example, a `time_trend` value of 2 means the time variable and its square will be included as predictors in the model. A value of 0 (the default) means time is not included as a predictor.
#' @param fixef a logical vector indicating whether unit fixed effects should be included as predictors. Default is `FALSE` to omit unit fixed effects.
#' @param identiy_only_log `logical`; whether to omit any models in which `log` is `TRUE` but the link in the `family` specification corresponds to something other than `"identity"`. Default is `TRUE`, and this should probably not be changed.
#' 
#' @returns
#' An `apm_models` object, which is a list containing the full cross (less any omitted combinations) of the model features specified in the arguments, with each combination a list. These have a `print()` method and can be combined using `c()`. Each model is named automatically, but these can be set manually using [names()] as well. Models can be removed by setting their value to `NULL`; see Examples.
#' 
#' @seealso [formula()], [family()]
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 1 baseline formulas,
#' # 2 families, 2 lags, 2 time trends
#' models1 <- apm_mod(crude_rate ~ 1,
#'                    family = list("gaussian", "quasipoisson"),
#'                    time_trend = 0:1,
#'                    lag = 0:1, fixef = TRUE)
#' models1
#' 
#' # Add a single other model with a square time trend
#' models2 <- apm_mod(crude_rate ~ 1,
#'                    family = "gaussian",
#'                    time_trend = 2,
#'                    fixef = FALSE)
#' models2
#' 
#' (models <- c(models1, models2))
#' 
#' # Remove a model
#' models[[4]] <- NULL
#' models


#' @export 
apm_mod <- function(formula_list, family = "gaussian", lag = 0L, diff_k = 0L,
                    log = FALSE, time_trend = 0L, fixef = FALSE, identiy_only_log = TRUE) {
  # Check arguments
  
  ## Check formula_list
  chk::chk_not_missing(formula_list, "`formula_list`")
  
  if (inherits(formula_list, "formula")) {
    formula_list <- list(formula_list)
  }
  else if (!is.list(formula_list) || !all(vapply(formula_list, inherits, logical(1L), "formula"))) {
    chk::err("`formula_list` must be a list of model formulas")
  }
  
  formula_list <- unique(formula_list)
  
  if (any(lengths(formula_list) < 3L)) {
    chk::err("all formulas in `formula_list` must have left-hand-side (outcome) variable")
  }
  
  ## Check family
  if (.okay_family(family)) {
    family <- list(family)
  }
  else if (!is.list(family) ||
           !all(vapply(family, .okay_family, logical(1L)))) {
    chk::err("`family` must be a list of model families")
  }
  
  family <- unique(family)
  
  family <- lapply(family, function(f) {
    if (is.character(f)) {
      if (f %in% c("negbin", "negative.binomial", "Negative Binomial")) {
        return(list(family = "Negative Binomial", link = "log"))
      }
      
      f <- get(f, mode = "function", envir = parent.frame(2L))
    }
    
    if (is.function(f)) {
      f <- f()
    }
    
    f
  })
  
  ## Check lag
  chk::chk_not_any_na(lag)
  chk::chk_whole_numeric(lag)
  chk::chk_gte(lag, 0)
  
  lag <- sort(unique(lag))
  
  ## Check diff_k
  chk::chk_not_any_na(diff_k)
  chk::chk_whole_numeric(diff_k)
  chk::chk_gte(diff_k, 0)
  
  diff_k <- sort(unique(diff_k))
  
  if (any(diff_k > 0)) {
    if (max(diff_k) <= min(lag))
      chk::wrn("`diff_k` will be ignored because all supplied values are less than the smallest value supplied to `lag`")
  }
  
  # Check log
  chk::chk_not_any_na(log)
  chk::chk_logical(log)
  
  log <- sort(unique(log))
  
  ## Check time_trend
  chk::chk_not_any_na(time_trend)
  chk::chk_whole_numeric(time_trend)
  chk::chk_gte(time_trend, 0)
  
  time_trend <- sort(unique(time_trend))
  
  # Check fixef
  chk::chk_not_any_na(fixef)
  chk::chk_logical(fixef)
  
  fixef <- sort(unique(fixef))
  
  grid <- expand.grid(formula = seq_along(formula_list),
                      family = seq_along(family),
                      lag = seq_along(lag),
                      diff_k = seq_along(diff_k),
                      log = seq_along(log),
                      time_trend = seq_along(time_trend),
                      fixef = seq_along(fixef))
  
  if (any(log)) {
    chk::chk_flag(identiy_only_log)
    
    ## Remove all combinations that involve log = TRUE with a non-identity link (probably invalid)
    if (identiy_only_log) {
      identity_links <- vapply(family, function(f) {
        f[["link"]] == "identity"
      }, logical(1L))
      
      if (!all(identity_links)) {
        grid <- grid[grid$family %in% which(identity_links) | !log[grid$log], , drop = FALSE]
      }
    }
  }
  
  if (any(diff_k > 0) && any(lag > 0)) {
    #Drop combinations where diff_k is less than lag because those are equivalent to having no diff_k
    grid <- grid[diff_k[grid$diff_k] == 0 | diff_k[grid$diff_k] > lag[grid$lag], , drop = FALSE]
  }
  
  out <- lapply(seq_len(nrow(grid)), function(i) {
    list(formula = formula_list[[grid$formula[i]]],
         family = family[[grid$family[i]]],
         lag = lag[[grid$lag[i]]],
         diff_k = diff_k[[grid$diff_k[i]]],
         log = log[[grid$log[i]]],
         time_trend = time_trend[[grid$time_trend[i]]],
         fixef = fixef[[grid$fixef[i]]])
  })
  
  for (i in which(duplicated(out))) {
    out[[i]] <- NULL
  }
  
  names(out) <- .name_mods(out)
  
  class(out) <- "apm_models"
  
  out
}

#' @exportS3Method print apm_models
print.apm_models <- function(x, ...) {
  for (i in seq_along(x)) {
    cat(sprintf("- Model %s: %s\n", i, names(x)[i]))
    cat(deparse1(x[[i]]$formula), "\n", sep = "")
    cat(sprintf("family: %s(link = %s)\n",
                x[[i]]$family$family,
                .add_quotes(x[[i]]$family$link, 2L)))
    cat(sprintf("outcome lag: %s\n", if (x[[i]]$lag == 0) "none" else toString(seq_len(x[[i]]$lag))))
    cat(sprintf("outcome diff: %s\n", if (x[[i]]$diff_k == 0) "none" else x[[i]]$diff_k))
    cat(sprintf("log outcome: %s\n", if (x[[i]]$log) "yes" else "no"))
    cat(sprintf("time trend: %s\n",
                if (x[[i]]$time_trend == 0) "none"
                else if (x[[i]]$time_trend == 1) "linear"
                else if (x[[i]]$time_trend == 2) "quadratic"
                else if (x[[i]]$time_trend == 3) "cubic"
                else sprintf("%s-degree polynomial", x[[i]]$time_trend)))
    cat(sprintf("unit fixed effects: %s\n", if (x[[i]]$fixef) "yes" else "no"))
    if (i != length(x)) cat("\n")
  }
  
  invisible(x)
}

#' @exportS3Method c apm_models
c.apm_models <- function(..., recursive = TRUE) {
  out <- NextMethod("c")
  
  for (i in which(duplicated(out))) {
    out[[i]] <- NULL
  }
  
  names(out) <- .name_mods(out)
  
  class(out) <- "apm_models"
  
  out
}

#' @exportS3Method `[` apm_models
`[.apm_models` <- function(..., recursive = TRUE) {
  out <- NextMethod("[")
  
  class(out) <- "apm_models"
  
  out
}

#' @exportS3Method `[[` apm_models
`[[.apm_models` <- function(..., recursive = TRUE) {
  NextMethod("[[")
}

.name_mods <- function(models) {
  ff <- unlist(lapply(models, function(m) {
    tt <- terms(m$formula)
    if (!is_null(attr(tt, "term.labels", TRUE))) {
      deparse1(tt)
    }
  }))
  
  id_models <- !is_null(ff) && length(unique(ff)) > 1L
  
  ff2 <- lapply(models, function(m) {
    fam <- m$family
    unlist(fam[c("family", "link")])
  })
  
  id_family <- length(unique(ff2)) > 1L
  
  if (id_family) {
    for (m in seq_along(models)) {
      if (ff2[[m]]["family"] == "Negative Binomial") {
        ff2[[m]]["family"] <- "NB"
        okay <- ff2[[m]]["link"] == "log"
      }
      else {
        famfun <- get0(ff2[[m]]["family"], mode = "function")
        okay <- {
          if (is_null(famfun)) FALSE
          else ff2[[m]]["link"] == eval(formals(famfun)[["link"]])
        }
      }
      
      if (okay || length(ff2[[m]]) == 1L) {
        ff2[[m]] <- ff2[[m]][1L]
      }
      else {
        ff2[[m]] <- sprintf("%s(link = %s)", ff2[[m]]["family"], ff2[[m]]["link"])
      }
    }
    
    id_family <- length(unique(ff2)) > 1L
  }
  
  vapply(seq_along(models), function(m) {
    model <- models[[m]]
    
    n <- character(0L)
    
    formula <- terms(model$formula)
    
    if (!is_null(attr(formula, "term.labels", TRUE))) {
      n <- {
        if (id_models) sprintf("baseline model %s", match(deparse1(formula), ff))
        else "baseline model"
      }
    }
    else if (model$time_trend == 0 && model$lag == 0 && !model$fixef) {
      n <- {
        if (attr(formula, "intercept", TRUE) == 0) "empty model"
        else "baseline mean"
      }
    }
    
    if (model$time_trend > 0) {
      n <- c(n, sprintf("%s trend", 
                        if (model$time_trend == 1) "linear"
                        else if (model$time_trend == 2) "quadratic"
                        else if (model$time_trend == 3) "cubic"
                        else sprintf("%s-degree polynomial", .ordinal(model$time_trend))))
    }
    
    if (model$lag > 0) {
      n <- c(n, sprintf("AR(%s)", model$lag))
    }
    
    if (model$fixef) {
      n <- c(n, "FE")
    }
    
    n2 <- character(0L)
    
    if (model$log) {
      n2 <- c(n2, "log")
    }
    
    if (model$diff_k > 0) {
      n2 <- c(n2, sprintf("%s diff", .ordinal(model$diff_k)))
    }
    
    if (id_family) {
      n2 <- c(n2, .firstup(ff2[[m]]))
    }
    
    n <- paste(n, collapse = " + ")
    
    if (!is_null(n2)) {
      n <- sprintf("%s (%s)", n, toString(n2))
    }
    
    n
  }, character(1L))
}
