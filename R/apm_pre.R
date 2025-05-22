#' Fit validation models to pre-treatment data
#' 
#' @description `apm_pre()` fits models to the pre-treatment data to compute the observed prediction errors for each model in each period and compute the Bayesian model averaging (BMA) weights eventually used in [apm_est()] to estimate the treatment effect.
#' 
#' @param models an `apm_models` object; the output of a call to [apm_mod()].
#' @param data a dataset containing all the variables named in the supplied models (i.e., the outcome and any predictors) as well as any variable named below.
#' @param weights an optional vector of weights (e.g., sampling weights) used to fit weighted regression models.
#' @param group_var string; the name of the treatment variable in `data` defining the "to be treated" and "not to be treated" groups. The corresponding variable should take on values of 0 and 1 only.
#' @param unit_var string; the name of the unit ID variable in `data`.
#' @param time_var string; the name of the variable in `data` containing the time variable.
#' @param val_times a numeric vector corresponding to the pre-treatment times that will be used as validation times when select the model with the optimal average expected prediction error.
#' @param nsim the number of simulated draws from the joint posterior of the fitted models to use to compute the BMA weights. Default is 1000. More is better but takes longer.
#' @param cl a cluster object created by [parallel::makeCluster()], or an integer to indicate number of child-processes (integer values are ignored on Windows) for parallel evaluations. It can also be `"future"` to use a future backend. `NULL` (default) refers to sequential evaluation. See the `cl` argument of [pbapply::pblapply()] for details.
#' @param verbose `logical`; whether to print information about the progress of the estimation, including a progress bar. Default is `TRUE`.
#' @param object an `apm_pre_fit` object; the output of a call to `apm_pre()`.
#' @param order how to order the summary; `NULL` (the default) is the same ordering as the models supplied to `apm_pre()`, `"weights"` orders the models by their computed BMA weights with the highest weights on top, and `"errors"` orders the models by their maximum absolute difference in prediction errors with the smallest errors on top.
#' @param \dots ignored.
#' 
#' @returns
#' `apm_pre()` returns an `apm_pre_fits` object, which is a list containing the models supplied to `models`, a grid of all fitted models, a list of all model fit objects, a list of all estimated coefficients, the joint covariance of the coefficients, the dataset supplied to `data`, and other components supplied to `apm_pre()`.
#' 
#' `summary()` produces a data frame containing the BMA weights and maximum absolute difference in mean prediction errors for each model, ordered according `order`. An asterisk appears next to the model with the smallest error.
#' 
#' @details
#' `apm_pre()` creates a grid of all models and all time points and fits all corresponding models. For each validation time supplied to `val_times`, each model is fit using all previous times. For example, for a validation time of 5, a model is fit with data only from periods 1-4.
#' 
#' [lm()], [glm()], or [MASS::glm.nb()] are used to fit the given models. The joint covariance matrix of all the coefficients is computed using the SUEST method described in Mize et al. (2019, p164), which is also used by the STATA command `suest`. This is equivalent to the covariance matrix computed by stacking the score equations for the models and fitting them using M-estimation and yields the equivalent of the HC0 covariance matrix for all within-model covariances. The covariance is clustered by `unit_id`.
#' 
#' To compute the BMA weights, random variate drawn from a multivariate normal distribution `nsim` times with mean vector equal to the concatenated model coefficients and covariance equal to the joint covariance matrix described above. For each iteration, the absolute average prediction errors are calculated for each model and validation period. A model is considered the "winner" if it its largest absolute average prediction error across validation periods is the smallest among all models. The BMA weight for each model is equal to the proportion of iterations in which that model was the "winner".
#' 
#' @seealso [lm()],[glm()], and [MASS::glm.nb()] for the functions used to fit the models; [apm_est()] to compute the ATT and its uncertainty; [plot.apm_pre_fits()] for plotting an `apm_pre_fits` object.
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- apm_mod(crude_rate ~ 1,
#'                   family = list("gaussian", "quasipoisson"),
#'                   time_trend = 0:1,
#'                   lag = 0:1)
#' models
#' 
#' # Fit the models to data
#' fits <- apm_pre(models, data = ptpdata,
#'                 group_var = "group",
#'                 time_var = "year",
#'                 val_times = 1999:2007,
#'                 unit_var = "state",
#'                 nsim = 200,
#'                 verbose = FALSE)
#' 
#' fits
#' 
#' summary(fits)
#' 
#' plot(fits, type = "weights")
#' 
#' plot(fits, type = "error")

#' @export 
apm_pre <- function(models, data, weights = NULL, group_var, time_var,
                    val_times, unit_var, nsim = 1000L, cl = NULL,
                    verbose = TRUE) {
  
  # Argument checks
  chk::chk_not_missing(models, "`models`")
  chk::chk_is(models, "apm_models")
  
  chk::chk_not_missing(data, "`data`")
  chk::chk_data(data)
  
  # Process and order dataset
  data <- as.data.frame(data)
  
  if (is_null(rownames(data))) {
    rownames(data) <- seq_len(nrow(data))
  }
  
  chk::chk_not_missing(group_var, "`group_var`")
  chk::chk_string(group_var)
  chk::chk_subset(group_var, names(data))
  if (length(unique(data[[group_var]])) != 2) {
    chk::err("the grouping variable must have exactly 2 unique values")
  }
  data[[group_var]] <- factor(data[[group_var]])
  levels(data[[group_var]]) <- group_levels <- c("0", "1")
  
  chk::chk_not_missing(time_var, "`time_var`")
  chk::chk_string(time_var)
  chk::chk_subset(time_var, names(data))
  
  chk::chk_not_missing(unit_var, "`unit_var`")
  chk::chk_string(unit_var)
  chk::chk_subset(unit_var, names(data))
  data[[unit_var]] <- factor(data[[unit_var]])
  
  #Order data by unit_var and time_var
  data <- data[order(data[[unit_var]], data[[time_var]]), , drop = FALSE]
  
  chk::chk_not_missing(val_times, "`val_times`")
  chk::chk_numeric(val_times)
  chk::chk_subset(val_times, data[[time_var]])
  
  if (is_null(weights)) {
    weights <- rep.int(1, nrow(data))
  }
  else {
    chk::chk_numeric(weights)
    if (length(weights) != nrow(data)) {
      chk::err(sprintf("`weights` must have length equal to the number of rows in `data` (=%s)", nrow(data)))
    }
  }
  
  chk::chk_count(nsim)
  chk::chk_gte(nsim, 10)
  
  chk::chk_flag(verbose)
  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))
  
  # Create grid of models to be fit
  grid <- as.data.frame(matrix(NA_integer_,
                               nrow = length(models) * length(val_times),
                               ncol = 2L,
                               dimnames = list(NULL, c("time_ind", "model"))))
  
  #Remove models with problematic lags
  first_time <- min(data[[time_var]])
  lags <- pmax(vapply(models, `[[`, numeric(1L), "lag"),
               vapply(models, `[[`, numeric(1L), "diff_k"))
  for (v in val_times) {
    lag_too_much <- v - lags <= first_time
    
    if (any(lag_too_much)) {
      chk::err("some models involve lags corresponding to a period prior to the earliest time. Decrease the `lag` or `diff_k` components of the supplied models or use later validation times")
    }
  }
  
  #Fit all estimates
  val_data <- val_weights <- val_fits <- val_coefs <- val_predict_prep <- val_groups <- vector("list", nrow(grid))
  
  #Get observed means at each time point
  times <- sort(unique(data[[time_var]]))
  times <- times[times <= max(val_times)]
  y <- .get_y(models, data)
  
  observed_val_means <- setNames(lapply(times, function(t) {
    setNames(
      vapply(group_levels, function(g) {
        mean(y[data[[time_var]] == t & data[[group_var]] == g])
      }, numeric(1L)),
      group_levels
    )
  }), times)
  
  pred_errors_array <- array(NA_real_,
                             dim = c(length(val_times), length(models), length(group_levels)),
                             dimnames = list(val_times, names(models), group_levels))
  
  if (verbose) {
    cat("Fitting models...")
  }
  
  f <- 1L
  for (i in seq_along(models)) {
    model <- models[[i]]
    
    mod <- .modify_formula_and_data(model, data, group_var,
                                    unit_var, time_var)
    
    d <- mod$data
    
    for (t in seq_along(val_times)) {
      
      val_time <- val_times[t]
      
      fit <- .fit_one_model(mod, weights = weights,
                            time_var = time_var,
                            val_time = val_time,
                            family = model$family)
      
      subset_i <- which(d[[time_var]] == val_time)
      
      val_data[[f]] <- d[subset_i, , drop = FALSE]
      val_weights[[f]] <- weights[subset_i]
      val_coefs[[f]] <- na.omit(coef(fit))
      
      val_predict_prep[[f]] <- .make_predict_prep(fit, val_data[[f]])
      
      val_fits[[f]] <- fit
      
      val_groups[[f]] <- setNames(lapply(group_levels, function(g) {
        which(val_data[[f]][[group_var]] == g)
      }), group_levels)
      
      #Compute pred error
      
      # Compute prediction errors for each model for each validation period using original coefs
      
      ##Generate predictions on validation data
      # p <- predict(fit, newdata = val_data[[f]], type = "response")
      p <- .predict_quick(val_coefs[[f]],
                          val_predict_prep[[f]],
                          val_fits[[f]]$family$linkinv)
      
      #Unlog if outcome is logged to keep on original scale
      if (model$log) {
        p <- exp(p)
      }
      
      predicted_val_means_i <- setNames(
        vapply(group_levels, function(g) {
          .wtd_mean(p, val_weights[[f]], val_groups[[f]][[g]])
        }, numeric(1L)),
        group_levels
      )
      
      val_time_c <- as.character(val_time)
      
      for (g in group_levels) {
        pred_errors_array[t, i, g] <- observed_val_means[[val_time_c]][g] - predicted_val_means_i[g]
      }
      
      grid[["time_ind"]][f] <- t
      grid[["model"]][f] <- i
      f <- f + 1L
    }
  }
  
  #Difference in average prediction errors
  pred_error_diffs_mat <- pred_errors_array[, , 2L, drop = FALSE] - pred_errors_array[, , 1L, drop = FALSE]
  
  #Simulate to get BMA weights
  
  ## Joint variance of all model coefficients, clustering for unit
  val_vcov <- vcovSUEST(val_fits, cluster = data[[unit_var]]) 
  
  if (verbose) {
    cat(" Done.\nSimulating to compute BMA weights...\n")
  }
  
  #Draw parameters
  sim_coefs <- MASS::mvrnorm(nsim,
                             mu = unlist(val_coefs),
                             Sigma = val_vcov)
  
  coefs_inds <- .list_ind(val_coefs)
  
  # Compute prediction errors for each model for each validation period for each simulation
  
  #out_mat: all prediction error diffs; length(times) x length(models) x nsim
  pred_error_diffs_sim <- simplify2array(pbapply::pblapply(seq_len(nsim), function(s) {
    
    mat <- matrix(NA_real_,
                  nrow = length(val_times),
                  ncol = length(models),
                  dimnames = list(val_times, names(models)))
    
    coefs <- sim_coefs[s, ]
    
    for (f in seq_len(nrow(grid))) {
      mi <- grid$model[f]
      ti <- grid$time_ind[f]
      
      val_time_c <- as.character(val_times[ti])
      
      #Compute pred error
      
      ##Generate predictions on validation data
      # fit[["coefficients"]] <- coefs[coefs_inds[[f]]]
      # p <- predict(fit, newdata = val_data[[f]], type = "response")
      p <- .predict_quick(coefs[coefs_inds[[f]]],
                          val_predict_prep[[f]],
                          val_fits[[f]]$family$linkinv)
      
      #Unlog if outcome is logged to keep on original scale
      if (models[[mi]]$log) {
        p <- exp(p)
      }
      
      predicted_val_means_s_i <- vapply(group_levels, function(g) {
        .wtd_mean(p, val_weights[[f]], val_groups[[f]][[g]])
      }, numeric(1L))
      
      mat[ti, mi] <- (observed_val_means[[val_time_c]][2L] - observed_val_means[[val_time_c]][1L]) -
        (predicted_val_means_s_i[2L] - predicted_val_means_s_i[1L])
    }
    
    mat
  }, cl = cl))
  
  optimal_models <- vapply(seq_len(nsim), function(s) {
    worst_pred_within_model <- .colMax(abs(pred_error_diffs_sim[, , s, drop = FALSE]))
    
    which.min(worst_pred_within_model)
  }, integer(1L))
  
  if (verbose) {
    cat("Done.\n")
  }
  
  observed_means <- do.call("rbind", observed_val_means)
  rownames(observed_means) <- names(observed_val_means)
  
  BMA_weights <- tabulate(optimal_models, nbins = length(models)) / nsim
  
  fits <- list(models = models,
               val_times = val_times,
               grid = grid,
               val_fits = val_fits,
               val_coefs = val_coefs,
               val_vcov = val_vcov,
               data = data,
               weights = weights,
               observed_means = observed_means,
               pred_errors = pred_errors_array,
               pred_error_diffs = pred_error_diffs_mat,
               BMA_weights = BMA_weights,
               nsim = nsim)
  
  attr(fits, "time_var") <- time_var
  attr(fits, "unit_var") <- unit_var
  attr(fits, "group_var") <- group_var
  
  class(fits) <- "apm_pre_fits"
  
  fits
}

#' @exportS3Method print apm_pre_fits
print.apm_pre_fits <- function(x, ...) {
  cat("An `apm_pre_fits` object\n")
  cat("\n")
  cat(sprintf(" - grouping variable: %s\n", attr(x, "group_var")))
  cat(sprintf(" - unit variable: %s\n", attr(x, "unit_var")))
  cat(sprintf(" - time variable: %s\n", attr(x, "time_var")))
  cat(sprintf("   - validation times: %s\n", toString(x[["val_times"]])))
  cat(sprintf(" - number of models compared: %s\n", length(x[["models"]])))
  cat(sprintf(" - number of simulation iterations: %s\n", x[["nsim"]]))
  cat("\n")
  cat("Use `summary()` or `plot()` to examine prediction errors and BMA weights.\n")
}

#' @rdname apm_pre
#' @exportS3Method summary apm_pre_fits
summary.apm_pre_fits <- function(object, order = NULL, ...) {
  out <- data.frame(bma = object[["BMA_weights"]],
                    err = .colMax(abs(object[["pred_error_diffs"]])),
                    row.names = names(object[["models"]]))
  
  names(out) <- c("BMA weights", "Max|errors|")
  
  if (!is_null(order)) {
    chk::chk_string(order)
    order <- .match_arg(order, c("weights", "errors"))
    
    if (order == "weights") {
      out <- out[order(out[[1L]], decreasing = TRUE), , drop = FALSE]
    }
    else {
      out <- out[order(out[[2L]]), , drop = FALSE]
    }
  }
  
  class(out) <- c("summary.apm_pre_fits", class(out))
  
  out
}

#' @exportS3Method print summary.apm_pre_fits
print.summary.apm_pre_fits <- function(x, digits = 3, ...) {
  out <- matrix(NA_character_, nrow = nrow(x),
                ncol = ncol(x) + 1L,
                dimnames = list(rownames(x),
                                c(colnames(x), "")))
  
  for (i in seq_along(x)) {
    out[, i] <- format(round(x[[i]], digits), digits = digits)
  }
  
  out[, 3L] <- ""
  out[which.min(x[[2L]]), 3L] <- "*"
  
  print.default(out, right = TRUE, quote = FALSE, ...)
  
  cat("\n")
  cat("Use `plot()` to plot prediction errors and BMA weights.\n")
  
  invisible(x)
}
