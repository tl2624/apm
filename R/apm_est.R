#' Estimate ATTS from models fits
#' 
#' @description `apm_est()` computes the ATTs from the models previously fit by [apm_pre()], choosing the optimal one by minimizing the largest absolute average prediction error across validation times. Optionally, this process can be simulated to arrive at a distribution of ATTs that accounts for the uncertainty in selecting the optimal model. `plot()` plots the resulting ATT(s).
#' 
#' @inheritParams apm_pre
#' @param fits an `apm_pre_fits` object; the output of a call to [apm_pre()].
#' @param post_time the value of the time variable considered post-treatment, for which the ATT is to be estimated.
#' @param M the sensitivity parameter for set identification. For `apm_est()`, the default is 0, i.e., under point identification. For `summary()`, this can be set to one or more positive values to produce uncertainty bounds for each value. Only allowed when not set to 0 in the call to `apm_est()`. See Details.
#' @param R the number of bootstrap iterations used to compute the sampling variance of the ATT. Default is 1000. More is better but takes longer.
#' @param all_models `logical`; whether to compute ATTs for all models (`TRUE`) or just those with BMA weights greater than 0 (`FALSE`, default). This will not effect the final estimates but leaving as `FALSE` can speed up computation when some models have BMA weights of 0.
#' @param x,object an `apm_est` object; the output of a call to `apm_est()`.
#' @param level the desired confidence level. Set to 0 to ignore sampling variation in computing the interval bounds. Default is .95.
#' @param label `logical`; whether to label the ATT estimates. Requires \pkg{ggrepel} to be installed. Default is `TRUE`.
#' @param size.weights `logicsl`; whether to size the points based on their BMA weights. Default is `TRUE`.
#' 
#' @returns
#' `apm_est()` returns an `apm_est` object, which contains the ATT estimates and their variance estimates. The following components are included:
#' \describe{
#' \item{BMA_att}{the BMA-weighted ATT}
#' \item{atts}{a matrix containing the ATT estimates from each model (when `all_models = FALSE`, only models with positive BMA weights are included)}
#' \item{BMA_var}{the total variance estimate for the BMA-weighted ATT incorporating the variance due to sampling and due to model selection}
#' \item{BMA_var_b}{the bootstrap-based component of the variance estimate for the BMA-weighted ATT due to sampling}
#' \item{BMA_var_m}{the component of the variance estimate for the BMA-weighted ATT due to model selection}
#' \item{M}{the value of the sensitivity parameter `M`}
#' \item{post_time}{the value supplied to `post_time`}
#' \item{pred_errors}{a matrix containing the difference in average prediction errors for each model and each pre-treatment validation period}
#' \item{BMA_weights}{the BMA weights computed by `apm_pre()` (when `all_models = FALSE`, only positive BMA weights are included)}
#' \item{boot_out}{an `fwb` object containing the bootstrap results}
#' }
#' 
#' `plot()` returns a `ggplot` object displaying the ATT for each model plotted against the maximum absolute difference in average prediction errors for that model. The model with the lowest maximum absolute difference in average prediction errors is displayed in red.
#' 
#' `summary()` produces a table with the BMA-weighted ATT, it's estimated standard error, and confidence interval limits. When `M` is greater than 0, additional rows for each value of `M` are included with the lower and upper bound. When `level` is greater than 0, these bounds include the uncertainty due to sampling and model selection; otherwise, they correspond to the set identification bounds for the ATT.
#' 
#' @details
#' `apm_est()` estimates the ATT from each model and combines them to form the BMA-weighted estimate of the ATT. Uncertainty for the BMA-weighted ATT is computed by combining two variance components, one that account for sampling and another that accounts for model selection. The component due to sampling is computed by bootstrapping the process of fitting the outcome model for the post-treatment outcome identified by `post_time` and computing the difference between the observed outcome mean difference and the model-predicted outcome mean difference. The fractional weighted bootstrap as implemented in [fwb::fwb()] is used to ensure no units are dropped from the analysis. In each bootstrap sample, the BMA-weighted ATT estimate is computed as the weighted average of the ATTs computed from the models using the fixed BMA weights computed by [apm_pre()], and the variance is computed as the empirical variance over the bootstrapped estimates. The variance component due to model selection is computed as the BMA-weighted variance of the original ATTs.
#' 
#' When `M` is greater than 0, bounds for set identification and their uncertainty are additionally computed. This involves bootstrapping the fitting of the pre-period models along with post-treatment models on order to compute the maximum absolute difference in average prediction errors for each model across validation periods. Each bootstrap sample produces a margin of error for each model computed as \eqn{M \times \delta_m} where \eqn{\delta_m} is the maximum absolute difference in average prediction errors for model \eqn{m}. Upper and lower bounds for the set-identified BMA-weighted ATT are computed as \eqn{\text{ATT}_m \pm M \times \delta_m}. The same procedure as above is then used to compute the variance of these bounds.
#' 
#' `summary()` displays the BMA-weighted ATT estimate, its standard error, and Wald confidence intervals. When `M` is greater than 0, bounds for the set-identified ATT are displayed in the confidence interval bound columns. The lower bound is computed as \eqn{\text{LB} - \sigma_{LB}Z_{l}} and the upper bound as \eqn{\text{UB} + \sigma_{UB}Z_{l}}, where \eqn{\text{LB}} and \eqn{\text{UB}} are the lower and upper bounds, \eqn{\sigma_{LB}} and \eqn{\sigma_{UB}} are their variances accounting for sampling and model selection, and \eqn{Z_{l}} is the critical Z-statistic for confidence level \eqn{l}. To display the set-identification bounds themselves, one should set `level = 0`.
#' 
#' @seealso [apm_pre()] for computing the BMA weights; [fwb::fwb()] for the fractional weighted bootstrap.
#' 
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 4 models: 2 time trends, 2 lags
#' models <- apm_mod(list(crude_rate ~ 1),
#'                    lag = 0:1,
#'                    time_trend = 0:1)
#' models
#' 
#' # Fit the models to data; unit_var must be supplied for
#' # fixed effects
#' fits <- apm_pre(models,
#'                  data = ptpdata,
#'                  group_var = "group",
#'                  time_var = "year",
#'                  val_times = 2004:2007,
#'                  unit_var = "state",
#'                  nsim = 100)
#' 
#' est <- apm_est(fits,
#'                 post_time = 2008,
#'                 M = 1,
#'                 R = 20)
#' 
#' est
#' 
#' # ATT estimate and bounds for M = 1
#' summary(est)
#' 
#' # Bounds for other values of M
#' summary(est, M = c(.5, 1, 1.5, 2))
#' 
#' # Set-ID bounds without uncertainty
#' summary(est, level = 0)
#' 
#' plot(est)

#' @export
apm_est <- function(fits, post_time, M = 0, R = 1000L, all_models = FALSE, cl = NULL, verbose = TRUE, ...) {
  chk::chk_is(fits, "apm_pre_fits")
  
  time_var <- attr(fits, "time_var")
  data <- fits$data
  
  chk::chk_not_missing(post_time, "`post_time`")
  chk::chk_number(post_time)
  chk::chk_subset(post_time, data[[time_var]])
  chk::chk_gt(post_time, fits$val_times)
  
  unit_var <- attr(fits, "unit_var")
  
  group_var <- attr(fits, "group_var")
  group_levels <- levels(data[[group_var]])
  
  chk::chk_number(M)
  chk::chk_gte(M, 0)
  
  chk::chk_count(R)
  chk::chk_gte(R, 10)
  
  weights <- fits$weights
  
  chk::chk_flag(all_models)
  chk::chk_flag(verbose)
  
  val_times <- fits$val_times
  models <- fits$models
  
  #Remove models that won't contribute
  if (!all_models) {
    models <- models[fits$BMA_weights > 0]
    BMA_weights <- fits$BMA_weights[fits$BMA_weights > 0]
  }
  else {
    BMA_weights <- fits$BMA_weights
  }
  
  mods <- lapply(models, function(model) {
    .modify_formula_and_data(model, data, group_var,
                             unit_var, time_var)
  })
  
  #FWB for ATTs
  .boot_fun <- function(.data, .weights, ...) {
    times <- post_time
    
    if (M > 0) {
      times <- c(times, val_times)
    }
    
    pred_error_mat <- matrix(NA_real_,
                             nrow = length(times),
                             ncol = length(models),
                             dimnames = list(times,
                                             names(models)))
    
    for (i in seq_along(models)) {
      model <- models[[i]]
      
      d <- mods[[i]]$data
      
      for (t in seq_along(times)) {
        
        time <- times[t]
        
        subset_i_v <- which(d[[time_var]] == time)
        
        .val_data_i_v <- d[subset_i_v,, drop = FALSE]
        .val_weights_i_v <- .weights[subset_i_v] * weights[subset_i_v]
        .val_groups_i_v <- setNames(lapply(group_levels, function(g) {
          which(.val_data_i_v[[group_var]] == g)
        }), group_levels)
        
        y <- model.response(model.frame(update(mods[[i]]$formula, . ~ 1),
                                        data = .val_data_i_v))
        
        if (model$log) {
          y <- exp(y)
        }
        
        observed_val_means_i <- setNames(
          vapply(group_levels, function(g) {
            .wtd_mean(y, .val_weights_i_v, .val_groups_i_v[[g]])
          }, numeric(1L)),
          group_levels
        )
        
        # Compute prediction errors for each model for each validation period
        fit <- .fit_one_model(mods[[i]],
                              weights = weights * .weights,
                              time_var = time_var,
                              val_time = time,
                              family = model$family)
        
        .val_predict_prep_i_v <- .make_predict_prep(fit, .val_data_i_v)
        
        ##Generate predictions on validation data
        # p <- predict(fit, newdata = .val_data_i_v, type = "response")
        p <- .predict_quick(na.omit(coef(fit)),
                            .val_predict_prep_i_v,
                            fit$family$linkinv)
        
        #Unlog if outcome is logged to keep on original scale
        if (model$log) {
          p <- exp(p)
        }
        
        predicted_val_means_i <- setNames(
          vapply(group_levels, function(g) {
            .wtd_mean(p, .val_weights_i_v, .val_groups_i_v[[g]])
          }, numeric(1L)),
          group_levels
        )
        
        pred_error_mat[t, i] <- (observed_val_means_i["1"] - observed_val_means_i["0"]) -
          (predicted_val_means_i["1"] - predicted_val_means_i["0"])
      }
    }
    
    #ATT is pred error for post_time
    .atts <- pred_error_mat[1,]
    
    if (M == 0) {
      return(.atts)
    }
    
    .max_abs_pred_errors <- apply(abs(pred_error_mat[-1, , drop = FALSE]), 2L, max)
    
    c(.atts, .max_abs_pred_errors)
  }
  
  boot_out <- fwb::fwb(data = data,
                       statistic = .boot_fun,
                       R = R,
                       cluster = data[[unit_var]],
                       verbose = verbose,
                       cl = cl, ...)
  
  att_inds <- seq_along(models)
  
  # ATT
  atts <- unname(boot_out[["t0"]][att_inds])
  atts_boot <- boot_out[["t"]][,att_inds, drop = FALSE]
  
  BMA_att <- sum(BMA_weights * atts)
  BMA_att_boot <- atts_boot %*% BMA_weights
  
  BMA_var_b <- var(BMA_att_boot)
  BMA_var_m <- sum(BMA_weights * (atts - BMA_att)^2)
  
  BMA_var <- BMA_var_b + BMA_var_m
  
  out <- list(BMA_att = c(ATT = unname(BMA_att)),
              atts = matrix(atts, ncol = 1L,
                            dimnames = list(names(models), "ATT")),
              BMA_var = c(ATT = unname(BMA_var)),
              BMA_var_b = c(ATT = unname(BMA_var_b)),
              BMA_var_m = c(ATT = unname(BMA_var_m)),
              M = M,
              post_time = post_time,
              observed_means = fits$observed_means,
              pred_errors = fits$pred_errors,
              pred_errors_diff = fits$pred_errors_diff,
              BMA_weights = BMA_weights,
              boot_out = boot_out)
  
  if (M > 0) {
    me <- M * unname(boot_out[["t0"]][-att_inds])
    me_boot <- M * boot_out[["t"]][,-att_inds, drop = FALSE]
    
    BMA_me <- sum(BMA_weights * me)
    BMA_me_boot <- me_boot %*% BMA_weights
    
    # ATT LB
    BMA_att_lb <- BMA_att - BMA_me
    
    BMA_lb_var_b <- var(BMA_att_boot - BMA_me_boot)
    BMA_lb_var_m <- sum(BMA_weights * ((atts - me) - BMA_att_lb)^2)
    
    BMA_lb_var <- BMA_lb_var_b + BMA_lb_var_m
    
    # ATT UB
    BMA_att_ub <- BMA_att + BMA_me
    
    BMA_ub_var_b <- var(BMA_att_boot + BMA_me_boot)
    BMA_ub_var_m <- sum(BMA_weights * ((atts + me) - BMA_att_ub)^2)
    
    BMA_ub_var <- BMA_ub_var_b + BMA_ub_var_m
    
    out[["BMA_att"]] <- c(out[["BMA_att"]],
                          LB = unname(BMA_att_lb),
                          UB = unname(BMA_att_ub))
    
    out[["atts"]] <- cbind(out[["atts"]],
                           LB = unname(atts - me),
                           UB = unname(atts + me))
    
    out[["BMA_var"]] <- c(out[["BMA_var"]],
                          LB = unname(BMA_lb_var),
                          UB = unname(BMA_ub_var))
    
    out[["BMA_var_b"]] <- c(out[["BMA_var_b"]],
                            LB = unname(BMA_lb_var_b),
                            UB = unname(BMA_ub_var_b))
    
    out[["BMA_var_m"]] <- c(out[["BMA_var_m"]],
                            LB = unname(BMA_lb_var_m),
                            UB = unname(BMA_ub_var_m))
  }
  
  attr(out, "time_var") <- time_var
  attr(out, "unit_var") <- unit_var
  attr(out, "group_var") <- group_var
  
  class(out) <- "apm_est"
  
  out
}

#' @exportS3Method print apm_est
print.apm_est <- function(x, ...) {
  cat("An `apm_est` object\n\n")
  cat(sprintf(" - grouping variable: %s\n", attr(x, "group_var")))
  cat(sprintf(" - unit variable: %s\n", attr(x, "unit_var")))
  cat(sprintf(" - time variable: %s\n", attr(x, "time_var")))
  cat(sprintf("   - validation times: %s\n", paste(x[["val_times"]], collapse = ", ")))
  cat(sprintf("   - post-treatment time: %s\n", paste(x[["post_time"]], collapse = ", ")))
  cat(sprintf(" - sensitivity parameter (M): %s\n", x[["M"]]))
  cat(sprintf(" - bootstrap replications: %s\n", nrow(x[["boot_out"]][["t"]])))
  cat("\n")
  cat("Use `summary()` or `plot()` to examine estimates and uncertainty bounds.\n")
  
  invisible(x)
}

#' @exportS3Method summary apm_est
#' @rdname apm_est
summary.apm_est <- function(object, level = .95, M = NULL, ...) {
  
  res <- data.frame(
    object[["BMA_att"]]["ATT"],
    sqrt(object[["BMA_var"]]["ATT"]))
  
  res <- cbind(res,
               res[[1L]] + res[[2L]] * qnorm((1 - level) / 2),
               res[[1L]] + res[[2L]] * qnorm(1 - (1 - level) / 2),
               res[[1L]] / res[[2L]],
               2 * pnorm(-abs(res[[1L]] / res[[2L]])))
  
  names(res) <- c("Estimate", "Std. Error", "CI low", "CI high", "z_value", "Pr(>|z|)")
  
  rownames(res) <- "ATT"
  
  if (is.null(M)) {
    M <- object[["M"]]
  }
  else {
    chk::chk_numeric(M)
    chk::chk_gte(M, 0)
    
    M <- sort(unique(M))
    M <- M[M > 0]
  }
  
  if (length(M) > 0) {
    if (object[["M"]] == 0) {
      chk::err("`M` cannot be nonzero when `M` was 0 in the call to `apm_est()`")
    }
    
    res2 <- do.call("rbind", lapply(M, function(m) {
      res_m <- as.data.frame(matrix(NA_real_, nrow = 1, ncol = ncol(res)))
      
      if (m == object[["M"]]) {
        BMA_att_lb <- object[["BMA_att"]]["LB"]
        BMA_att_ub <- object[["BMA_att"]]["UB"]
        
        BMA_lb_var <- object[["BMA_var"]]["LB"]
        BMA_ub_var <- object[["BMA_var"]]["UB"]
      }
      else {
        att_inds <- seq_along(object[["BMA_weights"]])
        
        atts <- unname(object[["boot_out"]][["t0"]][att_inds])
        atts_boot <- object[["boot_out"]][["t"]][,att_inds, drop = FALSE]
        
        BMA_att_boot <- atts_boot %*% object[["BMA_weights"]]
        
        me <- m * unname(object[["boot_out"]][["t0"]][-att_inds])
        me_boot <- m * object[["boot_out"]][["t"]][,-att_inds, drop = FALSE]
        
        BMA_me <- sum(object[["BMA_weights"]] * me)
        BMA_me_boot <- me_boot %*% object[["BMA_weights"]]
        
        # ATT LB
        BMA_att_lb <- object[["BMA_att"]]["ATT"] - BMA_me
        
        BMA_lb_var_b <- var(BMA_att_boot - BMA_me_boot)
        BMA_lb_var_m <- sum(object[["BMA_weights"]] * ((atts - me) - BMA_att_lb)^2)
        
        BMA_lb_var <- BMA_lb_var_b + BMA_lb_var_m
        
        # ATT UB
        BMA_att_ub <- object[["BMA_att"]]["ATT"] + BMA_me
        
        BMA_ub_var_b <- var(BMA_att_boot + BMA_me_boot)
        BMA_ub_var_m <- sum(object[["BMA_weights"]] * ((atts + me) - BMA_att_ub)^2)
        
        BMA_ub_var <- BMA_ub_var_b + BMA_ub_var_m
      }
      
      res_m[[3L]][1L] <- BMA_att_lb + sqrt(BMA_lb_var) * qnorm((1 - level) / 2)
      res_m[[4L]][1L] <- BMA_att_ub + sqrt(BMA_ub_var) * qnorm(1 - (1 - level) / 2)
      
      names(res_m) <- names(res)
      
      rownames(res_m) <- sprintf("M = %s", round(m, 2))
      
      res_m
    }))
    
    res <- rbind(res, res2)
  }
  
  class(res) <- c("summary.apm_est", class(res))
  
  res
}

#' @exportS3Method print summary.apm_est
print.summary.apm_est <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  printCoefmat(x, digits = digits, cs.ind = 1:4,
               tst.ind = 5,
               P.values = TRUE,
               has.Pvalue = TRUE,
               na.print = ".",
               zap.ind = 3:4)
  
  invisible(x)
}

#' @exportS3Method plot apm_est
#' @rdname apm_est
plot.apm_est <- function(x, label = TRUE, size.weights = TRUE, ...) {
  chk::chk_flag(label)

  max_abs_pred_error <- apply(abs(x[["pred_errors"]]), 2L, max)
  est <- x[["atts"]][,1L]
  
  labels <- {
    if (label) rownames(x[["atts"]])
    else NULL
  }
  
  plot_data <- data.frame(estimate = est,
                          pred_error = max_abs_pred_error[names(est)],
                          label = labels,
                          weights = x[["BMA_weights"]],
                          best = seq_along(est) == which.min(max_abs_pred_error[names(est)]))
  
  p <- ggplot(plot_data,
              aes(x = .data$pred_error,
                  y = .data$estimate)) +
    geom_point(aes(color = .data$best,
                   size = if (size.weights) .data$weights),
               alpha = .8) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    guides(color = "none", size = "none") +
    labs(x = "Maximum absolute pre-treatment difference in average prediction errors",
         y = "Estimate") +
    theme_bw()
  
  if (label) {
    p <- p + ggrepel::geom_text_repel(aes(label = .data$label),
                                      box.padding = 1,
                                      point.padding = 1,
                                      min.segment.length = .2)
  }
  
  p
}
