#' Compute the robustness changepoint
#' 
#' @description `robustness_bound()` computes the value of the sensitivity parameter M at which the robustness bounds change from excluding to including an ATT of 0.
#' 
#' @param object an `apm_est` object; the output of a call to [apm_est()]. `M` must have been set to a nonzero value to use `robustness_bound()`.
#' @param level the desired confidence level. Set to 0 to ignore sampling variation in computing the interval bounds. Default is .95.
#' 
#' @returns
#' A single number corresponding to the changepoint value of M. If there is no positive value of M for which the interval bounds cross 0, `NA` will be returned.
#' 
#' @seealso [summary.apm_est()] for examining the ATT and bounds for a given value of `M`; [uniroot()] for the function that finds the changepoint value of `M`.
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
#' #Changepoint value of M ignoring estimation uncertainty
#' (M <- robustness_bound(est, level = 0))
#' 
#' summary(est, level = 0, M = M)
#' 
#' #Changepoint value of M accounting for estimation uncertainty
#' (M <- robustness_bound(est, level = .95))
#' 
#' summary(est, level = .95, M = M)

#' @export
robustness_bound <- function(object, level = .95) {
  chk::chk_is(object, "apm_est")
  
  chk::chk_number(level)
  chk::chk_gte(level, 0)
  chk::chk_lt(level, 1)
  
  if (object[["M"]] == 0) {
    chk::err("`robustness_bound()` cannot be used when `M` was 0 in the call to `apm_est()`")
  }
  
  att_inds <- seq_along(object[["BMA_weights"]])
  
  atts <- unname(object[["boot_out"]][["t0"]][att_inds])
  atts_boot <- object[["boot_out"]][["t"]][,att_inds, drop = FALSE]
  
  BMA_att_boot <- atts_boot %*% object[["BMA_weights"]]
  
  .make_bounds <- function(m, b = "lower") {
    me <- m * unname(object[["boot_out"]][["t0"]][-att_inds])
    me_boot <- m * object[["boot_out"]][["t"]][,-att_inds, drop = FALSE]
    
    BMA_me <- sum(object[["BMA_weights"]] * me)
    BMA_me_boot <- me_boot %*% object[["BMA_weights"]]
    
    if (b == "lower") {
      # ATT LB
      BMA_att_lb <- object[["BMA_att"]]["ATT"] - BMA_me
      
      BMA_lb_var_b <- var(BMA_att_boot - BMA_me_boot)
      BMA_lb_var_m <- sum(object[["BMA_weights"]] * ((atts - me) - BMA_att_lb)^2)
      
      BMA_lb_var <- BMA_lb_var_b + BMA_lb_var_m
      
      drop(BMA_att_lb + sqrt(BMA_lb_var) * qnorm((1 - level) / 2))
    }
    else {
      # ATT UB
      BMA_att_ub <- object[["BMA_att"]]["ATT"] + BMA_me
      
      BMA_ub_var_b <- var(BMA_att_boot + BMA_me_boot)
      BMA_ub_var_m <- sum(object[["BMA_weights"]] * ((atts + me) - BMA_att_ub)^2)
      
      BMA_ub_var <- BMA_ub_var_b + BMA_ub_var_m
      
      drop(BMA_att_ub + sqrt(BMA_ub_var) * qnorm(1 - (1 - level) / 2))
    }
  }
  
  ATT <- object[["BMA_att"]]["ATT"]
  
  ss <- summary(object, level = level, M = 0)
  
  if (ss["ATT", "Estimate"] > 0) {
    if (ss["ATT", "CI low"] < 0) {
      return(NA_real_)
    }
    
    b <- "lower"
  }
  else {
    if (ss["ATT", "CI high"] > 0) {
      return(NA_real_)
    }
    
    b <- "upper"
  }

  u <- uniroot(.make_bounds, lower = 0, upper = 10, b = b,
               extendInt = switch(b, "upper" = "upX", "downX"),
               tol = 1e-8)
  
  u$root
}
