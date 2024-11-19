#' @title Plot outputs of `eepd_pre()`
#' 
#' @description `plot()` displays the Bayesian model averaging (BMA) weights for each model (computed by `eepd_fit()` as the posterior probability of selection) and the distribution of the difference in average prediction errors.
#' 
#' @param x an `eepd_pre_fits` object; the output of a call to [eepd_pre()].
#' @param type which values to plot: allowable options include `"weights"` to plot the BMA weights/posterior probabilities (default) and `"errors"` to plot the difference in average predictions errors for all models across validation periods. Abbreviations allowed.
#' @param abs `logical`; when `type = "errors"`, whether to plot the differences in average prediction errors in absolute value (`TRUE`, default) or not (`FALSE`). 
#' @param ncol when `type = "errors"`, the number of columns to use to display the plots. Default is 4.
#' @param clip_at when `type = "errors"`, the value (in robust z-score units) at which to clip the y-axis of the plot to prevent outliers from distorting it. Default is 15. Set to `Inf` to prevent clipping.
#' @param \dots ignored.
#' 
#' @returns
#' A `ggplot` object, which can be manipulated using `ggplot2` syntax (after loading `ggplot2`).
#' 
#' @details
#' When `type = "weights"`, `plot()` displays a bar plot with a bar for each model with height equal to the BMA weight/posterior probability of selection for that model. (Note that the plot margins can sometimes cut off the models names; use `theme(plot.margins =)` after loading `ggplot2` to extend the left margin of the plot to ensure all text is visible. Alternatively, the axis text can be rotated using `theme(axis.text.x =)`.)
#' 
#' When `type = "errors"`, `plot()` displays a lattice of bar plots, with a plot for each model displaying the difference in average prediction errors for each validation period. The period with the largest difference in average prediction errors will be shaded black. The model with the smallest maximum absolute difference in average prediction errors will have a gray label.
#' 
#' @seealso [eepd_pre()] to to compute the difference in average prediction errors and BMA weights; `ggplot2::geom_col()`, which is used to create the plots.
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- eepd_mod(crude_rate ~ 1,
#'                    family = "gaussian",
#'                    time_trend = 0:1,
#'                    lag = 0:1,
#'                    diff_k = 0:1)
#' models
#' 
#' # Fit the models to data
#' fits <- eepd_pre(models, data = ptpdata,
#'                  group_var = "group",
#'                  time_var = "year",
#'                  val_times = 1999:2007,
#'                  unit_var = "state",
#'                  nsim = 50)
#' fits
#' 
#' plot(fits, type = "weights")
#' 
#' plot(fits, type = "error", ncol = 2)

#' @exportS3Method plot eepd_pre_fits
plot.eepd_pre_fits <- function(x, type = "weights", abs = TRUE, ncol = 4, clip_at = 15, ...) {
  chk::chk_string(type)
  type <- match.arg(type, c("weights", "errors"))
  
  if (type == "weights") {
    
    models <- factor(names(x[["models"]]), levels = names(x[["models"]]))
    
    p <- ggplot() +
      geom_hline(yintercept = 0) +
      geom_col(aes(x = models, y = x[["BMA_weights"]])) +
      theme_bw() +
      labs(x = "Prediction Model", y = "Posterior Probability") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  else {
    chk::chk_flag(abs)
    chk::chk_count(ncol)
    chk::chk_gte(ncol, 1)
    chk::chk_number(clip_at)
    
    df <- x$grid
    
    if (nrow(df) != length(x[["pred_errors"]])) {
      chk::err("the number of models implied to have been fit by the input object's `grid` component does not equal the the number of average prediction errors calculated, indicating a malformed `eepd_pre_fit` object")
    }
    
    if (length(unique(x$grid$model)) != length(x$models)) {
      chk::err("the number of model specifications listed in the input object's `grid` component does not equal the the number of model specifications present, indicating a malformed `eepd_pre_fit` object")
    }
    
    df$pred_error <- as.vector(x[["pred_errors"]])
    
    if (abs) {
      df$pred_error <- abs(df$pred_error)
    }
    
    df$time <- factor(df$time_ind, levels = seq_along(x$val_times),
                      labels = x$val_times)
    df$model <- factor(df$model,
                       levels = seq_along(x$models),
                       labels = .firstup(names(x$models)))
    df$is_max <- factor(1, levels = 1:2, labels = c("no", "yes"))
    
    for (j in levels(df$model)) {
      df$is_max[df$model == j][which.max(df$pred_error[df$model == j])] <- "yes"
    }
    
    max_abs_pred_error <- apply(abs(x[["pred_errors"]]), 2, max)
    
    strip_cols <- rep(x = "white", times = length(models))
    strip_cols[which.min(max_abs_pred_error)] <- "gray"
    
    strip <- strip_themed(background_x = elem_list_rect(fill = strip_cols))
    
    #Adjust plot to accommodate extreme outliers
    p <- ggplot(df) +
      geom_col(aes(x = .data$time, y = .data$pred_error, fill = .data$is_max),
               width = .96) +
      geom_hline(yintercept = 0) +
      facet_wrap2(vars(.data$model), ncol = ncol, strip = strip,
                  labeller = label_wrap_gen(30)) +
      scale_fill_manual(values = c("yes" = "black", "no" = "gray70")) +
      labs(y = sprintf("%sDifference in Average Prediction Errors", if (abs) "Absolute " else ""),
           x = "Validation Period") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5))
    
    pe_std <- abs(df$pred_error - median(df$pred_err)) / mad(df$pred_error)
    if (any(pe_std > clip_at)) {
      ul <- max(df$pred_error[pe_std <= clip_at])
      ylim <- c(0, ul)
      p <- p + scale_y_continuous(expand = expansion()) +
        coord_cartesian(ylim = ylim)
    }
  }
  
  p
}
