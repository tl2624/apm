#' @title Plot outputs of `apm_pre()`
#' 
#' @description `plot()` displays the Bayesian model averaging (BMA) weights for each model (computed by `apm_fit()` as the posterior probability of selection) and the distribution of the difference in average prediction errors.
#' 
#' @param x an `apm_pre_fits` object; the output of a call to [apm_pre()].
#' @param type which values to plot: allowable options include `"weights"` to plot the BMA weights/posterior probabilities (default), `"errors"` to plot the difference in average predictions errors for all models across validation periods, `"predict"` to plot the time series and model predictions for each model, and `"corrected"` to plot the corrected predictions for the treated group for each model. Abbreviations allowed.
#' @param abs `logical`; when `type = "errors"`, whether to plot the differences in average prediction errors in absolute value (`TRUE`, default) or not (`FALSE`). 
#' @param ncol when `type` is `"errors"`, `"predict"`, or `"corrected"`, the number of columns to use to display the plots. Default is 4.
#' @param clip_at when `type = "errors"`, the value (in robust z-score units) at which to clip the y-axis of the plot to prevent outliers from distorting it. Default is 15. Set to `Inf` to prevent clipping.
#' @param model string; when `type = "predict"` or `type = "corrected"`, the model(s) to plot. Allowable values include `".optimal"` to plot the model with the smallest maximum absolute difference in average prediction errors, `".all"` to plot all models (excluding the BMA-weighted predictions), or the names of one or more specific models. Abbreviations allowed.
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
#' When `type = "predict"`, `plot()` displays a lattice of line plots, with a plot for each model displaying the observed and predicted outcomes for each validation period under each model. The observed outcomes are displayed as points, while the predicted outcomes are displayed as lines.
#' 
#' When `type = "corrected"`, `plot()` displays a lattice of line plots, with a plot for each model displaying the observed and corrected predictions for the treated group for each validation period under each model. The observed outcomes are displayed as points, while the corrected predictions are displayed as lines. Corrected predictions are computed as the observed outcome in the treated group minus the prediction error in the treated group plus the prediction error in the control group.
#' 
#' @seealso [apm_pre()] to to compute the difference in average prediction errors and BMA weights; `ggplot2::geom_col()`, which is used to create the plots.
#' 
#' @examples 
#' data("ptpdata")
#' 
#' # Combination of 8 models: 2 baseline formulas,
#' # 2 families, 2 lags
#' models <- apm_mod(crude_rate ~ 1,
#'                    family = "gaussian",
#'                    time_trend = 0:1,
#'                    lag = 0:1,
#'                    diff_k = 0:1)
#' models
#' 
#' # Fit the models to data
#' fits <- apm_pre(models, data = ptpdata,
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
#' 
#' plot(fits, type = "predict", model = ".optimal")
#' 
#' plot(fits, type = "corrected", model = ".optimal")

#' @exportS3Method plot apm_pre_fits
plot.apm_pre_fits <- function(x, type = "weights", abs = TRUE, ncol = 4L, clip_at = 15, model = ".optimal", ...) {
  chk::chk_string(type)
  type <- .match_arg(type, c("weights", "errors", "predict", "corrected"))
  
  if (length(x[["models"]]) != dim(x[["pred_errors"]])[2L] ||
      length(x[["models"]]) != dim(x[["pred_error_diffs"]])[2L] ||
      length(x[["models"]]) != length(x[["BMA_weights"]])) {
    chk::err("the `apm_pre_fit` object appears to be malformed")
  }
  
  dimnames(x[["pred_errors"]])[[2L]] <- names(x[["models"]])
  dimnames(x[["pred_error_diffs"]])[[2L]] <- names(x[["models"]])
  
  if (type == "weights") {
    
    d <- data.frame(
      models = factor(names(x[["models"]]), levels = names(x[["models"]])),
      BMA_weights = x[["BMA_weights"]]
    )
    
    p <- ggplot(d) +
      geom_hline(yintercept = 0) +
      geom_col(aes(x = .data[["models"]], y = .data[["BMA_weights"]])) +
      theme_bw() +
      labs(x = "Prediction Model", y = "Posterior Probability") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  else if (type == "errors") {
    chk::chk_flag(abs)
    chk::chk_count(ncol)
    chk::chk_gte(ncol, 1)
    chk::chk_number(clip_at)
    
    df <- x[["grid"]]
    
    if (nrow(df) != length(x[["pred_error_diffs"]])) {
      chk::err("the number of models implied to have been fit by the input object's `grid` component does not equal the the number of average prediction errors calculated, indicating a malformed `apm_pre_fit` object")
    }
    
    if (length(unique(x[["grid"]][["model"]])) != length(x[["models"]])) {
      chk::err("the number of model specifications listed in the input object's `grid` component does not equal the the number of model specifications present, indicating a malformed `apm_pre_fit` object")
    }
    
    df$pred_error_diffs <- as.vector(x[["pred_error_diffs"]])
    
    if (abs) {
      df$pred_error_diffs <- abs(df$pred_error_diffs)
    }
    
    df$time <- factor(df$time_ind, levels = seq_along(x$val_times),
                      labels = x[["val_times"]])
    
    df$model <- factor(df$model,
                       levels = seq_along(x$models),
                       labels = .firstup(names(x$models)))
    df$is_max <- factor(1, levels = 1:2, labels = c("no", "yes"))
    
    for (j in levels(df$model)) {
      df$is_max[df$model == j][which.max(df$pred_error_diffs[df$model == j])] <- "yes"
    }
    
    max_abs_pred_error <- .colMax(abs(x[["pred_error_diffs"]]))
    
    strip_cols <- rep.int("white", length(x$models))
    strip_cols[which.min(max_abs_pred_error)] <- "gray"
    
    strip <- strip_themed(background_x = elem_list_rect(fill = strip_cols))
    
    #Adjust plot to accommodate extreme outliers
    p <- ggplot(df) +
      geom_col(aes(x = as.numeric(.data$time),
                   y = .data$pred_error_diffs,
                   fill = .data$is_max),
               width = .96) +
      geom_hline(yintercept = 0) +
      facet_wrap2(vars(.data$model), ncol = ncol, strip = strip,
                  labeller = label_wrap_gen(30)) +
      scale_x_continuous(labels = levels(df$time),
                         breaks = seq_len(nlevels(df$time))) +
      scale_fill_manual(values = c("yes" = "black", "no" = "gray70")) +
      labs(y = sprintf("%sDifference in Average Prediction Errors",
                       if (abs) "Absolute " else ""),
           x = "Validation Period") +
      guides(fill = "none") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5))
    
    pe_std <- abs(df[["pred_error_diffs"]] - median(df[["pred_error_diffs"]])) / mad(df[["pred_error_diffs"]])
    if (any(pe_std > clip_at)) {
      ul <- max(df[["pred_error_diffs"]][pe_std <= clip_at])
      ylim <- c(0, ul)
      p <- p + scale_y_continuous(expand = expansion()) +
        coord_cartesian(ylim = ylim)
    }
  }
  else if (type == "predict") {
    chk::chk_character(model)
    
    chk::chk_count(ncol)
    chk::chk_gte(ncol, 1)
    
    model <- .match_arg(model, c(".optimal", ".all", names(x[["models"]])),
                        several.ok = TRUE)
    
    if (any(model == ".all")) {
      m_a <- which(model == ".all")[1L]
      model <- append(model, names(x[["models"]]), m_a - 1L)
      model <- model[model != ".all"]
    }
    
    if (any(model == ".optimal")) {
      model[model == ".optimal"] <- names(x[["models"]])[which.min(.colMax(abs(x[["pred_error_diffs"]])))]
    }
    
    model <- unique(model)
    
    pred_errors <- do.call("rbind", lapply(model, function(m) {
      pe <- x[["pred_errors"]][, m, ]
      
      merge(data.frame(time = rownames(x[["observed_means"]]),
                       model = m),
            as.data.frame(pe), by.x = "time", by.y = 0, all.x = TRUE)
    }))
    
    df <- reshape(as.data.frame(x[["observed_means"]]),
                  direction = "long", times = c("0", "1"),
                  varying = c("0", "1"), timevar = "group",
                  v.names = "observed", idvar = "time",
                  ids = rownames(x[["observed_means"]]))
    
    df_pred <- reshape(pred_errors,
                       direction = "long", times = c("0", "1"),
                       varying = c("0", "1"), timevar = "group",
                       v.names = "pred_error", idvar = c("time", "model"))
    
    df <- merge(df_pred, df,
                by = c("group", "time"),
                all = TRUE)
    
    df$pred <- df$observed - df$pred_error
    
    df$group <- factor(df$group, levels = c("0", "1"),
                       labels = c("Control", "Treated"))
    
    df$time <- as.factor(df$time)
    
    df$model <- factor(df$model, levels = model)
    
    p <- ggplot(df, aes(x = as.numeric(.data$time), color = .data$group)) +
      geom_point(aes(y = .data$observed)) +
      geom_line(aes(y = .data$pred, group = .data$group), na.rm = TRUE) +
      scale_x_continuous(labels = levels(df$time),
                         breaks = seq_len(nlevels(df$time))) +
      labs(x = "Period", y = "Outcome",
           title = "Treated and control groups' average predictions") +
      guides(color = guide_legend(title = element_blank(), position = "bottom")) +
      facet_wrap(vars(.data$model), ncol = ncol) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5))
  }
  else if (type == "corrected") {
    chk::chk_character(model)
    
    chk::chk_count(ncol)
    chk::chk_gte(ncol, 1)
    
    model <- .match_arg(model, c(".optimal", ".all", names(x[["models"]])),
                        several.ok = TRUE)
    
    if (any(model == ".all")) {
      m_a <- which(model == ".all")[1L]
      model <- append(model, names(x[["models"]]), m_a - 1L)
      model <- model[model != ".all"]
    }
    
    if (any(model == ".optimal")) {
      model[model == ".optimal"] <- names(x[["models"]])[which.min(.colMax(abs(x[["pred_error_diffs"]])))]
    }
    
    model <- unique(model)
    
    pred_errors <- do.call("rbind", lapply(model, function(m) {
      pe <- x[["pred_errors"]][, m, ]
      
      merge(data.frame(time = rownames(x[["observed_means"]]),
                       model = m),
            as.data.frame(pe), by.x = "time", by.y = 0, all.x = TRUE)
    }))
    
    names(pred_errors) <- c("time", "model", "pred_error.0", "pred_error.1")
    
    observed_means <- cbind(time = rownames(x[["observed_means"]]),
                            as.data.frame(x[["observed_means"]]))
    
    names(observed_means) <- c("time", "observed.0", "observed.1")
    
    df <- merge(pred_errors, observed_means, by = "time",
                all = TRUE)
    
    df$corrected.1 <- df[["observed.1"]] - df[["pred_error.1"]] + df[["pred_error.0"]]
    
    df$time <- as.factor(df$time)
    
    
    df$model <- factor(df$model, levels = model)
    
    p <- ggplot(df, aes(x = as.numeric(.data$time))) +
      geom_point(aes(y = .data$observed.1)) +
      geom_line(aes(y = .data$corrected.1), na.rm = TRUE) +
      scale_x_continuous(labels = levels(df$time),
                         breaks = seq_len(nlevels(df$time))) +
      labs(x = "Period", y = "Outcome",
           title = "Treated group's corrected predictions") +
      facet_wrap(vars(.data$model), ncol = ncol) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = .5))
  }
  
  p
}
