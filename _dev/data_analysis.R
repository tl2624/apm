# Setup -------------------------------------------------------------------
rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(scipen = 999)
#library(devtools)
#install_github("ngreifer/eepd")
library(parallel)
n_cores <- detectCores() - 1
library(tidyverse)
library(pbapply)
library(data.table)
library(eepd)

# Data analysis -----------------------------------------------------------

data("ptpdata")

years <- sort(unique(ptpdata$year))
treated_year <- 2008
val_years <- 1999:2007

## Same models as in Figure 2 of manuscript
models <- eepd_mod(list(crude_rate ~ 1,
                        crude_rate ~ 1 + year,
                        crude_rate ~ 1 + poly(year, 2)),
                   log = c(TRUE, FALSE),
                   lag = c(0, 1),
                   diff_k = c(0, 1),
                   fixef = TRUE)

model_names <- c("Baseline mean (log)", "Lin trend (log)", "Quad trend (log)", "LDV (log)", "Lin trend + LDV (log)",
                 "Quad trend + LDV (log)", "Baseline mean (log, first diff)", "Lin trend (log, first diff)",
                 "Quad trend (log, first diff)", "Baseline mean", "Lin trend", "Quad trend", "LDV",
                 "Lin trend + LDV", "Quad trend + LDV", "Baseline mean (first diff)", "Lin trend (first diff)",
                 "Quad trend (first diff)")

fits <- eepd_fit(models = models,
                 data = ptpdata,
                 group_var = "group",
                 time_var = "year",
                 val_times = val_years,
                 post_time = treated_year,
                 unit_var = "state")
plot(eepd_sim(fits = fits))

mod_point_ests <- sapply(X = 1:length(models),
                         FUN = function(x) { eepd_sim(fits = eepd_fit(models = eepd_mod(formula_list = models[[x]]$formula,
                                                                                        lag = models[[x]]$lag,
                                                                                        diff_k = models[[x]]$diff_k,
                                                                                        log = models[[x]]$log,
                                                                                        fixef = models[[x]]$fixef),
                                                                      data = ptpdata,
                                                                      group_var = "group",
                                                                      time_var = "year",
                                                                      val_times = val_years,
                                                                      post_time = treated_year,
                                                                      unit_var = "state"))$atts })

mod_max_abs_errors <- pblapply(X = 1:length(models),
                               FUN = function(x) { abs(eepd_sim(fits = eepd_fit(models = eepd_mod(formula_list = models[[x]]$formula,
                                                                                                  lag = models[[x]]$lag,
                                                                                                  diff_k = models[[x]]$diff_k,
                                                                                                  log = models[[x]]$log,
                                                                                                  fixef = models[[x]]$fixef),
                                                                                data = ptpdata,
                                                                                group_var = "group",
                                                                                time_var = "year",
                                                                                val_times = val_years,
                                                                                post_time = treated_year,
                                                                                unit_var = "state"))$sim_pred_errors[, , 1]) },
                               cl = n_cores)

## smallest absolute error in last pre-treatment period
which.min(sapply(X = 1:length(mod_max_abs_errors),
                 FUN = function(x) { mod_max_abs_errors[[x]]["2007"] }))

## smallest average absolute error
which.min(sapply(X = 1:length(mod_max_abs_errors),
                 FUN = function(x) { mean(mod_max_abs_errors[[x]]) }))
          
max_abs_errors <- sapply(X = 1:length(mod_max_abs_errors),
                         FUN = function(x) { max(mod_max_abs_errors[[x]]) })

model_err_plot_data <- data.frame(Model = rep(x = model_names,
                                              each = length(val_years)),
                                  Abs_error = unlist(mod_max_abs_errors),
                                  Year = as.numeric(names(unlist(mod_max_abs_errors))))
model_err_plot_data$Model <- factor(x = model_err_plot_data$Model,
                                    levels = model_names)

model_err_plot_data <- group_by(.data = model_err_plot_data,
                                Model) %>%
  mutate(max = ifelse(test = Abs_error == max(Abs_error),
                      yes = 1,
                      no = 0),
         max = as.factor(max))

library(ggh4x)
strip_cols <- rep(x = "white", times = length(models))
strip_cols[which.min(max_abs_errors)] <- "gray"
strip <- strip_themed(background_x = elem_list_rect(fill = strip_cols))

mod_pred_errors_barplot <- ggplot(data = model_err_plot_data,
                                  mapping = aes(x = Year,
                                                y = Abs_error,
                                                fill = ifelse(test = max == 1,
                                                              yes = "Highlighted",
                                                              no = "Normal"))) +
  geom_bar(stat = "identity") +
  facet_wrap2(facets = . ~ model_err_plot_data$Model,
              nrow = 3,
              ncol = 6,
              scales = "fixed",
              strip = strip) +
  scale_fill_manual(values = c("Highlighted" = "black",
                               "Normal" = "gray")) +
  theme_bw() +
  scale_x_continuous(breaks = sort(unique(model_err_plot_data$Year))) +
  # scale_y_continuous(limits = c(0, 1.25)) +
  theme(axis.text.x = element_text(angle = 90, size = 5),
        axis.text.y = element_text(size = 5),
        strip.text.x = element_text(size = 4.25),
        strip.background = element_rect(colour="black",
                                        fill="white"),
        panel.grid.major = element_blank(),
        legend.position = "none") +
  labs(x = "Validation Year",
       y = "Absolute Difference in Average Prediction Errors") 
mod_pred_errors_barplot
#ggsave(plot = mod_pred_errors_barplot,
#       file = "mod_pred_errors_barplot.pdf",
#       width = 6,
#       height = 4,
#       units = "in",
#       dpi = 600)

# Point estimates by abs pred error plot ----------------------------------

M <- 1

round(x = mod_point_ests[which.min(max_abs_errors)] - min(max_abs_errors), digits = 2)
round(x = mod_point_ests[which.min(max_abs_errors)] + min(max_abs_errors), digits = 2)

round(x = sqrt(mean((mod_point_ests - mean(mod_point_ests))^2)), digits = 2)

round(x = mod_point_ests[which.min(max_abs_errors)], digits = 2)
round(x = mod_point_ests[which.max(max_abs_errors)], digits = 2)


# Plot of maximum errors and point estimates with bounds ------------------

errors_ests_data <- data.frame(est = c(mod_point_ests,
                                       mod_point_ests - (M * max_abs_errors),
                                       mod_point_ests + (M * max_abs_errors)),
                               est_type = c(rep(x = "point", times = length(mod_point_ests)),
                                            rep(x = "lb", times = length(mod_point_ests)),
                                            rep(x = "ub", times = length(mod_point_ests))),
                               max_err = rep(x = max_abs_errors, times = 3),
                               model = "")

errors_ests_data$model[which.min(errors_ests_data$max_err[errors_ests_data$est_type == "point"])] <- "LDV (log)" #models[which.min(max_abs_errors)]
errors_ests_data$model[errors_ests_data$max_err == sort(max_abs_errors)[length(max_abs_errors)] & errors_ests_data$est_type == "point"] <- "Quad trend (first diff)" #models[order(max_abs_errors)[length(models)]]
errors_ests_data$model[errors_ests_data$max_err == sort(max_abs_errors)[13] & errors_ests_data$est_type == "point"] <- "Lin trend (first diff)" #models[order(max_abs_errors)[13]]
errors_ests_data$model[errors_ests_data$max_err == sort(max_abs_errors)[9] & errors_ests_data$est_type == "point"] <- "Lin trend (log)" #models[order(max_abs_errors)[9]]
seg_data <- data.frame(x = errors_ests_data$max_err,
                       xend = errors_ests_data$max_err, 
                       y = errors_ests_data$est[errors_ests_data$est_type == "lb"],
                       yend = errors_ests_data$est[errors_ests_data$est_type == "ub"],
                       model = errors_ests_data$model)

library(ggrepel)
ests_max_errs_plot <- ggplot(data = errors_ests_data,
                             mapping = aes(x = max_err,
                                           y = est,
                                           label = model)) +
  #geom_segment(data = seg_data,
  #             mapping = aes(x = x,
  #                           y = y,
  #                           xend = xend,
  #                           yend = yend),
  #             color = "gray",
  #             alpha = 1) +
  geom_point(data = filter(.data = errors_ests_data,
                           est_type == "point"),
             color = "black",
             size = 0.75) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  scale_x_continuous(limits = c(0.5, 2.1)) +
  scale_y_continuous(limits = c(-2.5, 4)) +
  labs(x = "Maximum absolute pre-treatment difference in average prediction errors",
       y = "Estimate") +
  geom_text_repel(size = 2.75)
ests_max_errs_plot
#ggsave(plot = ests_max_errs_plot,
#       file = "ests_max_errs_plot.pdf",
#       width = 6,
#       height = 4,
#       units = "in",
#       dpi = 600)


# BMA estimation and inference --------------------------------------------
sims <- 1000
set.seed(09291993)
res <- eepd_sim(fits = fits,
                nsim = sims,
                cl = n_cores,
                verbose = TRUE)

round(x = mean(res$optimal_models == which.min(max_abs_errors)), digits = 2)
round(x = mean(res$optimal_models == order(max_abs_errors)[2]), digits = 2)
models[order(max_abs_errors)[2]]

post_probs <- sapply(X = 1:length(models),
                     FUN = function(x) { mean(res$optimal_models == x) } )

## Plot posterior distribution
post_data <- data.frame(prob = post_probs,
                        model = "")
models[which(post_data$prob > 0)]
post_data$model[which(post_data$prob > 0)] <- model_names[which(post_data$prob > 0)]

post_dist_plot <- ggplot(data = filter(.data = post_data, model != ""),
                         mapping = aes(x = model,
                                       y = prob)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Prediction model",
       y = "Posterior probability") +
  theme(axis.text.x = element_text(angle = 90))
post_dist_plot
#ggsave(plot = post_dist_plot,
#       file = "post_dist_plot.pdf",
#       width = 6,
#       height = 4,
#       units = "in",
#       dpi = 600)

BMA_est <- sum(mod_point_ests * post_probs)
round(x = BMA_est, digits = 2)
round(x = mod_point_ests[which.min(max_abs_errors)], digits = 2)
mod_var <- sum((mod_point_ests - BMA_est)^2 * post_probs)

Ms <- seq(from = 0, to = 2, by = 0.01)
BMA_lb_ests <- pbsapply(X = Ms,
                        FUN = function(x) { sum((mod_point_ests - (x * max_abs_errors)) * post_probs) },
                        cl = n_cores)

changepoint_M <- Ms[which(BMA_lb_ests <= 0)[1]]

# Bootstrap inference -----------------------------------------------------

treated_states <- unique(ptpdata$state[ptpdata$group == 1])
control_states <- unique(ptpdata$state[ptpdata$group == 0])

n_1 <- length(treated_states)
n_0 <- length(control_states)

n_boots <- 1000
exp_est <- rep(x = NA, times = n_boots)
exp_lb_est <- rep(x = NA, times = n_boots)
exp_ub_est <- rep(x = NA, times = n_boots)

M <- 1.07

set.seed(09291993)
for(i in 1:n_boots){
  
  ## sample with replacement from N_treat treated states in population
  samp_treated_states = sample(x = treated_states, size = n_1, replace = TRUE)
  
  ## sample with replacement from control states in population
  samp_control_states = sample(x = control_states, size = n_0, replace = TRUE)
  
  samp_treat_data = do.call(what = "rbind",
                            args = pblapply(X = 1:n_1,
                                            FUN = function(x) { data.frame(id = x,
                                                                           ptpdata[which(ptpdata$state %in% samp_treated_states[x]),]) },
                                            cl = n_cores))
  
  samp_treat_data = setDT(samp_treat_data)

  samp_control_data = do.call(what = "rbind",
                              args = pblapply(X = 1:n_0,
                                              FUN = function(x) { data.frame(id = (n_1 + x),
                                                                             ptpdata[which(ptpdata$state %in% samp_control_states[x]),]) },
                                              cl = n_cores))
  
  samp_data = rbind(samp_treat_data, samp_control_data)
  
  boot_point_ests = pbsapply(X = 1:length(models),
                             FUN = function(x) { eepd_sim(fits = eepd_fit(models = eepd_mod(formula_list = models[[x]]$formula,
                                                                                            lag = models[[x]]$lag,
                                                                                            diff_k = models[[x]]$diff_k,
                                                                                            log = models[[x]]$log,
                                                                                            fixef = models[[x]]$fixef),
                                                                          data = samp_data,
                                                                          group_var = "group",
                                                                          time_var = "year",
                                                                          val_times = val_years,
                                                                          post_time = treated_year,
                                                                          unit_var = "id"))$atts },
                             cl = n_cores)
  
  boot_max_abs_errors = pbsapply(X = 1:length(models),
                               FUN = function(x) { max(abs(eepd_sim(fits = eepd_fit(models = eepd_mod(formula_list = models[[x]]$formula,
                                                                                                      lag = models[[x]]$lag,
                                                                                                      diff_k = models[[x]]$diff_k,
                                                                                                      log = models[[x]]$log,
                                                                                                      fixef = models[[x]]$fixef),
                                                                                    data = samp_data,
                                                                                    group_var = "group",
                                                                                    time_var = "year",
                                                                                    val_times = val_years,
                                                                                    post_time = treated_year,
                                                                                    unit_var = "id"))$sim_pred_errors[, , 1])) },
                               cl = n_cores)
  
  exp_est[i] = sum(boot_point_ests * post_probs)
  exp_lb_est[i] = sum((boot_point_ests - (M * boot_max_abs_errors)) * post_probs)
  exp_ub_est[i] = sum((boot_point_ests + (M * boot_max_abs_errors)) * post_probs)
  
}

alpha <- 0.05
boot_var <- mean((exp_est - mean(exp_est))^2)
round(x = sqrt(mod_var + boot_var), digits = 2)
c(round(x = BMA_est - qnorm(p = 1 - (alpha/2)) * sqrt(mod_var + boot_var), digits = 2),
  round(x = BMA_est + qnorm(p = 1 - (alpha/2)) * sqrt(mod_var + boot_var), digits = 2))


BMA_lb_est <- sum((mod_point_ests - (M * max_abs_errors)) * post_probs)
mod_var_lb_est <- sum(((mod_point_ests - (M * max_abs_errors)) - BMA_lb_est)^2 * post_probs)
boot_var_lb_est <- mean((exp_lb_est - mean(exp_lb_est))^2)

BMA_lb_est - (qnorm(p = 1 - (alpha/2)) * sqrt(mod_var_lb_est + boot_var_lb_est))

# Optimal model validation plot -------------------------------------------

ptpdata <- mutate(.data = ptpdata,
                  crude_rate_lag_1 = lag(x = crude_rate, n = 1L))
training_years <- 1995:1998
years <- c(training_years, val_years, treated_year)

## Optimal model is LDV on log scale with state fixed effects
fit_opt_mod <- eepd_fit(models = eepd_mod(formula_list = models[[which.min(max_abs_errors)]]$formula,
                                          lag = models[[which.min(max_abs_errors)]]$lag,
                                          diff_k = models[[which.min(max_abs_errors)]]$diff_k,
                                          log = models[[which.min(max_abs_errors)]]$log,
                                          fixef = models[[which.min(max_abs_errors)]]$fixef),
                        data = ptpdata,
                        group_var = "group",
                        time_var = "year",
                        val_times = val_years,
                        post_time = treated_year,
                        unit_var = "state")
## Generate model matrices for all years (including training years)
model_mats <- pblapply(X = years,
                       FUN = function(x) { model.matrix(object = lm(formula = log(crude_rate) ~ log(crude_rate_lag_1) + state,
                                                                    data = ptpdata,
                                                                    subset = year == x)) },
                       cl = n_cores)
names(model_mats) <- years

coefs <- fit_opt_mod$coefs
names(coefs) <- c(val_years, treated_year)

## Use estimated coefficients from first validation year (1999) to generate predictions for all training years (1994:1998)
opt_mod_preds <- pblapply(X = names(model_mats),
                          FUN = function(x) { 
                            
                            if(x %in% paste(training_years)){
                              
                              model_mats[[x]] %*% fit_opt_mod$coefs[[1]]
                              
                            } else{ model_mats[[x]] %*% coefs[[x]] }},
                          cl = n_cores)

ptpdata <- mutate(.data = ptpdata,
                  opt_mod_pred = NA)

for(i in 1:length(opt_mod_preds)){
  
  ptpdata[row.names(opt_mod_preds[[i]]), "opt_mod_pred"] = exp(opt_mod_preds[[i]])
  
}

opt_mod_val_data <- group_by(.data = ptpdata, year, group) %>%
  summarize(Outcome = mean(crude_rate),
            Pred = mean(opt_mod_pred),
            .groups = "keep")
library(reshape2)
opt_mod_val_data <- melt(data = opt_mod_val_data,
                         id.vars = c("year","group"),
                         value.name = "avg_y",
                         variable.name = "y_type")
opt_mod_val_data$comparison <- "Treated and control groups' average prediction errors"
opt_mod_val_data$model <- "LDV (log) + unit FEs"
opt_mod_val_data <- mutate(.data = opt_mod_val_data,
                           group = factor(x = group,
                                          levels = c(0, 1),
                                          labels = c("Control (neighboring states)",
                                                     "Treated (Missouri)")))

opt_mod_val_data_alt <- data.frame(year = c(years,
                                            years),
                                   group = "Treated (Missouri)",
                                   y_type = c(rep(x = "Outcome", times = length(years)),
                                              rep(x = "Pred", times = length(years))),
                                   avg_y = c(opt_mod_val_data$avg_y[which(opt_mod_val_data$group == "Treated (Missouri)" & opt_mod_val_data$y_type == "Outcome" & opt_mod_val_data$year %in% years)],
                                             opt_mod_val_data$avg_y[which(opt_mod_val_data$group == "Treated (Missouri)" & opt_mod_val_data$y_type == "Pred" & opt_mod_val_data$year %in% years)] +
                                               (opt_mod_val_data$avg_y[which(opt_mod_val_data$group == "Control (neighboring states)" & opt_mod_val_data$y_type == "Outcome" & opt_mod_val_data$year %in% years)] -
                                                  opt_mod_val_data$avg_y[which(opt_mod_val_data$group == "Control (neighboring states)" & opt_mod_val_data$y_type == "Pred" & opt_mod_val_data$year %in% years)])),
                                   comparison = "Treated group's corrected prediction errors",
                                   model = "LDV (log) + unit FEs")

opt_mod_val_plot_data <- rbind(opt_mod_val_data, opt_mod_val_data_alt)
opt_mod_val_plot_data <- mutate(.data = opt_mod_val_plot_data,
                                comparison = factor(x = comparison,
                                                    levels = c("Treated and control groups' average prediction errors",
                                                               "Treated group's corrected prediction errors"),
                                                    labels = c("Treated and control groups' average prediction errors",
                                                               "Treated group's corrected prediction errors")))
pre_treat_years <- 1994:2007
opt_model_val_plot <- ggplot(data = opt_mod_val_plot_data,
                             mapping = aes(x = year,
                                           y = avg_y,
                                           color = group)) +
  facet_wrap(facets = . ~ comparison,
             nrow = 1,
             ncol = 2,
             scales = "fixed") +
  geom_point(data = filter(.data = opt_mod_val_plot_data,
                           y_type == "Outcome"),
             size = 1.5) +
  geom_point(data = filter(.data = opt_mod_val_plot_data,
                           y_type == "Outcome"),
             size = 1.5) +
  geom_line(data = filter(.data = opt_mod_val_plot_data,
                          y_type == "Pred" & year %in% c(val_years, min(val_years)- 1)),
            linetype = "dashed") +
  geom_line(data = filter(.data = opt_mod_val_plot_data,
                          y_type == "Pred" & year %in% (2:(min(val_years) - 1))),
            linetype = "solid") +
  geom_line(data = filter(.data = opt_mod_val_plot_data,
                          y_type == "Pred" & year %in% c(val_years, min(val_years)- 1)),
            linetype = "dashed") +
  geom_line(data = filter(.data = opt_mod_val_plot_data,
                          y_type == "Pred" & year %in% (2:(min(val_years) - 1))),
            linetype = "solid") +
  scale_x_continuous(breaks = pre_treat_years,
                     limits = c(min(pre_treat_years), max(pre_treat_years)),
                     labels = pre_treat_years) +
  #scale_y_continuous(breaks = seq(from = 3, to = 7, by = 1),
  #                   limits = c(3, 7.5),
  #                   labels = c(3:7)) +
  theme_bw() +
  labs(x = "Year",
       y = "Average Gun Homicides (rate per 100k)") +
  scale_color_manual(values = c("#332288", "#117733")) +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 8))
opt_model_val_plot
#ggsave(plot = opt_model_val_plot,
#       file = "opt_model_val_plot.pdf",
#       width = 6,
#       height = 4,
#       units = "in",
#       dpi = 600)

## Check to make sure errors the same as produced in package
abs(opt_mod_val_data_alt$avg_y[opt_mod_val_data_alt$y_type == "Outcome"] -
      opt_mod_val_data_alt$avg_y[opt_mod_val_data_alt$y_type == "Pred"])
mod_max_abs_errors[[4]]

## Unpack more why prediction error is greatest in 2005
ptpdata <- mutate(.data = ptpdata,
                  pred_error = crude_rate - opt_mod_pred)

ptpdata$state[ptpdata$year == 2005][order(abs(ptpdata$pred_error[ptpdata$year == 2005]))]
coefs[["2005"]]

state_opt_mod_val_data <- group_by(.data = ptpdata, year, state) %>%
  summarize(Outcome = mean(crude_rate),
            Pred = mean(opt_mod_pred),
            .groups = "keep")
library(reshape2)
state_opt_mod_val_data <- melt(data = state_opt_mod_val_data,
                               id.vars = c("year","state"),
                               value.name = "avg_y",
                               variable.name = "y_type")

pred_errors_by_state <- ggplot(data = state_opt_mod_val_data,
                               mapping = aes(x = year,
                                             y = avg_y)) +
  facet_wrap(facets = . ~ state,
             nrow = 3,
             ncol = 3,
             scales = "fixed") +
  geom_point(data = filter(.data = state_opt_mod_val_data,
                           y_type == "Outcome"),
             size = 1.5) +
  geom_point(data = filter(.data = state_opt_mod_val_data,
                           y_type == "Outcome"),
             size = 1.5) +
  geom_line(data = filter(.data = state_opt_mod_val_data,
                          y_type == "Pred" & year %in% c(val_years, min(val_years)- 1)),
            linetype = "dashed") +
  geom_line(data = filter(.data = state_opt_mod_val_data,
                          y_type == "Pred" & year %in% (2:(min(val_years) - 1))),
            linetype = "solid") +
  geom_line(data = filter(.data = state_opt_mod_val_data,
                          y_type == "Pred" & year %in% c(val_years, min(val_years)- 1)),
            linetype = "dashed") +
  geom_line(data = filter(.data = state_opt_mod_val_data,
                          y_type == "Pred" & year %in% (2:(min(val_years) - 1))),
            linetype = "solid") +
  scale_x_continuous(breaks = pre_treat_years,
                     limits = c(min(pre_treat_years), max(pre_treat_years)),
                     labels = pre_treat_years) +
  #scale_y_continuous(breaks = seq(from = 3, to = 7, by = 1),
  #                   limits = c(3, 7.5),
  #                   labels = c(3:7)) +
  theme_bw() +
  labs(x = "Year",
       y = "Average Gun Homicides (rate per 100k)") +
  scale_color_manual(values = c("#332288", "#117733")) +
  theme(axis.text.x = element_text(angle = 90, size = 8),
        legend.title = element_blank(),
        legend.position = "bottom",
        strip.text.x = element_text(size = 8))
#ggsave(plot = pred_errors_by_state,
#       file = "pred_errors_by_state.pdf",
#       width = 6,
#       height = 4,
#       units = "in",
#       dpi = 600)