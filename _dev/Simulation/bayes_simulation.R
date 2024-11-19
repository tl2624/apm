# Setup -------------------------------------------------------------------
rm(list = ls())
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#library(devtools)
#install_github("ngreifer/eepd")
#set_github_pat()

# Parallelization for the loop over simulation iterations
library(doParallel)
library(foreach)

# Parallelization for the functions inside the loop
# these need to be passed to each core
library(parallel)
library(pbapply)
library(tidyverse)
library(data.table)
library(eepd)


n_cores <- 100 # this is the number of cores for the outer loop of the simulation (parallelized with foreach)
# n_cores * n_cores_inner must not exceed parallel::detectCores()
myCluster <- makeCluster(n_cores)
registerDoParallel(myCluster)

# Data cleaning ---------------------------------------------------------------
load(file = "schell_et_al_sim_data.RData")

## Calibrate years to Hasegawa et al (2019) data
years <- 1994:2008
## 2008 is only treatment year
treated_year <- 2008
val_years <- 1999:2007

## Recode data
data <- filter(.data = schell_et_al_sim_data, Year %in% years) %>%
  select(State, Year, Crude.Rate) %>%
  rename(state = State,
         year = Year,
         crude_rate = Crude.Rate) %>%
  mutate(state = as.character(state))


# Simulation setup --------------------------------------------------------

## Arbitrarily define treated and control populations of states
## analogous to the 1 treated and 8 control states in the
## Hasegawa et al (2019) data

## 5 treated states and 45 control states
N_treat <- 5
set.seed(11242017)
treated_states <- sample(x = unique(data$state), size = N_treat)
control_states <- setdiff(x = unique(data$state), treated_states)
data <- mutate(.data = data,
               group = ifelse(test = state %in% treated_states,
                              yes = 1,
                              no = 0),
               treat = ifelse(test = state %in% treated_states & year == treated_year,
                              yes = 1,
                              no = 0))
data <- setDT(x = data)

# Define candidate models ------------------------------------------------------

## Class of 6 candidate models
mods <- eepd_mod(list(crude_rate ~ 1,
                      crude_rate ~ 1 + year),
                 family = list("gaussian"),
                 lag = c(0, 1),
                 diff_k = c(0, 1),
                 fixef = FALSE)

mod_1 <- eepd_mod(list(crude_rate ~ 1), family = list("gaussian"), lag = 0, diff_k = 0, fixef = FALSE)
mod_2 <- eepd_mod(list(crude_rate ~ 1 + year), family = list("gaussian"), lag = 0, diff_k = 0, fixef = FALSE)
mod_3 <- eepd_mod(list(crude_rate ~ 1), family = list("gaussian"), lag = 1, diff_k = 0, fixef = FALSE)
mod_4 <- eepd_mod(list(crude_rate ~ 1 + year), family = list("gaussian"), lag = 1, diff_k = 0, fixef = FALSE)
mod_5 <- eepd_mod(list(crude_rate ~ 1), family = list("gaussian"), lag = 0, diff_k = 1, fixef = FALSE)
mod_6 <- eepd_mod(list(crude_rate ~ 1 + year), family = list("gaussian"), lag = 0, diff_k = 1, fixef = FALSE)

# Fit all models to population data --------------------------------------------

fits <- eepd_fit(models = mods,
                 data = data,
                 group_var = "group",
                 time_var = "year",
                 val_times = val_years,
                 post_time = treated_year,
                 unit_var = "state")

# Find optimal model for population distribution -------------------------------

pop_opt_mod <- eepd_mod(formula_list = list(crude_rate ~ 1),
                        family = list("gaussian"),
                        log = FALSE,
                        lag = 0,
                        diff_k = 0,
                        fixef = FALSE) # mods[[eepd_sim(fits = fits, nsim = 0)$optimal_models]]

pop_est_ATT <- eepd_sim(fits = fits,
                        nsim = 0)$atts

# Set parameters of simulation -------------------------------------------------

## Do simulations with 1, 3, 15 and 35 treated units
## and corresponding 8, 24, 120 and 280 control units
## (Ratio of treated to control units is fixed)

## Do simulations with M = 0, 1, 2

## Do simulations with population's optimal model and over all models
## Assess whether coverage differs between (a) known optimal model and (b) model selection

## number of simulations
n_sims <- 10^3

## sample size of treated population
n_treats <- c(1, 3, 15, 35)

## sample size of control population
n_conts <- c(8, 24, 120, 280)

## Number of draws from posterior over which model is optimal
n_draws <- 10^3

## Number of bootstrap replicates
n_boots <- 10^3

## List of vectors to store expected ATT estimate of each simulation for each sample size
#exp_ATT_ests <- replicate(n = length(n_treats),
#                          expr = rep(x = NA, times = n_sims),
#                          simplify = FALSE)

## List of vectors to store variance of ATT estimate of each simulation for each sample size
#var_ATT_ests <- replicate(n = length(n_treats),
#                          expr = rep(x = NA, times = n_sims),
#                          simplify = FALSE)

## List of matrices to store optimal model distribution of each simulation for each sample size
#opt_mods_dist <- replicate(n = length(n_treats),
#                           expr = matrix(data = NA,
#                                         nrow = n_sims,
#                                         ncol = length(mods),
#                                         dimnames = list(paste("sim", 1:n_sims, sep = "_"),
#                                                         paste("mod", 1:length(mods), sep = "_"))),
#                           simplify = FALSE)
#names(opt_mods_dist) <- paste("n_treat", n_treats, sep = "_")

## Double list of vectors to store bootstrapped expected ATT estimates of each simulation for each sample size
#boot_exp_ATT_ests <- replicate(n = length(n_treats),
#                               expr = replicate(n = n_sims,
#                                                expr = rep(x = NA, times = n_boots),
#                                                simplify = FALSE),
#                               simplify = FALSE)
#for(i in 1:length(boot_exp_ATT_ests)) { names(boot_exp_ATT_ests[[i]]) = paste("sim", 1:n_sims, sep = "_") }
#names(boot_exp_ATT_ests) <- paste("n_treat", n_treats, sep = "_")

## Function to calculate estimated ATT of each model for each bootstrap replicate
boot_ATT_est_fun <- function(treated_data,
                             control_data,
                             n_treat,
                             n_cont,
                             mod_names,
                             cl = (detectCores() - 1)){
  
  boot_samp_treats = sample(x = unique(treated_data[,id]),
                            size = n_treat,
                            replace = TRUE)
  
  boot_samp_controls = sample(x = unique(control_data[,id]),
                              size = n_cont,
                              replace = TRUE)
  
  ## bind sampled treated states
  boot_treat_data = do.call(what = "rbind",
                            args = pblapply(X = 1:n_treat,
                                            FUN = function(x) { treated_data[which(treated_data[,id] %in% boot_samp_treats[x]),] },
                                            cl = NULL))
  
  ## bind sampled control states
  boot_control_data = do.call(what = "rbind",
                              args = pblapply(X = 1:n_cont,
                                              FUN = function(x) { control_data[which(control_data[,id] %in% boot_samp_controls[x]),] },
                                              cl = NULL))
  
  ## bind sampled treated and control states
  boot_data = rbind(boot_treat_data, boot_control_data)
  
  # LAURA: Also fast, no need to parallelize
  return(list("mod_ATTs" = pbsapply(X = mod_names,
                                    FUN = function(x) { eepd_sim(fits = eepd_fit(models = get(x),
                                                                                 data = boot_data,
                                                                                 group_var = "group",
                                                                                 time_var = "year",
                                                                                 val_times = val_years,
                                                                                 post_time = treated_year,
                                                                                 unit_var = "id"),
                                                                 nsim = 0,
                                                                 cl = NULL)$atts },
                                    cl = NULL)))
  
  
}

# Do simulation ----------------------------------------------------------------

## Set seed for simulation
set.seed(09291992)
results <- list()

for(j in 1:length(n_treats)){
  
  # LAURA: parallelize this loop
  results[[j]] <- foreach(i = 1:n_sims, .packages=c('tidyverse','eepd','parallel','pbapply','data.table'),
                          .export=paste('mod',1:length(mods),sep="_") ,.inorder=FALSE) %dopar% {
                            
                            ## sample with replacement from N_treat treated states in population
                            samp_treated_states = sample(x = treated_states, size = n_treats[j], replace = TRUE)
                            
                            ## sample with replacement from control states in population
                            samp_control_states = sample(x = control_states, size = n_conts[j], replace = TRUE)
                            
                            ## bind sampled treated states
                            # LAURA: This step is cheap/fast, no need to parallelize
                            samp_treat_data = do.call(what = "rbind",
                                                      args = pblapply(X = 1:n_treats[j],
                                                                      FUN = function(x) { data.frame(id = x,
                                                                                                     data[which(data$state %in% samp_treated_states[x]),]) },
                                                                      cl = NULL))
                            
                            samp_treat_data = setDT(samp_treat_data)
                            
                            ## bind sampled control states
                            # LAURA: This step is cheap/fast, no need to parallelize
                            samp_control_data = do.call(what = "rbind",
                                                        args = pblapply(X = 1:n_conts[j],
                                                                        FUN = function(x) { data.frame(id = (n_treats[j] + x),
                                                                                                       data[which(data$state %in% samp_control_states[x]),]) },
                                                                        cl = NULL))
                            
                            samp_control_data = setDT(samp_control_data)
                            
                            ## bind sampled treated and control states
                            samp_data = rbind(samp_treat_data, samp_control_data)
                            
                            ## Get ATT estimate for each model
                            # LAURA: also cheap/fast; no need to parallelize
                            mod_ATT_ests = pbsapply(X = paste("mod", 1:6, sep = "_"),
                                                    FUN = function(x) { eepd_sim(fits = eepd_fit(models = get(x),
                                                                                                 data = samp_data,
                                                                                                 group_var = "group",
                                                                                                 time_var = "year",
                                                                                                 val_times = val_years,
                                                                                                 post_time = treated_year,
                                                                                                 unit_var = "id"),
                                                                                 nsim = 0,
                                                                                 cl = NULL)$atts },
                                                    cl = NULL)
                            
                            ## Simulation-based distribution of which model is optimal
                            # LAURA: this one is slow/expensive; its cost scales with nsim = n_draws
                            # Added optional cl argument to eepd_sim call to speed it up
                            optimal_models = eepd_sim(fits = eepd_fit(models = mods,
                                                                      data = samp_data,
                                                                      group_var = "group",
                                                                      time_var = "year",
                                                                      val_times = val_years,
                                                                      post_time = treated_year,
                                                                      unit_var = "id"),
                                                      nsim = n_draws, cl = NULL)$optimal_models
                            
                            ## Convert simulation-based distribution above to probability distribution
                            # LAURA: also cheap/fast; removed parallelization
                            opt_mods_dist = pbsapply(X = 1:length(mods),
                                                     FUN = function(x) { mean(x == optimal_models) },
                                                     cl = NULL)
                            
                            ## Calculate expected ATT
                            exp_ATT_ests = sum(mod_ATT_ests * opt_mods_dist)
                            
                            ## Calculate variance of ATT
                            var_ATT_ests = sum((mod_ATT_ests - exp_ATT_ests)^2 * opt_mods_dist)
                            
                            ## Estimate ATT for each model over n_boots bootstrap replicates
                            # LAURA: this is an expensive step that scales with n_boots
                            # Will benefit from parallelization, but on Windows,  integer number of cores is ignored
                            boot_ATT_ests = pbreplicate(n = n_boots,
                                                        expr = boot_ATT_est_fun(treated_data = samp_treat_data,
                                                                                control_data = samp_control_data,
                                                                                n_treat = n_treats[j],
                                                                                n_cont = n_conts[j],
                                                                                mod_names = paste("mod", 1:6, sep = "_"),
                                                                                cl = NULL)$mod_ATTs,
                                                        simplify = TRUE,
                                                        cl = NULL)
                            
                            ## Calculate expected ATT for each of the n_boots bootstrap replicates
                            # LAURA: fast/cheap, no need to parallelize
                            boot_exp_ATT_ests = pbsapply(X = 1:ncol(boot_ATT_ests),
                                                         FUN = function(x) { sum(boot_ATT_ests[,x] * opt_mods_dist) },
                                                         cl = NULL)
                            
                            ## save results for 100th, 200th, ... iterations
                            #if(i %in% seq(from = 100, to = n_sims, by = 100)){
                            
                            #save(samp_opt_mods, samp_ATT_ests, true_opt_mod_samp_ATT_ests,
                            #     samp_boot_opt_mods, samp_boot_ATT_ests, true_opt_mod_samp_boot_ATT_ests,
                            #     samp_LB_CI, samp_UB_CI, true_opt_mod_samp_LB_CI, true_opt_mod_samp_UB_CI,
                            #     file = "sim_res.RData")
                            
                            #} 
                            list('opt_mods_dist'=opt_mods_dist,'boot_exp_ATT_ests'=boot_exp_ATT_ests,
                                 'exp_ATT_ests'=exp_ATT_ests,'var_ATT_ests'=var_ATT_ests)
                          } # end of n_sims loop (i index)
  
} # end of n_treats loop (j index)

save(results, file = "bayes_sim_res.RData")

## Sim results stored as different list element for each sim
## Create new object of optimal model probability for each sim (under each sample size)
load(file = "bayes_sim_res.RData")
opt_mod_sim_res <- list()
for(i in 1:length(results)){
  
  opt_mod_sim_res[[i]] = sapply(X = 1:length(results[[i]]),
                                FUN = function(x) { results[[i]][[x]]$opt_mods_dist })
  
  }  

## Expected posterior probability of all 6 models (model 1 is optimal model)
apply(X = opt_mod_sim_res[[1]], MARGIN = 1, FUN = mean)
apply(X = opt_mod_sim_res[[2]], MARGIN = 1, FUN = mean)
apply(X = opt_mod_sim_res[[3]], MARGIN = 1, FUN = mean)
apply(X = opt_mod_sim_res[[4]], MARGIN = 1, FUN = mean)

## Expected ATT ests
exp_ATT_sim_res <- list()
for(i in 1:length(results)){
  
  exp_ATT_sim_res[[i]] = sapply(X = 1:length(results[[i]]),
                                FUN = function(x) { results[[i]][[x]]$exp_ATT_ests })
  
}  

## var ATT ests
var_ATT_sim_res <- list()
for(i in 1:length(results)){
  
  var_ATT_sim_res[[i]] = sapply(X = 1:length(results[[i]]),
                                FUN = function(x) { results[[i]][[x]]$var_ATT_ests })
  
}  

## Bootstrap vars
boot_var_sim_res <- list()
for(i in 1:length(results)){
  
  boot_var_sim_res[[i]] = sapply(X = 1:length(results[[i]]),
                                 FUN = function(x) { mean((results[[i]][[x]]$boot_exp_ATT_ests -
                                                             mean(results[[i]][[x]]$boot_exp_ATT_ests))^2) })
  
  }

## Bias
sapply(X = 1:length(exp_ATT_sim_res),
       FUN = function(x) { mean(exp_ATT_sim_res[[x]]) }) - pop_est_ATT

## Total estimated variance
total_var_sim_res <- lapply(X = 1:length(var_ATT_sim_res),
                            FUN = function(x) { var_ATT_sim_res[[x]] + boot_var_sim_res[[x]] })

## Coverage
CI_LB <- lapply(X = 1:length(exp_ATT_sim_res),
                FUN = function(x) { exp_ATT_sim_res[[x]] -
                    (qnorm(p = 0.975) * sqrt(total_var_sim_res[[x]])) })

CI_UB <- lapply(X = 1:length(exp_ATT_sim_res),
                FUN = function(x) { exp_ATT_sim_res[[x]] +
                    (qnorm(p = 0.975) * sqrt(total_var_sim_res[[x]])) })

coverage <- lapply(X = 1:length(CI_LB),
                   FUN = function(x) { mean(CI_LB[[x]] <= pop_est_ATT & pop_est_ATT <= CI_UB[[x]]) })

## Comparison of variances
true_vs_exp_vars <- sapply(X = 1:length(total_var_sim_res),
                           FUN = function(x) { cbind(mean(total_var_sim_res[[x]]),
                                                     mean((exp_ATT_sim_res[[x]] - mean(exp_ATT_sim_res[[x]]))^2)) })


mean((exp_ATT_sim_res[[1]] - (qnorm(p = 0.975) * sqrt(boot_var_sim_res[[1]]))) <= pop_est_ATT & 
       pop_est_ATT <= (exp_ATT_sim_res[[1]] + (qnorm(p = 0.975) * sqrt(boot_var_sim_res[[1]]))))
     




