# To do

## Major

* Use of OLS linear regression coefficient estimates and HC0 variance-covariance matrix by default; users can specify as an additional argument if they would like to use M-estimation. Much of the theory in the new manuscript is tailored to OLS regression and OLS is what most users do in practice anyway, so OLS makes more sense for the fitting of the prediction models as the default option. However, I do love having M-estimation in there as an alternative the user could specify.

NG: When I said "M-estimation" what I really meant was variance that accounts for covariance between parameters across models. We don't actually sue M-estimation; we use vcovSUEST, which is also what we would get from SUR. The manuscript mentions using the joint covariance matrix including covariances across models, so do you want that or no? Are you also saying to remove the capabilities for GLMs and just use lm() for all models? What are the two options that a user can provide? Setting all cross-model covariances to 0 vs estimating with vcovSUEST?

Only vcovSUEST; Want users to use any GLM.

* New variance estimation. This is on page 13 of the attached manuscript. This involves taking the sum of (1) variance of estimates due to model uncertainty holding the observed data fixed and (2) the variance of Bayesian model averaging estimates over bootstraps of the data holding the original posterior fixed. There is theoretical justification for this approach in this recent Biometrics article here. In simulations, we also show that it yields nominal coverage in samples that are at least moderately large. *Calculating part (1) of the variance consists of using draws from the quasi-posterior to generate simulation-based estimates of probability that each model is truly optimal model in the population, and then using this proportion of times each model is optimal to calculate the variance in the individual estimates under each model with respect to the simulation-based probabilities of each model’s being optimal*

NG: BMA variance is straightforward: weighted variance of ATTs from each model with weights equal to the proportion of times that model wins in simulation. Bootstrap component I don't understand as well: bootstrap ATTs for each model in post period, computed weighted average of ATTS in each bootstrap sample using BMA weights (fixed), compute sample variance? Why do we have to bootstrap when we we have Tau, Sigma, and p (vector of weights) so ATT = Tau p and SE = p Sigma p?

Add option to use analytic instead of bootstrap

* The overall estimated variance-covariance matrix should consist of only the validation periods. The post-treatment period should be excluded. The same for the estimated coefficients of each model. In other words, draws from the quasi-posterior should be only for the validation periods, not any other periods.

NG: This makes sense. Validation period estimation is just to get joiunt vcov for simulation, which is used to get the BMA weights (p). Totally separate from ATT.

* Adding M (the magnitude parameter) as an argument to the eepd_sim function; right now eepd_sim presumes M = 0. For M greater than 0, there should be two estimates of ATT bounds (lower and upper) and variance estimates for each of these bound estimates based on the approach in 2 above. *Note that with M > 0, draws from the quasi-posterior are used only to construct the probabilities that each model is the truly optimal one in the population. The differential prediction errors for each model and period are fixed in the sample data, i.e., not recalculated over each draw of coefficients from the quasi-posterior*

NG: I'll do this, but first want to be solid on M=0.

* Right now the use of unit fixed effects in predictions models seems to be such that, if we were to resample states from the data, then any state sampled more than once would be counted as its own unit with its own fixed effect. Yet we would like the unit fixed effects to be “state fixed effects” wherein if “Alabama” happens to show up twice in a bootstrap sample, there is only one state fixed effect that pertains to each of those two observations.

NG: I will look into this. As I understand, when we do a bootstrap, we do a weighted bootstrap, so units are weighted and fixed effects still refer to all observations with the same state. Even if we do a traditionak bootstrap, can be recast as a weighted bootstrap with integer weights.
 

## Minor

* Is it possible to adjust eepd_sim so that, if a user wants, it is easy to extract the point estimate in each time period under each model?
* Also, it would be helpful to be able to directly extract the maximum absolute prediction error for each model.
* It would also be helpful as well as the proportion of times each model is optimal over draws from the quasi-posterior.
* Can it also be made easy to extract the two parts going into the sum of the estimated variance in point 2 above?
* Finally, for the plot(est) output in the ReadMe file here, would it be possible to make the facet headings all shaded in white and only for the optimal model in the sample to be shaded in grey? I’m not sure if this is possible but it would help address a minor point one of the reviewers raised.