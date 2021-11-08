# SBWBayes
**new_post.R** contains the code that runs the TMB model, which is written in regniere_structured.cpp

**new_sims.R** contains the code that runs weather data simulations, and creates the figures in the main text.

**maximum_likelihood_estimation.R** contains the code that finds maximum likelihood estimates for the model parameters, estimates covariance matrices, draws from the resulting multivariate normal distribution and creates confidence intervals which are then compared to the Bayesian fits and the observed data.

The data is located in the repository **studensk/SBWPhenoData**; this must be installed using devtools prior to running the model.
