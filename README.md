# SBWBayes

(0) The raw data is located in the repository **studensk/SBWPhenoData**; this must be installed using devtools prior to running the model.

(1) **data_clean.R** reads and cleans phenology data for input into the model. (Outputs: all_days_df.csv, all_days_df_ind.csv)

(2a) **new_post.R** contains the code that runs the strictly parametric stage-structured TMB model, which is written in regniere_structured_parametric.cpp. (Outputs: model_results.csv)

(2b) **semipar_post.R** contains the code that runs the stage-structured TMB model with non-parametric adjustment parameters, which is written in regniere_structured.cpp. (Outputs: re_structured_ON_NEW.csv)

(3) **new_sims.R** contains the code that runs weather data simulations. (Outputs: ribbon_df.csv, dev_curves.csv, df_counts.csv)

(4) **sim_plots.R** generates figures related to full model runs.

(5a) **mle_cv.R** splits data and runs cross-validation for structured, strictly parametric model without Bayesian priors (Outputs: group_df.csv)

(5b) **cross_validation_parametric.R** splits data and runs cross-validation for structured, strictly parametric, Bayesian model (Outputs: data_grouped.csv, parametric_cv*.csv, parametric_diagnostics.csv)

(6a) **cv_mle_analysis.R** aggregates and summarizes MLE cross-validation results (Outputs: mle_parametric_cv_results.csv)

(6b) **cv_analysis_parametric.R** aggregates and summarizes Bayesian cross-validation results (Outputs: parametric_draws.csv, parametric_cv_results.csv)

(7) **cv_compare.R** plots cross-validation results




